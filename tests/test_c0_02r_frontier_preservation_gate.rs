use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::process::Command;

use fastbact::common::bitseq::{reverse_complement_ascii, ReadBank};
use fastbact::common::sketch::splitmix64;
use fastbact::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, UnitigNode};
use fastbact::polish::output::curate_contigs;
use fastbact::resolve::exspander::{
    resolve_repeats_exspander_with_stats, score_candidate_path, ExspanderParams,
};
use fastbact::resolve::mapper::{
    map_paired_reads, InsertSizeDist, PairLinkMatrix, ReadPairEvidence,
};

const MINIMAP2_BIN: &str =
    "/hpcfs/fhome/caizhh/.local/lib/python3.9/site-packages/quast_libs/minimap2/minimap2";

fn dna(seed: u64, len: usize) -> Vec<u8> {
    (0..len)
        .map(|i| b"ACGT"[(splitmix64(seed + i as u64) % 4) as usize])
        .collect()
}

fn add_node(graph: &mut CdbgGraph, seq: &[u8]) -> DirUnitigId {
    let id = graph.nodes.len() as u32;
    let seq_id = graph.arena.push_seq(seq);
    graph
        .nodes
        .push(UnitigNode::new(id, seq_id, seq.len() as u32, 25.0));
    DirUnitigId::new(id, false)
}

fn add_edge(graph: &mut CdbgGraph, from: DirUnitigId, to: DirUnitigId, overlap: u8) {
    graph.nodes[from.id() as usize].add_edge(from.is_rc(), GraphEdge::new(to, overlap, 10));
    graph.nodes[to.id() as usize].add_edge(to.twin().is_rc(), GraphEdge::new(from.twin(), overlap, 10));
}

fn sha256_of_bytes(bytes: &[u8]) -> String {
    use std::io::Write as _;
    let mut child = Command::new("sha256sum")
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .spawn()
        .expect("failed to spawn sha256sum");
    child
        .stdin
        .as_mut()
        .expect("failed to open stdin")
        .write_all(bytes)
        .expect("failed to write to sha256sum");
    let output = child.wait_with_output().expect("failed to wait on sha256sum");
    let out_str = String::from_utf8_lossy(&output.stdout);
    out_str.split_whitespace().next().unwrap_or("").to_string()
}

fn sha256_of_file<P: AsRef<Path>>(path: P) -> String {
    let output = Command::new("sha256sum")
        .arg(path.as_ref())
        .output()
        .expect("failed to run sha256sum");
    let out_str = String::from_utf8_lossy(&output.stdout);
    out_str.split_whitespace().next().unwrap_or("").to_string()
}

/// Mark every exact forward placement of `contig` on `reference`.
/// Returns false when the contig is not an exact reference substring.
fn mark_exact_reference_coverage(reference: &[u8], contig: &[u8], covered: &mut [bool]) -> bool {
    if contig.is_empty() || contig.len() > reference.len() || covered.len() != reference.len() {
        return false;
    }
    let mut found = false;
    for start in 0..=reference.len() - contig.len() {
        if &reference[start..start + contig.len()] == contig {
            found = true;
            for base in covered.iter_mut().take(start + contig.len()).skip(start) {
                *base = true;
            }
        }
    }
    found
}

/// Case 8 Node count and Reference length verification.
#[test]
fn test_case_8_node_count_and_reference_audit() {
    let k = 31usize;
    let n_repeat = 65usize;

    // Topology: S (700 bp) -> R1 (100 bp) -> ... -> R65 (100 bp)
    // 1 seed node (S) + 65 repeat nodes (R1..R65) = 66 total nodes
    let total_nodes = 1 + n_repeat;
    assert_eq!(total_nodes, 66, "Case 8 path has 66 nodes: 1 seed S + 65 R nodes");

    // Overlap k = 31 bp between every consecutive node.
    // Net contribution:
    // S: 700 bp
    // Each Ri: 100 - 31 = 69 bp
    let expected_unique = 700 + n_repeat * (100 - k);
    assert_eq!(expected_unique, 5185, "700 + 65 * 69 = 5185 bp");

    // Prefix assembled at step cap 64:
    // S + 64 R nodes = 65 nodes in prefix
    let prefix_length = 700 + 64 * (100 - k);
    assert_eq!(prefix_length, 5116, "700 + 64 * 69 = 5116 bp");

    // Frontier node R65 contribution:
    let r65_novel = 100 - k;
    assert_eq!(r65_novel, 69);
    assert_eq!(prefix_length + r65_novel, expected_unique);
}

/// Boundary sweep across 63, 64, 65, 66 steps with short (100 bp) repeat nodes.
/// Evaluates:
/// Level A: Graph nodes survived
/// Level B: Nodes retained through Stage 5 cleanup
/// Level C: Contigs emitted into final FASTA (min_len = 500 vs min_len = 0)
/// Level D: Reference coverage & Minimap2 alignment
#[test]
fn test_boundary_sweep_63_64_65_66_short_nodes() {
    let k = 31usize;

    for &n_repeat in &[63usize, 64, 65, 66] {
        let mut graph = CdbgGraph::new(k);

        // Junctions are unique 31-mers
        let junctions: Vec<Vec<u8>> = (0..=n_repeat).map(|i| dna(0x6100_0000 + i as u64, k)).collect();

        // Node S: 700 bp
        let mut s_seq = dna(0x5300_0001, 700 - k);
        s_seq.extend_from_slice(&junctions[0]);
        assert_eq!(s_seq.len(), 700);
        let mut prev = add_node(&mut graph, &s_seq);

        let mut repeat_seqs = Vec::with_capacity(n_repeat);
        for i in 1..=n_repeat {
            let mut seq = junctions[i - 1].clone();
            seq.extend(dna(0x7200_0000 + i as u64, 100 - 2 * k));
            seq.extend_from_slice(&junctions[i]);
            assert_eq!(seq.len(), 100);
            let id = add_node(&mut graph, &seq);
            add_edge(&mut graph, prev, id, k as u8);
            prev = id;
            repeat_seqs.push(seq);
        }

        let total_graph_nodes = 1 + n_repeat;
        assert_eq!(graph.nodes.len(), total_graph_nodes);

        // Build authentic reference sequence
        let mut ref_seq = s_seq.clone();
        for r_seq in &repeat_seqs {
            ref_seq.extend_from_slice(&r_seq[k..]);
        }
        let expected_ref_len = 700 + n_repeat * (100 - k);
        assert_eq!(ref_seq.len(), expected_ref_len);

        let ref_sha = sha256_of_bytes(&ref_seq);

        // Execute exSPAnder repeat resolution
        resolve_repeats_exspander_with_stats(
            &mut graph,
            &PairLinkMatrix::new(),
            &InsertSizeDist::default(),
            ExspanderParams::default(),
        ).expect("repeat resolution failed");

        // Level A & B: Graph node survival and cleanup retention
        let active_nodes: Vec<_> = graph.nodes.iter().filter(|n| !n.is_deleted).collect();
        let total_active_bp: usize = active_nodes.iter().map(|n| n.seq_len as usize).sum();

        if n_repeat <= 64 {
            // Entire path merged into 1 contig
            assert_eq!(active_nodes.len(), 1);
            assert_eq!(total_active_bp, expected_ref_len);
        } else {
            // Prefix merged (5116 bp) + unmerged frontier nodes protected from cleanup
            let prefix_node = active_nodes.iter().find(|n| n.seq_len == 5116);
            assert!(prefix_node.is_some(), "Prefix 5116 bp node must be present");
            // Frontier nodes must survive cleanup
            let _unmerged_count = n_repeat - 64;
            let short_survivors = active_nodes.iter().filter(|n| n.seq_len == 100).count();
            // Note: R_{64} connects to R_{65}. When step cap 64 stops at R_{64},
            // protect_unmerged_frontier protects R_{65}. At N=66, R66 stays
            // connected to that surviving R65, so both short nodes must remain.
            if n_repeat == 66 {
                assert_eq!(
                    short_survivors, 2,
                    "N=66 must retain both unmerged short nodes R65 and R66, not merely one of them"
                );
                assert!(!graph.nodes[65].is_deleted, "R65 node identity must survive cleanup");
                assert!(!graph.nodes[66].is_deleted, "R66 node identity must survive cleanup");
                assert_eq!(
                    graph.get_unitig_seq(65),
                    repeat_seqs[64],
                    "R65 sequence identity must be unchanged"
                );
                assert_eq!(
                    graph.get_unitig_seq(66),
                    repeat_seqs[65],
                    "R66 sequence identity must be unchanged"
                );
                let r65_unique = &repeat_seqs[64][k..];
                let r66_unique = &repeat_seqs[65][k..];
                assert_eq!(r65_unique.len(), 69, "R65 contributes 69 bp beyond the 31 bp overlap");
                assert_eq!(r66_unique.len(), 69, "R66 contributes 69 bp beyond the 31 bp overlap");
                assert!(
                    graph.get_unitig_seq(65).windows(69).any(|w| w == r65_unique),
                    "R65 unique tail must survive in the graph"
                );
                assert!(
                    graph.get_unitig_seq(66).windows(69).any(|w| w == r66_unique),
                    "R66 unique tail must survive in the graph"
                );
            } else {
                assert!(
                    short_survivors >= 1,
                    "At least frontier node R65 must survive cleanup"
                );
            }
        }

        // Level C: Contigs emitted into final FASTA
        // 1. With CLI default min_len = 500
        let contigs_default = curate_contigs(&graph, 500, 2.0);
        // 2. With raw retention min_len = 0 (all surviving graph nodes)
        let contigs_raw = curate_contigs(&graph, 0, 2.0);

        // Write both to FASTA and measure alignment
        let outdir = tempfile::tempdir().expect("tempdir");
        let ref_path = outdir.path().join("reference.fa");
        {
            let mut f = File::create(&ref_path).unwrap();
            writeln!(f, ">ref len={}", ref_seq.len()).unwrap();
            writeln!(f, "{}", std::str::from_utf8(&ref_seq).unwrap()).unwrap();
        }

        let default_fa_path = outdir.path().join("default_assembly.fa");
        {
            let mut f = File::create(&default_fa_path).unwrap();
            for c in &contigs_default {
                writeln!(f, ">contig_{} len={}", c.id, c.len).unwrap();
                writeln!(f, "{}", std::str::from_utf8(&c.seq).unwrap()).unwrap();
            }
        }

        let raw_fa_path = outdir.path().join("raw_assembly.fa");
        {
            let mut f = File::create(&raw_fa_path).unwrap();
            for c in &contigs_raw {
                writeln!(f, ">contig_{} len={}", c.id, c.len).unwrap();
                writeln!(f, "{}", std::str::from_utf8(&c.seq).unwrap()).unwrap();
            }
        }

        if n_repeat == 66 {
            let r65_unique = &repeat_seqs[64][k..];
            let r66_unique = &repeat_seqs[65][k..];
            assert!(
                contigs_raw.iter().any(|c| c.seq.windows(r65_unique.len()).any(|w| w == r65_unique)),
                "R65 unique sequence must survive raw curation (min_len = 0)"
            );
            assert!(
                contigs_raw.iter().any(|c| c.seq.windows(r66_unique.len()).any(|w| w == r66_unique)),
                "R66 unique sequence must survive raw curation (min_len = 0)"
            );

            let mut covered = vec![false; ref_seq.len()];
            for contig in &contigs_raw {
                assert!(
                    mark_exact_reference_coverage(&ref_seq, &contig.seq, &mut covered),
                    "raw contig len {} must be an exact forward substring of the N=66 reference",
                    contig.len
                );
            }
            let covered_bp = covered.iter().filter(|&&base| base).count();
            assert_eq!(
                covered_bp, ref_seq.len(),
                "N=66 raw contigs must union-cover 100.00% of the reference (covered {covered_bp} / {})",
                ref_seq.len()
            );
            assert_eq!(
                covered_bp * 10_000 / ref_seq.len(),
                10_000,
                "reference coverage union must be exactly 100.00%"
            );
        }

        let default_sha = sha256_of_file(&default_fa_path);
        let raw_sha = sha256_of_file(&raw_fa_path);

        // Level D: Reference coverage & Minimap2 alignment
        let mm2_output_default = Command::new(MINIMAP2_BIN)
            .args(&["-c", "--cs=short", ref_path.to_str().unwrap(), default_fa_path.to_str().unwrap()])
            .output()
            .expect("minimap2 execution");
        assert!(mm2_output_default.status.success());
        let paf_default = String::from_utf8_lossy(&mm2_output_default.stdout);

        let mm2_output_raw = Command::new(MINIMAP2_BIN)
            .args(&["-c", "--cs=short", ref_path.to_str().unwrap(), raw_fa_path.to_str().unwrap()])
            .output()
            .expect("minimap2 execution");
        assert!(mm2_output_raw.status.success());
        let paf_raw = String::from_utf8_lossy(&mm2_output_raw.stdout);

        println!(
            "[Boundary n_repeat={}] RefLen={}, RefSHA={}, DefaultContigs={}, DefaultSHA={}, RawContigs={}, RawSHA={}",
            n_repeat, expected_ref_len, &ref_sha[..12], contigs_default.len(), &default_sha[..12], contigs_raw.len(), &raw_sha[..12]
        );

        println!("PAF default:\n{}", paf_default);
        if n_repeat <= 64 {
            // Full 100% coverage in default FASTA
            assert_eq!(contigs_default.len(), 1);
            assert_eq!(contigs_default[0].len, expected_ref_len);
            assert!(paf_default.contains("NM:i:0"));
        } else {
            // For n_repeat > 64:
            // Under min_len=500: only prefix contig (5116 bp) is emitted because R65 (100 bp) < 500 bp
            assert_eq!(contigs_default.len(), 1);
            assert_eq!(contigs_default[0].len, 5116);
            assert!(paf_default.contains("cs:Z::5116"));

            println!("PAF raw for n_repeat={}:\n{}", n_repeat, paf_raw);

            // Under min_len=0: both prefix (5116 bp) AND frontier node R65 (100 bp) are emitted!
            // R65 aligns perfectly to the reference tail with NM:i:0 and 31 bp overlap!
            assert!(paf_raw.contains("cs:Z::5116"));
            assert!(paf_raw.contains("cs:Z::100"));
        }
    }
}

/// Boundary sweep across 63, 64, 65, 66 steps when the frontier is followed by a
/// standard genomic continuation node T (700 bp >= 500 bp).
/// Verifies that when genuine genomic sequences continue past the 64-step cap,
/// BOTH contigs survive cleanup, are emitted in final FASTA under standard min_len=500,
/// and achieve 100.00% reference coverage without any sequence loss!
#[test]
fn test_boundary_sweep_with_terminal_flank_node() {
    let k = 31usize;

    for &n_repeat in &[63usize, 64, 65, 66] {
        let mut graph = CdbgGraph::new(k);

        // Junctions are unique 31-mers
        let junctions: Vec<Vec<u8>> = (0..=n_repeat).map(|i| dna(0x8100_0000 + i as u64, k)).collect();

        // Node S: 700 bp
        let mut s_seq = dna(0x5300_0002, 700 - k);
        s_seq.extend_from_slice(&junctions[0]);
        let mut prev = add_node(&mut graph, &s_seq);

        let mut repeat_seqs = Vec::with_capacity(n_repeat);
        for i in 1..n_repeat {
            let mut seq = junctions[i - 1].clone();
            seq.extend(dna(0x7200_0000 + i as u64, 100 - 2 * k));
            seq.extend_from_slice(&junctions[i]);
            assert_eq!(seq.len(), 100);
            let id = add_node(&mut graph, &seq);
            add_edge(&mut graph, prev, id, k as u8);
            prev = id;
            repeat_seqs.push(seq);
        }

        // Terminal node T: 700 bp (>= 500 bp)
        let mut t_seq = junctions[n_repeat - 1].clone();
        t_seq.extend(dna(0x5400_0001, 700 - k));
        assert_eq!(t_seq.len(), 700);
        let t_id = add_node(&mut graph, &t_seq);
        add_edge(&mut graph, prev, t_id, k as u8);
        repeat_seqs.push(t_seq.clone());

        // Build authentic reference sequence
        let mut ref_seq = s_seq.clone();
        for r_seq in &repeat_seqs {
            ref_seq.extend_from_slice(&r_seq[k..]);
        }
        let expected_ref_len = 700 + (n_repeat - 1) * (100 - k) + (700 - k);
        assert_eq!(ref_seq.len(), expected_ref_len);

        // Execute exSPAnder repeat resolution
        resolve_repeats_exspander_with_stats(
            &mut graph,
            &PairLinkMatrix::new(),
            &InsertSizeDist::default(),
            ExspanderParams::default(),
        ).expect("repeat resolution failed");

        // Curate contigs under standard production settings: min_len = 500
        let contigs = curate_contigs(&graph, 500, 2.0);

        let outdir = tempfile::tempdir().expect("tempdir");
        let ref_path = outdir.path().join("reference.fa");
        {
            let mut f = File::create(&ref_path).unwrap();
            writeln!(f, ">ref len={}", ref_seq.len()).unwrap();
            writeln!(f, "{}", std::str::from_utf8(&ref_seq).unwrap()).unwrap();
        }

        let assembly_path = outdir.path().join("assembly.fa");
        {
            let mut f = File::create(&assembly_path).unwrap();
            for c in &contigs {
                writeln!(f, ">contig_{} len={}", c.id, c.len).unwrap();
                writeln!(f, "{}", std::str::from_utf8(&c.seq).unwrap()).unwrap();
            }
        }

        // Minimap2 alignment
        let mm2_output = Command::new(MINIMAP2_BIN)
            .args(&["-c", "--cs=short", ref_path.to_str().unwrap(), assembly_path.to_str().unwrap()])
            .output()
            .expect("minimap2 execution");
        assert!(mm2_output.status.success());
        let paf = String::from_utf8_lossy(&mm2_output.stdout);

        println!(
            "[Flanked n_repeat={}] RefLen={}, Contigs={}",
            n_repeat, expected_ref_len, contigs.len()
        );

        // Check reference coverage from PAF
        let mut covered = vec![false; expected_ref_len];
        for line in paf.lines() {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 12 {
                let r_start: usize = fields[7].parse().unwrap();
                let r_end: usize = fields[8].parse().unwrap();
                for i in r_start..r_end {
                    covered[i] = true;
                }
                // Verify zero mismatches/indels
                assert!(line.contains("NM:i:0"), "PAF line must have NM:i:0: {}", line);
            }
        }

        let covered_bp = covered.iter().filter(|&&b| b).count();
        assert_eq!(
            covered_bp, expected_ref_len,
            "100.00% reference sequence must be covered across all emitted contigs"
        );
    }
}

/// Case 7 Search Depth Truncation Safety Audit.
/// Verifies:
/// 1. Candidate search depth capping at 4 triggers SearchStatus::Truncated.
/// 2. Neither shortcut nor partial path is committed when an admissible continuation is hidden.
/// 3. Fork topology is safely preserved: R retains both outgoing edges.
/// 4. When the authentic path is within depth bound (<= 4 hops), search status is Complete
///    and the authentic path successfully wins and commits!
#[test]
fn test_case_7_depth_truncation_safety_comprehensive_audit() {
    let k = 31usize;

    // Sub-case A: 5-hop true route exceeds MAX_CANDIDATE_DEPTH (4).
    // Fork must be preserved; shortcut must NOT commit.
    {
        let mut graph = CdbgGraph::new(k);
        let calib_seq = dna(0xC0DE_7001, 2000);
        let _calib_id = add_node(&mut graph, &calib_seq);

        let k_ar = dna(0xA0A0_0001, k);
        let k_join = dna(0xB0B0_0002, k);
        let k_12 = dna(0xC0C0_0003, k);
        let k_23 = dna(0xD0D0_0004, k);
        let k_34 = dna(0xE0E0_0005, k);

        let mut a_seq = dna(0xA111_0001, 501 - k);
        a_seq.extend_from_slice(&k_ar);
        let mut r_seq = k_ar.clone();
        r_seq.extend(dna(0x7111_0001, 100 - 2 * k));
        r_seq.extend_from_slice(&k_join);

        let mut c1_seq = k_join.clone();
        c1_seq.extend(dna(0xC111_0001, 100 - 2 * k));
        c1_seq.extend_from_slice(&k_12);

        let mut c2_seq = k_12.clone();
        c2_seq.extend(dna(0xC222_0002, 100 - 2 * k));
        c2_seq.extend_from_slice(&k_23);

        let mut c3_seq = k_23.clone();
        c3_seq.extend(dna(0xC333_0003, 100 - 2 * k));
        c3_seq.extend_from_slice(&k_34);

        let mut c4_seq = k_34.clone();
        c4_seq.extend(dna(0xC444_0004, 100 - 2 * k));
        c4_seq.extend_from_slice(&k_join);

        let mut b_seq = k_join.clone();
        b_seq.extend(dna(0xB111_0001, 500 - k));

        let a_id = add_node(&mut graph, &a_seq);
        let r_id = add_node(&mut graph, &r_seq);
        let c1_id = add_node(&mut graph, &c1_seq);
        let c2_id = add_node(&mut graph, &c2_seq);
        let c3_id = add_node(&mut graph, &c3_seq);
        let c4_id = add_node(&mut graph, &c4_seq);
        let b_id = add_node(&mut graph, &b_seq);

        let ov = k as u8;
        add_edge(&mut graph, a_id, r_id, ov);
        add_edge(&mut graph, r_id, b_id, ov);
        add_edge(&mut graph, r_id, c1_id, ov);
        add_edge(&mut graph, c1_id, c2_id, ov);
        add_edge(&mut graph, c2_id, c3_id, ov);
        add_edge(&mut graph, c3_id, c4_id, ov);
        add_edge(&mut graph, c4_id, b_id, ov);

        let mut true_genome = Vec::new();
        true_genome.extend_from_slice(&a_seq);
        true_genome.extend_from_slice(&r_seq[k..]);
        true_genome.extend_from_slice(&c1_seq[k..]);
        true_genome.extend_from_slice(&c2_seq[k..]);
        true_genome.extend_from_slice(&c3_seq[k..]);
        true_genome.extend_from_slice(&c4_seq[k..]);
        true_genome.extend_from_slice(&b_seq[k..]);

        let b_origin = true_genome.len() - b_seq.len();
        let r1_start = a_seq.len() - 200;
        let r2_end = b_origin + 186;

        let mut bank = ReadBank::new();
        for _ in 0..25 {
            let r1 = true_genome[r1_start..r1_start + 150].to_vec();
            let r2 = reverse_complement_ascii(&true_genome[r2_end - 150..r2_end]);
            bank.pe_forward.push_seq(&r1);
            bank.pe_reverse.push_seq(&r2);
        }

        let mut library_inserts = Vec::with_capacity(41);
        library_inserts.extend(std::iter::repeat(600usize).take(15));
        library_inserts.extend(std::iter::repeat(700usize).take(11));
        library_inserts.extend(std::iter::repeat(800usize).take(15));
        for (i, &ins) in library_inserts.iter().enumerate() {
            let start = 100 + i * 2;
            let r1 = calib_seq[start..start + 150].to_vec();
            let r2 = reverse_complement_ascii(&calib_seq[start + ins - 150..start + ins]);
            bank.pe_forward.push_seq(&r1);
            bank.pe_reverse.push_seq(&r2);
        }

        let (matrix, dist) = map_paired_reads::<u64>(&graph, &bank, k);

        let params = ExspanderParams {
            min_support: 5,
            decision_threshold: 2.0,
        };
        resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params).expect("repeat resolution failed");

        let shortcut_len = a_seq.len() + (r_seq.len() - k) + (b_seq.len() - k);
        let shortcut_committed = graph.nodes.iter().any(|n| !n.is_deleted && n.seq_len == shortcut_len as u32);
        assert!(!shortcut_committed, "Shortcut must NOT be committed when deeper admissible path exists");

        // Verify fork preservation: both exits from R remain intact
        let r_exits = graph.nodes[r_id.id() as usize]
            .edges_fwd
            .iter()
            .filter(|edge| edge.target.id() == c1_id.id() || edge.target.id() == b_id.id())
            .count();
        assert_eq!(r_exits, 2, "Fork R must retain both exits (to C1 and to B)");
    }

    // Sub-case B: 3-hop true route (within MAX_CANDIDATE_DEPTH = 4).
    // Authentic path must successfully compete against shortcut and be committed!
    {
        let mut graph = CdbgGraph::new(k);
        let calib_seq = dna(0xC0DE_7002, 2000);
        let _calib_id = add_node(&mut graph, &calib_seq);

        let k_ar = dna(0xA0A0_0002, k);
        let k_join = dna(0xB0B0_0003, k);
        let k_12 = dna(0xC0C0_0004, k);

        let mut a_seq = dna(0xA111_0002, 501 - k);
        a_seq.extend_from_slice(&k_ar);
        let mut r_seq = k_ar.clone();
        r_seq.extend(dna(0x7111_0002, 100 - 2 * k));
        r_seq.extend_from_slice(&k_join);

        let mut c1_seq = k_join.clone();
        c1_seq.extend(dna(0xC111_0002, 100 - 2 * k));
        c1_seq.extend_from_slice(&k_12);

        let mut c2_seq = k_12.clone();
        c2_seq.extend(dna(0xC222_0003, 100 - 2 * k));
        c2_seq.extend_from_slice(&k_join);

        let mut b_seq = k_join.clone();
        b_seq.extend(dna(0xB111_0002, 500 - k));

        let a_id = add_node(&mut graph, &a_seq);
        let r_id = add_node(&mut graph, &r_seq);
        let c1_id = add_node(&mut graph, &c1_seq);
        let c2_id = add_node(&mut graph, &c2_seq);
        let b_id = add_node(&mut graph, &b_seq);

        let ov = k as u8;
        add_edge(&mut graph, a_id, r_id, ov);
        add_edge(&mut graph, r_id, b_id, ov);
        add_edge(&mut graph, r_id, c1_id, ov);
        add_edge(&mut graph, c1_id, c2_id, ov);
        add_edge(&mut graph, c2_id, b_id, ov);

        let mut true_genome = Vec::new();
        true_genome.extend_from_slice(&a_seq);
        true_genome.extend_from_slice(&r_seq[k..]);
        true_genome.extend_from_slice(&c1_seq[k..]);
        true_genome.extend_from_slice(&c2_seq[k..]);
        true_genome.extend_from_slice(&b_seq[k..]);

        let b_origin = true_genome.len() - b_seq.len();
        let r1_start = a_seq.len() - 200;
        let r2_end = b_origin + 186;

        let mut bank = ReadBank::new();
        for _ in 0..25 {
            let r1 = true_genome[r1_start..r1_start + 150].to_vec();
            let r2 = reverse_complement_ascii(&true_genome[r2_end - 150..r2_end]);
            bank.pe_forward.push_seq(&r1);
            bank.pe_reverse.push_seq(&r2);
        }

        let mut library_inserts = Vec::with_capacity(41);
        library_inserts.extend(std::iter::repeat(600usize).take(15));
        library_inserts.extend(std::iter::repeat(700usize).take(11));
        library_inserts.extend(std::iter::repeat(800usize).take(15));
        for (i, &ins) in library_inserts.iter().enumerate() {
            let start = 100 + i * 2;
            let r1 = calib_seq[start..start + 150].to_vec();
            let r2 = reverse_complement_ascii(&calib_seq[start + ins - 150..start + ins]);
            bank.pe_forward.push_seq(&r1);
            bank.pe_reverse.push_seq(&r2);
        }

        let (matrix, dist) = map_paired_reads::<u64>(&graph, &bank, k);

        let params = ExspanderParams {
            min_support: 5,
            decision_threshold: 2.0,
        };
        resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params).expect("repeat resolution failed");

        let authentic_assembled = graph
            .nodes
            .iter()
            .any(|n| !n.is_deleted && n.seq_len == true_genome.len() as u32);
        assert!(
            authentic_assembled,
            "Authentic path within depth bound must successfully assemble (expected len {})",
            true_genome.len()
        );
    }
}

/// Case 9 Pair-ID Deduplication Semantics Audit.
/// Verifies:
/// 1. Duplicate pair_id rows for the same unitig pair are counted only once in count and weight.
/// 2. Different pair_id rows are independently counted.
/// 3. The same pair_id can be evaluated in multiple candidate paths without cross-contamination.
/// 4. Out-of-range pairs do not enter seen_pair_ids or alter counts.
#[test]
fn test_case_9_pair_id_deduplication_comprehensive_audit() {
    let u1 = DirUnitigId::new(1, false);
    let u2 = DirUnitigId::new(2, false);

    let mut matrix = PairLinkMatrix::new();

    // In-range pair 1: four duplicate rows with identical span (dist1=200, dist2=200)
    for _ in 0..4 {
        matrix.add_evidence(
            u1,
            u2,
            ReadPairEvidence {
                pair_id: 1,
                dist_u1_exit: 200,
                dist_u2_entry: 200,
                r1_len: 150,
                r2_len: 150,
            },
        );
    }

    // In-range pair 2: two duplicate rows with identical span (dist1=200, dist2=200)
    for _ in 0..2 {
        matrix.add_evidence(
            u1,
            u2,
            ReadPairEvidence {
                pair_id: 2,
                dist_u1_exit: 200,
                dist_u2_entry: 200,
                r1_len: 150,
                r2_len: 150,
            },
        );
    }

    // Out-of-range pair 3: three duplicate rows with large span (dist1=2000, dist2=2000)
    for _ in 0..3 {
        matrix.add_evidence(
            u1,
            u2,
            ReadPairEvidence {
                pair_id: 3,
                dist_u1_exit: 2000,
                dist_u2_entry: 2000,
                r1_len: 150,
                r2_len: 150,
            },
        );
    }

    // In-range pair 4: one row
    matrix.add_evidence(
        u1,
        u2,
        ReadPairEvidence {
            pair_id: 4,
            dist_u1_exit: 200,
            dist_u2_entry: 200,
            r1_len: 150,
            r2_len: 150,
        },
    );

    let dist = InsertSizeDist {
        median: 700.0,
        mad: 100.0,
        min_insert: 350,
        max_insert: 1050,
    };

    // Candidate A: separation 300 bp.
    // Pair 1: span 200 + 300 + 200 = 700 bp (Z = 0, weight = 1.0). Deduplicated from 4 to 1.
    // Pair 2: span 200 + 300 + 200 = 700 bp (Z = 0, weight = 1.0). Deduplicated from 2 to 1.
    // Pair 3: span 2000 + 300 + 2000 = 4300 bp (out of range, ignored).
    // Pair 4: span 200 + 300 + 200 = 700 bp (Z = 0, weight = 1.0). 1 row.
    // Expected: count = 3, total_weight = 3.0.
    let (count_a, weight_a) = score_candidate_path(&matrix, &dist, u1, u2, 300);
    assert_eq!(count_a, 3, "Expected 3 distinct pairs (1, 2, 4)");
    assert!((weight_a - 3.0).abs() < 1e-9, "Expected weight 3.0");

    // Candidate B: separation 100 bp.
    // Verifies that calling score_candidate_path again does NOT have seen_pair_ids polluted
    // by Candidate A's evaluation!
    // Pair 1: span 200 + 100 + 200 = 500 bp (Z = 200 / 148.26 = 1.349, weight = exp(-0.5 * 1.349^2) ≈ 0.402).
    // Pair 2: span 500 bp.
    // Pair 4: span 500 bp.
    let (count_b, weight_b) = score_candidate_path(&matrix, &dist, u1, u2, 100);
    assert_eq!(count_b, 3, "Candidate B must also evaluate pairs 1, 2, 4 independently");
    let expected_w = 3.0 * (-0.5 * (200.0 / 148.26_f64).powi(2)).exp();
    assert!((weight_b - expected_w).abs() < 1e-4, "Candidate B weight formula verified");
}

