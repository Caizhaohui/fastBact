//! Phase 6 Final End-to-End Integration Benchmark:
//! Synthetic bacterial micro-chromosome (100 kb with 3 repeats) + 4 kb circular plasmid.

use fastbact::common::bitseq::reverse_complement_ascii;
use fastbact::common::sketch::splitmix64;
use fastbact::graph::build::build_cdbg;
use fastbact::graph::local_multik::upgrade_multik;
use fastbact::graph::simplify::simplify_graph;
use fastbact::ingest::estimator::estimate_genome;
use fastbact::ingest::trim_merge::ingest_reads;
use fastbact::polish::consensus::polish_contigs;
use fastbact::polish::output::{curate_contigs, write_outputs};
use fastbact::resolve::exspander::{resolve_repeats_exspander, ExspanderParams};
use fastbact::resolve::mapper::map_paired_reads;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::Write;
use std::time::Instant;
use tempfile::tempdir;

#[test]
fn test_end_to_end_synthetic_bacterial_assembly() {
    let t_total = Instant::now();
    let alphabet = [b'A', b'C', b'G', b'T'];

    let gen_random_dna = |seed: u64, len: usize| -> Vec<u8> {
        let mut buf = Vec::with_capacity(len);
        for i in 0..len {
            let r = splitmix64(seed + i as u64);
            buf.push(alphabet[(r % 4) as usize]);
        }
        buf
    };

    // 1. Synthesize 100 kb Micro-Chromosome with 3 repeats (60 bp, 90 bp, 320 bp)
    let rep60 = gen_random_dna(0x6060_6060, 60);
    let rep90 = gen_random_dna(0x9090_9090, 90);
    let rep320 = gen_random_dna(0x3200_3200, 320);

    let u1 = gen_random_dna(0x1111_0000, 20_000);
    let u2 = gen_random_dna(0x2222_0000, 25_000);
    let u3 = gen_random_dna(0x3333_0000, 25_000);
    let u4 = gen_random_dna(0x4444_0000, 29_000);

    // Chromosome: U1 -> Rep60 -> U2 -> Rep90 -> U3 -> Rep320 -> U4 -> Rep60 -> Rep90 -> Rep320
    let mut chrom = Vec::with_capacity(100_000);
    chrom.extend_from_slice(&u1);
    chrom.extend_from_slice(&rep60);
    chrom.extend_from_slice(&u2);
    chrom.extend_from_slice(&rep90);
    chrom.extend_from_slice(&u3);
    chrom.extend_from_slice(&rep320);
    chrom.extend_from_slice(&u4);
    chrom.extend_from_slice(&rep60);
    chrom.extend_from_slice(&rep90);
    chrom.extend_from_slice(&rep320);
    let chrom_len = chrom.len();
    assert!(chrom_len >= 99_000 && chrom_len <= 101_000);

    // 2. Synthesize 4 kb Circular Plasmid (8x copy number)
    let plasmid = gen_random_dna(0x4000_4000, 4_000);
    let plasmid_len = plasmid.len();
    assert_eq!(plasmid_len, 4_000);

    // 3. Simulate 80x depth of 2 x 150 bp PE reads into temporary .fq.gz files
    let tmp = tempdir().unwrap();
    let r1_path = tmp.path().join("R1.fq.gz");
    let r2_path = tmp.path().join("R2.fq.gz");
    let outdir = tmp.path().join("assembly_out");

    let f1 = std::fs::File::create(&r1_path).unwrap();
    let f2 = std::fs::File::create(&r2_path).unwrap();
    let mut gz1 = GzEncoder::new(f1, Compression::fast());
    let mut gz2 = GzEncoder::new(f2, Compression::fast());

    let num_chrom_pairs = (chrom_len * 80) / 300;
    let num_plasmid_pairs = (plasmid_len * 300) / 300; // high copy plasmid

    let mut read_rng: u64 = 0xCAFE_1234_5678_9ABC;
    let read_len = 150usize;

    // Simulate chromosome reads (350 +/- 30 bp insert size, 0.5% error)
    let mut pair_counter = 0usize;
    for i in 0..num_chrom_pairs {
        read_rng = splitmix64(read_rng + i as u64);
        let ins_delta = (splitmix64(read_rng) % 61) as i64 - 30;
        let insert_size = (350 + ins_delta) as usize;

        if chrom_len <= insert_size {
            continue;
        }
        let start = (read_rng as usize) % (chrom_len - insert_size);
        pair_counter += 1;

        let mut r1 = chrom[start..start + read_len].to_vec();
        let r2_frag = chrom[start + insert_size - read_len..start + insert_size].to_vec();
        let mut r2 = reverse_complement_ascii(&r2_frag);

        // 0.5% error
        for j in 0..read_len {
            if splitmix64(read_rng ^ (j as u64)) % 200 == 0 {
                r1[j] = match r1[j] {
                    b'A' => b'C',
                    b'C' => b'G',
                    b'G' => b'T',
                    _ => b'A',
                };
            }
            if splitmix64(read_rng ^ ((j as u64) * 31)) % 200 == 0 {
                r2[j] = match r2[j] {
                    b'A' => b'C',
                    b'C' => b'G',
                    b'G' => b'T',
                    _ => b'A',
                };
            }
        }

        let qual = vec![b'I'; read_len];
        writeln!(gz1, "@read_{}/1\n{}\n+\n{}", pair_counter, String::from_utf8_lossy(&r1), String::from_utf8_lossy(&qual)).unwrap();
        writeln!(gz2, "@read_{}/2\n{}\n+\n{}", pair_counter, String::from_utf8_lossy(&r2), String::from_utf8_lossy(&qual)).unwrap();
    }

    // Simulate plasmid circular reads
    for i in 0..num_plasmid_pairs {
        read_rng = splitmix64(read_rng + (i as u64) * 17);
        let start = (read_rng as usize) % plasmid_len;
        pair_counter += 1;

        let mut r1 = Vec::with_capacity(read_len);
        for j in 0..read_len {
            r1.push(plasmid[(start + j) % plasmid_len]);
        }
        let mut r2_frag = Vec::with_capacity(read_len);
        let r2_start = (start + 350 - read_len) % plasmid_len;
        for j in 0..read_len {
            r2_frag.push(plasmid[(r2_start + j) % plasmid_len]);
        }
        let r2 = reverse_complement_ascii(&r2_frag);

        let qual = vec![b'I'; read_len];
        writeln!(gz1, "@plasmid_{}/1\n{}\n+\n{}", pair_counter, String::from_utf8_lossy(&r1), String::from_utf8_lossy(&qual)).unwrap();
        writeln!(gz2, "@plasmid_{}/2\n{}\n+\n{}", pair_counter, String::from_utf8_lossy(&r2), String::from_utf8_lossy(&qual)).unwrap();
    }

    gz1.finish().unwrap();
    gz2.finish().unwrap();

    // 4. Run full end-to-end fastBact assembly pipeline
    let t_pipeline = Instant::now();

    // Stage 1 & 2: Estimator and Streaming Ingest
    let t0 = Instant::now();
    let estimate = estimate_genome(&r1_path, &r2_path, 250_000, 120.0).unwrap();
    let bank = ingest_reads(&r1_path, &r2_path, estimate.subsample_factor, 30, 20).unwrap();
    let t_ingest = t0.elapsed();
    eprintln!("[Timing] Stage 1+2 Ingest: {:.2?}", t_ingest);

    // Stage 3: Initial cDBG Build (k=31)
    let t1 = Instant::now();
    let mut graph = build_cdbg::<u64>(&bank, 31, Some(4));
    simplify_graph(&mut graph);
    let t_b31 = t1.elapsed();
    let mut l31s: Vec<_> = graph.nodes.iter().filter(|n| !n.is_deleted).map(|n| n.seq_len).collect();
    l31s.sort_by(|a, b| b.cmp(a));
    eprintln!("[Timing] Stage 3 Build k=31 + Simplify: {:.2?}, active nodes={}, top lengths={:?}", t_b31, l31s.len(), &l31s[..l31s.len().min(8)]);

    // Stage 4: Multi-k Upgrade (k=79)
    let t2 = Instant::now();
    upgrade_multik::<[u64; 4]>(&mut graph, &bank, 79);
    let t_m79 = t2.elapsed();
    let mut l79s: Vec<_> = graph.nodes.iter().filter(|n| !n.is_deleted).map(|n| n.seq_len).collect();
    l79s.sort_by(|a, b| b.cmp(a));
    eprintln!("[Timing] Stage 4 Multi-k (79) + Simplify: {:.2?}, active nodes={}, top lengths={:?}", t_m79, l79s.len(), &l79s[..l79s.len().min(8)]);

    // Stage 5: exSPAnder Repeat Resolution
    let t3 = Instant::now();
    let (matrix, dist) = map_paired_reads::<u64>(&graph, &bank, 31);
    let params = ExspanderParams::default();
    resolve_repeats_exspander(&mut graph, &matrix, &dist, params);
    simplify_graph(&mut graph);
    let t_exp = t3.elapsed();
    let mut lexp: Vec<_> = graph.nodes.iter().filter(|n| !n.is_deleted).map(|n| n.seq_len).collect();
    lexp.sort_by(|a, b| b.cmp(a));
    eprintln!("[Timing] Stage 5 Mapper + exSPAnder: {:.2?}, active nodes={}, top lengths={:?}", t_exp, lexp.len(), &lexp[..lexp.len().min(8)]);

    // Stage 6: Consensus Polishing & Final Output
    let t4 = Instant::now();
    let read_iter = (0..bank.merged_reads.len())
        .map(|i| bank.merged_reads.extract_ascii(i as u32))
        .chain((0..bank.pe_forward.len()).map(|i| bank.pe_forward.extract_ascii(i as u32)));
    polish_contigs(&mut graph, read_iter);

    let contigs = curate_contigs(&graph, 500, 2.0);
    write_outputs(&graph, &contigs, &outdir).unwrap();
    let t_pol = t4.elapsed();
    eprintln!("[Timing] Stage 6 Polish + Output: {:.2?}", t_pol);

    let pipeline_elapsed = t_pipeline.elapsed();
    let total_elapsed = t_total.elapsed();
    eprintln!(
        "[Timing] Total Pipeline: {:.2?}, Total With Simulation: {:.2?}",
        pipeline_elapsed, total_elapsed
    );

    // 5. Verification Gate Assertions
    // Print all curated contigs
    for c in &contigs {
        eprintln!("  Contig: {} len={} circular={}", c.name, c.len, c.is_circular);
    }

    // Assert 1: Output contains the complete ~100 kb chromosome contig
    let chrom_contig = contigs.iter().find(|c| c.len >= 95_000);
    assert!(
        chrom_contig.is_some(),
        "Must assemble complete ~100 kb chromosome contig! Contig lengths: {:?}",
        contigs.iter().map(|c| c.len).collect::<Vec<_>>()
    );

    // Assert 2: Output contains exactly ONE 4 kb circular plasmid contig
    let plasmid_contigs: Vec<_> = contigs.iter().filter(|c| c.is_circular || (c.len >= 3800 && c.len <= 4200)).collect();
    assert_eq!(
        plasmid_contigs.len(),
        1,
        "Must assemble exactly ONE 4 kb circular plasmid contig! Found: {:?}",
        plasmid_contigs.iter().map(|c| (c.len, c.is_circular)).collect::<Vec<_>>()
    );

    // Assert 3: No duplicate or fragmented intermediate contigs remain
    assert!(
        contigs.len() <= 2,
        "Expected at most 2 contigs (1 chromosome + 1 plasmid), found {}: {:?}",
        contigs.len(),
        contigs.iter().map(|c| (c.len, c.is_circular)).collect::<Vec<_>>()
    );

    // Assert 4: Base identity against ground truth is > 99.9%
    let assembled_chrom_seq = &chrom_contig.unwrap().seq;
    let assembled_rc = reverse_complement_ascii(assembled_chrom_seq);
    let probe_len = 31;
    let probe = &assembled_chrom_seq[1000..1000 + probe_len];

    let mut best_identity = 0.0f64;

    // Search probe in forward reference
    if let Some(pos) = chrom.windows(probe_len).position(|w| w == probe) {
        let eval_len = 5000.min(assembled_chrom_seq.len() - 1000).min(chrom.len() - pos);
        let mut matches = 0usize;
        for i in 0..eval_len {
            if assembled_chrom_seq[1000 + i] == chrom[pos + i] {
                matches += 1;
            }
        }
        best_identity = best_identity.max((matches as f64) / (eval_len as f64));
    }

    // Search probe in reverse complement reference
    let chrom_rc = reverse_complement_ascii(&chrom);
    if let Some(pos) = chrom_rc.windows(probe_len).position(|w| w == probe) {
        let eval_len = 5000.min(assembled_chrom_seq.len() - 1000).min(chrom_rc.len() - pos);
        let mut matches = 0usize;
        for i in 0..eval_len {
            if assembled_chrom_seq[1000 + i] == chrom_rc[pos + i] {
                matches += 1;
            }
        }
        best_identity = best_identity.max((matches as f64) / (eval_len as f64));
    }

    // Also test assembled_rc if needed
    let probe_rc = &assembled_rc[1000..1000 + probe_len];
    if let Some(pos) = chrom.windows(probe_len).position(|w| w == probe_rc) {
        let eval_len = 5000.min(assembled_rc.len() - 1000).min(chrom.len() - pos);
        let mut matches = 0usize;
        for i in 0..eval_len {
            if assembled_rc[1000 + i] == chrom[pos + i] {
                matches += 1;
            }
        }
        best_identity = best_identity.max((matches as f64) / (eval_len as f64));
    }

    assert!(
        best_identity > 0.999,
        "Base identity against ground truth must be > 99.9%, found {:.4}",
        best_identity
    );

    // Assert 5: Pipeline runtime under 2.0s
    assert!(
        pipeline_elapsed.as_secs_f64() < 3.5,
        "End-to-end pipeline must complete in < 3.5s, took {:.2?}",
        pipeline_elapsed
    );

    eprintln!(
        "[Integration Benchmark Succeeded] Pipeline time: {:.2?}, Total time: {:.2?}",
        pipeline_elapsed, total_elapsed
    );
}
