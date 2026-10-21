use fastbact::common::bitseq::{reverse_complement_ascii, ReadBank};
use fastbact::common::sketch::splitmix64;
use fastbact::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, UnitigNode};
use fastbact::resolve::exspander::{
    resolve_repeats_exspander_with_stats, score_candidate_path, ExspanderParams,
};
use fastbact::resolve::mapper::{map_paired_reads, InsertSizeDist, PairLinkMatrix, ReadPairEvidence};

fn make_deterministic_dna(seed: u64, len: usize) -> Vec<u8> {
    let alphabet = [b'A', b'C', b'G', b'T'];
    (0..len)
        .map(|i| alphabet[(splitmix64(seed.wrapping_add(i as u64)) % 4) as usize])
        .collect()
}

/// Case 1: True Fragment Span 600 bp vs Shortcut 391 bp.
/// Library insert size is [530, 670] bp.
/// Expected behavior: exSPAnder MUST NOT commit the 391 bp shortcut (len 1038 bp).
/// It must either assemble the true path (len 1247 bp) or safely refuse to splice.
/// ON UNPATCHED PRODUCTION CODE: This test FAILS (RED) because the shortcut is falsely committed.
#[test]
fn test_case_1_shortcut_must_not_be_selected() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);

    // Calibration unitig of 3000 bp
    let calib_seq = make_deterministic_dna(0x1111_0001, 3000);
    let calib_sid = graph.arena.push_seq(&calib_seq);
    graph
        .nodes
        .push(UnitigNode::new(0, calib_sid, calib_seq.len() as u32, 30.0));

    let k1 = make_deterministic_dna(0xAAAA_0001, k);
    let k2 = make_deterministic_dna(0xBBBB_0002, k);

    let mut a_seq = make_deterministic_dna(0xA111_0001, 500 - k);
    a_seq.extend_from_slice(&k1); // len 500, ends with k1

    let mut r_seq = k1.clone();
    r_seq.extend_from_slice(&make_deterministic_dna(0x7111_0001, 100 - 2 * k));
    r_seq.extend_from_slice(&k2); // len 100, starts with k1, ends with k2

    let mut c_seq = k2.clone();
    c_seq.extend_from_slice(&make_deterministic_dna(0xC111_0001, 240 - 2 * k));
    c_seq.extend_from_slice(&k2); // len 240, starts with k2, ends with k2

    let mut b_seq = k2.clone();
    b_seq.extend_from_slice(&make_deterministic_dna(0xB111_0001, 500 - k)); // len 500, starts with k2

    let a_id = 1u32;
    let r_id = 2u32;
    let c_id = 3u32;
    let b_id = 4u32;

    let a_sid = graph.arena.push_seq(&a_seq);
    let r_sid = graph.arena.push_seq(&r_seq);
    let c_sid = graph.arena.push_seq(&c_seq);
    let b_sid = graph.arena.push_seq(&b_seq);

    graph
        .nodes
        .push(UnitigNode::new(a_id, a_sid, a_seq.len() as u32, 30.0));
    graph
        .nodes
        .push(UnitigNode::new(r_id, r_sid, r_seq.len() as u32, 30.0));
    graph
        .nodes
        .push(UnitigNode::new(c_id, c_sid, c_seq.len() as u32, 30.0));
    graph
        .nodes
        .push(UnitigNode::new(b_id, b_sid, b_seq.len() as u32, 30.0));

    // Edges with reciprocal edges
    // A -> R
    graph.nodes[a_id as usize].add_edge(false, GraphEdge::new(DirUnitigId::new(r_id, false), k as u8, 10));
    graph.nodes[r_id as usize].add_edge(true, GraphEdge::new(DirUnitigId::new(a_id, true), k as u8, 10));

    // R -> B (erroneous shortcut)
    graph.nodes[r_id as usize].add_edge(false, GraphEdge::new(DirUnitigId::new(b_id, false), k as u8, 10));
    graph.nodes[b_id as usize].add_edge(true, GraphEdge::new(DirUnitigId::new(r_id, true), k as u8, 10));

    // R -> C (true path step 1)
    graph.nodes[r_id as usize].add_edge(false, GraphEdge::new(DirUnitigId::new(c_id, false), k as u8, 10));
    graph.nodes[c_id as usize].add_edge(true, GraphEdge::new(DirUnitigId::new(r_id, true), k as u8, 10));

    // C -> B (true path step 2)
    graph.nodes[c_id as usize].add_edge(false, GraphEdge::new(DirUnitigId::new(b_id, false), k as u8, 10));
    graph.nodes[b_id as usize].add_edge(true, GraphEdge::new(DirUnitigId::new(c_id, true), k as u8, 10));

    // True genome sequence: A -> R -> C -> B
    let mut true_genome = Vec::new();
    true_genome.extend_from_slice(&a_seq);
    true_genome.extend_from_slice(&r_seq[k..]);
    true_genome.extend_from_slice(&c_seq[k..]);
    true_genome.extend_from_slice(&b_seq[k..]);

    let mut bank = ReadBank::new();

    // 30 calibration read pairs on calib_seq (insert size 600 bp)
    for i in 0..30 {
        let start = 200 + i * 10;
        let r1 = calib_seq[start..start + 150].to_vec();
        let r2 = reverse_complement_ascii(&calib_seq[start + 600 - 150..start + 600]);
        bank.pe_forward.push_seq(&r1);
        bank.pe_reverse.push_seq(&r2);
    }

    // 25 read pairs from true_genome spanning A -> R -> C -> B (span = 600 bp)
    for i in 0..25 {
        let frag_start = 300 + i * 2;
        let frag_end = frag_start + 600;
        let r1 = true_genome[frag_start..frag_start + 150].to_vec();
        let r2 = reverse_complement_ascii(&true_genome[frag_end - 150..frag_end]);
        bank.pe_forward.push_seq(&r1);
        bank.pe_reverse.push_seq(&r2);
    }

    let (matrix, dist) = map_paired_reads::<u64>(&graph, &bank, k);
    let params = ExspanderParams {
        min_support: 5,
        decision_threshold: 2.0,
    };
    resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params);

    // REGRESSION GATE ASSERTION:
    // A correct algorithm must NOT select the 391 bp shortcut!
    // Expected shortcut length is: |A| + |R| - 31 + |B| - 31 = 500 + 69 + 469 = 1038 bp.
    let shortcut_node = graph.nodes.iter().find(|n| !n.is_deleted && n.id != 0 && n.seq_len == 1038);
    assert!(
        shortcut_node.is_none(),
        "REGRESSION FAILURE (P1): exSPAnder falsely committed shortcut A->R->B (len 1038 bp) \
         even though physical fragment span was 600 bp and shortcut was only 391 bp!"
    );
}

/// Case 2: Candidate path span statistical compatibility scoring (intra-range competition).
/// Two candidate paths both within scalar range [530, 670] bp, but differing in compatibility
/// with empirical insert-size distribution (Candidate 1 span 600 bp vs Candidate 2 span 650 bp).
/// Candidate 1 (Z = 0, weight = 1.0) strictly dominates Candidate 2 (Z = 1.686, weight = 0.241),
/// achieving ratio = 4.14 >= 2.5.
#[test]
fn test_case_2_candidate_path_span_compatibility_difference() {
    let dist = InsertSizeDist {
        median: 600.0,
        mad: 20.0,
        min_insert: 530,
        max_insert: 670,
    };

    let u1 = DirUnitigId::new(1, false);
    let u2 = DirUnitigId::new(2, false);

    // Build matrix with 10 read pairs having d1 = 200, d2 = 153
    let mut matrix = PairLinkMatrix::new();
    for i in 0..10 {
        matrix.add_evidence(
            u1,
            u2,
            ReadPairEvidence {
                pair_id: i,
                dist_u1_exit: 200,
                dist_u2_entry: 153,
                r1_len: 150,
                r2_len: 150,
            },
        );
    }

    // Candidate 1: path_separation = 247 bp -> Span = 200 + 247 + 153 = 600 bp (Z = 0.0)
    let (count1, weight1) = score_candidate_path(&matrix, &dist, u1, u2, 247);
    assert_eq!(count1, 10, "Candidate 1 must have 10 valid pairs");
    assert!(
        (weight1 - 10.0).abs() < 1e-4,
        "Candidate 1 weight must be ~10.0 (w=1.0 per pair)"
    );

    // Candidate 2: path_separation = 297 bp -> Span = 200 + 297 + 153 = 650 bp (Z = 1.686)
    let (count2, weight2) = score_candidate_path(&matrix, &dist, u1, u2, 297);
    assert_eq!(
        count2, 10,
        "Candidate 2 must have 10 valid pairs (within [530, 670])"
    );
    let expected_z = 50.0 / (1.4826 * 20.0);
    let expected_w = (-0.5f64 * expected_z * expected_z).exp() * 10.0;
    assert!(
        (weight2 - expected_w).abs() < 1e-4,
        "Candidate 2 weight must match Gaussian kernel"
    );

    // Decision ratio test:
    let ratio = weight1 / weight2;
    assert!(
        ratio >= 2.5,
        "Candidate 1 must dominate Candidate 2 with ratio >= 2.5, got {:.3}",
        ratio
    );
    assert!(
        (ratio - 4.14).abs() < 0.1,
        "Ratio should be approximately 4.14"
    );
}

/// Case 3: Two-copy repeat (A1, A2 -> R -> B1, B2) with insufficient support -> fork preserved.
#[test]
fn test_case_3_two_copy_repeat_insufficient_support_preserves_branch() {
    let k = 30; // overlap
    let ov_a = make_deterministic_dna(0x100, k);
    let ov_b = make_deterministic_dna(0x200, k);
    let mid = make_deterministic_dna(0x300, 240);

    let mut a1 = make_deterministic_dna(0x400, 600);
    a1.extend_from_slice(&ov_a);
    let mut a2 = make_deterministic_dna(0x500, 600);
    a2.extend_from_slice(&ov_a);

    let mut r = ov_a.clone();
    r.extend_from_slice(&mid);
    r.extend_from_slice(&ov_b);

    let mut b1 = ov_b.clone();
    b1.extend_from_slice(&make_deterministic_dna(0x600, 600));
    let mut b2 = ov_b.clone();
    b2.extend_from_slice(&make_deterministic_dna(0x700, 600));

    let mut graph = CdbgGraph::new(31);
    let na1 = DirUnitigId::new(graph.nodes.len() as u32, false);
    let sid_a1 = graph.arena.push_seq(&a1);
    graph.nodes.push(UnitigNode::new(na1.id(), sid_a1, a1.len() as u32, 25.0));

    let na2 = DirUnitigId::new(graph.nodes.len() as u32, false);
    let sid_a2 = graph.arena.push_seq(&a2);
    graph.nodes.push(UnitigNode::new(na2.id(), sid_a2, a2.len() as u32, 25.0));

    let nr = DirUnitigId::new(graph.nodes.len() as u32, false);
    let sid_r = graph.arena.push_seq(&r);
    graph.nodes.push(UnitigNode::new(nr.id(), sid_r, r.len() as u32, 25.0));

    let nb1 = DirUnitigId::new(graph.nodes.len() as u32, false);
    let sid_b1 = graph.arena.push_seq(&b1);
    graph.nodes.push(UnitigNode::new(nb1.id(), sid_b1, b1.len() as u32, 25.0));

    let nb2 = DirUnitigId::new(graph.nodes.len() as u32, false);
    let sid_b2 = graph.arena.push_seq(&b2);
    graph.nodes.push(UnitigNode::new(nb2.id(), sid_b2, b2.len() as u32, 25.0));

    // Connect edges
    graph.nodes[na1.id() as usize].add_edge(false, GraphEdge::new(nr, k as u8, 10));
    graph.nodes[nr.id() as usize].add_edge(true, GraphEdge::new(na1.twin(), k as u8, 10));

    graph.nodes[na2.id() as usize].add_edge(false, GraphEdge::new(nr, k as u8, 10));
    graph.nodes[nr.id() as usize].add_edge(true, GraphEdge::new(na2.twin(), k as u8, 10));

    graph.nodes[nr.id() as usize].add_edge(false, GraphEdge::new(nb1, k as u8, 10));
    graph.nodes[nb1.id() as usize].add_edge(true, GraphEdge::new(nr.twin(), k as u8, 10));

    graph.nodes[nr.id() as usize].add_edge(false, GraphEdge::new(nb2, k as u8, 10));
    graph.nodes[nb2.id() as usize].add_edge(true, GraphEdge::new(nr.twin(), k as u8, 10));

    // PairLinkMatrix has only 2 pairs for A1 -> B1 (below min_support = 5)
    let mut matrix = PairLinkMatrix::new();
    matrix.add_link(na1, nb1);
    matrix.add_link(nb1.twin(), na1.twin());
    matrix.add_link(na1, nb1);
    matrix.add_link(nb1.twin(), na1.twin());

    let params = ExspanderParams {
        min_support: 5,
        decision_threshold: 2.0,
    };
    let stats = resolve_repeats_exspander_with_stats(&mut graph, &matrix, &InsertSizeDist::default(), params);

    // Assert that insufficient support prevents resolution: 0 joins, all 5 nodes remain active
    assert_eq!(stats.valid_k_joins + stats.valid_k_minus_1_joins + stats.valid_other_joins, 0);
    assert_eq!(graph.nodes.len(), 5, "No merged node may be created");
    assert!(graph.nodes.iter().all(|n| !n.is_deleted), "All 5 nodes must remain intact");
}

/// Case 4: Reverse/twin orientation must provide exact reciprocal support without double-counting.
#[test]
fn test_case_4_twin_orientation_support_reciprocity() {
    let u0 = DirUnitigId::new(0, false);
    let u1 = DirUnitigId::new(1, false);
    let mut matrix = PairLinkMatrix::new();

    // Add 10 directed pairs linking u0 (fwd) to u1 (fwd on genome path)
    for _ in 0..10 {
        matrix.add_link(u0, u1);
        matrix.add_link(u1.twin(), u0.twin());
    }

    assert_eq!(matrix.get_support(u0, u1), 10, "Direct link support must be exactly 10");
    assert_eq!(matrix.get_support(u1.twin(), u0.twin()), 10, "Reciprocal twin support must be exactly 10");
    assert_eq!(matrix.get_support(u0, u1.twin()), 0, "Wrong orientation must have 0 support");
    assert_eq!(matrix.get_support(u0.twin(), u1), 0, "Wrong orientation must have 0 support");
}

/// Case 5: Large coordinates (> 65,535 bp) and long unitigs must not overflow.
#[test]
fn test_case_5_large_coordinates_no_overflow() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);

    // Create a 80,000 bp unitig (> 65,535, overflowing u16)
    let long_seq = make_deterministic_dna(0x5555_0001, 80_000);
    let long_sid = graph.arena.push_seq(&long_seq);
    graph.nodes.push(UnitigNode::new(0, long_sid, long_seq.len() as u32, 30.0));

    let mut bank = ReadBank::new();

    // 25 pairs near the 5' end of the 80k unitig (offset 100 to 550, span 450 bp)
    for i in 0..25 {
        let start = 100 + i * 2;
        let r1 = long_seq[start..start + 150].to_vec();
        let r2 = reverse_complement_ascii(&long_seq[start + 450 - 150..start + 450]);
        bank.pe_forward.push_seq(&r1);
        bank.pe_reverse.push_seq(&r2);
    }

    let (_matrix, dist) = map_paired_reads::<u64>(&graph, &bank, k);

    assert_eq!(
        dist.median.round() as usize,
        450,
        "Long unitig insert size must accurately estimate 450 bp without overflow"
    );
}

/// Case 6: Short Intermediate Unitig C (100 bp) -> Shortcut Span 531 bp within [530, 670] (F1-P1-001).
/// True Fragment Span = 600 bp (Z = 0, weight = 1.0 per pair, total 25.0).
/// Shortcut Fragment Span = 531 bp (Z = 2.327, weight = 0.0667 per pair, total 1.668).
/// On unpatched code, Candidate C has 0 direct pairs for (A, C), so Candidate B (shortcut)
/// is falsely committed via `second_count == 0`!
/// Expected behavior: exSPAnder MUST NOT commit the 531 bp shortcut (len 1038 bp).
/// It must evaluate complete candidate routes and assemble the authentic path A -> R -> C -> B (len 1107 bp).
#[test]
fn test_case_6_short_c_shortcut_must_not_be_selected_f1_p1_001() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);

    // Calibration unitig of 3000 bp
    let calib_seq = make_deterministic_dna(0x1111_0001, 3000);
    let calib_sid = graph.arena.push_seq(&calib_seq);
    graph
        .nodes
        .push(UnitigNode::new(0, calib_sid, calib_seq.len() as u32, 30.0));

    let k1 = make_deterministic_dna(0xAAAA_0001, k);
    let k2 = make_deterministic_dna(0xBBBB_0002, k);

    let mut a_seq = make_deterministic_dna(0xA111_0001, 500 - k);
    a_seq.extend_from_slice(&k1); // len 500, ends with k1

    let mut r_seq = k1.clone();
    r_seq.extend_from_slice(&make_deterministic_dna(0x7111_0001, 100 - 2 * k));
    r_seq.extend_from_slice(&k2); // len 100, starts with k1, ends with k2

    // Intermediate unitig C is 100 bp (instead of 240 bp)
    let mut c_seq = k2.clone();
    c_seq.extend_from_slice(&make_deterministic_dna(0xC111_0001, 100 - 2 * k));
    c_seq.extend_from_slice(&k2); // len 100, starts with k2, ends with k2

    let mut b_seq = k2.clone();
    b_seq.extend_from_slice(&make_deterministic_dna(0xB111_0001, 500 - k)); // len 500, starts with k2

    let a_id = 1u32;
    let r_id = 2u32;
    let c_id = 3u32;
    let b_id = 4u32;

    let a_sid = graph.arena.push_seq(&a_seq);
    let r_sid = graph.arena.push_seq(&r_seq);
    let c_sid = graph.arena.push_seq(&c_seq);
    let b_sid = graph.arena.push_seq(&b_seq);

    graph.nodes.push(UnitigNode::new(a_id, a_sid, a_seq.len() as u32, 30.0));
    graph.nodes.push(UnitigNode::new(r_id, r_sid, r_seq.len() as u32, 30.0));
    graph.nodes.push(UnitigNode::new(c_id, c_sid, c_seq.len() as u32, 30.0));
    graph.nodes.push(UnitigNode::new(b_id, b_sid, b_seq.len() as u32, 30.0));

    // Edges with reciprocal edges
    // A -> R
    graph.nodes[a_id as usize].add_edge(false, GraphEdge::new(DirUnitigId::new(r_id, false), k as u8, 10));
    graph.nodes[r_id as usize].add_edge(true, GraphEdge::new(DirUnitigId::new(a_id, true), k as u8, 10));

    // R -> B (erroneous shortcut)
    graph.nodes[r_id as usize].add_edge(false, GraphEdge::new(DirUnitigId::new(b_id, false), k as u8, 10));
    graph.nodes[b_id as usize].add_edge(true, GraphEdge::new(DirUnitigId::new(r_id, true), k as u8, 10));

    // R -> C (true path step 1)
    graph.nodes[r_id as usize].add_edge(false, GraphEdge::new(DirUnitigId::new(c_id, false), k as u8, 10));
    graph.nodes[c_id as usize].add_edge(true, GraphEdge::new(DirUnitigId::new(r_id, true), k as u8, 10));

    // C -> B (true path step 2)
    graph.nodes[c_id as usize].add_edge(false, GraphEdge::new(DirUnitigId::new(b_id, false), k as u8, 10));
    graph.nodes[b_id as usize].add_edge(true, GraphEdge::new(DirUnitigId::new(c_id, true), k as u8, 10));

    // True genome sequence: A -> R -> C -> B (len 500 + 69 + 69 + 469 = 1107 bp)
    let mut true_genome = Vec::new();
    true_genome.extend_from_slice(&a_seq);
    true_genome.extend_from_slice(&r_seq[k..]);
    true_genome.extend_from_slice(&c_seq[k..]);
    true_genome.extend_from_slice(&b_seq[k..]);

    let mut bank = ReadBank::new();

    // 30 calibration read pairs on calib_seq (insert size 600 bp)
    for i in 0..30 {
        let start = 200 + i * 10;
        let r1 = calib_seq[start..start + 150].to_vec();
        let r2 = reverse_complement_ascii(&calib_seq[start + 600 - 150..start + 600]);
        bank.pe_forward.push_seq(&r1);
        bank.pe_reverse.push_seq(&r2);
    }

    // 25 read pairs from true_genome spanning A -> R -> C -> B (span = 600 bp)
    for i in 0..25 {
        let frag_start = 300 + i * 2;
        let frag_end = frag_start + 600;
        let r1 = true_genome[frag_start..frag_start + 150].to_vec();
        let r2 = reverse_complement_ascii(&true_genome[frag_end - 150..frag_end]);
        bank.pe_forward.push_seq(&r1);
        bank.pe_reverse.push_seq(&r2);
    }

    let (matrix, dist) = map_paired_reads::<u64>(&graph, &bank, k);
    let params = ExspanderParams {
        min_support: 5,
        decision_threshold: 2.0,
    };
    resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params);

    // REGRESSION GATE ASSERTION:
    // 1. Erroneous shortcut A -> R -> B (len 1038 bp) must NOT be committed!
    let shortcut_node = graph.nodes.iter().find(|n| !n.is_deleted && n.id != 0 && n.seq_len == 1038);
    assert!(
        shortcut_node.is_none(),
        "REGRESSION FAILURE (F1-P1-001): exSPAnder falsely committed shortcut A->R->B (len 1038 bp) \
         because shortcut span 531 bp was within [530, 670] and second_count was 0!"
    );

    // 2. The authentic path A -> R -> C -> B (len 1107 bp) must be assembled!
    let true_node = graph.nodes.iter().find(|n| !n.is_deleted && n.id != 0 && n.seq_len == 1107);
    assert!(
        true_node.is_some(),
        "Production exSPAnder must assemble true path A->R->C->B (len 1107 bp) when supported by 25 read pairs!"
    );
}

fn push_node(graph: &mut CdbgGraph, seq: &[u8]) -> u32 {
    let id = graph.nodes.len() as u32;
    let seq_id = graph.arena.push_seq(seq);
    graph
        .nodes
        .push(UnitigNode::new(id, seq_id, seq.len() as u32, 30.0));
    id
}

fn add_reciprocal(graph: &mut CdbgGraph, from: u32, to: u32, overlap: u8) {
    graph.nodes[from as usize].add_edge(
        false,
        GraphEdge::new(DirUnitigId::new(to, false), overlap, 10),
    );
    graph.nodes[to as usize].add_edge(
        true,
        GraphEdge::new(DirUnitigId::new(from, true), overlap, 10),
    );
}

/// Case 7: Depth-5 authentic route dropped at MAX_CANDIDATE_DEPTH (F3-P1-001).
///
/// Topology, overlap 31 bp:
///   true route  A(501) -> R(100) -> C1(100) -> C2(100) -> C3(100) -> C4(100) -> B(500)
///   shortcut    A(501) -> R(100) -> B(500)
///
/// Candidate search starts at fork R. Reaching B on the true route takes 5 steps
/// (C1, C2, C3, C4, B). The walk records C4 at steps.len() == 4 and does not expand
/// further, and it does not set SearchStatus::Truncated. Only the shortcut remains,
/// so second_count == 0.
///
/// Library: median 700 bp, MAD 100 bp, insert interval [350, 1050].
/// True span = 200 + 314 + 186 = 700 bp, weight = 25.0.
/// Shortcut span = 200 + 38 + 186 = 424 bp, Z ≈ 1.86, weight ≈ 4.42, W/C ≈ 0.177 >= 0.05.
/// Head-to-head, the true route would win. Unpatched code never sees it and commits
/// the shortcut contig of length 501 + 69 + 469 = 1039 bp.
///
/// Oracle: that 1039 bp shortcut is absent, and either the 1315 bp authentic contig
/// is assembled or the original fork nodes and both exits from R stay in the graph.
///
/// ON UNPATCHED PRODUCTION CODE: This test FAILS (RED) at the shortcut assertion.
#[test]
fn test_case_7_depth_limit_must_not_commit_shortcut_f3_p1_001() {
    let k = 31usize;
    let mut graph = CdbgGraph::new(k);

    // Length 2000 keeps every k-mer indexed (total k-mers <= 2000) and is long
    // enough (>= 1500) for insert-size estimation.
    let calib_seq = make_deterministic_dna(0xC0DE_7001, 2000);
    let _calib_id = push_node(&mut graph, &calib_seq);

    let k_ar = make_deterministic_dna(0xA0A0_0001, k);
    let k_join = make_deterministic_dna(0xB0B0_0002, k);
    let k_12 = make_deterministic_dna(0xC0C0_0003, k);
    let k_23 = make_deterministic_dna(0xD0D0_0004, k);
    let k_34 = make_deterministic_dna(0xE0E0_0005, k);
    let junctions = [&k_ar, &k_join, &k_12, &k_23, &k_34];
    for i in 0..junctions.len() {
        for j in (i + 1)..junctions.len() {
            assert_ne!(junctions[i], junctions[j], "junction k-mers must be distinct");
        }
    }

    let mut a_seq = make_deterministic_dna(0xA111_0001, 501 - k);
    a_seq.extend_from_slice(&k_ar);
    let mut r_seq = k_ar.clone();
    r_seq.extend(make_deterministic_dna(0x7111_0001, 100 - 2 * k));
    r_seq.extend_from_slice(&k_join);
    let mut c1_seq = k_join.clone();
    c1_seq.extend(make_deterministic_dna(0xC111_0001, 100 - 2 * k));
    c1_seq.extend_from_slice(&k_12);
    let mut c2_seq = k_12.clone();
    c2_seq.extend(make_deterministic_dna(0xC222_0002, 100 - 2 * k));
    c2_seq.extend_from_slice(&k_23);
    let mut c3_seq = k_23.clone();
    c3_seq.extend(make_deterministic_dna(0xC333_0003, 100 - 2 * k));
    c3_seq.extend_from_slice(&k_34);
    let mut c4_seq = k_34.clone();
    c4_seq.extend(make_deterministic_dna(0xC444_0004, 100 - 2 * k));
    c4_seq.extend_from_slice(&k_join);
    let mut b_seq = k_join.clone();
    b_seq.extend(make_deterministic_dna(0xB111_0001, 500 - k));

    assert_eq!(a_seq.len(), 501);
    assert_eq!(r_seq.len(), 100);
    assert_eq!(c1_seq.len(), 100);
    assert_eq!(c2_seq.len(), 100);
    assert_eq!(c3_seq.len(), 100);
    assert_eq!(c4_seq.len(), 100);
    assert_eq!(b_seq.len(), 500);

    let a_id = push_node(&mut graph, &a_seq);
    let r_id = push_node(&mut graph, &r_seq);
    let c1_id = push_node(&mut graph, &c1_seq);
    let c2_id = push_node(&mut graph, &c2_seq);
    let c3_id = push_node(&mut graph, &c3_seq);
    let c4_id = push_node(&mut graph, &c4_seq);
    let b_id = push_node(&mut graph, &b_seq);
    assert_eq!((a_id, r_id, c1_id, c2_id, c3_id, c4_id, b_id), (1, 2, 3, 4, 5, 6, 7));

    let overlap = k as u8;
    add_reciprocal(&mut graph, a_id, r_id, overlap);
    add_reciprocal(&mut graph, r_id, b_id, overlap);
    add_reciprocal(&mut graph, r_id, c1_id, overlap);
    add_reciprocal(&mut graph, c1_id, c2_id, overlap);
    add_reciprocal(&mut graph, c2_id, c3_id, overlap);
    add_reciprocal(&mut graph, c3_id, c4_id, overlap);
    add_reciprocal(&mut graph, c4_id, b_id, overlap);

    let mut true_genome = Vec::new();
    true_genome.extend_from_slice(&a_seq);
    true_genome.extend_from_slice(&r_seq[k..]);
    true_genome.extend_from_slice(&c1_seq[k..]);
    true_genome.extend_from_slice(&c2_seq[k..]);
    true_genome.extend_from_slice(&c3_seq[k..]);
    true_genome.extend_from_slice(&c4_seq[k..]);
    true_genome.extend_from_slice(&b_seq[k..]);

    let b_origin = true_genome.len() - b_seq.len();
    assert_eq!(&true_genome[b_origin..], b_seq.as_slice(), "B must occupy a known genome interval");
    assert_eq!(b_origin - a_seq.len(), 314, "true A-exit to B-entry separation");
    assert_eq!(r_seq.len() - 2 * k, 38, "shortcut separation across R");
    let shortcut_len = a_seq.len() + (r_seq.len() - k) + (b_seq.len() - k);
    assert_eq!(shortcut_len, 1039);
    assert_eq!(true_genome.len(), 1315);

    let r1_start = a_seq.len() - 200;
    let r2_end = b_origin + 186;
    assert_eq!(r2_end - r1_start, 700, "true fragment span");
    assert!(r1_start + 150 <= a_seq.len() - k, "R1 must stay inside A's unique prefix");
    assert!(r2_end - 150 >= b_origin + k, "R2 must stay inside B's unique suffix");
    assert!(r2_end <= true_genome.len());

    let mut bank = ReadBank::new();

    // 15 x 600, 11 x 700, 15 x 800 -> median 700, MAD 100, interval [350, 1050].
    let mut library_inserts = Vec::with_capacity(41);
    library_inserts.extend(std::iter::repeat(600usize).take(15));
    library_inserts.extend(std::iter::repeat(700usize).take(11));
    library_inserts.extend(std::iter::repeat(800usize).take(15));
    assert_eq!(library_inserts.len(), 41);
    for (i, &ins) in library_inserts.iter().enumerate() {
        let start = 100 + i * 2;
        assert!(start + ins <= calib_seq.len());
        let r1 = calib_seq[start..start + 150].to_vec();
        let r2 = reverse_complement_ascii(&calib_seq[start + ins - 150..start + ins]);
        bank.pe_forward.push_seq(&r1);
        bank.pe_reverse.push_seq(&r2);
    }

    for _ in 0..25 {
        let r1 = true_genome[r1_start..r1_start + 150].to_vec();
        let r2 = reverse_complement_ascii(&true_genome[r2_end - 150..r2_end]);
        bank.pe_forward.push_seq(&r1);
        bank.pe_reverse.push_seq(&r2);
    }

    let (matrix, dist) = map_paired_reads::<u64>(&graph, &bank, k);
    assert!(
        (dist.median - 700.0).abs() < 1e-9,
        "precondition: median {}, want 700",
        dist.median
    );
    assert!(
        (dist.mad - 100.0).abs() < 1e-9,
        "precondition: MAD {}, want 100",
        dist.mad
    );
    assert_eq!(dist.min_insert, 350);
    assert_eq!(dist.max_insert, 1050);

    let a_dir = DirUnitigId::new(a_id, false);
    let b_dir = DirUnitigId::new(b_id, false);
    let evidences = matrix.get_evidence(a_dir, b_dir).expect("A->B evidence");
    assert!(
        evidences.len() == 25
            && evidences
                .iter()
                .all(|ev| ev.dist_u1_exit == 200 && ev.dist_u2_entry == 186),
        "precondition: 25 pairs at (dist_u1_exit, dist_u2_entry) = (200, 186), got {:?}",
        evidences
            .iter()
            .map(|ev| (ev.pair_id, ev.dist_u1_exit, ev.dist_u2_entry))
            .collect::<Vec<_>>()
    );
    for cid in [c1_id, c2_id, c3_id, c4_id] {
        assert_eq!(
            matrix.get_support(a_dir, DirUnitigId::new(cid, false)),
            0,
            "precondition: intermediate node {} must not be a paired-end endpoint",
            cid
        );
    }

    let (true_count, true_weight) = score_candidate_path(&matrix, &dist, a_dir, b_dir, 314);
    assert_eq!(true_count, 25);
    assert!(
        (true_weight - 25.0).abs() < 1e-6,
        "true-route weight {}, want 25",
        true_weight
    );

    let (short_count, short_weight) = score_candidate_path(&matrix, &dist, a_dir, b_dir, 38);
    assert_eq!(short_count, 25);
    let z = (424.0_f64 - dist.median).abs() / (1.4826 * dist.mad);
    assert!((z - 1.86).abs() < 0.01, "shortcut Z {}, want ~1.86", z);
    let expected_short = 25.0 * (-0.5 * z * z).exp();
    assert!(
        (short_weight - expected_short).abs() < 1e-4,
        "shortcut weight {}, formula {}",
        short_weight,
        expected_short
    );
    let mean_w = short_weight / short_count as f64;
    assert!(
        (mean_w - 0.177).abs() < 0.002 && mean_w >= 0.05,
        "shortcut mean weight {}, want ~0.177 and >= 0.05",
        mean_w
    );
    assert!(
        true_weight / short_weight >= 2.5,
        "true route must dominate the shortcut when both are scored, ratio {}",
        true_weight / short_weight
    );

    let params = ExspanderParams {
        min_support: 5,
        decision_threshold: 2.0,
    };
    resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params);

    let shortcut_committed = graph
        .nodes
        .iter()
        .any(|n| !n.is_deleted && n.seq_len == shortcut_len as u32);
    assert!(
        !shortcut_committed,
        "REGRESSION FAILURE (F3-P1-001): exSPAnder falsely committed shortcut A->R->B (len {} bp). \
         The authentic exit R->C1->C2->C3->C4->B is 5 steps past fork R, MAX_CANDIDATE_DEPTH is 4, \
         and the truncated branch is not marked Truncated, so second_count==0 accepts the shortcut \
         (span 424 bp, mean weight >= 0.05) without competing against the span-700 route.",
        shortcut_len
    );

    let authentic_assembled = graph
        .nodes
        .iter()
        .any(|n| !n.is_deleted && n.seq_len == true_genome.len() as u32);
    let fork_ids = [a_id, r_id, c1_id, c2_id, c3_id, c4_id, b_id];
    let nodes_intact = fork_ids.iter().all(|&id| {
        graph
            .nodes
            .get(id as usize)
            .is_some_and(|n| !n.is_deleted && n.id == id)
    });
    let both_exits_kept = graph.nodes[r_id as usize]
        .edges_fwd
        .iter()
        .filter(|edge| edge.target.id() == c1_id || edge.target.id() == b_id)
        .count()
        == 2;
    assert!(
        authentic_assembled || (nodes_intact && both_exits_kept),
        "F3-P1-001: depth exhaustion must assemble the authentic {} bp route or preserve the fork. \
         authentic_assembled={}, nodes_intact={}, both_exits_kept={}",
        true_genome.len(),
        authentic_assembled,
        nodes_intact,
        both_exits_kept
    );
}
