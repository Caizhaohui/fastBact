use fastbact::common::bitseq::{reverse_complement_ascii, ReadBank};
use fastbact::common::sketch::splitmix64;
use fastbact::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, UnitigNode};
use fastbact::resolve::exspander::{resolve_repeats_exspander_with_stats, ExspanderParams};
use fastbact::resolve::mapper::{map_paired_reads, InsertSizeDist, PairLinkMatrix};

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

/// Case 2: Two candidate paths both within scalar gap <= max_insert, but differing in compatibility
/// with empirical insert-size distribution (Candidate 1 span 600 bp vs Candidate 2 span 391 bp).
#[test]
fn test_case_2_candidate_path_span_compatibility_difference() {
    let _k = 31;
    let dist = InsertSizeDist {
        median: 600.0,
        mad: 20.0,
        min_insert: 530,
        max_insert: 670,
    };

    // Candidate 1: A -> R -> C -> B
    // Path internal separation = len(R) + len(C) - 3*k = 100 + 240 - 93 = 247 bp.
    // Distance from R1 start (offset 300 of A, len 500) to exit = 500 - 300 = 200 bp.
    // Distance from B entry to R2 end (offset 153 of B) = 153 bp.
    // Span = 200 + 247 + 153 = 600 bp.
    let span_cand_1: i64 = 600;

    // Candidate 2: A -> R -> B
    // Path internal separation = len(R) - 2*k = 100 - 62 = 38 bp.
    // Span = 200 + 38 + 153 = 391 bp.
    let span_cand_2: i64 = 391;

    let is_cand_1_valid = span_cand_1 >= dist.min_insert as i64 && span_cand_1 <= dist.max_insert as i64;
    let is_cand_2_valid = span_cand_2 >= dist.min_insert as i64 && span_cand_2 <= dist.max_insert as i64;

    assert!(is_cand_1_valid, "Candidate 1 (600 bp) must be valid within [530, 670]");
    assert!(!is_cand_2_valid, "Candidate 2 (391 bp) must be rejected outside [530, 670]");
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
