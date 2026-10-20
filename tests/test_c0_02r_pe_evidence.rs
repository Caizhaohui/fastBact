use fastbact::common::sketch::splitmix64;
use fastbact::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, UnitigNode};
use fastbact::resolve::exspander::{resolve_repeats_exspander_with_stats, ExspanderParams};
use fastbact::resolve::mapper::{InsertSizeDist, PairLinkMatrix};

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

fn link_both_strands(matrix: &mut PairLinkMatrix, from: DirUnitigId, to: DirUnitigId, n: usize) {
    for _ in 0..n {
        matrix.add_link(from, to);
        matrix.add_link(to.twin(), from.twin());
    }
}

/// Two-copy repeat R: A1,A2 -> R -> B1,B2. Only the A1/B1 pair has read-pair support.
fn two_copy_repeat(mid_len: usize) -> (CdbgGraph, PairLinkMatrix, Vec<u8>) {
    let ov_a = dna(0x100, 30);
    let ov_b = dna(0x200, 30);
    let mid = dna(0x300, mid_len);

    let pa1 = dna(0x400, 670);
    let mut a1 = pa1.clone();
    a1.extend_from_slice(&ov_a);
    let mut a2 = dna(0x500, 620);
    a2.extend_from_slice(&ov_a);
    let mut r = ov_a.clone();
    r.extend_from_slice(&mid);
    r.extend_from_slice(&ov_b);
    let mut b1 = ov_b.clone();
    b1.extend_from_slice(&dna(0x600, 570));
    let mut b2 = ov_b.clone();
    b2.extend_from_slice(&dna(0x700, 520));

    let mut expected = pa1;
    expected.extend_from_slice(&ov_a);
    expected.extend_from_slice(&mid);
    expected.extend_from_slice(&b1);

    let mut graph = CdbgGraph::new(31);
    let na1 = add_node(&mut graph, &a1);
    let na2 = add_node(&mut graph, &a2);
    let nr = add_node(&mut graph, &r);
    let nb1 = add_node(&mut graph, &b1);
    let nb2 = add_node(&mut graph, &b2);
    add_edge(&mut graph, na1, nr, 30);
    add_edge(&mut graph, na2, nr, 30);
    add_edge(&mut graph, nr, nb1, 30);
    add_edge(&mut graph, nr, nb2, 30);

    let mut matrix = PairLinkMatrix::new();
    link_both_strands(&mut matrix, na1, nb1, 10);
    (graph, matrix, expected)
}

#[test]
fn link_support_is_exact_per_directed_pair() {
    let u0 = DirUnitigId::new(0, false);
    let u1 = DirUnitigId::new(1, false);
    let mut matrix = PairLinkMatrix::new();
    for _ in 0..10 {
        matrix.add_link(u0, u1);
    }
    assert_eq!(matrix.get_support(u0, u1), 10);
    assert_eq!(matrix.get_support(u0, u1.twin()), 0);
    assert_eq!(matrix.get_support(u0.twin(), u1), 0);
    assert_eq!(matrix.get_support(u0.twin(), u1.twin()), 0);
}

#[test]
fn repeat_wider_than_insert_size_is_not_resolved_by_pair_links() {
    let (mut graph, matrix, _) = two_copy_repeat(1440);
    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    );
    assert_eq!(stats.valid_k_joins + stats.valid_k_minus_1_joins + stats.valid_other_joins, 0);
    assert_eq!(graph.nodes.len(), 5, "no merged node may be created");
    assert!(graph.nodes.iter().all(|n| !n.is_deleted));
}

#[test]
fn repeat_within_insert_size_is_resolved_with_exact_sequence() {
    let (mut graph, matrix, expected) = two_copy_repeat(240);
    resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    );
    let merged: Vec<_> = graph
        .nodes
        .iter()
        .filter(|n| !n.is_deleted && n.seq_len as usize == expected.len())
        .collect();
    assert_eq!(merged.len(), 1);
    assert_eq!(graph.get_unitig_seq(merged[0].id), expected);
}

#[test]
fn test_short_repeat_cycle_terminates_without_looping_f1_p2_003() {
    let mut graph = CdbgGraph::new(31);
    let ov = dna(0x999, 30);

    // S is a long seed unitig (600 bp >= 500 bp)
    let mut s_seq = dna(0x101, 570);
    s_seq.extend_from_slice(&ov);
    let s_id = add_node(&mut graph, &s_seq);

    // R1 is a short repeat (100 bp < 500 bp)
    let mut r1_seq = ov.clone();
    r1_seq.extend_from_slice(&dna(0x202, 40));
    r1_seq.extend_from_slice(&ov);
    let r1_id = add_node(&mut graph, &r1_seq);

    // R2 is a short repeat (100 bp < 500 bp)
    let mut r2_seq = ov.clone();
    r2_seq.extend_from_slice(&dna(0x303, 40));
    r2_seq.extend_from_slice(&ov);
    let r2_id = add_node(&mut graph, &r2_seq);

    // Edges forming S -> R1 -> R2 -> R1 (unbranched cycle of short nodes)
    add_edge(&mut graph, s_id, r1_id, 30);
    add_edge(&mut graph, r1_id, r2_id, 30);
    add_edge(&mut graph, r2_id, r1_id, 30);

    let matrix = PairLinkMatrix::new();
    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    );

    // The repeat cycle must terminate deterministically without infinite loop!
    // Path validation may reject splicing or terminate cleanly.
    // Graph must remain in a valid state.
    assert!(
        stats.valid_k_joins + stats.valid_k_minus_1_joins + stats.valid_other_joins <= 2,
        "Short repeat cycle must not produce infinite joins"
    );
}
