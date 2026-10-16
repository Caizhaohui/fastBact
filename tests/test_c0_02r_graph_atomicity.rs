//! FB-P0-02R / C0-02R: Graph Mutation Safety and Atomicity Tests (G1–G6).
//!
//! Validates:
//! G1: Verification failure leaves full graph snapshot identical (nodes, edges, seqs, deleted flags)
//! G2: Path rejection does not pollute used-state for other legal paths
//! G3: Consecutive candidate paths: failure of one does not affect another legal path
//! G4: Short nodes belonging to rejected paths are protected from cleanup deletion
//! G5: Isolated short nodes are deleted according to documented policy
//! G6: Successful commit has no abnormal deletion, duplicate creation, or lost nodes

use fastbact::common::bitseq::reverse_complement_ascii;
use fastbact::common::sketch::splitmix64;
use fastbact::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, UnitigNode};
use fastbact::resolve::exspander::{
    resolve_repeats_exspander_with_stats, ExspanderParams,
};
use fastbact::resolve::mapper::{InsertSizeDist, PairLinkMatrix};

fn make_deterministic_dna(seed: u64, len: usize) -> Vec<u8> {
    let alphabet = [b'A', b'C', b'G', b'T'];
    let mut buf = Vec::with_capacity(len);
    for i in 0..len {
        let r = splitmix64(seed + i as u64);
        buf.push(alphabet[(r % 4) as usize]);
    }
    assert_ne!(buf, reverse_complement_ascii(&buf));
    buf
}

fn create_unitig(graph: &mut CdbgGraph, id: u32, seq: &[u8], cov: f32) -> DirUnitigId {
    let seq_id = graph.arena.push_seq(seq);
    let node = UnitigNode::new(id, seq_id, seq.len() as u32, cov);
    if id as usize == graph.nodes.len() {
        graph.nodes.push(node);
    } else if (id as usize) < graph.nodes.len() {
        graph.nodes[id as usize] = node;
    } else {
        panic!("Non-contiguous node id");
    }
    DirUnitigId::new(id, false)
}

fn add_bidirected_edge(
    graph: &mut CdbgGraph,
    from: DirUnitigId,
    to: DirUnitigId,
    overlap: u8,
    support: u16,
) {
    graph.nodes[from.id() as usize].add_edge(from.is_rc(), GraphEdge::new(to, overlap, support));
    graph.nodes[to.id() as usize].add_edge(
        to.twin().is_rc(),
        GraphEdge::new(from.twin(), overlap, support),
    );
}

#[derive(Clone, Debug, PartialEq)]
struct GraphSnapshot {
    nodes_count: usize,
    node_data: Vec<(u32, u32, f32, bool, bool, Vec<GraphEdge>, Vec<GraphEdge>, Vec<u8>)>,
}

fn capture_graph_snapshot(graph: &CdbgGraph) -> GraphSnapshot {
    let node_data = graph
        .nodes
        .iter()
        .map(|n| {
            (
                n.id,
                n.seq_len,
                n.kmer_cov,
                n.is_deleted,
                n.is_circular,
                n.edges_fwd.to_vec(),
                n.edges_rev.to_vec(),
                graph.get_unitig_seq(n.id),
            )
        })
        .collect();
    GraphSnapshot {
        nodes_count: graph.nodes.len(),
        node_data,
    }
}

/// G1: Verification failure leaves full graph snapshot identical
#[test]
fn test_g1_verification_failure_leaves_graph_snapshot_identical() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);

    let seq0 = make_deterministic_dna(0x6101, 600);
    let seq1 = make_deterministic_dna(0x6102, 600);

    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);

    // Edge claims overlap 31, but sequences do not match!
    add_bidirected_edge(&mut graph, u0, u1, 31, 10);

    let mut matrix = PairLinkMatrix::new();
    for _ in 0..10 {
        matrix.add_link(u0, u1);
    }

    let snapshot_before = capture_graph_snapshot(&graph);

    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams {
            min_support: 5,
            decision_threshold: 2.5,
        },
    );

    assert_eq!(stats.valid_k_joins, 0);
    assert_eq!(stats.valid_k_minus_1_joins, 0);
    assert!(stats.rejected_invalid_joins >= 1);

    let snapshot_after = capture_graph_snapshot(&graph);

    assert_eq!(
        snapshot_before, snapshot_after,
        "Full graph snapshot must be 100% identical after verification failure"
    );
}

/// G2: Path rejection does not pollute used-state for other legal paths
#[test]
fn test_g2_path_rejection_does_not_pollute_used_state_for_other_paths() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);

    // u0 (len 700): longer seed, connects to u1 but with invalid sequence
    // u2 (len 600): second seed, connects to u1 with valid sequence
    // u1 (len 600): shared candidate target
    let shared_ov = make_deterministic_dna(0x6201, 31);
    let p0 = make_deterministic_dna(0x6202, 669);
    let p2 = make_deterministic_dna(0x6203, 569);
    let s1 = make_deterministic_dna(0x6204, 569);

    let seq0 = p0; // len 669, does not end with shared_ov
    let mut seq2 = p2;
    seq2.extend_from_slice(&shared_ov); // len 600, ends with shared_ov

    let mut seq1 = shared_ov.clone();
    seq1.extend_from_slice(&s1); // len 600, starts with shared_ov

    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0); // 669 bp (first seed)
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0); // 600 bp
    let u2 = create_unitig(&mut graph, 2, &seq2, 25.0); // 600 bp (second seed)

    // u0 -> u1 bogus edge (seq0 does not match seq1)
    add_bidirected_edge(&mut graph, u0, u1, 31, 10);
    // u2 -> u1 valid edge (seq2 matches seq1 at 31 bp)
    add_bidirected_edge(&mut graph, u2, u1, 31, 10);

    let mut matrix = PairLinkMatrix::new();
    for _ in 0..10 {
        matrix.add_link(u0, u1);
        matrix.add_link(u2, u1);
    }

    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams {
            min_support: 5,
            decision_threshold: 2.5,
        },
    );

    // First seed u0 failed, but second seed u2 MUST succeed because u1 was not polluted by used_in_path!
    assert_eq!(stats.valid_k_joins, 1, "Seed u2 must legally splice u1");
    assert!(stats.rejected_invalid_joins >= 1, "Seed u0 must be rejected");

    let active_nodes: Vec<_> = graph.nodes.iter().filter(|n| !n.is_deleted).collect();
    assert_eq!(active_nodes.len(), 2, "u0 remains intact, u2+u1 merged into 1 contig");
    assert!(!graph.nodes[0].is_deleted, "u0 must remain active");
    assert!(graph.nodes[1].is_deleted, "u1 was merged by u2");
    assert!(graph.nodes[2].is_deleted, "u2 was merged with u1");
}

/// G3: Consecutive candidate paths: failure of one does not affect another legal path
#[test]
fn test_g3_consecutive_candidate_paths_isolated_failure() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);

    // Pair A: u0 (600 bp) and u1 (600 bp), invalid edge (fails)
    let seq0 = make_deterministic_dna(0x6301, 600);
    let seq1 = make_deterministic_dna(0x6302, 600);
    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);
    add_bidirected_edge(&mut graph, u0, u1, 31, 10);

    // Pair B: u2 (700 bp) and u3 (700 bp), valid 31 bp overlap (succeeds)
    let ov_seq = make_deterministic_dna(0x6303, 31);
    let p2 = make_deterministic_dna(0x6304, 669);
    let s3 = make_deterministic_dna(0x6305, 669);
    let mut seq2 = p2;
    seq2.extend_from_slice(&ov_seq);
    let mut seq3 = ov_seq;
    seq3.extend_from_slice(&s3);
    let u2 = create_unitig(&mut graph, 2, &seq2, 25.0);
    let u3 = create_unitig(&mut graph, 3, &seq3, 25.0);
    add_bidirected_edge(&mut graph, u2, u3, 31, 10);

    let mut matrix = PairLinkMatrix::new();
    for _ in 0..10 {
        matrix.add_link(u0, u1);
        matrix.add_link(u2, u3);
    }

    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams {
            min_support: 5,
            decision_threshold: 2.5,
        },
    );

    assert_eq!(stats.valid_k_joins, 1);
    assert!(stats.rejected_invalid_joins >= 1);

    // u0 and u1 intact; u2 and u3 merged
    assert!(!graph.nodes[0].is_deleted);
    assert!(!graph.nodes[1].is_deleted);
    assert!(graph.nodes[2].is_deleted);
    assert!(graph.nodes[3].is_deleted);
    let active_nodes: Vec<_> = graph.nodes.iter().filter(|n| !n.is_deleted).collect();
    assert_eq!(active_nodes.len(), 3, "Active nodes: u0, u1, merged u2+u3");
}

/// G4: Short nodes belonging to rejected paths are protected from cleanup deletion
#[test]
fn test_g4_short_nodes_in_rejected_paths_preserved_by_cleanup() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);

    let seq0 = make_deterministic_dna(0x6401, 600);
    let seq1 = make_deterministic_dna(0x6402, 300); // Short node (< 500 bp)

    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);

    // Bogus edge triggering validation failure
    add_bidirected_edge(&mut graph, u0, u1, 31, 10);

    let mut matrix = PairLinkMatrix::new();
    for _ in 0..10 {
        matrix.add_link(u0, u1);
    }

    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams {
            min_support: 5,
            decision_threshold: 2.5,
        },
    );

    assert_eq!(stats.valid_k_joins, 0);
    assert!(stats.rejected_invalid_joins >= 1);

    // Contract 2: u1 is < 500 bp, but it was involved in a rejected candidate path!
    // It MUST NOT be deleted by post-exspander cleanup!
    assert!(
        !graph.nodes[1].is_deleted,
        "Short unitig u1 (<500 bp) from rejected path must be preserved by cleanup protection"
    );
    assert!(!graph.nodes[0].is_deleted, "Seed u0 must remain active");
}

/// G5: Isolated short nodes are deleted according to documented policy
#[test]
fn test_g5_unrelated_isolated_short_nodes_deleted_by_policy() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);

    let seq0 = make_deterministic_dna(0x6501, 600);
    let seq1 = make_deterministic_dna(0x6502, 300); // Short node in rejected path
    let seq2 = make_deterministic_dna(0x6503, 200); // Completely isolated short node (<500 bp)

    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);
    let _u2 = create_unitig(&mut graph, 2, &seq2, 25.0); // No edges to u0 or u1!

    add_bidirected_edge(&mut graph, u0, u1, 31, 10);

    let mut matrix = PairLinkMatrix::new();
    for _ in 0..10 {
        matrix.add_link(u0, u1);
    }

    resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams {
            min_support: 5,
            decision_threshold: 2.5,
        },
    );

    // Contract 2: u1 was part of rejected path -> protected
    assert!(
        !graph.nodes[1].is_deleted,
        "u1 in rejected path must be preserved"
    );

    // Contract 3: u2 is an unrelated isolated short node (<500 bp) -> deleted by cleanup policy
    assert!(
        graph.nodes[2].is_deleted,
        "Unrelated isolated short node u2 (<500 bp) must be deleted by cleanup policy"
    );
}

/// G6: Successful commit has no abnormal deletion, duplicate creation, or lost nodes
#[test]
fn test_g6_successful_commit_no_abnormal_deletion_or_lost_nodes() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);
    let ov_len = 31;

    let p0 = make_deterministic_dna(0x6601, 569);
    let shared = make_deterministic_dna(0x6602, ov_len);
    let s1 = make_deterministic_dna(0x6603, 569);
    let seq2 = make_deterministic_dna(0x6604, 800); // Independent bystander node

    let mut seq0 = p0;
    seq0.extend_from_slice(&shared); // len = 600

    let mut seq1 = shared;
    seq1.extend_from_slice(&s1); // len = 600

    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);
    let _u2 = create_unitig(&mut graph, 2, &seq2, 25.0);

    add_bidirected_edge(&mut graph, u0, u1, ov_len as u8, 10);

    let mut matrix = PairLinkMatrix::new();
    for _ in 0..10 {
        matrix.add_link(u0, u1);
    }

    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams {
            min_support: 5,
            decision_threshold: 2.5,
        },
    );

    assert_eq!(stats.valid_k_joins, 1);

    // Total nodes in array: 0, 1, 2 (initial) + 3 (new spliced) = 4
    assert_eq!(graph.nodes.len(), 4);

    // u0 and u1 must be marked deleted
    assert!(graph.nodes[0].is_deleted);
    assert!(graph.nodes[1].is_deleted);

    // Bystander u2 must remain active and unmodified
    assert!(!graph.nodes[2].is_deleted);
    assert_eq!(graph.get_unitig_seq(2), seq2);

    // Newly created node u3 must be active with exact merged length
    assert!(!graph.nodes[3].is_deleted);
    assert_eq!(graph.nodes[3].seq_len, 1169);

    // Active nodes count is exactly 2: u2 and u3
    let active_nodes: Vec<_> = graph.nodes.iter().filter(|n| !n.is_deleted).collect();
    assert_eq!(active_nodes.len(), 2);
    let active_ids: Vec<u32> = active_nodes.iter().map(|n| n.id).collect();
    assert_eq!(active_ids, vec![2, 3]);

    // Total active bases = 800 + 1169 = 1969 bp
    let total_active_bases: usize = active_nodes.iter().map(|n| n.seq_len as usize).sum();
    assert_eq!(total_active_bases, 1969);
}
