//! FB-P0-01R / C0-02R: Directed Edge Identity and Verification Contract Tests (E1–E8).
//!
//! Validates:
//! E1: Non-existent real edge, but suffix/prefix match
//! E2: incoming_edge=None on non-first step
//! E3: incoming_edge.target incorrect
//! E4: incoming_edge.overlap_k inconsistent with step.overlap
//! E5: Forged support field or copied edge metadata
//! E6: Backward path with missing reciprocal edge
//! E7: Multiple overlap candidates have identical support (ambiguity rejection)
//! E8: Normal path with correct k and k-1 overlaps still passes

use fastbact::common::bitseq::reverse_complement_ascii;
use fastbact::common::sketch::splitmix64;
use fastbact::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, UnitigNode};
use fastbact::resolve::exspander::{
    resolve_edge_between, resolve_repeats_exspander_with_stats, validate_and_splice_path,
    EdgeResolveResult, ExspanderParams, PathStep,
};
use fastbact::resolve::mapper::PairLinkMatrix;

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

/// E1: Non-existent real edge, but suffix/prefix match
#[test]
fn test_e1_nonexistent_edge_matching_sequence_rejected() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);
    let ov_len = 31;

    let p0 = make_deterministic_dna(0xE101, 569);
    let shared = make_deterministic_dna(0xE102, ov_len);
    let s1 = make_deterministic_dna(0xE103, 569);

    let mut seq0 = p0;
    seq0.extend_from_slice(&shared);

    let mut seq1 = shared.clone();
    seq1.extend_from_slice(&s1);

    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);

    // Notice: NO edge added in graph between u0 and u1!
    // Try to validate a forged path
    let path = vec![
        PathStep {
            dir: u0,
            incoming_edge: None,
            overlap: 0,
            support: 0,
        },
        PathStep {
            dir: u1,
            incoming_edge: Some(GraphEdge::new(u1, ov_len as u8, 10)),
            overlap: ov_len,
            support: 10,
        },
    ];

    let res = validate_and_splice_path(&graph, &path);
    assert!(res.is_err());
    assert_eq!(
        res.unwrap_err(),
        "Incoming edge does not exist in predecessor outgoing edges"
    );
}

/// E2: incoming_edge=None on non-first step
#[test]
fn test_e2_none_incoming_edge_on_non_first_step_rejected() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);
    let ov_len = 31;

    let p0 = make_deterministic_dna(0xE201, 569);
    let shared = make_deterministic_dna(0xE202, ov_len);
    let s1 = make_deterministic_dna(0xE203, 569);

    let mut seq0 = p0;
    seq0.extend_from_slice(&shared);

    let mut seq1 = shared;
    seq1.extend_from_slice(&s1);

    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);
    add_bidirected_edge(&mut graph, u0, u1, ov_len as u8, 10);

    // Step 1 illegally lacks incoming_edge
    let path = vec![
        PathStep {
            dir: u0,
            incoming_edge: None,
            overlap: 0,
            support: 0,
        },
        PathStep {
            dir: u1,
            incoming_edge: None, // Illegal on non-first step!
            overlap: ov_len,
            support: 10,
        },
    ];

    let res = validate_and_splice_path(&graph, &path);
    assert!(res.is_err());
    assert_eq!(res.unwrap_err(), "Non-initial step missing incoming edge");
}

/// E3: incoming_edge.target incorrect
#[test]
fn test_e3_incoming_edge_target_mismatch_rejected() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);
    let ov_len = 31;

    let seq0 = make_deterministic_dna(0xE301, 600);
    let seq1 = make_deterministic_dna(0xE302, 600);
    let seq2 = make_deterministic_dna(0xE303, 600);

    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);
    let u2 = create_unitig(&mut graph, 2, &seq2, 25.0);

    add_bidirected_edge(&mut graph, u0, u1, ov_len as u8, 10);
    add_bidirected_edge(&mut graph, u0, u2, ov_len as u8, 10);

    // Step 1 dir is u1, but incoming_edge claims target is u2
    let path = vec![
        PathStep {
            dir: u0,
            incoming_edge: None,
            overlap: 0,
            support: 0,
        },
        PathStep {
            dir: u1,
            incoming_edge: Some(GraphEdge::new(u2, ov_len as u8, 10)), // Target mismatch!
            overlap: ov_len,
            support: 10,
        },
    ];

    let res = validate_and_splice_path(&graph, &path);
    assert!(res.is_err());
    assert_eq!(
        res.unwrap_err(),
        "Incoming edge target does not match step dir"
    );
}

/// E4: incoming_edge.overlap_k inconsistent with step.overlap
#[test]
fn test_e4_incoming_edge_overlap_mismatch_rejected() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);

    let seq0 = make_deterministic_dna(0xE401, 600);
    let seq1 = make_deterministic_dna(0xE402, 600);

    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);

    add_bidirected_edge(&mut graph, u0, u1, 30, 10);

    // incoming_edge has overlap 30, but step.overlap claims 31
    let path = vec![
        PathStep {
            dir: u0,
            incoming_edge: None,
            overlap: 0,
            support: 0,
        },
        PathStep {
            dir: u1,
            incoming_edge: Some(GraphEdge::new(u1, 30, 10)),
            overlap: 31, // Inconsistent with incoming_edge.overlap_k!
            support: 10,
        },
    ];

    let res = validate_and_splice_path(&graph, &path);
    assert!(res.is_err());
    assert_eq!(
        res.unwrap_err(),
        "Incoming edge overlap_k does not match step overlap"
    );
}

/// E5: Forged support field or copied edge metadata
#[test]
fn test_e5_forged_support_metadata_rejected() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);
    let ov_len = 31;

    let p0 = make_deterministic_dna(0xE501, 569);
    let shared = make_deterministic_dna(0xE502, ov_len);
    let s1 = make_deterministic_dna(0xE503, 569);

    let mut seq0 = p0;
    seq0.extend_from_slice(&shared);

    let mut seq1 = shared;
    seq1.extend_from_slice(&s1);

    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);

    // Graph edge has authentic support = 10
    add_bidirected_edge(&mut graph, u0, u1, ov_len as u8, 10);

    // PathStep carries forged support = 999
    let path = vec![
        PathStep {
            dir: u0,
            incoming_edge: None,
            overlap: 0,
            support: 0,
        },
        PathStep {
            dir: u1,
            incoming_edge: Some(GraphEdge::new(u1, ov_len as u8, 999)), // Forged support!
            overlap: ov_len,
            support: 999,
        },
    ];

    let res = validate_and_splice_path(&graph, &path);
    assert!(res.is_err());
    assert_eq!(
        res.unwrap_err(),
        "Incoming edge support does not match graph edge support"
    );
}

/// E6: Backward path with missing reciprocal edge
#[test]
fn test_e6_backward_path_missing_reciprocal_edge_rejected() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);
    let ov_len = 31;

    let p0 = make_deterministic_dna(0xE601, 569);
    let shared = make_deterministic_dna(0xE602, ov_len);
    let s1 = make_deterministic_dna(0xE603, 619); // u1 longer (650 > 600) to become seed

    let mut seq0 = p0;
    seq0.extend_from_slice(&shared); // len = 600

    let mut seq1 = shared;
    seq1.extend_from_slice(&s1); // len = 650

    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);

    // Add only unidirectional edge from u1.twin() to u0.twin(), but DO NOT add reciprocal u0 -> u1
    graph.nodes[u1.id() as usize].add_edge(u1.twin().is_rc(), GraphEdge::new(u0.twin(), ov_len as u8, 10));

    let mut matrix = PairLinkMatrix::new();
    for _ in 0..10 {
        matrix.add_link(u1.twin(), u0.twin());
    }

    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &fastbact::resolve::mapper::InsertSizeDist::default(),
        ExspanderParams {
            min_support: 5,
            decision_threshold: 2.5,
        },
    );

    // Path must be conservatively rejected due to missing reciprocal edge
    assert_eq!(stats.valid_k_joins, 0);
    assert!(stats.rejected_invalid_joins >= 1);

    // Graph state: both nodes remain intact, no merged node
    let active_nodes: Vec<_> = graph.nodes.iter().filter(|n| !n.is_deleted).collect();
    assert_eq!(active_nodes.len(), 2, "Both nodes must be preserved when reciprocal edge is missing");
}

/// E7: Multiple overlap candidates have identical support (tied support)
#[test]
fn test_e7_multiple_overlap_candidates_tied_support_ambiguous_rejected() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);

    let prefix0 = make_deterministic_dna(0xE701, 569);
    let homopolymer = vec![b'A'; 31];
    let suffix1 = make_deterministic_dna(0xE702, 569);

    let mut seq0 = prefix0;
    seq0.extend_from_slice(&homopolymer);

    let mut seq1 = homopolymer;
    seq1.extend_from_slice(&suffix1);

    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);

    // Both overlap 30 and 31 have identical support = 15
    add_bidirected_edge(&mut graph, u0, u1, 30, 15);
    add_bidirected_edge(&mut graph, u0, u1, 31, 15);

    let res = resolve_edge_between(&graph, u0, u1);
    assert_eq!(
        res,
        EdgeResolveResult::Ambiguous,
        "Equal support must resolve to Ambiguous"
    );

    let mut matrix = PairLinkMatrix::new();
    for _ in 0..10 {
        matrix.add_link(u0, u1);
    }

    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &fastbact::resolve::mapper::InsertSizeDist::default(),
        ExspanderParams {
            min_support: 5,
            decision_threshold: 2.5,
        },
    );

    assert_eq!(stats.valid_k_joins, 0);
    assert_eq!(stats.valid_k_minus_1_joins, 0);
    assert!(stats.rejected_ambiguous_joins >= 1);

    let active_nodes: Vec<_> = graph.nodes.iter().filter(|n| !n.is_deleted).collect();
    assert_eq!(active_nodes.len(), 2, "Ambiguous candidates must not be spliced");
}

/// E8: Normal path with correct k and k-1 overlaps still passes
#[test]
fn test_e8_valid_mixed_k_and_k_minus_1_path_succeeds() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);

    let ov_k = 31;
    let ov_km1 = 30;

    let p0 = make_deterministic_dna(0xE801, 569);
    let s01 = make_deterministic_dna(0xE802, ov_k);
    let m1 = make_deterministic_dna(0xE803, 539);
    let s12 = make_deterministic_dna(0xE804, ov_km1);
    let suff2 = make_deterministic_dna(0xE805, 570);

    let mut seq0 = p0.clone();
    seq0.extend_from_slice(&s01);

    let mut seq1 = s01.clone();
    seq1.extend_from_slice(&m1);
    seq1.extend_from_slice(&s12);

    let mut seq2 = s12.clone();
    seq2.extend_from_slice(&suff2);

    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);
    let u2 = create_unitig(&mut graph, 2, &seq2, 25.0);

    add_bidirected_edge(&mut graph, u0, u1, ov_k as u8, 10);
    add_bidirected_edge(&mut graph, u1, u2, ov_km1 as u8, 10);

    let mut expected_seq = p0;
    expected_seq.extend_from_slice(&s01);
    expected_seq.extend_from_slice(&m1);
    expected_seq.extend_from_slice(&s12);
    expected_seq.extend_from_slice(&suff2);
    assert_eq!(expected_seq.len(), 1739);

    let edge01 = graph.nodes[0].find_edge(false, u1, ov_k as u8).copied().unwrap();
    let edge12 = graph.nodes[1].find_edge(false, u2, ov_km1 as u8).copied().unwrap();

    let path = vec![
        PathStep {
            dir: u0,
            incoming_edge: None,
            overlap: 0,
            support: 0,
        },
        PathStep {
            dir: u1,
            incoming_edge: Some(edge01),
            overlap: ov_k,
            support: 10,
        },
        PathStep {
            dir: u2,
            incoming_edge: Some(edge12),
            overlap: ov_km1,
            support: 10,
        },
    ];

    let spliced_seq = validate_and_splice_path(&graph, &path).expect("Valid path must splice");
    assert_eq!(
        spliced_seq, expected_seq,
        "Spliced sequence must exactly match expected 1739 bp sequence"
    );

    let mut matrix = PairLinkMatrix::new();
    for _ in 0..10 {
        matrix.add_link(u0, u1);
        matrix.add_link(u1, u2);
        matrix.add_link(u0, u2);
    }

    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &fastbact::resolve::mapper::InsertSizeDist::default(),
        ExspanderParams {
            min_support: 5,
            decision_threshold: 2.5,
        },
    );

    assert_eq!(stats.valid_k_joins, 1);
    assert_eq!(stats.valid_k_minus_1_joins, 1);

    let active_nodes: Vec<_> = graph.nodes.iter().filter(|n| !n.is_deleted).collect();
    assert_eq!(active_nodes.len(), 1);
    assert_eq!(graph.get_unitig_seq(active_nodes[0].id), expected_seq);
}
