use fastbact::common::bitseq::reverse_complement_ascii;
use fastbact::common::sketch::splitmix64;
use fastbact::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, UnitigNode};
use fastbact::resolve::exspander::{resolve_repeats_exspander, ExspanderParams};
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

fn run_exspander(graph: &mut CdbgGraph, matrix: &PairLinkMatrix) {
    let dist = InsertSizeDist::default();
    let params = ExspanderParams {
        min_support: 5,
        decision_threshold: 2.5,
    };
    resolve_repeats_exspander(graph, matrix, &dist, params);
}

/// T1: Two unitigs, valid overlap = k - 1 (30 bp)
#[test]
fn test_t1_two_unitigs_valid_overlap_k_minus_1() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);
    let ov_len = k - 1; // 30 bp

    let prefix0 = make_deterministic_dna(0x1001, 570);
    let shared = make_deterministic_dna(0x1002, ov_len);
    let suffix1 = make_deterministic_dna(0x1003, 570);

    let mut seq0 = prefix0;
    seq0.extend_from_slice(&shared);

    let mut seq1 = shared.clone();
    seq1.extend_from_slice(&suffix1);

    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);
    add_bidirected_edge(&mut graph, u0, u1, ov_len as u8, 10);

    let mut matrix = PairLinkMatrix::new();
    for _ in 0..10 {
        matrix.add_link(u0, u1);
    }

    let mut expected_seq = seq0.clone();
    expected_seq.extend_from_slice(&suffix1); // length = 600 + 600 - 30 = 1170

    run_exspander(&mut graph, &matrix);

    let active_nodes: Vec<_> = graph.nodes.iter().filter(|n| !n.is_deleted).collect();
    assert_eq!(active_nodes.len(), 1, "Should merge into exactly 1 active unitig");
    let merged_seq = graph.get_unitig_seq(active_nodes[0].id);
    assert_eq!(merged_seq, expected_seq, "Spliced DNA sequence must match exactly");
}

/// T2: Two unitigs, valid overlap = k (31 bp)
#[test]
fn test_t2_two_unitigs_valid_overlap_k() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);
    let ov_len = k; // 31 bp

    let prefix0 = make_deterministic_dna(0x2001, 569);
    let shared = make_deterministic_dna(0x2002, ov_len);
    let suffix1 = make_deterministic_dna(0x2003, 569);

    let mut seq0 = prefix0;
    seq0.extend_from_slice(&shared);

    let mut seq1 = shared.clone();
    seq1.extend_from_slice(&suffix1);

    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);
    add_bidirected_edge(&mut graph, u0, u1, ov_len as u8, 10);

    let mut matrix = PairLinkMatrix::new();
    for _ in 0..10 {
        matrix.add_link(u0, u1);
    }

    let mut expected_seq = seq0.clone();
    expected_seq.extend_from_slice(&suffix1); // length = 600 + 600 - 31 = 1169

    run_exspander(&mut graph, &matrix);

    let active_nodes: Vec<_> = graph.nodes.iter().filter(|n| !n.is_deleted).collect();
    assert_eq!(active_nodes.len(), 1, "Should merge into exactly 1 active unitig");
    let merged_seq = graph.get_unitig_seq(active_nodes[0].id);
    assert_eq!(merged_seq.len(), expected_seq.len(), "Length must be 1169 bp, not 1170 bp");
    assert_eq!(merged_seq, expected_seq, "Spliced DNA sequence must match exactly without duplicated base");
}

/// T3: Same directed node pair simultaneously has k and k-1 overlap edges
#[test]
fn test_t3_same_directed_pair_both_k_and_k_minus_1_edges() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);

    // Overlap k=31: u0 ends with 31 As, u1 starts with 31 As
    // Both 31-bp and 30-bp overlaps match at sequence level
    let prefix0 = make_deterministic_dna(0x3001, 569);
    let homopolymer = vec![b'A'; 31];
    let suffix1 = make_deterministic_dna(0x3002, 569);

    let mut seq0 = prefix0;
    seq0.extend_from_slice(&homopolymer);

    let mut seq1 = homopolymer.clone();
    seq1.extend_from_slice(&suffix1);

    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);

    // Add both edges: overlap 31 has higher support (20), overlap 30 has lower support (10)
    add_bidirected_edge(&mut graph, u0, u1, 31, 20);
    add_bidirected_edge(&mut graph, u0, u1, 30, 10);

    let mut matrix = PairLinkMatrix::new();
    for _ in 0..15 {
        matrix.add_link(u0, u1);
    }

    // Expected: edge with higher support (overlap 31) is unambiguously chosen
    let mut expected_seq = seq0.clone();
    expected_seq.extend_from_slice(&suffix1); // length = 1169

    run_exspander(&mut graph, &matrix);

    let active_nodes: Vec<_> = graph.nodes.iter().filter(|n| !n.is_deleted).collect();
    assert_eq!(active_nodes.len(), 1, "Should merge into exactly 1 active unitig");
    let merged_seq = graph.get_unitig_seq(active_nodes[0].id);
    assert_eq!(merged_seq, expected_seq, "Edge with higher support (k=31) must be selected");
}

/// T4: Reverse-complement orientation
#[test]
fn test_t4_reverse_complement_orientation() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);
    let ov_len = k; // 31 bp

    let prefix0 = make_deterministic_dna(0x4001, 569);
    let shared = make_deterministic_dna(0x4002, ov_len);
    let suffix1_rc = make_deterministic_dna(0x4003, 569);

    let mut seq0 = prefix0;
    seq0.extend_from_slice(&shared);

    // u1 in reverse-complement orientation starts with shared
    let mut u1_stranded_seq = shared.clone();
    u1_stranded_seq.extend_from_slice(&suffix1_rc);
    let seq1 = reverse_complement_ascii(&u1_stranded_seq);

    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);
    let u1_rc = u1.twin();

    add_bidirected_edge(&mut graph, u0, u1_rc, ov_len as u8, 10);

    let mut matrix = PairLinkMatrix::new();
    for _ in 0..10 {
        matrix.add_link(u0, u1_rc);
    }

    let mut expected_seq = seq0.clone();
    expected_seq.extend_from_slice(&suffix1_rc); // length = 1169

    run_exspander(&mut graph, &matrix);

    let active_nodes: Vec<_> = graph.nodes.iter().filter(|n| !n.is_deleted).collect();
    assert_eq!(active_nodes.len(), 1, "Should merge into exactly 1 active unitig");
    let merged_seq = graph.get_unitig_seq(active_nodes[0].id);
    assert_eq!(merged_seq, expected_seq, "Reverse-complement orientation must splice correctly");
}

/// T5: Three unitigs, path mixing k and k-1 overlap
#[test]
fn test_t5_three_unitigs_mixing_k_and_k_minus_1() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);

    let ov_k = k; // 31 bp
    let ov_k_minus_1 = k - 1; // 30 bp

    let p0 = make_deterministic_dna(0x5001, 569);
    let s01 = make_deterministic_dna(0x5002, ov_k);
    let m1 = make_deterministic_dna(0x5003, 539);
    let s12 = make_deterministic_dna(0x5004, ov_k_minus_1);
    let suff2 = make_deterministic_dna(0x5005, 570);

    // u0: p0 + s01 (len 600)
    let mut seq0 = p0.clone();
    seq0.extend_from_slice(&s01);

    // u1: s01 + m1 + s12 (len 31 + 539 + 30 = 600)
    let mut seq1 = s01.clone();
    seq1.extend_from_slice(&m1);
    seq1.extend_from_slice(&s12);

    // u2: s12 + suff2 (len 30 + 570 = 600)
    let mut seq2 = s12.clone();
    seq2.extend_from_slice(&suff2);

    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);
    let u2 = create_unitig(&mut graph, 2, &seq2, 25.0);

    add_bidirected_edge(&mut graph, u0, u1, ov_k as u8, 10);
    add_bidirected_edge(&mut graph, u1, u2, ov_k_minus_1 as u8, 10);

    let mut matrix = PairLinkMatrix::new();
    for _ in 0..10 {
        matrix.add_link(u0, u1);
        matrix.add_link(u1, u2);
        matrix.add_link(u0, u2);
    }

    // Expected: p0 + s01 + m1 + s12 + suff2
    // Length = 569 + 31 + 539 + 30 + 570 = 1739 bp
    let mut expected_seq = p0;
    expected_seq.extend_from_slice(&s01);
    expected_seq.extend_from_slice(&m1);
    expected_seq.extend_from_slice(&s12);
    expected_seq.extend_from_slice(&suff2);
    assert_eq!(expected_seq.len(), 1739);

    run_exspander(&mut graph, &matrix);

    let active_nodes: Vec<_> = graph.nodes.iter().filter(|n| !n.is_deleted).collect();
    assert_eq!(active_nodes.len(), 1, "Should merge 3 unitigs into 1 contig");
    let merged_seq = graph.get_unitig_seq(active_nodes[0].id);
    assert_eq!(merged_seq, expected_seq, "Mixed k and k-1 overlaps must assemble base-accurate sequence");
}

/// T6: Edge recorded overlap does not match actual suffix/prefix
#[test]
fn test_t6_edge_overlap_mismatches_actual_sequence() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);

    let seq0 = make_deterministic_dna(0x6001, 600);
    let seq1 = make_deterministic_dna(0x6002, 600);

    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);

    // Add edge claiming overlap 31, but actual sequences are random and do NOT match!
    add_bidirected_edge(&mut graph, u0, u1, 31, 10);

    let mut matrix = PairLinkMatrix::new();
    for _ in 0..10 {
        matrix.add_link(u0, u1);
    }

    run_exspander(&mut graph, &matrix);

    // Must REJECT splicing because sequence overlap is invalid!
    let active_nodes: Vec<_> = graph.nodes.iter().filter(|n| !n.is_deleted).collect();
    assert_eq!(
        active_nodes.len(),
        2,
        "Splicing must be rejected when edge overlap does not match actual sequence"
    );
}

/// T7: Candidate path without valid connection, must reject splicing
#[test]
fn test_t7_no_valid_connection_rejects_splicing() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);

    let seq0 = make_deterministic_dna(0x7001, 600);
    let seq1 = make_deterministic_dna(0x7002, 600);

    let _u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let _u1 = create_unitig(&mut graph, 1, &seq1, 25.0);

    // No edge at all!
    let matrix = PairLinkMatrix::new();

    run_exspander(&mut graph, &matrix);

    let active_nodes: Vec<_> = graph.nodes.iter().filter(|n| !n.is_deleted).collect();
    assert_eq!(active_nodes.len(), 2, "No connection: both unitigs remain intact");
}

/// T8: After splicing failure, original graph nodes, edges, sequences, and deleted state remain unchanged
#[test]
fn test_t8_failure_preserves_graph_state_atomically() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);

    let seq0 = make_deterministic_dna(0x8001, 600);
    let seq1 = make_deterministic_dna(0x8002, 300); // smaller node (< 500)

    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);

    // Bogus edge to trigger rejection
    add_bidirected_edge(&mut graph, u0, u1, 31, 10);

    let initial_node_count = graph.nodes.len();
    let initial_node0_deleted = graph.nodes[0].is_deleted;
    let initial_node1_deleted = graph.nodes[1].is_deleted;
    let initial_edges_count: usize = graph
        .nodes
        .iter()
        .map(|n| n.edges_fwd.len() + n.edges_rev.len())
        .sum();

    let matrix = PairLinkMatrix::new();

    run_exspander(&mut graph, &matrix);

    assert_eq!(graph.nodes.len(), initial_node_count, "Node count unchanged");
    assert_eq!(
        graph.nodes[0].is_deleted, initial_node0_deleted,
        "Node 0 deletion state unchanged"
    );
    assert_eq!(
        graph.nodes[1].is_deleted, initial_node1_deleted,
        "Node 1 (<500 bp) must NOT be deleted when splicing fails"
    );
    assert_eq!(graph.get_unitig_seq(0), seq0, "Node 0 sequence unchanged");
    assert_eq!(graph.get_unitig_seq(1), seq1, "Node 1 sequence unchanged");
    let current_edges_count: usize = graph
        .nodes
        .iter()
        .map(|n| n.edges_fwd.len() + n.edges_rev.len())
        .sum();
    assert_eq!(
        current_edges_count, initial_edges_count,
        "Edge count unchanged"
    );
}

/// T9: Same target has multiple different overlap edges, selection basis explicit, not arbitrary traversal order
#[test]
fn test_t9_multiple_overlap_edges_explicit_selection_and_ambiguity_rejection() {
    let k = 31;

    // Subtest 9A: Higher support on k-1 (30 bp) selects k-1
    {
        let mut graph = CdbgGraph::new(k);
        let prefix0 = make_deterministic_dna(0x9001, 569);
        let homopolymer = vec![b'A'; 31];
        let suffix1 = make_deterministic_dna(0x9002, 569);

        let mut seq0 = prefix0.clone();
        seq0.extend_from_slice(&homopolymer);

        let mut seq1 = homopolymer.clone();
        seq1.extend_from_slice(&suffix1);

        let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
        let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);

        // Edge with overlap 30 has support 20; overlap 31 has support 10
        add_bidirected_edge(&mut graph, u0, u1, 30, 20);
        add_bidirected_edge(&mut graph, u0, u1, 31, 10);

        let mut matrix = PairLinkMatrix::new();
        for _ in 0..15 {
            matrix.add_link(u0, u1);
        }

        // Expected: selects overlap 30. Length = 600 + 600 - 30 = 1170
        let mut expected_seq = seq0;
        expected_seq.push(b'A'); // extra A because overlap is 30 instead of 31
        expected_seq.extend_from_slice(&suffix1);

        run_exspander(&mut graph, &matrix);

        let active_nodes: Vec<_> = graph.nodes.iter().filter(|n| !n.is_deleted).collect();
        assert_eq!(active_nodes.len(), 1);
        let merged_seq = graph.get_unitig_seq(active_nodes[0].id);
        assert_eq!(merged_seq, expected_seq, "Higher support on k-1 edge must select k-1");
    }

    // Subtest 9B: Equal support on both edges -> Ambiguous join, must reject
    {
        let mut graph = CdbgGraph::new(k);
        let prefix0 = make_deterministic_dna(0x9003, 569);
        let homopolymer = vec![b'A'; 31];
        let suffix1 = make_deterministic_dna(0x9004, 569);

        let mut seq0 = prefix0;
        seq0.extend_from_slice(&homopolymer);

        let mut seq1 = homopolymer;
        seq1.extend_from_slice(&suffix1);

        let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
        let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);

        // Both edges have equal support 15
        add_bidirected_edge(&mut graph, u0, u1, 30, 15);
        add_bidirected_edge(&mut graph, u0, u1, 31, 15);

        let mut matrix = PairLinkMatrix::new();
        for _ in 0..15 {
            matrix.add_link(u0, u1);
        }

        run_exspander(&mut graph, &matrix);

        let active_nodes: Vec<_> = graph.nodes.iter().filter(|n| !n.is_deleted).collect();
        assert_eq!(
            active_nodes.len(),
            2,
            "Ambiguous join with tied support must be rejected"
        );
    }
}

/// T10: Backward path reversal and twin conversion still references correct edge and overlap
#[test]
fn test_t10_backward_path_reversal_references_correct_edge_and_overlap() {
    let k = 31;
    let mut graph = CdbgGraph::new(k);
    let ov_len = k; // 31 bp

    let prefix0 = make_deterministic_dna(0xA001, 569);
    let shared = make_deterministic_dna(0xA002, ov_len);
    let suffix1 = make_deterministic_dna(0xA003, 569);

    let mut seq0 = prefix0;
    seq0.extend_from_slice(&shared);

    let mut seq1 = shared.clone();
    seq1.extend_from_slice(&suffix1);

    // Make u1 longer so it is selected as the first seed
    // u1 length = 650, u0 length = 600
    let extra = make_deterministic_dna(0xA004, 50);
    seq1.extend_from_slice(&extra);

    let u0 = create_unitig(&mut graph, 0, &seq0, 25.0);
    let u1 = create_unitig(&mut graph, 1, &seq1, 25.0);

    add_bidirected_edge(&mut graph, u0, u1, ov_len as u8, 10);

    let mut matrix = PairLinkMatrix::new();
    for _ in 0..15 {
        // PE links connecting u0 and u1
        matrix.add_link(u0, u1);
        matrix.add_link(u1.twin(), u0.twin());
    }

    // Since u1 is longer (650 > 600), exspander seeds from u1.
    // Backward extension from u1 extends to u0.
    // Resulting contig must be u0 + u1[31..]
    let mut expected_seq = seq0.clone();
    expected_seq.extend_from_slice(&seq1[ov_len..]);

    run_exspander(&mut graph, &matrix);

    let active_nodes: Vec<_> = graph.nodes.iter().filter(|n| !n.is_deleted).collect();
    assert_eq!(active_nodes.len(), 1, "Should merge via backward extension");
    let merged_seq = graph.get_unitig_seq(active_nodes[0].id);
    assert_eq!(
        merged_seq, expected_seq,
        "Backward extension reversal and twin conversion must preserve correct k overlap"
    );
}
