use fastbact::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, UnitigNode};
use fastbact::resolve::exspander::{
    resolve_repeats_exspander_with_stats, validate_and_splice_path, ExspanderParams, PathStep,
};
use fastbact::resolve::mapper::{InsertSizeDist, PairLinkMatrix};

fn fixture(supports: &[u16]) -> (CdbgGraph, Vec<PathStep>) {
    let mut graph = CdbgGraph::new(31);
    for (id, len) in [(0, 600), (1, 650)] {
        let seq = graph.arena.push_seq(&vec![b'A'; len]);
        graph
            .nodes
            .push(UnitigNode::new(id, seq, u32::try_from(len).unwrap(), 25.0));
    }
    let from = DirUnitigId::new(0, false);
    let to = DirUnitigId::new(1, false);
    for (&support, overlap) in supports.iter().zip([30, 31]) {
        graph.nodes[0].add_edge(false, GraphEdge::new(to, overlap, support));
        graph.nodes[1].add_edge(true, GraphEdge::new(from.twin(), overlap, support));
    }
    let edge = graph.nodes[0].edges(false)[0];
    let path = vec![
        PathStep {
            dir: from,
            incoming_edge: None,
            overlap: 0,
            support: 0,
        },
        PathStep {
            dir: to,
            incoming_edge: Some(edge),
            overlap: usize::from(edge.overlap_k),
            support: u32::from(edge.support),
        },
    ];
    (graph, path)
}

#[test]
fn e7_direct_validator_rejects_tied_overlaps() {
    let (graph, path) = fixture(&[10, 10]);
    assert!(validate_and_splice_path(&graph, &path).is_err());
}

#[test]
fn e7_direct_validator_rejects_lower_support_edge() {
    let (graph, path) = fixture(&[10, 20]);
    assert!(validate_and_splice_path(&graph, &path).is_err());
}

#[test]
fn e7_direct_validator_accepts_selected_edge() {
    let (graph, mut path) = fixture(&[20, 10]);
    // PathStep.support may contain paired-read evidence rather than GraphEdge support.
    path[1].support = 99;
    assert_eq!(
        validate_and_splice_path(&graph, &path).unwrap(),
        vec![b'A'; 1220]
    );
}

#[test]
fn e6_backward_reciprocal_support_must_match_selected_edge() {
    let (mut graph, _) = fixture(&[10]);
    graph.nodes[1].edges_rev[0].support = 11;
    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &PairLinkMatrix::new(),
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    );
    assert_eq!(stats.valid_k_minus_1_joins, 0);
    assert!(stats.rejected_invalid_joins > 0);
    assert_eq!(graph.nodes.len(), 2);
    assert!(graph.nodes.iter().all(|node| !node.is_deleted));
}

#[test]
fn e1_matching_sequence_on_wrong_strand_is_rejected() {
    let (mut graph, path) = fixture(&[10]);
    let edge = graph.nodes[0].edges_fwd.remove(0);
    graph.nodes[0].edges_rev.push(edge);
    assert!(validate_and_splice_path(&graph, &path).is_err());
}
