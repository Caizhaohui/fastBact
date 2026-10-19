use fastbact::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, UnitigNode};
use fastbact::resolve::exspander::{
    resolve_repeats_exspander_with_stats, validate_and_splice_path, ExspanderParams, PathStep,
};
use fastbact::resolve::mapper::{InsertSizeDist, PairLinkMatrix};

fn node(graph: &mut CdbgGraph, len: u32) -> DirUnitigId {
    let id = u32::try_from(graph.nodes.len()).unwrap();
    let seq = graph
        .arena
        .push_seq(&vec![b'A'; usize::try_from(len).unwrap()]);
    graph.nodes.push(UnitigNode::new(id, seq, len, 25.0));
    DirUnitigId::new(id, false)
}

fn edge(graph: &mut CdbgGraph, from: DirUnitigId, to: DirUnitigId) {
    graph.nodes[usize::try_from(from.id()).unwrap()].add_edge(false, GraphEdge::new(to, 31, 10));
    graph.nodes[usize::try_from(to.id()).unwrap()]
        .add_edge(true, GraphEdge::new(from.twin(), 31, 10));
}

fn resolve(graph: &mut CdbgGraph, links: &PairLinkMatrix) {
    resolve_repeats_exspander_with_stats(
        graph,
        links,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    );
}

#[test]
fn g5_connected_short_nodes_are_not_isolated_cleanup_artifacts() {
    let mut graph = CdbgGraph::new(31);
    let a = node(&mut graph, 200);
    let b = node(&mut graph, 300);
    edge(&mut graph, a, b);
    let before = format!("{graph:?}");
    resolve(&mut graph, &PairLinkMatrix::new());
    assert_eq!(format!("{graph:?}"), before);
}

#[test]
fn g5_incoming_only_short_node_is_connected() {
    let mut graph = CdbgGraph::new(31);
    let a = node(&mut graph, 700);
    let b = node(&mut graph, 200);
    graph.nodes[0].is_circular = true;
    graph.nodes[0].add_edge(false, GraphEdge::new(b, 31, 10));
    let before = format!("{graph:?}");
    resolve(&mut graph, &PairLinkMatrix::new());
    assert_eq!(
        format!("{graph:?}"),
        before,
        "incoming edge from {a:?} must protect short target"
    );
}

#[test]
fn g1_failed_full_validation_preserves_arena_and_entire_graph() {
    let mut graph = CdbgGraph::new(31);
    let a = node(&mut graph, 700);
    let b = node(&mut graph, 200);
    edge(&mut graph, a, b);
    graph.nodes[1].edges_rev[0].support = 11;
    let before = format!("{graph:?}");
    resolve(&mut graph, &PairLinkMatrix::new());
    assert_eq!(format!("{graph:?}"), before);
    assert_eq!(graph.get_unitig_seq(0), vec![b'A'; 700]);
    assert_eq!(graph.get_unitig_seq(1), vec![b'A'; 200]);
}

#[test]
fn g2_validation_failure_keeps_shared_short_path_available() {
    let mut graph = CdbgGraph::new(31);
    let bad = node(&mut graph, 750);
    let repeat = node(&mut graph, 100);
    let target = node(&mut graph, 600);
    let alternative = node(&mut graph, 650);
    edge(&mut graph, bad, repeat);
    edge(&mut graph, repeat, target);
    edge(&mut graph, alternative, repeat);
    graph.nodes[1].edges_rev[0].support = 11;
    let bad_before = format!("{:?}", graph.nodes[0]);
    resolve(&mut graph, &PairLinkMatrix::new());
    assert_eq!(format!("{:?}", graph.nodes[0]), bad_before);
    assert_eq!(graph.nodes.len(), 5);
    assert!(!graph.nodes[4].is_deleted);
    assert_eq!(graph.get_unitig_seq(4), vec![b'A'; 1288]);
    assert!(graph.nodes[1..4].iter().all(|n| n.is_deleted));
}

#[test]
fn g6_legal_short_repeat_reuse_commits_exact_sequence_once() {
    let mut graph = CdbgGraph::new(31);
    let a = node(&mut graph, 700);
    let repeat = node(&mut graph, 100);
    let b = node(&mut graph, 650);
    let c = node(&mut graph, 600);
    edge(&mut graph, a, repeat);
    edge(&mut graph, repeat, b);
    edge(&mut graph, b, repeat);
    edge(&mut graph, repeat, c);
    let mut links = PairLinkMatrix::new();
    for _ in 0..10 {
        links.add_link(a, b);
        links.add_link(b, c);
    }
    resolve(&mut graph, &links);
    assert_eq!(graph.nodes.len(), 5);
    assert!(graph.nodes[..4].iter().all(|n| n.is_deleted));
    assert!(!graph.nodes[4].is_deleted);
    assert_eq!(graph.get_unitig_seq(4), vec![b'A'; 2026]);
}

#[test]
fn g6_long_unique_node_reuse_remains_rejected() {
    let mut graph = CdbgGraph::new(31);
    let a = node(&mut graph, 700);
    let b = node(&mut graph, 650);
    edge(&mut graph, a, b);
    edge(&mut graph, b, a);
    let path = vec![
        PathStep {
            dir: a,
            incoming_edge: None,
            overlap: 0,
            support: 0,
        },
        PathStep {
            dir: b,
            incoming_edge: Some(graph.nodes[0].edges_fwd[0]),
            overlap: 31,
            support: 10,
        },
        PathStep {
            dir: a,
            incoming_edge: Some(graph.nodes[1].edges_fwd[0]),
            overlap: 31,
            support: 10,
        },
    ];
    let before = format!("{graph:?}");
    assert!(validate_and_splice_path(&graph, &path).is_err());
    assert_eq!(format!("{graph:?}"), before);
}
