//! Temporary F23 probe (not part of the permanent suite): dump active-node
//! composition for the unique-anchor fixture to explain the S+T carrier.

use fastbact::common::sketch::splitmix64;
use fastbact::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, UnitigNode};
use fastbact::resolve::exspander::{resolve_repeats_exspander_with_stats, ExspanderParams};
use fastbact::resolve::mapper::{InsertSizeDist, PairLinkMatrix, ReadPairEvidence};

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
    graph.nodes[to.id() as usize].add_edge(
        to.twin().is_rc(),
        GraphEdge::new(from.twin(), overlap, 10),
    );
}

#[test]
fn probe_topology_scenario1() {
    let sl = |id: u64| id * 4096;
    let junction: Vec<u8> = dna(sl(1000), 30);
    let junction_b: Vec<u8> = dna(sl(1001), 30);

    let mut s_seq = dna(sl(0), 570);
    s_seq.extend_from_slice(&junction);
    let mut c_seqs = Vec::new();
    for i in 0..3u64 {
        let mut seq = junction.clone();
        seq.extend(dna(sl(1 + i), 140));
        seq.extend_from_slice(&junction);
        c_seqs.push(seq);
    }
    let mut f_seq = junction.clone();
    f_seq.extend(dna(sl(4), 140));
    f_seq.extend_from_slice(&junction_b);
    let mut w_seq = junction_b.clone();
    w_seq.extend(dna(sl(5), 570));
    let mut t_seq = junction_b.clone();
    t_seq.extend(dna(sl(6), 570));

    let mut graph = CdbgGraph::new(31);
    let s = add_node(&mut graph, &s_seq);
    let mut prev = s;
    for seq in &c_seqs {
        let id = add_node(&mut graph, seq);
        add_edge(&mut graph, prev, id, 30);
        prev = id;
    }
    let f = add_node(&mut graph, &f_seq);
    add_edge(&mut graph, prev, f, 30);
    let w = add_node(&mut graph, &w_seq);
    let t = add_node(&mut graph, &t_seq);
    add_edge(&mut graph, f, w, 30);
    add_edge(&mut graph, f, t, 30);

    let s_tail = s_seq[100..220].to_vec();
    let w_tail = w_seq[200..320].to_vec();
    let t_tail = t_seq[200..320].to_vec();

    let mut matrix = PairLinkMatrix::new();
    let ev = |pair_id, d1, d2| ReadPairEvidence { pair_id, dist_u1_exit: d1, dist_u2_entry: d2, r1_len: 150, r2_len: 150 };
    for i in 0..16u32 {
        matrix.add_evidence(DirUnitigId::new(3, false), w, ev(i, 120, 120));
    }
    for i in 16..22u32 {
        matrix.add_evidence(DirUnitigId::new(3, false), t, ev(i, 120, 120));
    }
    for i in 22..26u32 {
        matrix.add_evidence(s, w, ev(i, 120, 120));
    }
    for i in 26..30u32 {
        matrix.add_evidence(s, t, ev(i, 120, 120));
    }

    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("resolution failed");

    println!("stats: k={} k-1={} other={} amb={} inv={}",
        stats.valid_k_joins, stats.valid_k_minus_1_joins, stats.valid_other_joins,
        stats.rejected_ambiguous_joins, stats.rejected_invalid_joins);

    for node in &graph.nodes {
        if node.is_deleted {
            println!("node {:>2}: DELETED", node.id);
            continue;
        }
        let seq = graph.get_unitig_seq(node.id);
        let has = |tail: &[u8]| seq.windows(tail.len()).any(|x| x == tail);
        println!(
            "node {:>2}: len {:>5} s_tail={} w_tail={} t_tail={} head={:?}",
            node.id,
            node.seq_len,
            has(&s_tail),
            has(&w_tail),
            has(&t_tail),
            std::str::from_utf8(&seq[..20.min(seq.len())]).unwrap_or("?"),
        );
    }
}
