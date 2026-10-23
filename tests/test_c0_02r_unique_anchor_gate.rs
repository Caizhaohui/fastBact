//! Phase F23 horizon guard — repeat-join withholding gate.
//!
//! Forensic target: R1-001 `contig_16` chimeric relocation (runtime decision
//! `DEC_000022`, evidence chain F20 -> F21-A -> F21-B -> F22-D). While
//! extending through an unbranched chain of short (< 500 bp) repeat nodes, the
//! search anchor degrades onto collapsed shared repeat sequence. When the walk
//! stands beyond the last long (>= 500 bp) anchor's insert horizon at a
//! branching repeat, unique-anchored evidence is arithmetically impossible and
//! the only remaining links are repeat-internal ones, which cannot distinguish
//! which genomic copy they originate from: ectopic pairs from the homologous
//! copy then out-vote the authentic exit and a chimeric join across megabases
//! is committed.
//!
//! Contract (F23 horizon guard, narrowed per Codex root-cause review
//! `codex_f23_root_cause_analysis.md`):
//!   1. A branching fork whose overlap-adjusted anchor-to-fork distance
//!      exceeds `max_insert` is withheld (`withheld_no_unique_anchor`), not
//!      decided on repeat-internal evidence.
//!   2. Within the horizon, forks keep the pre-F23 `curr_dir`-anchored
//!      semantics unchanged — the accepted T26/T27/T30/T34 contracts. The
//!      known hazard of a genuinely multi-copy short anchor inside the
//!      horizon is a documented scope boundary, not covered by this guard.
//!
//! Topology (both scenarios, overlap 30, k = 31), node ids in comment:
//!
//!   S(600) -> C1 -> C2 -> C3 -> F -> { W(600), T(600) }
//!   0        1     2     3     4     5        6
//!
//! W is the ectopic exit whose link evidence originates from the homologous
//! copy; T the authentic exit. The repeat-internal rows are keyed by C3 —
//! where the unguarded code's scoring anchor stands — and cannot distinguish
//! genomic copies.
//!
//! ON UNPATCHED PRODUCTION CODE the horizon scenario FAILS (RED, jobs
//! 3010464/3010541): the ectopic exit wins 16 : 6 with weighted ratio
//! 2.67 >= decision_threshold 2.5 and the chimeric join S..F..W is committed.

use fastbact::common::sketch::splitmix64;
use fastbact::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, UnitigNode};
use fastbact::resolve::exspander::{resolve_repeats_exspander_with_stats, ExspanderParams};
use fastbact::resolve::mapper::{InsertSizeDist, PairLinkMatrix, ReadPairEvidence};

fn dna(seed: u64, len: usize) -> Vec<u8> {
    (0..len)
        .map(|i| b"ACGT"[(splitmix64(seed + i as u64) % 4) as usize])
        .collect()
}

/// Per-node input slots for `dna`, 4096 apart. `dna(s, len)` draws
/// `splitmix64(s + i)`, so sequences whose seed ranges overlap are shifts of
/// one another; every fixture sequence gets a disjoint 4096-wide slot
/// (all lengths are < 4096) to keep tails genuinely unique.
const fn slot(id: u64) -> u64 {
    id * 4096
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

fn add_rows(
    matrix: &mut PairLinkMatrix,
    u1: DirUnitigId,
    u2: DirUnitigId,
    pair_ids: std::ops::Range<u32>,
    dist_u1_exit: u32,
    dist_u2_entry: u32,
) {
    for pair_id in pair_ids {
        matrix.add_evidence(
            u1,
            u2,
            ReadPairEvidence {
                pair_id,
                dist_u1_exit,
                dist_u2_entry,
                r1_len: 150,
                r2_len: 150,
            },
        );
    }
}

fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    haystack.windows(needle.len()).any(|w| w == needle)
}

/// Node ids of the fixture, by creation order.
const C3: u32 = 3;
const W: u32 = 5;
const T: u32 = 6;

struct Fixture {
    graph: CdbgGraph,
    matrix: PairLinkMatrix,
    s_tail: Vec<u8>,
    w_tail: Vec<u8>,
    t_tail: Vec<u8>,
}

/// Builds S(600) -> C1 -> C2 -> C3 -> F -> {W(600), T(600)} with 30 bp
/// overlaps; `chain_len` / `fork_len` control the walked chain geometry.
fn build_fixture(chain_len: usize, fork_len: usize) -> Fixture {
    // Slots: 1000/1001 junctions, 0 = S, 1..=3 chain, 4 = F, 5 = W, 6 = T.
    let junction: Vec<u8> = dna(slot(1000), 30);
    let junction_b: Vec<u8> = dna(slot(1001), 30);

    let mut s_seq = dna(slot(0), 600 - 30);
    s_seq.extend_from_slice(&junction);

    let mut c_seqs: Vec<Vec<u8>> = Vec::new();
    for i in 0..3u64 {
        let mut seq = junction.clone();
        seq.extend(dna(slot(1 + i), chain_len - 60));
        seq.extend_from_slice(&junction);
        c_seqs.push(seq);
    }

    let mut f_seq = junction.clone();
    f_seq.extend(dna(slot(4), fork_len - 60));
    f_seq.extend_from_slice(&junction_b);

    let mut w_seq = junction_b.clone();
    w_seq.extend(dna(slot(5), 600 - 30));
    let mut t_seq = junction_b.clone();
    t_seq.extend(dna(slot(6), 600 - 30));

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

    Fixture {
        graph,
        matrix: PairLinkMatrix::new(),
        s_tail: s_seq[100..220].to_vec(),
        w_tail: w_seq[200..320].to_vec(),
        t_tail: t_seq[200..320].to_vec(),
    }
}

fn active_sequences(graph: &CdbgGraph) -> Vec<Vec<u8>> {
    graph
        .nodes
        .iter()
        .filter(|node| !node.is_deleted)
        .map(|node| graph.get_unitig_seq(node.id))
        .collect()
}

fn u(id: u32) -> DirUnitigId {
    DirUnitigId::new(id, false)
}

/// Horizon scenario — the DEC_000022 geometry: the anchor-to-fork distance is
/// 3 x (300 - 30) = 810 bp on the walked chain, plus the fork hop 200 - 30 =
/// 170, i.e. 980 bp beyond the 700 bp insert horizon. Repeat-internal rows
/// anchor the ectopic exit W with 16 pairs against the authentic exit T with
/// 6; unguarded code commits W (chimera S..F..W). The guard must withhold the
/// fork and conserve both exits.
#[test]
fn fork_beyond_unique_anchor_horizon_withholds_instead_of_chimera() {
    let mut fx = build_fixture(300, 200);

    // Repeat-internal rows (read1 inside the shared chain node C3): exactly
    // the evidence class that carries no locus specificity. 16 : 6 with
    // weighted ratio 2.67 would commit W unguarded.
    add_rows(&mut fx.matrix, u(C3), u(W), 0..16, 120, 120);
    add_rows(&mut fx.matrix, u(C3), u(T), 16..22, 120, 120);

    let stats = resolve_repeats_exspander_with_stats(
        &mut fx.graph,
        &fx.matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    assert_eq!(
        stats.withheld_no_unique_anchor, 1,
        "the beyond-horizon fork must be withheld by the horizon guard"
    );

    let seqs = active_sequences(&fx.graph);
    let chimera_sw = seqs
        .iter()
        .any(|seq| contains(seq, &fx.s_tail) && contains(seq, &fx.w_tail));
    let chimera_st = seqs
        .iter()
        .any(|seq| contains(seq, &fx.s_tail) && contains(seq, &fx.t_tail));
    let t_retained = seqs.iter().any(|seq| contains(seq, &fx.t_tail));
    let w_retained = seqs.iter().any(|seq| contains(seq, &fx.w_tail));

    assert!(
        !chimera_sw && !chimera_st,
        "REGRESSION FAILURE (F23 horizon guard): a contig splices the unique \
         anchor region onto a fork exit across the shared repeat chain (S+W: \
         {chimera_sw}, S+T: {chimera_st}). Unguarded code out-votes the \
         authentic exit with ectopic repeat-internal pairs (16:6, weighted \
         ratio 2.67 >= 2.5) — the DEC_000022 chimeric join mechanism. A fork \
         beyond the unique anchor's insert horizon must be withheld.",
    );
    assert!(
        t_retained && w_retained,
        "withholding must conserve both exits: T retained = {t_retained}, \
         W retained = {w_retained}",
    );
}

/// Within-horizon passthrough control — the guard must be inert when the fork
/// stays inside the anchor's horizon: chain 3 x (100 - 30) = 210 plus the
/// fork hop 150 - 30 = 120 gives 330 <= 700. The accepted pre-F23 contract
/// (T26/T27/T30/T34 family) requires the fork to be decided on repeat-
/// internal, curr_dir-anchored evidence in that regime; the documented
/// multi-copy-anchor hazard inside the horizon is a scope boundary of this
/// guard, not a regression.
#[test]
fn fork_within_unique_anchor_horizon_keeps_pre_f23_semantics() {
    let mut fx = build_fixture(100, 150);

    add_rows(&mut fx.matrix, u(C3), u(W), 0..16, 100, 60);
    add_rows(&mut fx.matrix, u(C3), u(T), 16..22, 100, 60);

    let stats = resolve_repeats_exspander_with_stats(
        &mut fx.graph,
        &fx.matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    assert_eq!(
        stats.withheld_no_unique_anchor, 0,
        "the within-horizon fork must not be intercepted by the horizon guard"
    );

    let seqs = active_sequences(&fx.graph);
    let bridged = seqs
        .iter()
        .any(|seq| contains(seq, &fx.s_tail) && contains(seq, &fx.w_tail));
    let t_retained = seqs.iter().any(|seq| contains(seq, &fx.t_tail));

    assert!(
        bridged,
        "REGRESSION FAILURE (F23 horizon guard over-reach): a within-horizon \
         fork was withheld or fragmented. Pre-F23 semantics commit the \
         repeat-internal 16:6 bridge at fork distance 330 <= 700; the guard \
         must stay inert inside the horizon (accepted T-series contract).",
    );
    assert!(
        t_retained,
        "the non-selected exit must be conserved by the commit path"
    );
}
