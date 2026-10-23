//! Phase F23 Option A — Unique Anchor Invariant gate.
//!
//! Forensic target: R1-001 `contig_16` chimeric relocation (runtime decision
//! `DEC_000022`, evidence chain F20 -> F21-A -> F21-B -> F22-D). While
//! extending through an unbranched chain of short (< 500 bp) repeat nodes,
//! the search anchor degrades onto collapsed shared repeat sequence. Candidate
//! exits are then scored with paired-end links anchored at that repeat node,
//! whose reads carry no locus specificity: ectopic pairs from the homologous
//! genomic copy out-vote the authentic exit and a chimeric join across
//! megabases is committed.
//!
//! Fix contract (`docs/development/evidence/C0_02R/phase_f22/
//! counterfactual_strategy_matrix.md`, Option A):
//!   Rule 1: fork candidate scoring anchors at the most recent long
//!           (>= 500 bp) node of the walk.
//!   Rule 2: repeat-internal (short-to-short) PE links never decide a fork.
//!   Rule 3: a fork beyond the unique anchor's insert-size horizon is
//!           withheld (`WITHHELD_NO_UNIQUE_ANCHOR`), not spliced on
//!           repeat-internal evidence.
//!
//! Topology (both scenarios, overlap 30, k = 31), node ids in comment:
//!
//!   S(600) -> C1 -> C2 -> C3 -> F -> { W(600), T(600) }
//!   0        1     2     3     4     5        6
//!
//! C1..C3 and F are short shared-repeat nodes; W is the ectopic exit whose
//! link evidence originates from the homologous copy, T the authentic exit.
//! The repeat-internal rows are keyed by C3 — the scoring anchor the
//! unpatched code uses — which is exactly the defect: those rows cannot
//! distinguish which genomic copy of C3 they originate from.
//!
//! ON UNPATCHED PRODUCTION CODE both tests FAIL (RED): the ectopic exit wins
//! 16 : 6 with weighted ratio 2.67 >= decision_threshold 2.5 and the chimeric
//! join S..F..W is committed.

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

    let _ = (f, w, t);

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

/// Scenario 1 — the DEC_000022 geometry: the unique anchor's horizon is
/// exhausted before the fork (anchor gap 3 x (200 - 30) = 510 bp, exit
/// separation 140 + 510 = 650 bp), so no unique-anchored pair can be
/// in-window. Repeat-internal rows anchor the ectopic exit W with 16 pairs
/// against the authentic exit T with 6. Unpatched code commits W: a chimeric
/// contig containing both S-unique and W-unique sequence. The contract
/// requires a clean withhold: no join through F at all.
#[test]
fn fork_beyond_unique_anchor_horizon_withholds_instead_of_chimera() {
    let mut fx = build_fixture(200, 200);

    // Repeat-internal rows (read1 inside the shared chain node C3): exactly
    // the evidence class that carries no locus specificity.
    add_rows(&mut fx.matrix, u(C3), u(W), 0..16, 120, 120);
    add_rows(&mut fx.matrix, u(C3), u(T), 16..22, 120, 120);
    // Unique-anchored rows exist physically but span 120 + 650 + 120 = 890 bp,
    // outside [250, 700]: the window check, not row absence, must exclude
    // them (Rule 3 via arithmetic, not via missing data).
    add_rows(&mut fx.matrix, u(0), u(W), 22..26, 120, 120);
    add_rows(&mut fx.matrix, u(0), u(T), 26..30, 120, 120);

    resolve_repeats_exspander_with_stats(
        &mut fx.graph,
        &fx.matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

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
        "REGRESSION FAILURE (F23 unique-anchor gate, horizon case): a contig \
         splices the unique anchor region onto a fork exit across the shared \
         repeat chain (S+W chimera: {chimera_sw}, S+T: {chimera_st}). \
         Unpatched code out-votes the authentic exit with ectopic \
         repeat-internal pairs (16:6, weighted ratio 2.67 >= 2.5) — the \
         DEC_000022 chimeric join mechanism. A fork beyond the unique \
         anchor's insert horizon must be withheld, not spliced.",
    );
    assert!(
        t_retained && w_retained,
        "withholding must conserve both exits: T retained = {t_retained}, \
         W retained = {w_retained}",
    );
}

/// Scenario 2 — anchor within horizon: anchor gap 3 x (100 - 30) = 210 bp,
/// exit separation 90 + 210 = 300 bp, so unique-anchored evidence for the
/// authentic exit T (8 pairs, span 100 + 300 + 60 = 460 bp) is admissible
/// and must decide the fork. The ectopic repeat-internal rows (16 vs 6)
/// must be ignored entirely (Rule 2).
#[test]
fn fork_within_unique_anchor_horizon_resolves_on_anchored_evidence() {
    let mut fx = build_fixture(100, 150);

    // Repeat-internal rows: ectopic W 16, authentic T 6 — the unpatched code
    // decides this fork on exactly these.
    add_rows(&mut fx.matrix, u(C3), u(W), 0..16, 100, 60);
    add_rows(&mut fx.matrix, u(C3), u(T), 16..22, 100, 60);
    // Unique-anchored rows: in-window evidence for T only. W has none: a
    // unique anchor never pairs with the other locus's exit.
    add_rows(&mut fx.matrix, u(0), u(T), 22..30, 100, 60);

    resolve_repeats_exspander_with_stats(
        &mut fx.graph,
        &fx.matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    let seqs = active_sequences(&fx.graph);
    let resolved = seqs
        .iter()
        .any(|seq| contains(seq, &fx.s_tail) && contains(seq, &fx.t_tail));
    let chimera = seqs
        .iter()
        .any(|seq| contains(seq, &fx.s_tail) && contains(seq, &fx.w_tail));

    assert!(
        !chimera,
        "REGRESSION FAILURE (F23 unique-anchor gate, horizon-inside case): a \
         contig splices the unique anchor onto the ectopic exit W. Unpatched \
         code decides the fork on repeat-internal pairs (16:6) even though \
         unique-anchored in-window evidence for T exists.",
    );
    assert!(
        resolved,
        "the authentic exit T must be resolved onto the anchor through the \
         short chain when unique-anchored in-window evidence exists; got no \
         such contig (fragmentation or withhold instead of resolution)",
    );
}
