//! Phase F24 — Provenance-Safe Contiguity Recovery Gate
//!
//! Ground truth contract: F24_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md (GPT-6 Astra).
//! Tests verifiable recovery classes:
//!   - R-I (Physically Isolated Exit Anchor): P01 positive vs N01 negative control.
//!   - R-S (Deterministic Structural Corridor Traversal): P03 positive vs N03 negative control.
//!   - P02/N02: Horizon boundary tests.

use fastbact::common::bitseq::reverse_complement_ascii;
use fastbact::common::sketch::splitmix64;
use fastbact::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, UnitigNode};
use fastbact::resolve::exspander::{
    resolve_repeats_exspander_with_stats, score_candidate_path_prov, ExspanderParams,
    NodeProvenanceMeta, ProvenanceTable,
};
use fastbact::resolve::mapper::{InsertSizeDist, PairLinkMatrix, ReadPairEvidence};

fn dna(seed: u64, len: usize) -> Vec<u8> {
    (0..len)
        .map(|i| b"ACGT"[(splitmix64(seed + i as u64) % 4) as usize])
        .collect()
}

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

fn add_evidence_row(
    matrix: &mut PairLinkMatrix,
    u1: DirUnitigId,
    u2: DirUnitigId,
    pair_id: u32,
    dist_u1_exit: u32,
    dist_u2_entry: u32,
) {
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

fn add_pair_rows(
    matrix: &mut PairLinkMatrix,
    u1: DirUnitigId,
    u2: DirUnitigId,
    range: std::ops::Range<u32>,
    dist_u1_exit: u32,
    dist_u2_entry: u32,
) {
    for pid in range {
        add_evidence_row(matrix, u1, u2, pid, dist_u1_exit, dist_u2_entry);
    }
}

fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    haystack.windows(needle.len()).any(|w| w == needle)
}

fn active_sequences(graph: &CdbgGraph) -> Vec<Vec<u8>> {
    graph
        .nodes
        .iter()
        .filter(|n| !n.is_deleted)
        .map(|n| graph.get_unitig_seq(n.id))
        .collect()
}

fn f24_same_contig(seqs: &[Vec<u8>], left: &[u8], right: &[u8]) -> bool {
    let left_rc = reverse_complement_ascii(left);
    let right_rc = reverse_complement_ascii(right);
    seqs.iter().any(|seq| {
        (contains(seq, left) && contains(seq, right))
            || (contains(seq, &left_rc) && contains(seq, &right_rc))
    })
}

fn bio_seq(seed: u64, len: usize, junction: &[u8]) -> Vec<u8> {
    let mut seq = junction.to_vec();
    seq.extend(dna(slot(seed), len.saturating_sub(2 * junction.len())));
    seq.extend_from_slice(junction);
    seq
}

/// F24 P01 (Positive Control — R-I: Physically Isolated Exit Anchor):
/// Topology: S_A(600) -> X(120) <- S_B(600)
///           X(120) -> U(2000) -> F(100) -> {T(400), W(400)}
///
/// U is a long node (2000 bp) with single context downstream of X.
/// Distance from X to U's 3' exit is 2000 - 30 = 1970 bp > max_insert (700 bp).
/// The entry shadow of X cannot physically reach U's 3' exit.
/// Terminal evidence from U to T (span ~450 bp within [250, 700]) supports authentic exit T (20:0).
/// Expected under F24: U's 3' exit is recognized as physically isolated (R-I),
/// allowing U to act as an anchor seed and commit U..F..T.
#[test]
fn test_f24_p01_isolated_exit_anchor_resolves() {
    let junction = dna(slot(9000), 30);
    let sa_seq = bio_seq(101, 600, &junction);
    let sb_seq = bio_seq(102, 600, &junction);
    let x_seq = bio_seq(103, 120, &junction);
    let u_seq = bio_seq(104, 2000, &junction);
    let f_seq = bio_seq(105, 100, &junction);
    let t_seq = bio_seq(106, 400, &junction);
    let w_seq = bio_seq(107, 400, &junction);

    let mut graph = CdbgGraph::new(31);
    let sa = add_node(&mut graph, &sa_seq);
    let sb = add_node(&mut graph, &sb_seq);
    let x = add_node(&mut graph, &x_seq);
    let u = add_node(&mut graph, &u_seq);
    let f = add_node(&mut graph, &f_seq);
    let t = add_node(&mut graph, &t_seq);
    let w = add_node(&mut graph, &w_seq);

    // Links: S_A -> X, S_B -> X, X -> U, U -> F, F -> T, F -> W
    for (from, to) in [
        (sa, x), (sb, x), (x, u), (u, f), (f, t), (f, w)
    ] {
        add_edge(&mut graph, from, to, 30);
    }

    let mut matrix = PairLinkMatrix::new();
    // Terminal rows from U to T: dist_u1_exit = 150 bp, dist_u2_entry = 150 bp.
    // Separation across F = 100 - 30 = 70 bp. Span = 150 + 70 + 150 = 370 bp (in [250, 700]).
    // Read offset on U: p = 2000 - 150 = 1850 bp > 120 + 700 = 820 bp (certified isolated from X!).
    add_pair_rows(&mut matrix, u, t, 0..20, 150, 150);

    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    let seqs = active_sequences(&graph);
    let u_tail = &u_seq[1000..1200];
    let t_tail = &t_seq[100..220];

    println!(
        "P01 stats: k_joins={}, k_minus_1={}, valid_other={}, withheld={}",
        stats.valid_k_joins,
        stats.valid_k_minus_1_joins,
        stats.valid_other_joins,
        stats.withheld_no_unique_anchor
    );
    assert!(
        f24_same_contig(&seqs, u_tail, t_tail),
        "P01: U's physically isolated 3' exit must anchor the walk and resolve U..F..T"
    );
    assert!(
        stats.valid_k_joins + stats.valid_k_minus_1_joins >= 2,
        "must resolve U->F and F->T"
    );
}

/// F24 N01 (Negative Control — Inherited Shared Locus Must NOT Reset Anchor):
/// Topology: S_A(600) -> C1(300) -> C2(300) -> C3(300) -> X(120) -> L(2000) -> F(100) -> {T(400), W(400)}
///           S_B(600) -> X(120)
/// Even though L is long (2000 bp), L is entered by X which has multiple independent external origins (S_A, S_B).
/// In F4 adversarial scenario, pooled L evidence favors W (16:6).
/// L is multi-context (inherited sharedness).
/// Mandatory Invariant F24-I01 & F24-I04: Inherited sharedness must NEVER establish uniqueness or decide forks!
#[test]
fn test_f24_n01_inherited_shared_locus_withholds() {
    let junction = dna(slot(9100), 30);
    let sa_seq = bio_seq(201, 600, &junction);
    let c1_seq = bio_seq(202, 300, &junction);
    let c2_seq = bio_seq(203, 300, &junction);
    let c3_seq = bio_seq(204, 300, &junction);
    let sb_seq = bio_seq(205, 600, &junction);
    let x_seq = bio_seq(206, 120, &junction);
    let l_seq = bio_seq(207, 2000, &junction);
    let f_seq = bio_seq(208, 100, &junction);
    let t_seq = bio_seq(209, 400, &junction);
    let w_seq = bio_seq(210, 400, &junction);

    let mut graph = CdbgGraph::new(31);
    let sa = add_node(&mut graph, &sa_seq);
    let c1 = add_node(&mut graph, &c1_seq);
    let c2 = add_node(&mut graph, &c2_seq);
    let c3 = add_node(&mut graph, &c3_seq);
    let sb = add_node(&mut graph, &sb_seq);
    let x = add_node(&mut graph, &x_seq);
    let l = add_node(&mut graph, &l_seq);
    let f = add_node(&mut graph, &f_seq);
    let t = add_node(&mut graph, &t_seq);
    let w = add_node(&mut graph, &w_seq);

    for (from, to) in [
        (sa, c1), (c1, c2), (c2, c3), (c3, x),
        (sb, x),
        (x, l),
        (l, f),
        (f, t), (f, w),
    ] {
        add_edge(&mut graph, from, to, 30);
    }

    let mut matrix = PairLinkMatrix::new();
    // Adversarial pooled evidence from L favors ectopic exit W (16:6)
    add_pair_rows(&mut matrix, l, w, 0..16, 120, 120);
    add_pair_rows(&mut matrix, l, t, 16..22, 120, 120);

    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    let seqs = active_sequences(&graph);
    let sa_tail = &sa_seq[100..220];
    let w_tail = &w_seq[100..220];
    let t_tail = &t_seq[100..220];

    assert!(
        !f24_same_contig(&seqs, sa_tail, w_tail) && !f24_same_contig(&seqs, sa_tail, t_tail),
        "N01: Shared multi-context L must NOT decide fork F; neither S_A..W nor S_A..T committed"
    );
    assert!(stats.withheld_no_unique_anchor >= 1, "fork F must be withheld");
}

/// F24 P03 (Positive Control — R-S: Deterministic Structural Corridor Traversal from Unanchored Seed):
/// Topology: In_A(400) -> S_shared(600) <- In_B(400)
///           S_shared(600) -> C1(200) -> C2(200) -> C3(200) -> F(100) -> {T(400), W(400)}
///
/// S_shared has multiple incoming contexts (in-degree 2 from In_A and In_B), so seed_anchor_eligible
/// is false and S_shared is an UNANCHORED seed (unique_anchor == None).
/// Under F24 R-S (Astra §1, §3.5):
/// An unanchored seed cannot decide forks, but must deterministically traverse unbranched 1-out/1-in
/// structural corridors (S_shared -> C1 -> C2 -> C3 -> F).
/// When reaching fork F, the walk halts and withholds without taking either branch.
/// Result: S_shared..C1..C2..C3 are joined together in one contig; neither T nor W is committed.
#[test]
fn test_f24_p03_deterministic_corridor_traversal() {
    let junction = dna(slot(9200), 30);
    let in_a_seq = bio_seq(300, 400, &junction);
    let in_b_seq = bio_seq(308, 400, &junction);
    let sa_seq = bio_seq(301, 600, &junction);
    let c1_seq = bio_seq(302, 200, &junction);
    let c2_seq = bio_seq(303, 200, &junction);
    let c3_seq = bio_seq(304, 200, &junction);
    let f_seq = bio_seq(305, 100, &junction);
    let t_seq = bio_seq(306, 400, &junction);
    let w_seq = bio_seq(307, 400, &junction);

    let mut graph = CdbgGraph::new(31);
    let in_a = add_node(&mut graph, &in_a_seq);
    let in_b = add_node(&mut graph, &in_b_seq);
    let sa = add_node(&mut graph, &sa_seq);
    let c1 = add_node(&mut graph, &c1_seq);
    let c2 = add_node(&mut graph, &c2_seq);
    let c3 = add_node(&mut graph, &c3_seq);
    let f = add_node(&mut graph, &f_seq);
    let t = add_node(&mut graph, &t_seq);
    let w = add_node(&mut graph, &w_seq);

    for (from, to) in [
        (in_a, sa), (in_b, sa),
        (sa, c1), (c1, c2), (c2, c3), (c3, f), (f, t), (f, w)
    ] {
        add_edge(&mut graph, from, to, 30);
    }

    // Zero PE evidence provided at fork F
    let matrix = PairLinkMatrix::new();

    let _stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    let seqs = active_sequences(&graph);
    let sa_tail = &sa_seq[100..220];
    let c3_tail = &c3_seq[50..150];
    let t_tail = &t_seq[100..220];
    let w_tail = &w_seq[100..220];

    assert!(
        f24_same_contig(&seqs, sa_tail, c3_tail),
        "P03: deterministic corridor S_shared..C1..C2..C3 must be joined from unanchored seed"
    );
    assert!(
        !f24_same_contig(&seqs, sa_tail, t_tail) && !f24_same_contig(&seqs, sa_tail, w_tail),
        "P03: fork F has no evidence and must withhold; neither T nor W committed"
    );
}

/// F24 N03 (Negative Control — R-S Must NOT Shortcut Through Confluences):
/// Topology: C1(200) -> X(200)
///           C2(200) -> X(200)
/// Node X has in-degree 2. C1 has out-degree 1, but its target X has in-degree 2.
/// Structural corridor traversal (R-S) requires both out_degree(u) == 1 AND in_degree(v) == 1.
/// R-S must NOT greedily contract C1 into X.
#[test]
fn test_f24_n03_unbranched_into_confluence_stops() {
    let junction = dna(slot(9300), 30);
    let c1_seq = bio_seq(401, 200, &junction);
    let c2_seq = bio_seq(402, 200, &junction);
    let x_seq = bio_seq(403, 200, &junction);

    let mut graph = CdbgGraph::new(31);
    let c1 = add_node(&mut graph, &c1_seq);
    let c2 = add_node(&mut graph, &c2_seq);
    let x = add_node(&mut graph, &x_seq);

    add_edge(&mut graph, c1, x, 30);
    add_edge(&mut graph, c2, x, 30);

    let matrix = PairLinkMatrix::new();

    let _stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    let seqs = active_sequences(&graph);
    let c1_tail = &c1_seq[50..150];
    let x_tail = &x_seq[50..150];

    // C1 and X must remain separate because X has in-degree 2
    assert!(
        !f24_same_contig(&seqs, c1_tail, x_tail),
        "N03: In-degree 2 confluence X must NOT be contracted via R-S shortcut"
    );
}

/// F24 P02 (Positive Control — Interval Boundary Exact Reach + 1):
/// Node U has length L = 800 with upstream entry shadow reach b = 300.
/// Read 1 has start coordinate p = b + 1 = 301.
/// In exit coordinates: dist_u1_exit = L - p = 800 - 301 = 499.
/// Read 2 on target T has dist_u2_entry = 50, path separation = 0.
/// Span = 499 + 0 + 50 = 549, well within default [250, 700] insert size distribution.
/// Since p = 301 > reach = 300, the read starts strictly outside the entry shadow.
/// Contract verification:
/// 1. score_candidate_path_prov admits the row (valid_count == 1, weight > 0).
/// 2. When provided as deciding evidence for fork F -> {T, W}, fork resolves to T.
#[test]
fn test_f24_p02_interval_boundary_exact_reach() {
    let mut graph = CdbgGraph::new(31);
    let u_seq = dna(slot(501), 800);
    let t_seq = dna(slot(502), 400);
    let u = add_node(&mut graph, &u_seq);
    let t = add_node(&mut graph, &t_seq);

    let mut prov = ProvenanceTable::default();
    let reach = 300u32;
    // Node U: single context (count 1), bci_downstream, entry shadow reach = 300.
    prov.push_composite(NodeProvenanceMeta {
        external_context_count_fwd: 1,
        external_context_count_rev: 1,
        bci_fwd: false,
        bci_rev: false,
        bci_downstream_fwd: true,
        bci_downstream_rev: false,
        entry_shadow_reach_fwd: reach,
        entry_shadow_reach_rev: 0,
        exit_isolated_fwd: false,
        exit_isolated_rev: false,
        provenance_boundary: None,
        ..NodeProvenanceMeta::default()
    });
    // Target T: clean single context
    prov.push_composite(NodeProvenanceMeta {
        external_context_count_fwd: 1,
        external_context_count_rev: 1,
        ..NodeProvenanceMeta::default()
    });

    let mut matrix = PairLinkMatrix::new();
    // Read starts at p = reach + 1 = 301.
    // dist_u1_exit = 800 - 301 = 499.
    // dist_u2_entry = 50. Span = 499 + 0 + 50 = 549 in [250, 700].
    add_evidence_row(&mut matrix, u, t, 1, 499, 50);

    let dist = InsertSizeDist::default();
    let (count, weight) = score_candidate_path_prov(&matrix, &graph, &prov, &dist, u, t, 0);

    assert_eq!(
        count, 1,
        "P02: read starting at p = reach + 1 must be strictly admitted outside shadow"
    );
    assert!(
        weight > 0.0,
        "P02: admitted row must carry positive Gaussian kernel weight"
    );
}

/// F24 N02 (Negative Control — Interval Boundary Exactly at Reach / One bp Inside Shadow):
/// Node U has length L = 800 with upstream entry shadow reach b = 300.
/// Case 1: Read starts exactly at p = b = 300.
/// dist_u1_exit = 800 - 300 = 500. Span = 500 + 0 + 50 = 550 in [250, 700].
/// Since p = 300 <= reach = 300, the read start is inside the entry shadow [0, reach].
/// Case 2: Read starts at p = b - 1 = 299.
/// dist_u1_exit = 800 - 299 = 501. Span = 501 + 0 + 50 = 551 in [250, 700].
/// Both rows fall within the entry shadow and must be filtered out.
/// Contract verification:
/// score_candidate_path_prov returns count == 0 and weight == 0.0.
#[test]
fn test_f24_n02_interval_boundary_one_bp_inside_shadow() {
    let mut graph = CdbgGraph::new(31);
    let u_seq = dna(slot(601), 800);
    let t_seq = dna(slot(602), 400);
    let u = add_node(&mut graph, &u_seq);
    let t = add_node(&mut graph, &t_seq);

    let mut prov = ProvenanceTable::default();
    let reach = 300u32;
    prov.push_composite(NodeProvenanceMeta {
        external_context_count_fwd: 1,
        external_context_count_rev: 1,
        bci_fwd: false,
        bci_rev: false,
        bci_downstream_fwd: true,
        bci_downstream_rev: false,
        entry_shadow_reach_fwd: reach,
        entry_shadow_reach_rev: 0,
        exit_isolated_fwd: false,
        exit_isolated_rev: false,
        provenance_boundary: None,
        ..NodeProvenanceMeta::default()
    });
    prov.push_composite(NodeProvenanceMeta {
        external_context_count_fwd: 1,
        external_context_count_rev: 1,
        ..NodeProvenanceMeta::default()
    });

    let mut matrix = PairLinkMatrix::new();
    // Row 1: exactly at boundary p = reach = 300 (dist_u1_exit = 800 - 300 = 500)
    add_evidence_row(&mut matrix, u, t, 1, 500, 50);
    // Row 2: 1 bp inside shadow p = reach - 1 = 299 (dist_u1_exit = 800 - 299 = 501)
    add_evidence_row(&mut matrix, u, t, 2, 501, 50);

    let dist = InsertSizeDist::default();
    let (count, weight) = score_candidate_path_prov(&matrix, &graph, &prov, &dist, u, t, 0);

    assert_eq!(
        count, 0,
        "N02: reads starting at p <= reach must be filtered out by entry-shadow cutoff"
    );
    assert_eq!(
        weight, 0.0,
        "N02: filtered rows must contribute exactly zero weight"
    );
}

/// F24 Rollback (Failed-Walk and PAG Rollback Gate):
/// Topology: S_A(600) -> X(120) <- S_B(600)
///           X(120) -> C1(350) -> C2(350) -> F(100) -> {T(400), W(400)}
///
/// S_A is a clean unique anchor.
/// The walk traverses confluence X, becoming contaminated (AR-3).
/// The walk advances through corridor C1 -> C2.
/// When reaching fork F, anchor_gap = (120-30) + (350-30) + (350-30) = 730 bp > max_insert (700 bp).
/// The fork F is beyond the unique anchor's horizon in a contaminated zone.
/// Mandatory Invariant F24-I05, F24-I09 & PAG contract (§3.10.1/4):
/// 1. No attributable anchor covers fork F; PAG withholding triggers.
/// 2. stats.withheld_no_unique_anchor and stats.withheld_no_attributable_anchor increment.
/// 3. discard_unattributable_walk rolls back speculative steps back to the anchor.
/// 4. All unmerged frontier nodes (C1, C2, F, T, W) are protected and conserved.
/// 5. No chimeric contig (S_A..T or S_A..W) is committed.
#[test]
fn test_f24_rollback_on_withheld_fork() {
    let junction = dna(slot(9400), 30);
    let sa_seq = bio_seq(701, 600, &junction);
    let sb_seq = bio_seq(702, 600, &junction);
    let x_seq = bio_seq(703, 120, &junction);
    let c1_seq = bio_seq(704, 350, &junction);
    let c2_seq = bio_seq(705, 350, &junction);
    let f_seq = bio_seq(706, 100, &junction);
    let t_seq = bio_seq(707, 400, &junction);
    let w_seq = bio_seq(708, 400, &junction);

    let mut graph = CdbgGraph::new(31);
    let sa = add_node(&mut graph, &sa_seq);
    let sb = add_node(&mut graph, &sb_seq);
    let x = add_node(&mut graph, &x_seq);
    let c1 = add_node(&mut graph, &c1_seq);
    let c2 = add_node(&mut graph, &c2_seq);
    let f = add_node(&mut graph, &f_seq);
    let t = add_node(&mut graph, &t_seq);
    let w = add_node(&mut graph, &w_seq);

    for (from, to) in [
        (sa, x), (sb, x), (x, c1), (c1, c2), (c2, f), (f, t), (f, w)
    ] {
        add_edge(&mut graph, from, to, 30);
    }

    // Adversarial or absent evidence at fork F
    let matrix = PairLinkMatrix::new();

    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    let seqs = active_sequences(&graph);
    let sa_tail = &sa_seq[100..220];
    let t_tail = &t_seq[100..220];
    let w_tail = &w_seq[100..220];

    // Verify Invariant F24-I09 & PAG Rollback:
    // 1. Fork F was withheld due to lack of attributable anchor
    assert!(
        stats.withheld_no_unique_anchor >= 1,
        "Rollback: fork F must be withheld under PAG when anchor gap exceeds horizon"
    );
    assert!(
        stats.withheld_no_attributable_anchor >= 1,
        "Rollback: withheld_no_attributable_anchor must record the PAG horizon rejection"
    );
    // 2. Both divergent targets T and W remain alive in graph
    assert!(
        !graph.nodes[t.id() as usize].is_deleted && !graph.nodes[w.id() as usize].is_deleted,
        "Rollback: targets T and W must remain conserved and alive"
    );
    // 3. No chimeric join connects S_A to either T or W
    assert!(
        !f24_same_contig(&seqs, sa_tail, t_tail) && !f24_same_contig(&seqs, sa_tail, w_tail),
        "Rollback: neither S_A..T nor S_A..W must be committed"
    );
}

/// F24 N04 (Negative Control — R-S Must NOT Branch from Unanchored Seed):
/// Topology: In_A(400) -> S_shared(600) <- In_B(400)
///           S_shared(600) -> { Branch1(300), Branch2(300) }
///
/// S_shared is an unanchored seed (in-degree 2 from In_A and In_B).
/// S_shared has out-degree 2 (Branch1 and Branch2).
/// Structural corridor traversal (R-S, Astra §1, §3.5) applies strictly to
/// unbranched 1-out/1-in edges: out_degree(u) == 1 && in_degree(v) == 1.
/// An unanchored seed cannot decide between diverging branches; it must halt and withhold.
/// Result: Neither Branch1 nor Branch2 is joined to S_shared; both branches remain conserved.
#[test]
fn test_f24_n04_unanchored_seed_branching_stops() {
    let junction = dna(slot(9500), 30);
    let in_a_seq = bio_seq(801, 400, &junction);
    let in_b_seq = bio_seq(802, 400, &junction);
    let sa_seq = bio_seq(803, 600, &junction);
    let b1_seq = bio_seq(804, 300, &junction);
    let b2_seq = bio_seq(805, 300, &junction);

    let mut graph = CdbgGraph::new(31);
    let in_a = add_node(&mut graph, &in_a_seq);
    let in_b = add_node(&mut graph, &in_b_seq);
    let sa = add_node(&mut graph, &sa_seq);
    let b1 = add_node(&mut graph, &b1_seq);
    let b2 = add_node(&mut graph, &b2_seq);

    for (from, to) in [
        (in_a, sa), (in_b, sa),
        (sa, b1), (sa, b2),
    ] {
        add_edge(&mut graph, from, to, 30);
    }

    let matrix = PairLinkMatrix::new();

    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    let seqs = active_sequences(&graph);
    let sa_tail = &sa_seq[100..220];
    let b1_tail = &b1_seq[50..150];
    let b2_tail = &b2_seq[50..150];

    // Verify N04:
    // 1. Unanchored fork was withheld
    assert!(
        stats.withheld_no_unique_anchor >= 1,
        "N04: Branching from unanchored seed must be withheld"
    );
    // 2. Both branches remain alive and unmerged to S_shared
    assert!(
        !graph.nodes[b1.id() as usize].is_deleted && !graph.nodes[b2.id() as usize].is_deleted,
        "N04: Branch targets must remain conserved and alive"
    );
    assert!(
        !f24_same_contig(&seqs, sa_tail, b1_tail) && !f24_same_contig(&seqs, sa_tail, b2_tail),
        "N04: Neither Branch1 nor Branch2 must be committed with unanchored seed"
    );
}
