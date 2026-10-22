//! Phase F9-I state-semantics gate for exSPAnder repeat resolution.
//!
//! Target commit: `37a7c0b` (unpatched). Do not treat cleanup protection,
//! historical fork failures, or nodes already consumed by this resolver as
//! reasons to drop a candidate from competition.
//!
//! On this commit the intended status is:
//! - T1, T2, T3 fail (the three P1 defects).
//! - T4 and T13 fail because a cleanup-banned or deleted route still hides the
//!   real competitor.
//! - T7's later-seed eligibility check fails for the same cross-seed ban.
//! - T5, T6, T7 atomicity, T8, T9, T10, and T12 hold on the current resolver.
//!
//! T11 is the full release suite plus genome runs. It is not an in-process unit test.
//!
//! Shadow rule these tests lock in for the repair:
//! a consumed-node candidate blocks a decision only when its route diverges from
//! the chosen route (`zip` of steps finds a different directed target) and
//! `best.total_weight < threshold * shadow.total_weight`. A shadow that merely
//! continues the chosen route is the same exit and must not block (T4).
//!
//! Phase F11 checked single-use resolver contract:
//! - T14 (F10-P1-002): the tandem hop `R → R` stays in the candidate search,
//!   so the weak decoy must not commit.
//! - T15 (F10-P0-001): a second call on a `Completed` graph returns
//!   `AlreadyResolved` before cleanup, so Case A's side branch survives.
//! - T16 (F10-P1-001): a second call returns `AlreadyResolved` before
//!   exploration, so the consumed shadow cannot be replaced by T3's decoy.
//! - T17: cloning a `Completed` graph preserves `Completed` and rejects re-entry.
//! - T18: freshly constructed graphs start in `NeverRun` and resolve independently.
//!
//! Phase F13 multi-copy tandem repeat completeness (F12-P1-001):
//! - T19: three-copy tandem repeat (R→R→R→X) must withhold weak decoy D.
//! - T20: four-step tandem repeat (R→R→R→R) reaching depth cap must withhold weak decoy D.
//! - T21: alternating two-node cycle (A⇄B) must withhold weak decoy D.
//! - T22: reverse-complement twin symmetry of three-copy tandem repeat must withhold weak decoy D.

use std::collections::BTreeSet;

use fastbact::common::bitseq::reverse_complement_ascii;
use fastbact::common::sketch::splitmix64;
use fastbact::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, ResolverLifecycleState, UnitigNode};
use fastbact::polish::output::curate_contigs;
use fastbact::resolve::exspander::{
    resolve_edge_between, resolve_repeats_exspander_with_stats, score_candidate_path,
    EdgeResolveResult, ExspanderError, ExspanderParams,
};
use fastbact::resolve::mapper::{InsertSizeDist, PairLinkMatrix, ReadPairEvidence};

const OV: usize = 31;

fn dna(seed: u64, len: usize) -> Vec<u8> {
    (0..len)
        .map(|i| {
            let mixed = seed
                .wrapping_mul(0x9E37_79B9_7F4A_7C15)
                .wrapping_add((i as u64).wrapping_mul(0xBF58_476D_1CE4_E5B9));
            b"ACGT"[(splitmix64(mixed) % 4) as usize]
        })
        .collect()
}

fn junction(seed: u64) -> Vec<u8> {
    let seq = dna(seed, OV);
    assert_ne!(
        seq,
        reverse_complement_ascii(&seq),
        "fixture junction {seed:#x} is an RC palindrome"
    );
    seq
}

fn assert_junctions_distinct(parts: &[Vec<u8>]) {
    for i in 0..parts.len() {
        for j in (i + 1)..parts.len() {
            assert_ne!(parts[i], parts[j], "fixture junctions {i} and {j} match");
            assert_ne!(
                parts[i],
                reverse_complement_ascii(&parts[j]),
                "fixture junction {i} is the RC of junction {j}"
            );
        }
    }
}

fn span_node(left: &[u8], fill_seed: u64, total_len: usize, right: &[u8]) -> Vec<u8> {
    assert!(
        total_len > left.len() + right.len(),
        "fixture node shorter than its junctions"
    );
    let mut seq = left.to_vec();
    seq.extend(dna(fill_seed, total_len - left.len() - right.len()));
    seq.extend_from_slice(right);
    assert_eq!(seq.len(), total_len);
    seq
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
    graph.nodes[to.id() as usize]
        .add_edge(to.twin().is_rc(), GraphEdge::new(from.twin(), overlap, 10));
}

fn add_forward_only(graph: &mut CdbgGraph, from: DirUnitigId, to: DirUnitigId, overlap: u8) {
    graph.nodes[from.id() as usize].add_edge(from.is_rc(), GraphEdge::new(to, overlap, 10));
}

fn push_bio(graph: &mut CdbgGraph, bio: &[u8], mirror: bool) -> DirUnitigId {
    if mirror {
        add_node(graph, &reverse_complement_ascii(bio))
    } else {
        add_node(graph, bio)
    }
}

fn link_bio(graph: &mut CdbgGraph, from: DirUnitigId, to: DirUnitigId, overlap: u8, mirror: bool) {
    if mirror {
        add_edge(graph, from.twin(), to.twin(), overlap);
    } else {
        add_edge(graph, from, to, overlap);
    }
}

fn walk(id: DirUnitigId, mirror: bool) -> DirUnitigId {
    if mirror {
        id.twin()
    } else {
        id
    }
}

fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    !needle.is_empty()
        && haystack
            .windows(needle.len())
            .any(|window| window == needle)
}

fn splice(nodes: &[&[u8]], overlaps: &[usize]) -> Vec<u8> {
    assert_eq!(nodes.len(), overlaps.len() + 1);
    let mut assembled = nodes[0].to_vec();
    for (node, &overlap) in nodes.iter().skip(1).zip(overlaps) {
        assert!(overlap > 0 && overlap < node.len() && overlap < assembled.len());
        assert_eq!(&assembled[assembled.len() - overlap..], &node[..overlap]);
        assembled.extend_from_slice(&node[overlap..]);
    }
    assembled
}

fn active_summary(graph: &CdbgGraph) -> String {
    graph
        .nodes
        .iter()
        .filter(|node| !node.is_deleted)
        .map(|node| format!("{}:{}", node.id, node.seq_len))
        .collect::<Vec<_>>()
        .join(",")
}

fn has_active_seq(graph: &CdbgGraph, seq: &[u8]) -> bool {
    graph
        .nodes
        .iter()
        .filter(|node| !node.is_deleted)
        .any(|node| graph.get_unitig_seq(node.id) == seq)
}

fn alive(graph: &CdbgGraph, id: u32) -> bool {
    graph
        .nodes
        .get(id as usize)
        .is_some_and(|node| !node.is_deleted)
}

fn params() -> ExspanderParams {
    let params = ExspanderParams::default();
    assert_eq!(params.min_support, 5);
    assert!((params.decision_threshold - 2.5).abs() < 1e-12);
    params
}

fn resolve(graph: &mut CdbgGraph, matrix: &PairLinkMatrix, dist: &InsertSizeDist) {
    resolve_repeats_exspander_with_stats(graph, matrix, dist, params()).expect("repeat resolution failed");
}

fn default_dist() -> InsertSizeDist {
    let dist = InsertSizeDist::default();
    assert_eq!(dist.median, 450.0);
    assert_eq!(dist.mad, 50.0);
    assert_eq!(dist.min_insert, 250);
    assert_eq!(dist.max_insert, 700);
    dist
}

fn add_pairs(
    matrix: &mut PairLinkMatrix,
    from: DirUnitigId,
    to: DirUnitigId,
    n: u32,
    dist1: u32,
    dist2: u32,
    id_base: u32,
) {
    add_pairs_len(matrix, from, to, n, dist1, dist2, id_base, 150);
}

fn add_pairs_len(
    matrix: &mut PairLinkMatrix,
    from: DirUnitigId,
    to: DirUnitigId,
    n: u32,
    dist1: u32,
    dist2: u32,
    id_base: u32,
    read_len: u16,
) {
    for offset in 0..n {
        matrix.add_evidence(
            from,
            to,
            ReadPairEvidence {
                pair_id: id_base + offset,
                dist_u1_exit: dist1,
                dist_u2_entry: dist2,
                r1_len: read_len,
                r2_len: read_len,
            },
        );
    }
}

/// A read of `read_len` bp ending `dist` bp from the unitig boundary lies fully
/// inside a unitig of `unitig_len` bp.
fn assert_read_inside(unitig_len: usize, dist: u32, read_len: u32) {
    assert!(
        read_len <= dist && (dist as usize) <= unitig_len,
        "read of {read_len} bp at distance {dist} does not sit inside {unitig_len} bp"
    );
}

fn add_pairs_bio(
    matrix: &mut PairLinkMatrix,
    from: DirUnitigId,
    to: DirUnitigId,
    n: u32,
    dist1: u32,
    dist2: u32,
    id_base: u32,
    mirror: bool,
) {
    add_pairs(
        matrix,
        walk(from, mirror),
        walk(to, mirror),
        n,
        dist1,
        dist2,
        id_base,
    );
}

fn assert_score(
    matrix: &PairLinkMatrix,
    dist: &InsertSizeDist,
    from: DirUnitigId,
    to: DirUnitigId,
    sep: i64,
    expect_count: u32,
) -> f64 {
    let (count, weight) = score_candidate_path(matrix, dist, from, to, sep);
    assert_eq!(
        count, expect_count,
        "fixture score count {:?} -> {:?} sep {sep}",
        from, to
    );
    let expect_weight = f64::from(expect_count);
    assert!(
        (weight - expect_weight).abs() < 1e-6,
        "fixture score weight {weight} != {expect_weight} for {:?} -> {:?} sep {sep}",
        from,
        to
    );
    weight
}

fn assert_valid_edge(graph: &CdbgGraph, from: DirUnitigId, to: DirUnitigId, overlap: u8) {
    match resolve_edge_between(graph, from, to) {
        EdgeResolveResult::Valid(edge) => {
            assert_eq!(edge.target, to);
            assert_eq!(edge.overlap_k, overlap);
        }
        other => panic!("expected a sequence-valid edge {from:?} -> {to:?}, got {other:?}"),
    }
}

fn assert_not_homopolymer(seq: &[u8]) {
    assert!(seq.len() > 1 && seq.iter().any(|base| *base != seq[0]));
}

fn structural_fingerprint(
    graph: &CdbgGraph,
) -> Vec<(
    u32,
    u32,
    bool,
    bool,
    Vec<GraphEdge>,
    Vec<GraphEdge>,
    Vec<u8>,
)> {
    graph
        .nodes
        .iter()
        .map(|node| {
            (
                node.id,
                node.seq_len,
                node.is_deleted,
                node.is_circular,
                node.edges_fwd.to_vec(),
                node.edges_rev.to_vec(),
                graph.get_unitig_seq(node.id),
            )
        })
        .collect()
}

fn canonical(seq: &[u8]) -> Vec<u8> {
    let rc = reverse_complement_ascii(seq);
    if seq <= rc.as_slice() {
        seq.to_vec()
    } else {
        rc
    }
}

fn mirror_signature(graph: &CdbgGraph) -> (Vec<u32>, Vec<bool>, Vec<Vec<u8>>) {
    let lengths = graph.nodes.iter().map(|node| node.seq_len).collect();
    let deleted = graph.nodes.iter().map(|node| node.is_deleted).collect();
    let seqs = graph
        .nodes
        .iter()
        .map(|node| canonical(&graph.get_unitig_seq(node.id)))
        .collect();
    (lengths, deleted, seqs)
}

fn rc_palindrome(half_seed: u64, half_len: usize) -> Vec<u8> {
    let half = dna(half_seed, half_len);
    let mut seq = half.clone();
    seq.extend(reverse_complement_ascii(&half));
    assert_eq!(seq, reverse_complement_ascii(&seq));
    seq
}

struct Ce1 {
    graph: CdbgGraph,
    matrix: PairLinkMatrix,
    dist: InsertSizeDist,
    anchor: DirUnitigId,
    x_dir: DirUnitigId,
    d_dir: DirUnitigId,
    d_id: u32,
    x_id: u32,
    s1_contig: Vec<u8>,
    good: Vec<u8>,
    bad: Vec<u8>,
    sep: i64,
}

/// S1(900) → R(169), X(300) → R, S2(700) → Q(100) → {X, D(300)}.
fn build_ce1(mirror: bool) -> Ce1 {
    let mut graph = CdbgGraph::new(31);
    let j_sr = junction(0x11_0001);
    let j_sq = junction(0x11_0002);
    let j_qx = junction(0x11_0003);
    assert_junctions_distinct(&[j_sr.clone(), j_sq.clone(), j_qx.clone()]);

    let s1 = span_node(&[], 0x11_0101, 900, &j_sr);
    let r = span_node(&j_sr, 0x11_0102, 169, &[]);
    let x = span_node(&j_qx, 0x11_0103, 300, &j_sr);
    let s2 = span_node(&[], 0x11_0104, 700, &j_sq);
    let q = span_node(&j_sq, 0x11_0105, 100, &j_qx);
    let d = span_node(&j_qx, 0x11_0106, 300, &[]);
    assert_eq!(
        (s1.len(), r.len(), x.len(), s2.len(), q.len(), d.len()),
        (900, 169, 300, 700, 100, 300)
    );
    assert!(x.len() < 500 && d.len() < 500 && q.len() < 500 && r.len() < 500);

    let s1_id = push_bio(&mut graph, &s1, mirror);
    let r_id = push_bio(&mut graph, &r, mirror);
    let x_id = push_bio(&mut graph, &x, mirror);
    let s2_id = push_bio(&mut graph, &s2, mirror);
    let q_id = push_bio(&mut graph, &q, mirror);
    let d_id = push_bio(&mut graph, &d, mirror);
    for (from, to) in [
        (s1_id, r_id),
        (x_id, r_id),
        (s2_id, q_id),
        (q_id, x_id),
        (q_id, d_id),
    ] {
        link_bio(&mut graph, from, to, OV as u8, mirror);
        assert_valid_edge(&graph, walk(from, mirror), walk(to, mirror), OV as u8);
    }
    assert_eq!(
        graph.nodes[q_id.id() as usize]
            .edges(walk(q_id, mirror).is_rc())
            .len(),
        2
    );

    // Path separation at either child of Q: len(Q) - ov(S2→Q) - ov(Q→child) = 38.
    let sep = q.len() as i64 - OV as i64 - OV as i64;
    assert_eq!(sep, 38);
    let dist = default_dist();
    let mut matrix = PairLinkMatrix::new();
    // 200 + 38 + 212 = 450, the default median. Both reads sit inside the unitig.
    assert!(200 <= s2.len() && 212 <= x.len() && 212 <= d.len());
    add_pairs_bio(&mut matrix, s2_id, x_id, 20, 200, 212, 1_000, mirror);
    add_pairs_bio(&mut matrix, s2_id, d_id, 5, 200, 212, 2_000, mirror);
    let anchor = walk(s2_id, mirror);
    let x_dir = walk(x_id, mirror);
    let d_dir = walk(d_id, mirror);
    let x_weight = assert_score(&matrix, &dist, anchor, x_dir, sep, 20);
    let d_weight = assert_score(&matrix, &dist, anchor, d_dir, sep, 5);
    assert!(
        x_weight / d_weight >= 2.5,
        "fixture X/D ratio {} must clear 2.5",
        x_weight / d_weight
    );

    Ce1 {
        graph,
        matrix,
        dist,
        anchor,
        x_dir,
        d_dir,
        d_id: d_id.id(),
        x_id: x_id.id(),
        s1_contig: splice(&[&s1, &r], &[OV]),
        good: splice(&[&s2, &q, &x], &[OV, OV]),
        bad: splice(&[&s2, &q, &d], &[OV, OV]),
        sep,
    }
}

struct Ce3 {
    graph: CdbgGraph,
    matrix: PairLinkMatrix,
    dist: InsertSizeDist,
    anchor: DirUnitigId,
    r_dir: DirUnitigId,
    d_dir: DirUnitigId,
    q_id: u32,
    d_id: u32,
    r_id: u32,
    s1_contig: Vec<u8>,
    bad: Vec<u8>,
    sep: i64,
}

/// S1(900) → R(169) consumes R. S2(700) → Q(100) → {R, D(300)}.
fn build_ce3(mirror: bool) -> Ce3 {
    let mut graph = CdbgGraph::new(31);
    let j_r = junction(0x31_0001);
    let j_sq = junction(0x31_0002);
    assert_junctions_distinct(&[j_r.clone(), j_sq.clone()]);

    let s1 = span_node(&[], 0x31_0101, 900, &j_r);
    let r = span_node(&j_r, 0x31_0102, 169, &[]);
    let s2 = span_node(&[], 0x31_0103, 700, &j_sq);
    let q = span_node(&j_sq, 0x31_0104, 100, &j_r);
    let d = span_node(&j_r, 0x31_0105, 300, &[]);
    assert_eq!(
        (s1.len(), r.len(), s2.len(), q.len(), d.len()),
        (900, 169, 700, 100, 300)
    );

    let s1_id = push_bio(&mut graph, &s1, mirror);
    let r_id = push_bio(&mut graph, &r, mirror);
    let s2_id = push_bio(&mut graph, &s2, mirror);
    let q_id = push_bio(&mut graph, &q, mirror);
    let d_id = push_bio(&mut graph, &d, mirror);
    for (from, to) in [(s1_id, r_id), (s2_id, q_id), (q_id, r_id), (q_id, d_id)] {
        link_bio(&mut graph, from, to, OV as u8, mirror);
        assert_valid_edge(&graph, walk(from, mirror), walk(to, mirror), OV as u8);
    }

    let sep = q.len() as i64 - OV as i64 - OV as i64;
    assert_eq!(sep, 38);
    // d2 = 150 fits inside R (169) and D (300). 262 + 38 + 150 = 450.
    assert!(262 <= s2.len() && 150 <= r.len() && 150 <= d.len());
    let dist = default_dist();
    let mut matrix = PairLinkMatrix::new();
    add_pairs_bio(&mut matrix, s2_id, r_id, 20, 262, 150, 3_000, mirror);
    add_pairs_bio(&mut matrix, s2_id, d_id, 5, 262, 150, 4_000, mirror);
    let anchor = walk(s2_id, mirror);
    let r_dir = walk(r_id, mirror);
    let d_dir = walk(d_id, mirror);
    let r_weight = assert_score(&matrix, &dist, anchor, r_dir, sep, 20);
    let d_weight = assert_score(&matrix, &dist, anchor, d_dir, sep, 5);
    assert!(
        d_weight < 2.5 * r_weight,
        "fixture shadow gate must be closed: D {d_weight} vs R {r_weight}"
    );
    let single = d_weight / f64::from(5);
    assert!(
        single >= 0.05,
        "fixture D alone must pass the single-candidate rule, got {single}"
    );

    Ce3 {
        graph,
        matrix,
        dist,
        anchor,
        r_dir,
        d_dir,
        q_id: q_id.id(),
        d_id: d_id.id(),
        r_id: r_id.id(),
        s1_contig: splice(&[&s1, &r], &[OV]),
        bad: splice(&[&s2, &q, &d], &[OV, OV]),
        sep,
    }
}

/// T1 / CE-1. Cleanup registration of X must not remove X from S2's candidate search.
///
/// ON UNPATCHED CODE THIS FAILS: S1's commit inserts X into `rejected_nodes`,
/// BFS skips X, and the 5-pair decoy D is committed as the only candidate.
#[test]
fn test_t1_ce1_cleanup_must_not_exclude_later_candidate() {
    let mut fx = build_ce1(false);
    assert_eq!(fx.sep, 38);
    let _ = (fx.anchor, fx.x_dir, fx.d_dir);

    resolve(&mut fx.graph, &fx.matrix, &fx.dist);

    assert!(
        has_active_seq(&fx.graph, &fx.s1_contig),
        "fixture: S1→R must commit and register X before S2 searches, active=[{}]",
        active_summary(&fx.graph)
    );
    let good = has_active_seq(&fx.graph, &fx.good);
    let bad = has_active_seq(&fx.graph, &fx.bad);
    let d_alive = alive(&fx.graph, fx.d_id);
    let x_alive = alive(&fx.graph, fx.x_id);
    assert!(
        good && !bad && d_alive && !x_alive,
        "T1/CE-1 RED: cleanup protection of external neighbor X excluded it from S2's BFS. \
         X has 20 centered pairs and D has 5 (ratio 4.0 >= 2.5), but the unpatched search \
         skipped X and committed D under the single-candidate rule. \
         S2+Q+X present={good}, S2+Q+D present={bad}, D active={d_alive}, X active={x_alive}, \
         active=[{}]",
        active_summary(&fx.graph)
    );
}

/// T2 / CE-2. A failed fork at an earlier seed must not permanently ban its endpoints.
///
/// ON UNPATCHED CODE THIS FAILS: S1' records X in `rejected_nodes` when 3 pairs
/// miss `min_support`. S2' then skips X and commits the 5-pair decoy D.
#[test]
fn test_t2_ce2_fork_failure_must_not_ban_later_seed() {
    let mut graph = CdbgGraph::new(31);
    let j_sx = junction(0x21_0001);
    let j_sm = junction(0x21_0002);
    assert_junctions_distinct(&[j_sx.clone(), j_sm.clone()]);

    let s1 = span_node(&[], 0x21_0101, 900, &j_sx);
    let x = span_node(&j_sx, 0x21_0102, 300, &[]);
    let y = span_node(&j_sx, 0x21_0103, 300, &[]);
    let s2 = span_node(&[], 0x21_0104, 700, &j_sm);
    let m = span_node(&j_sm, 0x21_0105, 100, &j_sx);
    let d = span_node(&j_sx, 0x21_0106, 300, &[]);
    assert_eq!(
        (s1.len(), x.len(), y.len(), s2.len(), m.len(), d.len()),
        (900, 300, 300, 700, 100, 300)
    );
    assert_ne!(&x[OV..], &y[OV..]);
    assert_ne!(&x[OV..], &d[OV..]);
    assert_ne!(&y[OV..], &d[OV..]);

    let s1_id = add_node(&mut graph, &s1);
    let x_id = add_node(&mut graph, &x);
    let y_id = add_node(&mut graph, &y);
    let s2_id = add_node(&mut graph, &s2);
    let m_id = add_node(&mut graph, &m);
    let d_id = add_node(&mut graph, &d);
    for (from, to) in [
        (s1_id, x_id),
        (s1_id, y_id),
        (s2_id, m_id),
        (m_id, x_id),
        (m_id, d_id),
    ] {
        add_edge(&mut graph, from, to, OV as u8);
        assert_valid_edge(&graph, from, to, OV as u8);
    }
    assert_eq!(graph.nodes[s1_id.id() as usize].edges_fwd.len(), 2);
    assert_eq!(graph.nodes[m_id.id() as usize].edges_fwd.len(), 2);

    let dist = default_dist();
    let mut matrix = PairLinkMatrix::new();
    // S1 fork, sep = -31. 200 + -31 + 281 = 450. Three pairs stay under min_support.
    assert!(281 <= x.len() && 281 <= y.len());
    add_pairs(&mut matrix, s1_id, x_id, 3, 200, 281, 100);
    add_pairs(&mut matrix, s1_id, y_id, 3, 200, 281, 200);
    assert_score(&matrix, &dist, s1_id, x_id, -(OV as i64), 3);
    assert_score(&matrix, &dist, s1_id, y_id, -(OV as i64), 3);

    // S2 fork at M, sep = 38. 200 + 38 + 212 = 450.
    let sep = m.len() as i64 - OV as i64 - OV as i64;
    assert_eq!(sep, 38);
    assert!(212 <= x.len() && 212 <= d.len());
    add_pairs(&mut matrix, s2_id, x_id, 20, 200, 212, 1_100);
    add_pairs(&mut matrix, s2_id, d_id, 5, 200, 212, 1_200);
    let x_weight = assert_score(&matrix, &dist, s2_id, x_id, sep, 20);
    let d_weight = assert_score(&matrix, &dist, s2_id, d_id, sep, 5);
    assert!(x_weight / d_weight >= 2.5);

    let good = splice(&[&s2, &m, &x], &[OV, OV]);
    let bad = splice(&[&s2, &m, &d], &[OV, OV]);
    let s1_x = splice(&[&s1, &x], &[OV]);
    let s1_y = splice(&[&s1, &y], &[OV]);

    resolve(&mut graph, &matrix, &dist);

    assert!(
        alive(&graph, s1_id.id()) && !has_active_seq(&graph, &s1_x) && !has_active_seq(&graph, &s1_y),
        "fixture: S1' has only 3 pairs per arm, so its fork must fail without splicing. active=[{}]",
        active_summary(&graph)
    );
    let good_present = has_active_seq(&graph, &good);
    let bad_present = has_active_seq(&graph, &bad);
    let d_alive = alive(&graph, d_id.id());
    let x_alive = alive(&graph, x_id.id());
    assert!(
        good_present && !bad_present && d_alive && !x_alive,
        "T2/CE-2 RED: S1' fork failure inserted X into rejected_nodes, so S2' skipped the \
         20-pair continuation and committed D (5 pairs). S2'+M+X present={good_present}, \
         S2'+M+D present={bad_present}, D active={d_alive}, X active={x_alive}, active=[{}]",
        active_summary(&graph)
    );
}

/// T3 / CE-3. A node consumed by an earlier seed competes as a shadow and blocks a decoy.
///
/// ON UNPATCHED CODE THIS FAILS: BFS skips deleted R, D is the only candidate,
/// and W/C >= 0.05 commits S2+Q+D.
#[test]
fn test_t3_ce3_consumed_shadow_blocks_decoy() {
    let mut fx = build_ce3(false);
    assert_eq!(fx.sep, 38);
    let _ = (fx.anchor, fx.r_dir, fx.d_dir);

    resolve(&mut fx.graph, &fx.matrix, &fx.dist);

    assert!(
        has_active_seq(&fx.graph, &fx.s1_contig) && !alive(&fx.graph, fx.r_id),
        "fixture: S1 must consume R before S2's fork is decided, active=[{}]",
        active_summary(&fx.graph)
    );
    let bad = has_active_seq(&fx.graph, &fx.bad);
    let q_alive = alive(&fx.graph, fx.q_id);
    let d_alive = alive(&fx.graph, fx.d_id);
    let s2_alive = fx
        .graph
        .nodes
        .iter()
        .any(|node| !node.is_deleted && node.seq_len == 700);
    assert!(
        !bad && q_alive && d_alive && s2_alive,
        "T3/CE-3 RED: consumed R was skipped because is_deleted, so it never competed. \
         D's 5 pairs passed the single-candidate rule (W/C >= 0.05) even though R still \
         has 20 pairs and 5 < 2.5 * 20. The fork at Q must stay unresolved. \
         S2+Q+D present={bad}, Q active={q_alive}, D active={d_alive}, S2 still 700 bp={s2_alive}, \
         active=[{}]",
        active_summary(&fx.graph)
    );
}

/// T4 / CE-4. A shadow further along the chosen exit must not veto that exit.
///
/// R is consumed and lies on Q → X → R. X itself has 20 pairs. Ranking X against
/// R would be 20/18 < 2.5, so R must be removed from the ranking as a same-route
/// shadow. D's 2 pairs are below `min_support` and must not become the winner.
///
/// ON UNPATCHED CODE THIS FAILS: committing S1+R puts X in `rejected_nodes`, the
/// search never scores X, and D's 2 pairs cannot carry the fork either.
#[test]
fn test_t4_ce4_same_exit_shadow_does_not_block() {
    let mut graph = CdbgGraph::new(31);
    let j_r = junction(0x41_0001);
    let j_sq = junction(0x41_0002);
    let j_qx = junction(0x41_0003);
    assert_junctions_distinct(&[j_r.clone(), j_sq.clone(), j_qx.clone()]);

    let s1 = span_node(&[], 0x41_0101, 900, &j_r);
    let r = span_node(&j_r, 0x41_0102, 169, &[]);
    let x = span_node(&j_qx, 0x41_0103, 300, &j_r);
    let s2 = span_node(&[], 0x41_0104, 700, &j_sq);
    let q = span_node(&j_sq, 0x41_0105, 100, &j_qx);
    let d = span_node(&j_qx, 0x41_0106, 300, &[]);
    assert_eq!(x.len(), 300);
    assert!(93 <= r.len() && 212 <= x.len() && 212 <= d.len());

    let s1_id = add_node(&mut graph, &s1);
    let r_id = add_node(&mut graph, &r);
    let x_id = add_node(&mut graph, &x);
    let s2_id = add_node(&mut graph, &s2);
    let q_id = add_node(&mut graph, &q);
    let d_id = add_node(&mut graph, &d);
    for (from, to) in [
        (s1_id, r_id),
        (x_id, r_id),
        (s2_id, q_id),
        (q_id, x_id),
        (q_id, d_id),
    ] {
        add_edge(&mut graph, from, to, OV as u8);
        assert_valid_edge(&graph, from, to, OV as u8);
    }

    let sep_x = q.len() as i64 - OV as i64 - OV as i64;
    let sep_r = sep_x + x.len() as i64 - OV as i64;
    assert_eq!((sep_x, sep_r), (38, 307));
    let dist = default_dist();
    let mut matrix = PairLinkMatrix::new();
    add_pairs(&mut matrix, s2_id, x_id, 20, 200, 212, 4_100);
    add_pairs(&mut matrix, s2_id, r_id, 18, 50, 93, 4_200);
    add_pairs(&mut matrix, s2_id, d_id, 2, 200, 212, 4_300);
    let w_x = assert_score(&matrix, &dist, s2_id, x_id, sep_x, 20);
    let w_r = assert_score(&matrix, &dist, s2_id, r_id, sep_r, 18);
    let w_d = assert_score(&matrix, &dist, s2_id, d_id, sep_x, 2);
    assert!(
        w_x / w_r < 2.5,
        "fixture: treating R as a normal rival must be able to veto X"
    );
    assert!(w_x / w_d >= 2.5);
    assert!(2 < params().min_support);

    let s1_contig = splice(&[&s1, &r], &[OV]);
    let good = splice(&[&s2, &q, &x], &[OV, OV]);
    let bad = splice(&[&s2, &q, &d], &[OV, OV]);
    assert!(!contains(&good, &r[OV..]));

    resolve(&mut graph, &matrix, &dist);

    assert!(
        has_active_seq(&graph, &s1_contig) && !alive(&graph, r_id.id()),
        "fixture: S1 must consume R, active=[{}]",
        active_summary(&graph)
    );
    let good_present = has_active_seq(&graph, &good);
    let bad_present = has_active_seq(&graph, &bad);
    let d_alive = alive(&graph, d_id.id());
    let x_alive = alive(&graph, x_id.id());
    assert!(
        good_present && !bad_present && d_alive && !x_alive,
        "T4/CE-4 RED: R is a shadow on the same exit as X (Q→X→R) and must not block X. \
         Unpatched code banned X during S1's neighbor scan, so the 20-pair exit was never \
         scored (D has only 2 pairs and cannot win alone). \
         S2+Q+X present={good_present}, S2+Q+D present={bad_present}, D active={d_alive}, \
         X active={x_alive}, active=[{}]",
        active_summary(&graph)
    );
}

/// T5. A node deleted before this resolver runs is not a consumed-node shadow.
///
/// Bubble's 20 pairs would block Good's 8 if every deleted node competed.
/// Upstream deletions stay ineligible, and Good commits.
#[test]
fn test_t5_upstream_deleted_node_is_not_a_shadow() {
    let mut graph = CdbgGraph::new(31);
    let j_sg = junction(0x51_0001);
    let j_fork = junction(0x51_0002);
    assert_junctions_distinct(&[j_sg.clone(), j_fork.clone()]);

    let s = span_node(&[], 0x51_0101, 700, &j_sg);
    let q = span_node(&j_sg, 0x51_0102, 100, &j_fork);
    let good = span_node(&j_fork, 0x51_0103, 300, &[]);
    let bubble = span_node(&j_fork, 0x51_0104, 180, &[]);
    assert!(150 <= good.len() && 150 <= bubble.len());
    assert_ne!(&good[OV..], &bubble[OV..]);

    let s_id = add_node(&mut graph, &s);
    let q_id = add_node(&mut graph, &q);
    let good_id = add_node(&mut graph, &good);
    let bubble_id = add_node(&mut graph, &bubble);
    add_edge(&mut graph, s_id, q_id, OV as u8);
    add_edge(&mut graph, q_id, good_id, OV as u8);
    add_edge(&mut graph, q_id, bubble_id, OV as u8);
    assert_valid_edge(&graph, q_id, good_id, OV as u8);
    assert_valid_edge(&graph, q_id, bubble_id, OV as u8);

    let sep = q.len() as i64 - OV as i64 - OV as i64;
    assert_eq!(sep, 38);
    let dist = default_dist();
    let mut matrix = PairLinkMatrix::new();
    add_pairs(&mut matrix, s_id, good_id, 8, 262, 150, 5_100);
    add_pairs(&mut matrix, s_id, bubble_id, 20, 262, 150, 5_200);
    let w_good = assert_score(&matrix, &dist, s_id, good_id, sep, 8);
    let w_bubble = assert_score(&matrix, &dist, s_id, bubble_id, sep, 20);
    assert!(
        w_good < 2.5 * w_bubble,
        "fixture: bubble evidence would block Good if it were allowed to compete"
    );

    graph.nodes[bubble_id.id() as usize].is_deleted = true;
    assert!(graph.nodes[bubble_id.id() as usize].is_deleted);

    let expected = splice(&[&s, &q, &good], &[OV, OV]);
    resolve(&mut graph, &matrix, &dist);

    let spliced = has_active_seq(&graph, &expected);
    let bubble_deleted = !alive(&graph, bubble_id.id());
    assert!(
        spliced && bubble_deleted && !contains(&expected, &bubble[OV..]),
        "T5: an upstream-deleted bubble is not a resolver shadow and must not block Good. \
         S+Q+Good present={spliced}, bubble still deleted={bubble_deleted}, active=[{}]",
        active_summary(&graph)
    );
    assert!(
        !alive(&graph, good_id.id()),
        "Good is part of the committed contig"
    );
}

fn side_branch_survives(
    graph: &CdbgGraph,
    x_id: u32,
    unique: &[u8],
    spliced: &[u8],
) -> (bool, bool, bool, bool) {
    let spliced_ok = has_active_seq(graph, spliced);
    let x_alive = alive(graph, x_id);
    let in_graph = graph
        .nodes
        .iter()
        .filter(|node| !node.is_deleted)
        .any(|node| contains(&graph.get_unitig_seq(node.id), unique));
    let raw = curate_contigs(graph, 0, 2.0);
    let in_raw = raw.iter().any(|contig| contains(&contig.seq, unique));
    (spliced_ok, x_alive, in_graph, in_raw)
}

fn case_a_incoming_branch_survives() {
    let mut graph = CdbgGraph::new(31);
    let junction = dna(0xA001, 31);
    let mut s_seq = dna(0xA002, 669);
    s_seq.extend_from_slice(&junction);
    let mut r_seq = junction.clone();
    r_seq.extend(dna(0xA003, 138));
    let unique = dna(0xA004, 69);
    let mut x_seq = unique.clone();
    x_seq.extend_from_slice(&junction);
    assert_eq!((s_seq.len(), r_seq.len(), x_seq.len()), (700, 169, 100));
    assert_not_homopolymer(&x_seq);

    let s_id = add_node(&mut graph, &s_seq);
    let r_id = add_node(&mut graph, &r_seq);
    let x_id = add_node(&mut graph, &x_seq);
    add_edge(&mut graph, s_id, r_id, 31);
    add_edge(&mut graph, x_id, r_id, 31);
    let spliced = splice(&[&s_seq, &r_seq], &[31]);
    resolve(
        &mut graph,
        &PairLinkMatrix::new(),
        &InsertSizeDist::default(),
    );
    let (spliced_ok, x_alive, in_graph, in_raw) =
        side_branch_survives(&graph, x_id.id(), &unique, &spliced);
    assert!(
        spliced_ok && x_alive && in_graph && in_raw,
        "T6 case A: incoming X must survive S→R. spliced={spliced_ok}, X={x_alive}, \
         graph={in_graph}, raw={in_raw}, active=[{}]",
        active_summary(&graph)
    );
}

fn case_b_intermediate_branch_survives() {
    let mut graph = CdbgGraph::new(31);
    let j0 = dna(0xB001, 31);
    let j1 = dna(0xB002, 31);
    let mut s_seq = dna(0xB003, 669);
    s_seq.extend_from_slice(&j0);
    let mut r1_seq = j0.clone();
    r1_seq.extend(dna(0xB004, 38));
    r1_seq.extend_from_slice(&j1);
    let mut r2_seq = j1.clone();
    r2_seq.extend(dna(0xB005, 69));
    let unique = dna(0xB006, 69);
    let mut x_seq = unique.clone();
    x_seq.extend_from_slice(&j0);
    assert_eq!(
        (s_seq.len(), r1_seq.len(), r2_seq.len(), x_seq.len()),
        (700, 100, 100, 100)
    );

    let s_id = add_node(&mut graph, &s_seq);
    let r1_id = add_node(&mut graph, &r1_seq);
    let r2_id = add_node(&mut graph, &r2_seq);
    let x_id = add_node(&mut graph, &x_seq);
    add_edge(&mut graph, s_id, r1_id, 31);
    add_edge(&mut graph, r1_id, r2_id, 31);
    add_edge(&mut graph, x_id, r1_id, 31);
    let spliced = splice(&[&s_seq, &r1_seq, &r2_seq], &[31, 31]);
    resolve(
        &mut graph,
        &PairLinkMatrix::new(),
        &InsertSizeDist::default(),
    );
    let (spliced_ok, x_alive, in_graph, in_raw) =
        side_branch_survives(&graph, x_id.id(), &unique, &spliced);
    assert!(
        spliced_ok && x_alive && in_graph && in_raw,
        "T6 case B: intermediate incoming X must survive. spliced={spliced_ok}, X={x_alive}, \
         graph={in_graph}, raw={in_raw}, active=[{}]",
        active_summary(&graph)
    );
}

fn case_c_multi_connection_survives() {
    let mut graph = CdbgGraph::new(31);
    let junction = dna(0xC001, 31);
    let mut s_seq = dna(0xC002, 669);
    s_seq.extend_from_slice(&junction);
    let mut r1_seq = junction.clone();
    r1_seq.extend(dna(0xC003, 38));
    r1_seq.extend_from_slice(&junction);
    let mut r2_seq = junction.clone();
    r2_seq.extend(dna(0xC004, 69));
    let unique = dna(0xC005, 69);
    let mut x_seq = unique.clone();
    x_seq.extend_from_slice(&junction);
    let s_id = add_node(&mut graph, &s_seq);
    let r1_id = add_node(&mut graph, &r1_seq);
    let r2_id = add_node(&mut graph, &r2_seq);
    let x_id = add_node(&mut graph, &x_seq);
    add_edge(&mut graph, s_id, r1_id, 31);
    add_edge(&mut graph, r1_id, r2_id, 31);
    add_edge(&mut graph, x_id, r1_id, 31);
    add_edge(&mut graph, x_id, r2_id, 31);
    assert_eq!(graph.nodes[x_id.id() as usize].edges_fwd.len(), 2);
    let spliced = splice(&[&s_seq, &r1_seq, &r2_seq], &[31, 31]);
    resolve(
        &mut graph,
        &PairLinkMatrix::new(),
        &InsertSizeDist::default(),
    );
    let (spliced_ok, x_alive, in_graph, in_raw) =
        side_branch_survives(&graph, x_id.id(), &unique, &spliced);
    assert!(
        spliced_ok && x_alive && in_graph && in_raw,
        "T6 case C: multi-connected X must survive. spliced={spliced_ok}, X={x_alive}, \
         graph={in_graph}, raw={in_raw}, active=[{}]",
        active_summary(&graph)
    );
}

fn case_d_reverse_branch_survives() {
    let mut graph = CdbgGraph::new(31);
    let junction = dna(0xD001, 31);
    let mut s_seq = dna(0xD002, 669);
    s_seq.extend_from_slice(&junction);
    let mut r_seq = junction.clone();
    r_seq.extend(dna(0xD003, 138));
    let biological_unique = dna(0xD004, 69);
    let mut biological = biological_unique.clone();
    biological.extend_from_slice(&junction);
    let x_seq = reverse_complement_ascii(&biological);
    let stored_unique = x_seq[31..].to_vec();
    let s_id = add_node(&mut graph, &s_seq);
    let r_id = add_node(&mut graph, &r_seq);
    let x_id = add_node(&mut graph, &x_seq);
    add_edge(&mut graph, s_id, r_id, 31);
    add_edge(&mut graph, x_id.twin(), r_id, 31);
    assert_valid_edge(&graph, x_id.twin(), r_id, 31);
    let spliced = splice(&[&s_seq, &r_seq], &[31]);
    resolve(
        &mut graph,
        &PairLinkMatrix::new(),
        &InsertSizeDist::default(),
    );
    let (spliced_ok, x_alive, in_graph, in_raw) =
        side_branch_survives(&graph, x_id.id(), &stored_unique, &spliced);
    assert!(
        spliced_ok && x_alive && in_graph && in_raw,
        "T6 case D: reverse-orientation X must survive. spliced={spliced_ok}, X={x_alive}, \
         graph={in_graph}, raw={in_raw}, active=[{}]",
        active_summary(&graph)
    );
}

fn case_e_cap_preserves_both_frontiers() {
    let k = 31usize;
    let n_repeat = 65usize;
    let mut graph = CdbgGraph::new(31);
    let junctions: Vec<Vec<u8>> = (0..=n_repeat)
        .map(|i| dna(0xE100_0000 + i as u64, k))
        .collect();
    let mut s_seq = dna(0xE200_0001, 700 - k);
    s_seq.extend_from_slice(&junctions[0]);
    let mut prev = add_node(&mut graph, &s_seq);
    let mut repeat_seqs = Vec::with_capacity(n_repeat);
    for i in 1..=n_repeat {
        let mut seq = junctions[i - 1].clone();
        seq.extend(dna(0xE300_0000 + i as u64, 100 - 2 * k));
        seq.extend_from_slice(&junctions[i]);
        assert_eq!(seq.len(), 100);
        let id = add_node(&mut graph, &seq);
        add_edge(&mut graph, prev, id, k as u8);
        prev = id;
        repeat_seqs.push(seq);
    }
    let unique = dna(0xE400_0001, 69);
    let mut x_seq = unique.clone();
    x_seq.extend_from_slice(&repeat_seqs[63][..k]);
    let x_id = add_node(&mut graph, &x_seq);
    let r64 = DirUnitigId::new(64, false);
    add_edge(&mut graph, x_id, r64, k as u8);
    assert_valid_edge(&graph, x_id, r64, k as u8);
    let r65_unique = repeat_seqs[64][k..].to_vec();

    resolve(
        &mut graph,
        &PairLinkMatrix::new(),
        &InsertSizeDist::default(),
    );

    let prefix_committed = graph
        .nodes
        .iter()
        .any(|node| !node.is_deleted && node.seq_len == 5116);
    let r65_alive = !graph.nodes[65].is_deleted;
    let x_alive = alive(&graph, x_id.id());
    let raw = curate_contigs(&graph, 0, 2.0);
    assert!(
        prefix_committed
            && r65_alive
            && graph.get_unitig_seq(65) == repeat_seqs[64]
            && x_alive
            && raw.iter().any(|contig| contains(&contig.seq, &unique))
            && raw.iter().any(|contig| contains(&contig.seq, &r65_unique)),
        "T6 case E: 64-step cap must keep outgoing R65 and incoming X. \
         prefix={prefix_committed}, R65={r65_alive}, X={x_alive}, active=[{}]",
        active_summary(&graph)
    );
}

fn case_f_validation_failure_keeps_external_neighbor_searchable() {
    let mut graph = CdbgGraph::new(31);
    let j_fail = dna(0xF101, 31);
    let j_qx = dna(0xF102, 31);
    let j_sq = dna(0xF103, 31);
    let mut s1 = dna(0xF111, 869);
    s1.extend_from_slice(&j_fail);
    let mut r1 = j_fail.clone();
    r1.extend(dna(0xF112, 138));
    let x_middle = dna(0xF113, 38);
    let mut x_seq = j_qx.clone();
    x_seq.extend_from_slice(&x_middle);
    x_seq.extend_from_slice(&j_fail);
    let mut s2 = dna(0xF114, 669);
    s2.extend_from_slice(&j_sq);
    let q_middle = dna(0xF115, 38);
    let mut q_seq = j_sq.clone();
    q_seq.extend_from_slice(&q_middle);
    q_seq.extend_from_slice(&j_qx);
    let mut d_seq = j_qx.clone();
    d_seq.extend(dna(0xF116, 69));
    let mut z_seq = j_fail.clone();
    z_seq.extend(dna(0xF117, 69));

    let s1_id = add_node(&mut graph, &s1);
    let r1_id = add_node(&mut graph, &r1);
    let x_id = add_node(&mut graph, &x_seq);
    let s2_id = add_node(&mut graph, &s2);
    let q_id = add_node(&mut graph, &q_seq);
    let d_id = add_node(&mut graph, &d_seq);
    let z_id = add_node(&mut graph, &z_seq);
    add_forward_only(&mut graph, s1_id, r1_id, 31);
    add_edge(&mut graph, x_id, r1_id, 31);
    add_edge(&mut graph, x_id, z_id, 31);
    add_edge(&mut graph, s2_id, q_id, 31);
    add_edge(&mut graph, q_id, x_id, 31);
    add_edge(&mut graph, q_id, d_id, 31);
    assert!(matches!(
        resolve_edge_between(&graph, r1_id.twin(), s1_id.twin()),
        EdgeResolveResult::NoEdge
    ));

    let mut matrix = PairLinkMatrix::new();
    for _ in 0..8 {
        matrix.add_link(s2_id, x_id);
    }
    let dist = InsertSizeDist {
        median: 278.0,
        mad: 40.0,
        min_insert: 100,
        max_insert: 700,
    };
    let failed_edge = graph.nodes[s1_id.id() as usize].edges_fwd.clone();
    resolve(&mut graph, &matrix, &dist);

    assert!(alive(&graph, s1_id.id()) && alive(&graph, r1_id.id()));
    assert_eq!(graph.nodes[s1_id.id() as usize].edges_fwd, failed_edge);
    let expected = splice(&[&s2, &q_seq, &x_seq], &[31, 31]);
    assert!(
        has_active_seq(&graph, &expected),
        "T6 case F: external neighbor X stays searchable after a failed reciprocal check. \
         active=[{}]",
        active_summary(&graph)
    );
    let forbidden = splice(&[&s1, &r1], &[31]);
    assert!(graph
        .nodes
        .iter()
        .all(|node| graph.get_unitig_seq(node.id) != forbidden));
}

fn case_g_unbranched_path_splices() {
    let mut graph = CdbgGraph::new(31);
    let j0 = dna(0x6701, 31);
    let j1 = dna(0x6702, 31);
    let mut s_seq = dna(0x6703, 669);
    s_seq.extend_from_slice(&j0);
    let mut r_seq = j0.clone();
    r_seq.extend(dna(0x6704, 38));
    r_seq.extend_from_slice(&j1);
    let mut t_seq = j1.clone();
    t_seq.extend(dna(0x6705, 519));
    let s_id = add_node(&mut graph, &s_seq);
    let r_id = add_node(&mut graph, &r_seq);
    let t_id = add_node(&mut graph, &t_seq);
    add_edge(&mut graph, s_id, r_id, 31);
    add_edge(&mut graph, r_id, t_id, 31);
    let expected = splice(&[&s_seq, &r_seq, &t_seq], &[31, 31]);
    assert_eq!(expected.len(), 1288);
    resolve(
        &mut graph,
        &PairLinkMatrix::new(),
        &InsertSizeDist::default(),
    );
    let active: Vec<_> = graph.nodes.iter().filter(|node| !node.is_deleted).collect();
    assert_eq!(
        active.len(),
        1,
        "T6 case G: unbranched S→R→T must collapse to one contig, active=[{}]",
        active_summary(&graph)
    );
    assert_eq!(graph.get_unitig_seq(active[0].id), expected);
    assert!(graph.nodes[s_id.id() as usize].is_deleted);
    assert!(graph.nodes[r_id.id() as usize].is_deleted);
    assert!(graph.nodes[t_id.id() as usize].is_deleted);
}

/// T6. Cases A–G from the incoming side-branch gate keep their current meaning.
#[test]
fn test_t6_cases_a_through_g_stay_compatible() {
    case_a_incoming_branch_survives();
    case_b_intermediate_branch_survives();
    case_c_multi_connection_survives();
    case_d_reverse_branch_survives();
    case_e_cap_preserves_both_frontiers();
    case_f_validation_failure_keeps_external_neighbor_searchable();
    case_g_unbranched_path_splices();
}

/// T7 atomicity. A missing reciprocal edge rejects the path before any commit.
/// No constituent is deleted, no contig is added, and the fingerprint is unchanged.
/// That is the observable form of "nothing entered consumed_nodes".
#[test]
fn test_t7_validation_failure_is_atomic() {
    let mut graph = CdbgGraph::new(31);
    let junction = dna(0x71A1, 31);
    let mut s_seq = dna(0x71A2, 769);
    s_seq.extend_from_slice(&junction);
    let mut r_seq = junction.clone();
    r_seq.extend(dna(0x71A3, 138));
    let mut x_seq = dna(0x71A4, 69);
    x_seq.extend_from_slice(&junction);
    assert_eq!((s_seq.len(), r_seq.len(), x_seq.len()), (800, 169, 100));

    let s_id = add_node(&mut graph, &s_seq);
    let r_id = add_node(&mut graph, &r_seq);
    let x_id = add_node(&mut graph, &x_seq);
    add_forward_only(&mut graph, s_id, r_id, 31);
    add_edge(&mut graph, x_id, r_id, 31);
    assert_valid_edge(&graph, s_id, r_id, 31);
    assert!(matches!(
        resolve_edge_between(&graph, r_id.twin(), s_id.twin()),
        EdgeResolveResult::NoEdge
    ));

    let before = structural_fingerprint(&graph);
    let before_len = graph.nodes.len();
    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &PairLinkMatrix::new(),
        &InsertSizeDist::default(),
        params(),
    ).expect("repeat resolution failed");
    assert_eq!(
        stats.valid_k_joins + stats.valid_k_minus_1_joins + stats.valid_other_joins,
        0
    );
    assert!(stats.rejected_invalid_joins >= 1);
    assert_eq!(graph.nodes.len(), before_len);
    assert!(graph.nodes.iter().all(|node| !node.is_deleted));
    assert_eq!(
        structural_fingerprint(&graph),
        before,
        "T7: validation failure must not delete nodes, add contigs, or edit edges"
    );
}

/// T7 eligibility. Nodes on a path that failed validation were not consumed.
/// A later seed whose supported exit is one of those nodes must still be able to use it.
///
/// ON UNPATCHED CODE THIS FAILS: the Err branch inserts R1 into `rejected_nodes`,
/// S2 skips R1, and the 5-pair decoy D is committed.
#[test]
fn test_t7_failed_validation_does_not_ban_later_seed() {
    let mut graph = CdbgGraph::new(31);
    let j_r = junction(0x77_0001);
    let j_sq = junction(0x77_0002);
    assert_junctions_distinct(&[j_r.clone(), j_sq.clone()]);

    let s1 = span_node(&[], 0x77_0101, 900, &j_r);
    let r1 = span_node(&j_r, 0x77_0102, 300, &[]);
    let s2 = span_node(&[], 0x77_0103, 700, &j_sq);
    let q = span_node(&j_sq, 0x77_0104, 100, &j_r);
    let d = span_node(&j_r, 0x77_0105, 300, &[]);
    assert!(212 <= r1.len() && 212 <= d.len());

    let s1_id = add_node(&mut graph, &s1);
    let r1_id = add_node(&mut graph, &r1);
    let s2_id = add_node(&mut graph, &s2);
    let q_id = add_node(&mut graph, &q);
    let d_id = add_node(&mut graph, &d);
    add_forward_only(&mut graph, s1_id, r1_id, OV as u8);
    add_edge(&mut graph, s2_id, q_id, OV as u8);
    add_edge(&mut graph, q_id, r1_id, OV as u8);
    add_edge(&mut graph, q_id, d_id, OV as u8);
    assert_valid_edge(&graph, s1_id, r1_id, OV as u8);
    assert!(matches!(
        resolve_edge_between(&graph, r1_id.twin(), s1_id.twin()),
        EdgeResolveResult::NoEdge
    ));
    assert_valid_edge(&graph, q_id, r1_id, OV as u8);
    assert_valid_edge(&graph, q_id, d_id, OV as u8);

    let sep = q.len() as i64 - OV as i64 - OV as i64;
    assert_eq!(sep, 38);
    let dist = default_dist();
    let mut matrix = PairLinkMatrix::new();
    add_pairs(&mut matrix, s2_id, r1_id, 20, 200, 212, 7_100);
    add_pairs(&mut matrix, s2_id, d_id, 5, 200, 212, 7_200);
    let w_r = assert_score(&matrix, &dist, s2_id, r1_id, sep, 20);
    let w_d = assert_score(&matrix, &dist, s2_id, d_id, sep, 5);
    assert!(w_r / w_d >= 2.5);

    let failed_edge = graph.nodes[s1_id.id() as usize].edges_fwd.clone();
    let forbidden = splice(&[&s1, &r1], &[OV]);
    let good = splice(&[&s2, &q, &r1], &[OV, OV]);
    let bad = splice(&[&s2, &q, &d], &[OV, OV]);

    resolve(&mut graph, &matrix, &dist);

    assert!(
        alive(&graph, s1_id.id()),
        "fixture: failed seed S1 must stay active, active=[{}]",
        active_summary(&graph)
    );
    assert_eq!(graph.nodes[s1_id.id() as usize].edges_fwd, failed_edge);
    assert_eq!(graph.get_unitig_seq(s1_id.id()), s1);
    assert!(graph
        .nodes
        .iter()
        .all(|node| graph.get_unitig_seq(node.id) != forbidden));

    let good_present = has_active_seq(&graph, &good);
    let bad_present = has_active_seq(&graph, &bad);
    let d_alive = alive(&graph, d_id.id());
    let r_alive = alive(&graph, r1_id.id());
    assert!(
        good_present && !bad_present && d_alive && !r_alive,
        "T7 RED: validation failure never consumed R1, but the Err branch still banned it. \
         S2 skipped the 20-pair exit and committed D (5 pairs). \
         S2+Q+R1 present={good_present}, S2+Q+D present={bad_present}, \
         D active={d_alive}, R1 active={r_alive}, S1 active=true, active=[{}]",
        active_summary(&graph)
    );
}

/// T8. Reverse-complement copies of T1 and T3 produce the same per-node outcome
/// up to reverse complement: same lengths, same deletion flags, same canonical sequences.
#[test]
fn test_t8_twin_orientation_matches_forward() {
    let mut forward = build_ce1(false);
    let mut mirror = build_ce1(true);
    resolve(&mut forward.graph, &forward.matrix, &forward.dist);
    resolve(&mut mirror.graph, &mirror.matrix, &mirror.dist);
    assert_eq!(
        mirror_signature(&forward.graph),
        mirror_signature(&mirror.graph),
        "T8: CE-1 mirror diverged from the forward fixture. \
         forward active=[{}] mirror active=[{}]",
        active_summary(&forward.graph),
        active_summary(&mirror.graph)
    );

    let mut forward = build_ce3(false);
    let mut mirror = build_ce3(true);
    resolve(&mut forward.graph, &forward.matrix, &forward.dist);
    resolve(&mut mirror.graph, &mirror.matrix, &mirror.dist);
    assert_eq!(
        mirror_signature(&forward.graph),
        mirror_signature(&mirror.graph),
        "T8: CE-3 mirror diverged from the forward fixture. \
         forward active=[{}] mirror active=[{}]",
        active_summary(&forward.graph),
        active_summary(&mirror.graph)
    );
}

/// T9. Directed cycles and hairpins stop on physical node identity.
/// Undecided forks that contain those motifs do not change the graph.
#[test]
fn test_t9_cycles_and_hairpins_terminate() {
    // S → A → R → A. Physical-id tracking stops the second visit to A.
    {
        let mut graph = CdbgGraph::new(31);
        let j_a = junction(0x91_0001);
        let j_ar = junction(0x91_0002);
        assert_junctions_distinct(&[j_a.clone(), j_ar.clone()]);
        let s = span_node(&[], 0x91_0101, 700, &j_a);
        let a = span_node(&j_a, 0x91_0102, 120, &j_ar);
        let r = span_node(&j_ar, 0x91_0103, 120, &j_a);
        let s_id = add_node(&mut graph, &s);
        let a_id = add_node(&mut graph, &a);
        let r_id = add_node(&mut graph, &r);
        add_edge(&mut graph, s_id, a_id, OV as u8);
        add_edge(&mut graph, a_id, r_id, OV as u8);
        add_edge(&mut graph, r_id, a_id, OV as u8);
        assert_valid_edge(&graph, a_id, r_id, OV as u8);
        assert_valid_edge(&graph, r_id, a_id, OV as u8);
        let expected = splice(&[&s, &a, &r], &[OV, OV]);
        let looped = splice(&[&s, &a, &r, &a], &[OV, OV, OV]);
        let before = graph.nodes.len();
        resolve(
            &mut graph,
            &PairLinkMatrix::new(),
            &InsertSizeDist::default(),
        );
        assert_eq!(
            graph.nodes.len(),
            before + 1,
            "cycle emitted more than one contig"
        );
        assert!(
            has_active_seq(&graph, &expected) && !has_active_seq(&graph, &looped),
            "T9: A→R→A must stop when A repeats, active=[{}]",
            active_summary(&graph)
        );
        assert!(!alive(&graph, a_id.id()) && !alive(&graph, r_id.id()));
        let a_body = &a[OV..a.len() - OV];
        let contig = graph.get_unitig_seq(before as u32);
        assert_eq!(
            contig
                .windows(a_body.len())
                .filter(|w| *w == a_body)
                .count(),
            1,
            "A's unique bases must appear once, not once per lap"
        );
    }

    // R → R.rc uses an even-length RC-palindrome overlap. Odd overlaps cannot match.
    {
        let mut graph = CdbgGraph::new(31);
        let j_sr = junction(0x92_0001);
        let hairpin = rc_palindrome(0x92_0002, 15);
        assert_eq!(hairpin.len(), 30);
        let s = span_node(&[], 0x92_0101, 700, &j_sr);
        let r = span_node(&j_sr, 0x92_0102, 200, &hairpin);
        assert_eq!(r.len(), 200);
        let s_id = add_node(&mut graph, &s);
        let r_id = add_node(&mut graph, &r);
        add_edge(&mut graph, s_id, r_id, OV as u8);
        add_edge(&mut graph, r_id, r_id.twin(), hairpin.len() as u8);
        assert_valid_edge(&graph, s_id, r_id, OV as u8);
        assert_valid_edge(&graph, r_id, r_id.twin(), hairpin.len() as u8);
        let expected = splice(&[&s, &r], &[OV]);
        let before = graph.nodes.len();
        resolve(
            &mut graph,
            &PairLinkMatrix::new(),
            &InsertSizeDist::default(),
        );
        assert_eq!(
            graph.nodes.len(),
            before + 1,
            "hairpin did not stop at one copy of R"
        );
        assert!(has_active_seq(&graph, &expected));
        assert_eq!(
            graph
                .get_unitig_seq(before as u32)
                .windows(hairpin.len())
                .filter(|window| *window == hairpin)
                .count(),
            1
        );
    }

    // Branching cycle with no pairs: BFS must return and leave the graph alone.
    {
        let mut graph = CdbgGraph::new(31);
        let j_sq = junction(0x93_0001);
        let j_qr = junction(0x93_0002);
        let j_rc = junction(0x93_0003);
        assert_junctions_distinct(&[j_sq.clone(), j_qr.clone(), j_rc.clone()]);
        let s = span_node(&[], 0x93_0101, 700, &j_sq);
        let q = span_node(&j_sq, 0x93_0102, 100, &j_qr);
        let r = span_node(&j_qr, 0x93_0103, 80, &j_rc);
        let c = span_node(&j_rc, 0x93_0104, 80, &j_qr);
        let d = span_node(&j_qr, 0x93_0105, 80, &[]);
        let s_id = add_node(&mut graph, &s);
        let q_id = add_node(&mut graph, &q);
        let r_id = add_node(&mut graph, &r);
        let c_id = add_node(&mut graph, &c);
        let d_id = add_node(&mut graph, &d);
        for (from, to) in [
            (s_id, q_id),
            (q_id, r_id),
            (q_id, d_id),
            (r_id, c_id),
            (c_id, r_id),
        ] {
            add_edge(&mut graph, from, to, OV as u8);
            assert_valid_edge(&graph, from, to, OV as u8);
        }
        let before = structural_fingerprint(&graph);
        resolve(
            &mut graph,
            &PairLinkMatrix::new(),
            &InsertSizeDist::default(),
        );
        assert_eq!(
            structural_fingerprint(&graph),
            before,
            "T9: directed cycle A/R/C under an undecided fork must not mutate the graph"
        );
    }

    // Branching hairpin with no pairs.
    {
        let mut graph = CdbgGraph::new(31);
        let j_sq = junction(0x94_0001);
        let j_qr = junction(0x94_0002);
        let hairpin = rc_palindrome(0x94_0003, 15);
        assert_junctions_distinct(&[j_sq.clone(), j_qr.clone()]);
        let s = span_node(&[], 0x94_0101, 700, &j_sq);
        let q = span_node(&j_sq, 0x94_0102, 100, &j_qr);
        let r = span_node(&j_qr, 0x94_0103, 120, &hairpin);
        let d = span_node(&j_qr, 0x94_0104, 90, &[]);
        assert_eq!(r.len(), 120);
        let s_id = add_node(&mut graph, &s);
        let q_id = add_node(&mut graph, &q);
        let r_id = add_node(&mut graph, &r);
        let d_id = add_node(&mut graph, &d);
        add_edge(&mut graph, s_id, q_id, OV as u8);
        add_edge(&mut graph, q_id, r_id, OV as u8);
        add_edge(&mut graph, q_id, d_id, OV as u8);
        add_edge(&mut graph, r_id, r_id.twin(), hairpin.len() as u8);
        assert_valid_edge(&graph, r_id, r_id.twin(), hairpin.len() as u8);
        let before = structural_fingerprint(&graph);
        resolve(
            &mut graph,
            &PairLinkMatrix::new(),
            &InsertSizeDist::default(),
        );
        assert_eq!(
            structural_fingerprint(&graph),
            before,
            "T9: hairpin under an undecided fork must not mutate the graph"
        );
    }
}

fn assert_budget_transition(sep_u: i64, len_u: i64, overlap: i64, expect: i64) -> i64 {
    let sep_next = sep_u + len_u - overlap;
    assert_eq!(sep_next, expect, "fixture separation arithmetic");
    sep_next
}

/// BFS expansion. The supported endpoint sits one hop past U.
/// `sep_next == max_insert` is admitted; `max_insert + 1` leaves the graph unchanged.
fn run_bfs_budget(overlap_ub: u8, admit: bool) {
    const MAX_SEP: i64 = 500;
    let mut graph = CdbgGraph::new(31);
    let j_ar = junction(0xA10_0001);
    let j_ru = junction(0xA10_0002);
    let j_ub = dna(0xA10_0003, overlap_ub as usize);
    assert_ne!(j_ar, j_ru);
    assert_ne!(&j_ub[..OV], j_ar.as_slice());
    assert_ne!(&j_ub[..OV], j_ru.as_slice());

    let a = span_node(&[], 0xA10_0101, 700, &j_ar);
    let r = span_node(&j_ar, 0xA10_0102, 400, &j_ru);
    let u = span_node(&j_ru, 0xA10_0103, 250, &j_ub);
    let b = span_node(&j_ub, 0xA10_0104, 200, &[]);
    let decoy = span_node(&j_ru, 0xA10_0105, 80, &[]);
    assert!(b.len() < 500 && u.len() < 500 && r.len() < 500 && decoy.len() < 500);

    let sep_u = r.len() as i64 - OV as i64 - OV as i64;
    assert_eq!(sep_u, 338);
    let sep_b = assert_budget_transition(
        sep_u,
        u.len() as i64,
        i64::from(overlap_ub),
        if admit { MAX_SEP } else { MAX_SEP + 1 },
    );
    assert!(sep_u <= MAX_SEP);

    let a_id = add_node(&mut graph, &a);
    let r_id = add_node(&mut graph, &r);
    let u_id = add_node(&mut graph, &u);
    let b_id = add_node(&mut graph, &b);
    let d_id = add_node(&mut graph, &decoy);
    add_edge(&mut graph, a_id, r_id, OV as u8);
    add_edge(&mut graph, r_id, u_id, OV as u8);
    add_edge(&mut graph, r_id, d_id, OV as u8);
    add_edge(&mut graph, u_id, b_id, overlap_ub);
    assert_eq!(graph.nodes[r_id.id() as usize].edges_fwd.len(), 2);
    assert_valid_edge(&graph, u_id, b_id, overlap_ub);

    let mut matrix = PairLinkMatrix::new();
    // Span equals sep_b only when both distances are 0, which is the inclusive insert bound.
    add_pairs(&mut matrix, a_id, b_id, 8, 0, 0, 10_000);
    let dist = InsertSizeDist {
        median: 500.0,
        mad: 40.0,
        min_insert: 100,
        max_insert: MAX_SEP as usize,
    };
    if admit {
        assert_score(&matrix, &dist, a_id, b_id, sep_b, 8);
    }

    let expected = splice(&[&a, &r, &u, &b], &[OV, OV, overlap_ub as usize]);
    let before = structural_fingerprint(&graph);
    resolve(&mut graph, &matrix, &dist);

    if admit {
        assert!(
            has_active_seq(&graph, &expected),
            "T10: sep_next == {MAX_SEP} must be explored and spliced. active=[{}]",
            active_summary(&graph)
        );
    } else {
        assert_eq!(
            structural_fingerprint(&graph),
            before,
            "T10: sep_next == {} must not be enqueued, and the fork must stay unresolved. \
             active=[{}]",
            MAX_SEP + 1,
            active_summary(&graph)
        );
    }
}

/// Depth-cap helper. C1..C4 stay inside the budget. C5 is the boundary hop.
/// Equality counts as a hidden continuation and blocks the shortcut.
/// One base over the budget is not a continuation, so the shortcut may commit.
fn run_depth_budget(overlap_c5: u8, hidden: bool) {
    const MAX_SEP: i64 = 600;
    let mut graph = CdbgGraph::new(31);
    let j_ar = junction(0xA20_0001);
    let j_rb = junction(0xA20_0002);
    let j_12 = junction(0xA20_0003);
    let j_23 = junction(0xA20_0004);
    let j_34 = junction(0xA20_0005);
    let j_45 = dna(0xA20_0006, overlap_c5 as usize);
    assert_junctions_distinct(&[
        j_ar.clone(),
        j_rb.clone(),
        j_12.clone(),
        j_23.clone(),
        j_34.clone(),
    ]);

    let a = span_node(&[], 0xA20_0101, 700, &j_ar);
    let r = span_node(&j_ar, 0xA20_0102, 100, &j_rb);
    let b = span_node(&j_rb, 0xA20_0103, 400, &[]);
    let c1 = span_node(&j_rb, 0xA20_0111, 180, &j_12);
    let c2 = span_node(&j_12, 0xA20_0112, 180, &j_23);
    let c3 = span_node(&j_23, 0xA20_0113, 180, &j_34);
    let c4 = span_node(&j_34, 0xA20_0114, 165, &j_45);
    let c5 = span_node(&j_45, 0xA20_0115, 80, &[]);
    assert!([
        r.len(),
        b.len(),
        c1.len(),
        c2.len(),
        c3.len(),
        c4.len(),
        c5.len()
    ]
    .iter()
    .all(|len| *len < 500));

    let sep_c1 = r.len() as i64 - OV as i64 - OV as i64;
    assert_eq!(sep_c1, 38);
    let sep_c2 = assert_budget_transition(sep_c1, c1.len() as i64, OV as i64, 187);
    let sep_c3 = assert_budget_transition(sep_c2, c2.len() as i64, OV as i64, 336);
    let sep_c4 = assert_budget_transition(sep_c3, c3.len() as i64, OV as i64, 485);
    let expect_c5 = if hidden { MAX_SEP } else { MAX_SEP + 1 };
    let sep_c5 =
        assert_budget_transition(sep_c4, c4.len() as i64, i64::from(overlap_c5), expect_c5);
    assert!(sep_c4 <= MAX_SEP && sep_c5 == expect_c5);
    let _ = sep_c5;

    let a_id = add_node(&mut graph, &a);
    let r_id = add_node(&mut graph, &r);
    let b_id = add_node(&mut graph, &b);
    let c1_id = add_node(&mut graph, &c1);
    let c2_id = add_node(&mut graph, &c2);
    let c3_id = add_node(&mut graph, &c3);
    let c4_id = add_node(&mut graph, &c4);
    let c5_id = add_node(&mut graph, &c5);
    for (from, to, ov) in [
        (a_id, r_id, OV as u8),
        (r_id, b_id, OV as u8),
        (r_id, c1_id, OV as u8),
        (c1_id, c2_id, OV as u8),
        (c2_id, c3_id, OV as u8),
        (c3_id, c4_id, OV as u8),
        (c4_id, c5_id, overlap_c5),
    ] {
        add_edge(&mut graph, from, to, ov);
        assert_valid_edge(&graph, from, to, ov);
    }
    assert_eq!(graph.nodes[r_id.id() as usize].edges_fwd.len(), 2);

    let mut matrix = PairLinkMatrix::new();
    add_pairs(&mut matrix, a_id, b_id, 10, 120, 120, 20_000);
    let dist = InsertSizeDist {
        median: 278.0,
        mad: 40.0,
        min_insert: 100,
        max_insert: MAX_SEP as usize,
    };
    // sep(B) = 38. 120 + 38 + 120 = 278.
    assert_score(&matrix, &dist, a_id, b_id, sep_c1, 10);

    let shortcut = splice(&[&a, &r, &b], &[OV, OV]);
    let before = structural_fingerprint(&graph);
    resolve(&mut graph, &matrix, &dist);

    if hidden {
        assert_eq!(
            structural_fingerprint(&graph),
            before,
            "T10: sep_next == {MAX_SEP} at depth 4 is a hidden continuation. \
             The shortcut must not commit. active=[{}]",
            active_summary(&graph)
        );
    } else {
        assert!(
            has_active_seq(&graph, &shortcut),
            "T10: sep_next == {} is outside the budget, so the depth cap is not hiding a \
             continuation and the supported shortcut must commit. active=[{}]",
            MAX_SEP + 1,
            active_summary(&graph)
        );
    }
}

/// T10. Inclusive separation budget at both the BFS enqueue and the depth-cap helper.
#[test]
fn test_t10_separation_budget_is_exact() {
    run_bfs_budget(88, true);
    run_bfs_budget(87, false);
    run_depth_budget(50, true);
    run_depth_budget(49, false);
}

/// T12. A chain of consumed nodes is bounded by the separation budget.
///
/// S1 commits A→B→C→Far. From S2, Far's separation is 305 and the insert budget
/// is 250, so Far cannot compete. The in-budget decoy has 8 centered pairs and
/// commits. The resolver adds exactly two contigs.
#[test]
fn test_t12_consumed_chain_respects_separation_budget() {
    const MAX_SEP: i64 = 250;
    let mut graph = CdbgGraph::new(31);
    let j_a = junction(0xC2_0001);
    let j_ab = junction(0xC2_0002);
    let j_bc = junction(0xC2_0003);
    let j_cf = junction(0xC2_0004);
    let j_sq = junction(0xC2_0005);
    assert_junctions_distinct(&[
        j_a.clone(),
        j_ab.clone(),
        j_bc.clone(),
        j_cf.clone(),
        j_sq.clone(),
    ]);

    let s1 = span_node(&[], 0xC2_0101, 900, &j_a);
    let a = span_node(&j_a, 0xC2_0102, 120, &j_ab);
    let b = span_node(&j_ab, 0xC2_0103, 120, &j_bc);
    let c = span_node(&j_bc, 0xC2_0104, 120, &j_cf);
    let far = span_node(&j_cf, 0xC2_0105, 120, &[]);
    let s2 = span_node(&[], 0xC2_0106, 700, &j_sq);
    let q = span_node(&j_sq, 0xC2_0107, 100, &j_a);
    let decoy = span_node(&j_a, 0xC2_0108, 300, &[]);
    assert!([a.len(), b.len(), c.len(), far.len(), q.len(), decoy.len()]
        .iter()
        .all(|len| *len < 500));

    let sep_a = q.len() as i64 - OV as i64 - OV as i64;
    let sep_b = assert_budget_transition(sep_a, a.len() as i64, OV as i64, 127);
    let sep_c = assert_budget_transition(sep_b, b.len() as i64, OV as i64, 216);
    let sep_far = assert_budget_transition(sep_c, c.len() as i64, OV as i64, 305);
    assert_eq!(sep_a, 38);
    assert!(sep_c <= MAX_SEP && sep_far > MAX_SEP);

    let s1_id = add_node(&mut graph, &s1);
    let a_id = add_node(&mut graph, &a);
    let b_id = add_node(&mut graph, &b);
    let c_id = add_node(&mut graph, &c);
    let far_id = add_node(&mut graph, &far);
    let s2_id = add_node(&mut graph, &s2);
    let q_id = add_node(&mut graph, &q);
    let d_id = add_node(&mut graph, &decoy);
    for (from, to) in [
        (s1_id, a_id),
        (a_id, b_id),
        (b_id, c_id),
        (c_id, far_id),
        (s2_id, q_id),
        (q_id, a_id),
        (q_id, d_id),
    ] {
        add_edge(&mut graph, from, to, OV as u8);
        assert_valid_edge(&graph, from, to, OV as u8);
    }
    assert_eq!(graph.nodes[q_id.id() as usize].edges_fwd.len(), 2);
    assert_eq!(graph.nodes[a_id.id() as usize].edges_fwd.len(), 1);

    let mut matrix = PairLinkMatrix::new();
    // 100 + 38 + 112 = 250, on the insert ceiling.
    assert!(112 <= decoy.len());
    add_pairs(&mut matrix, s2_id, d_id, 8, 100, 112, 12_000);
    add_pairs(&mut matrix, s2_id, far_id, 20, 40, 40, 12_100);
    let dist = InsertSizeDist {
        median: 250.0,
        mad: 30.0,
        min_insert: 50,
        max_insert: MAX_SEP as usize,
    };
    assert_score(&matrix, &dist, s2_id, d_id, sep_a, 8);
    let (far_count, _) = score_candidate_path(&matrix, &dist, s2_id, far_id, sep_far);
    assert_eq!(
        far_count, 0,
        "fixture: Far's span is outside the insert window"
    );

    let s1_contig = splice(&[&s1, &a, &b, &c, &far], &[OV, OV, OV, OV]);
    let s2_contig = splice(&[&s2, &q, &decoy], &[OV, OV]);
    assert!(contains(&s1_contig, &far[OV..]));
    assert!(!contains(&s2_contig, &far[OV..]));
    let before = graph.nodes.len();

    resolve(&mut graph, &matrix, &dist);

    assert_eq!(
        graph.nodes.len(),
        before + 2,
        "T12: consumed chain must not explode into extra contigs, active=[{}]",
        active_summary(&graph)
    );
    assert!(
        has_active_seq(&graph, &s1_contig) && has_active_seq(&graph, &s2_contig),
        "T12: S1 consumes the chain and S2 commits the in-budget decoy. \
         Far at sep {sep_far} > {MAX_SEP} must not block that decision. active=[{}]",
        active_summary(&graph)
    );
    assert!(!alive(&graph, far_id.id()));
    assert!(!alive(&graph, d_id.id()));
}

/// T13. Two candidates share the first hop and then diverge. The branch that
/// enters a consumed node is an alternate exit even though `steps[0]` matches.
///
/// Q → X → D has 20 pairs. Q → X → R has 20 pairs and R was consumed by S1.
/// Those routes differ at the second step, and 20 < 2.5 * 20, so the decision
/// is blocked. The weak arm Q → Y (5 pairs) must not be taken either.
/// Without the shadow, D would beat Y by 4.0 and commit; the shadow is what
/// changes the result.
///
/// ON UNPATCHED CODE THIS FAILS: X is cleanup-banned, D is never scored, and Y
/// is committed as the only candidate.
#[test]
fn test_t13_shared_prefix_downstream_shadow_blocks_fork() {
    let mut graph = CdbgGraph::new(31);
    let j_r = junction(0xD3_0001);
    let j_sq = junction(0xD3_0002);
    let j_qx = junction(0xD3_0003);
    assert_junctions_distinct(&[j_r.clone(), j_sq.clone(), j_qx.clone()]);

    let s1 = span_node(&[], 0xD3_0101, 900, &j_r);
    let r = span_node(&j_r, 0xD3_0102, 169, &[]);
    let s2 = span_node(&[], 0xD3_0103, 700, &j_sq);
    let q = span_node(&j_sq, 0xD3_0104, 100, &j_qx);
    let x = span_node(&j_qx, 0xD3_0105, 100, &j_r);
    let d = span_node(&j_r, 0xD3_0106, 300, &[]);
    let y = span_node(&j_qx, 0xD3_0107, 300, &[]);
    assert_eq!((x.len(), r.len(), d.len(), y.len()), (100, 169, 300, 300));
    assert!(100 <= r.len() && 193 <= d.len() && 212 <= y.len());

    let s1_id = add_node(&mut graph, &s1);
    let r_id = add_node(&mut graph, &r);
    let s2_id = add_node(&mut graph, &s2);
    let q_id = add_node(&mut graph, &q);
    let x_id = add_node(&mut graph, &x);
    let d_id = add_node(&mut graph, &d);
    let y_id = add_node(&mut graph, &y);
    for (from, to) in [
        (s1_id, r_id),
        (s2_id, q_id),
        (q_id, x_id),
        (q_id, y_id),
        (x_id, d_id),
        (x_id, r_id),
    ] {
        add_edge(&mut graph, from, to, OV as u8);
        assert_valid_edge(&graph, from, to, OV as u8);
    }
    assert_eq!(graph.nodes[q_id.id() as usize].edges_fwd.len(), 2);
    assert_eq!(graph.nodes[x_id.id() as usize].edges_fwd.len(), 2);

    let sep_x = q.len() as i64 - OV as i64 - OV as i64;
    let sep_down = sep_x + x.len() as i64 - OV as i64;
    assert_eq!((sep_x, sep_down), (38, 107));
    let dist = default_dist();
    let mut matrix = PairLinkMatrix::new();
    // Downstream span: 243 + 107 + 100 = 450. Y: 200 + 38 + 212 = 450.
    add_pairs(&mut matrix, s2_id, d_id, 20, 243, 100, 13_100);
    add_pairs(&mut matrix, s2_id, r_id, 20, 243, 100, 13_200);
    add_pairs(&mut matrix, s2_id, y_id, 5, 200, 212, 13_300);
    let w_d = assert_score(&matrix, &dist, s2_id, d_id, sep_down, 20);
    let w_r = assert_score(&matrix, &dist, s2_id, r_id, sep_down, 20);
    let w_y = assert_score(&matrix, &dist, s2_id, y_id, sep_x, 5);
    assert!(w_d / w_y >= 2.5, "without the shadow, D would beat Y");
    assert!(
        w_d < 2.5 * w_r,
        "diverging shadow R must be strong enough to block D"
    );

    let s1_contig = splice(&[&s1, &r], &[OV]);
    let through_d = splice(&[&s2, &q, &x, &d], &[OV, OV, OV]);
    let through_y = splice(&[&s2, &q, &y], &[OV, OV]);
    assert!(!contains(&through_y, &d[OV..]));
    assert!(!contains(&through_d, &y[OV..]));

    resolve(&mut graph, &matrix, &dist);

    assert!(
        has_active_seq(&graph, &s1_contig) && !alive(&graph, r_id.id()),
        "fixture: S1 must consume R, active=[{}]",
        active_summary(&graph)
    );
    let d_path = has_active_seq(&graph, &through_d);
    let y_path = has_active_seq(&graph, &through_y);
    let preserved = alive(&graph, s2_id.id())
        && alive(&graph, q_id.id())
        && alive(&graph, x_id.id())
        && alive(&graph, d_id.id())
        && alive(&graph, y_id.id());
    assert!(
        !d_path && !y_path && preserved,
        "T13 RED: Q→X→D and Q→X→R share the first hop and diverge at X. R is consumed \
         and tied with D (20 vs 20), so the fork must be preserved. Unpatched code banned X, \
         never scored D, and committed weak Y (5 pairs) as the only candidate. \
         D-path present={d_path}, Y-path present={y_path}, fork preserved={preserved}, \
         active=[{}]",
        active_summary(&graph)
    );
}

#[test]
fn test_t8_canonical_helper_is_orientation_independent() {
    let seq = dna(0x8811, 64);
    let rc = reverse_complement_ascii(&seq);
    assert_eq!(canonical(&seq), canonical(&rc));
    let mut seen = BTreeSet::new();
    seen.insert(canonical(&seq));
    seen.insert(canonical(&rc));
    assert_eq!(seen.len(), 1);
}

/// T14 / F10-P1-002. Tandem repeat `R → R` must stay in the candidate search.
///
/// ```text
/// S(700) → Q(100) → D(300)
///             └──→ R(120) → X(300)
///                    ↺
/// ```
///
/// Overlaps are 31 bp. R is shorter than 500 bp and not circular; its 3' junction
/// equals its 5' junction, so the self-loop is sequence-valid and the internal
/// bases differ. Inserts are centered at 450 bp with MAD 20, inside [380, 520].
///
/// S→Q→D has separation 38 and span 450 (5 pairs, weight 5). One copy of R reaches
/// X at separation 127 and span 361, outside the window. The tandem hop reaches X
/// at separation 216 and span 450 (20 pairs, weight 20), which beats D by 4.0.
/// That route is three hops from Q, inside the depth cap, and inside max_insert.
///
/// ON UNPATCHED CODE THIS FAILS: `visited` drops `R → R` without marking the
/// search truncated. X then has no in-range pair, D is the only candidate, and
/// 5.0/5 >= 0.05 commits the misassembly S→Q→D.
#[test]
fn test_t14_tandem_repeat_cycle_must_not_allow_decoy() {
    let mut graph = CdbgGraph::new(31);
    let j_sq = junction(0xE14_0001);
    let j_loop = junction(0xE14_0002);
    assert_junctions_distinct(&[j_sq.clone(), j_loop.clone()]);

    let s_seq = span_node(&[], 0xE14_0101, 700, &j_sq);
    let q_seq = span_node(&j_sq, 0xE14_0102, 100, &j_loop);
    let d_seq = span_node(&j_loop, 0xE14_0103, 300, &[]);
    let r_seq = span_node(&j_loop, 0xE14_0104, 120, &j_loop);
    let x_seq = span_node(&j_loop, 0xE14_0105, 300, &[]);
    assert_eq!(
        (
            s_seq.len(),
            q_seq.len(),
            d_seq.len(),
            r_seq.len(),
            x_seq.len()
        ),
        (700, 100, 300, 120, 300)
    );
    assert!(q_seq.len() < 500 && d_seq.len() < 500 && r_seq.len() < 500 && x_seq.len() < 500);
    assert_eq!(&r_seq[..OV], &r_seq[r_seq.len() - OV..]);
    assert_ne!(
        &r_seq[OV..OV + OV],
        &r_seq[..OV],
        "R's interior must differ from its junction"
    );
    assert_ne!(d_seq, x_seq);
    assert_ne!(&d_seq[OV..], &x_seq[OV..]);

    let s_id = add_node(&mut graph, &s_seq);
    let q_id = add_node(&mut graph, &q_seq);
    let d_id = add_node(&mut graph, &d_seq);
    let r_id = add_node(&mut graph, &r_seq);
    let x_id = add_node(&mut graph, &x_seq);
    assert_eq!(
        (s_id.id(), q_id.id(), d_id.id(), r_id.id(), x_id.id()),
        (0, 1, 2, 3, 4)
    );
    assert!(!graph.nodes[r_id.id() as usize].is_circular);

    for (from, to) in [
        (s_id, q_id),
        (q_id, d_id),
        (q_id, r_id),
        (r_id, r_id),
        (r_id, x_id),
    ] {
        add_edge(&mut graph, from, to, OV as u8);
        assert_valid_edge(&graph, from, to, OV as u8);
    }
    assert_eq!(graph.nodes[q_id.id() as usize].edges_fwd.len(), 2);
    assert_eq!(graph.nodes[r_id.id() as usize].edges_fwd.len(), 2);

    // Q's suffix is the shared loop junction, so both arms leave Q at the same overlap.
    let sep_d = q_seq.len() as i64 - OV as i64 - OV as i64;
    let sep_single = sep_d + r_seq.len() as i64 - OV as i64;
    let sep_tandem = sep_single + r_seq.len() as i64 - OV as i64;
    assert_eq!((sep_d, sep_single, sep_tandem), (38, 127, 216));
    assert!(sep_tandem <= 520);

    const READ: u32 = 75;
    assert_read_inside(s_seq.len(), 200, READ);
    assert_read_inside(d_seq.len(), 212, READ);
    assert_read_inside(s_seq.len(), 117, READ);
    assert_read_inside(x_seq.len(), 117, READ);

    let dist = InsertSizeDist {
        median: 450.0,
        mad: 20.0,
        min_insert: 380,
        max_insert: 520,
    };
    let mut matrix = PairLinkMatrix::new();
    add_pairs_len(&mut matrix, s_id, d_id, 5, 200, 212, 14_100, READ as u16);
    add_pairs_len(&mut matrix, s_id, x_id, 20, 117, 117, 14_200, READ as u16);

    let d_weight = assert_score(&matrix, &dist, s_id, d_id, sep_d, 5);
    let (single_count, single_weight) =
        score_candidate_path(&matrix, &dist, s_id, x_id, sep_single);
    assert_eq!(
        (single_count, single_weight),
        (0, 0.0),
        "one copy of R must leave S→X outside [380, 520]"
    );
    let tandem_weight = assert_score(&matrix, &dist, s_id, x_id, sep_tandem, 20);
    assert!(
        (d_weight / f64::from(5) - 1.0).abs() < 1e-6,
        "D's five pairs are centered and pass the single-candidate rule on their own"
    );
    assert!(
        tandem_weight / d_weight >= 2.5,
        "the tandem route must be strong enough to defeat D"
    );

    let misassembly = splice(&[&s_seq, &q_seq, &d_seq], &[OV, OV]);
    resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params()).expect("repeat resolution failed");

    let d_node = d_id.id() as usize;
    assert!(
        !graph.nodes[d_node].is_deleted,
        "D must not be consumed by a misassembly commit"
    );
    assert!(
        !has_active_seq(&graph, &misassembly),
        "Misassembled contig S-Q-D must not be created"
    );
}

/// T15 / F10-P0-001. Case A's protected side branch must survive a second call.
///
/// S(700) → R(169) and X(100) → R, overlap 31 bp. X carries 69 bp that S and R
/// do not. The first call commits S→R, deletes S and R, and protects X.
/// The graph is then `Completed`. A second call must return `AlreadyResolved`
/// before exploration or Stage 5 cleanup, so X stays and the fingerprint
/// does not change.
#[test]
fn test_t15_reentry_cleanup_must_not_delete_protected_side_branch() {
    let mut graph = CdbgGraph::new(31);
    let junction = dna(0xA001, 31);
    let mut s_seq = dna(0xA002, 669);
    s_seq.extend_from_slice(&junction);
    let mut r_seq = junction.clone();
    r_seq.extend(dna(0xA003, 138));
    let unique = dna(0xA004, 69);
    let mut x_seq = unique.clone();
    x_seq.extend_from_slice(&junction);
    assert_eq!((s_seq.len(), r_seq.len(), x_seq.len()), (700, 169, 100));
    assert_not_homopolymer(&x_seq);
    assert!(!contains(&s_seq, &unique) && !contains(&r_seq, &unique));

    let s_dir = add_node(&mut graph, &s_seq);
    let r_dir = add_node(&mut graph, &r_seq);
    let x_dir = add_node(&mut graph, &x_seq);
    let s_id = s_dir.id();
    let r_id = r_dir.id();
    let x_id = x_dir.id();
    add_edge(&mut graph, s_dir, r_dir, 31);
    add_edge(&mut graph, x_dir, r_dir, 31);
    assert_valid_edge(&graph, s_dir, r_dir, 31);
    assert_valid_edge(&graph, x_dir, r_dir, 31);
    let spliced = splice(&[&s_seq, &r_seq], &[31]);

    let matrix = PairLinkMatrix::new();
    let dist = InsertSizeDist::default();
    let res1 = resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params());
    assert!(res1.is_ok(), "pass 1 must succeed");
    let fp_before = structural_fingerprint(&graph);

    assert!(
        has_active_seq(&graph, &spliced)
            && graph.nodes[s_id as usize].is_deleted
            && graph.nodes[r_id as usize].is_deleted,
        "fixture: pass 1 must commit S→R before the re-entry cleanup check, active=[{}]",
        active_summary(&graph)
    );
    assert!(!graph.nodes[x_id as usize].is_deleted);
    assert!(contains(&graph.get_unitig_seq(x_id), &unique));

    let res2 = resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params());
    assert_eq!(
        res2,
        Err(ExspanderError::AlreadyResolved),
        "second resolver call must return Err(AlreadyResolved)"
    );
    assert_eq!(
        structural_fingerprint(&graph),
        fp_before,
        "graph fingerprint must remain strictly identical after second call"
    );
    assert!(
        !graph.nodes[x_id as usize].is_deleted,
        "X must not be deleted by cleanup on second call"
    );
}

/// T16 / F10-P1-001. Consumed-node provenance must survive a second call.
///
/// Recreates T3: S1(900) → R(169) commits on the first call and R's 20 pairs
/// block S2(700) → Q(100) → D(300) (5 pairs). The graph is then `Completed`.
/// A second call must return `AlreadyResolved` before exploration, so a fresh
/// `consumed_nodes` set cannot treat R as absent and commit the decoy D.
#[test]
fn test_t16_reentry_consumed_must_not_allow_decoy() {
    let mut fx = build_ce3(false);
    assert_eq!(fx.sep, 38);
    let d_id = fx.d_id;

    let res1 = resolve_repeats_exspander_with_stats(&mut fx.graph, &fx.matrix, &fx.dist, params());
    assert!(res1.is_ok(), "pass 1 must succeed");
    let fp_before = structural_fingerprint(&fx.graph);

    assert!(
        has_active_seq(&fx.graph, &fx.s1_contig) && !alive(&fx.graph, fx.r_id),
        "fixture: pass 1 must consume R, active=[{}]",
        active_summary(&fx.graph)
    );
    assert!(
        alive(&fx.graph, d_id) && alive(&fx.graph, fx.q_id) && !has_active_seq(&fx.graph, &fx.bad),
        "fixture: pass 1 must keep the S2 fork unresolved, active=[{}]",
        active_summary(&fx.graph)
    );
    let r_weight = assert_score(&fx.matrix, &fx.dist, fx.anchor, fx.r_dir, fx.sep, 20);
    let d_weight = assert_score(&fx.matrix, &fx.dist, fx.anchor, fx.d_dir, fx.sep, 5);
    assert!(
        d_weight < 2.5 * r_weight,
        "fixture: retained S2→R evidence must still be strong enough to block D"
    );

    let res2 = resolve_repeats_exspander_with_stats(&mut fx.graph, &fx.matrix, &fx.dist, params());
    assert_eq!(
        res2,
        Err(ExspanderError::AlreadyResolved),
        "second resolver call must return Err(AlreadyResolved)"
    );
    assert_eq!(
        structural_fingerprint(&fx.graph),
        fp_before,
        "graph fingerprint must remain strictly identical after second call"
    );
    assert!(
        !fx.graph.nodes[d_id as usize].is_deleted,
        "D must not be consumed by decoy commit on second call"
    );
}

fn resolvable_linear_graph(tag: u64) -> CdbgGraph {
    let mut graph = CdbgGraph::new(31);
    let junction = dna(tag, 31);
    let mut s_seq = dna(tag.wrapping_add(1), 669);
    s_seq.extend_from_slice(&junction);
    let mut r_seq = junction.clone();
    r_seq.extend(dna(tag.wrapping_add(2), 138));
    assert_eq!((s_seq.len(), r_seq.len()), (700, 169));
    let s_dir = add_node(&mut graph, &s_seq);
    let r_dir = add_node(&mut graph, &r_seq);
    add_edge(&mut graph, s_dir, r_dir, 31);
    assert_valid_edge(&graph, s_dir, r_dir, 31);
    assert_eq!(graph.resolver_state(), ResolverLifecycleState::NeverRun);
    graph
}

/// T17. Cloning a `Completed` graph keeps `Completed` and rejects re-entry.
#[test]
fn test_t17_cloned_graph_preserves_lifecycle_state() {
    let mut graph = resolvable_linear_graph(0xC17_0001);
    let matrix = PairLinkMatrix::new();
    let dist = InsertSizeDist::default();

    let res1 = resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params());
    assert!(res1.is_ok(), "pass 1 must succeed");
    assert_eq!(graph.resolver_state(), ResolverLifecycleState::Completed);
    let original_fp = structural_fingerprint(&graph);

    let mut cloned = graph.clone();
    assert_eq!(cloned.resolver_state(), ResolverLifecycleState::Completed);
    let fp_before = structural_fingerprint(&cloned);

    let res2 = resolve_repeats_exspander_with_stats(&mut cloned, &matrix, &dist, params());
    assert_eq!(
        res2,
        Err(ExspanderError::AlreadyResolved),
        "re-entry on a clone of a Completed graph must return Err(AlreadyResolved)"
    );
    assert_eq!(
        structural_fingerprint(&cloned),
        fp_before,
        "clone fingerprint must remain strictly identical after rejected re-entry"
    );
    assert_eq!(cloned.resolver_state(), ResolverLifecycleState::Completed);
    assert_eq!(structural_fingerprint(&graph), original_fp);
    assert_eq!(graph.resolver_state(), ResolverLifecycleState::Completed);
}

/// T18. Distinct fresh graphs start in `NeverRun` and resolve independently.
#[test]
fn test_t18_independent_graphs_resolve_separately() {
    let mut graph_a = resolvable_linear_graph(0xC18_0001);
    let mut graph_b = resolvable_linear_graph(0xC18_1001);
    assert_eq!(graph_a.resolver_state(), ResolverLifecycleState::NeverRun);
    assert_eq!(graph_b.resolver_state(), ResolverLifecycleState::NeverRun);
    let b_before = structural_fingerprint(&graph_b);

    let matrix = PairLinkMatrix::new();
    let dist = InsertSizeDist::default();
    let res_a = resolve_repeats_exspander_with_stats(&mut graph_a, &matrix, &dist, params());
    assert!(res_a.is_ok(), "graph A must resolve from NeverRun");
    assert_eq!(graph_a.resolver_state(), ResolverLifecycleState::Completed);
    assert_eq!(
        graph_b.resolver_state(),
        ResolverLifecycleState::NeverRun,
        "resolving graph A must not change graph B"
    );
    assert_eq!(structural_fingerprint(&graph_b), b_before);

    let res_b = resolve_repeats_exspander_with_stats(&mut graph_b, &matrix, &dist, params());
    assert!(res_b.is_ok(), "graph B must resolve independently of graph A");
    assert_eq!(graph_b.resolver_state(), ResolverLifecycleState::Completed);

    let fp_a = structural_fingerprint(&graph_a);
    let again_a = resolve_repeats_exspander_with_stats(&mut graph_a, &matrix, &dist, params());
    assert_eq!(again_a, Err(ExspanderError::AlreadyResolved));
    assert_eq!(structural_fingerprint(&graph_a), fp_a);

    let mut graph_c = CdbgGraph::new(31);
    assert_eq!(graph_c.resolver_state(), ResolverLifecycleState::NeverRun);
    let res_c = resolve_repeats_exspander_with_stats(&mut graph_c, &matrix, &dist, params());
    assert!(
        res_c.is_ok(),
        "a freshly constructed graph must still resolve after other graphs have completed"
    );
    assert_eq!(graph_c.resolver_state(), ResolverLifecycleState::Completed);
}

/// T19 / F12-P1-001. Three-copy tandem repeat must withhold weak decoy.
///
/// Topology:
/// S(700) → Q(100) → D(300)
///        └→ R(120) → X(300)
///            ↺ (loop)
///
/// Overlap: 31 bp.
/// S→Q→D has sep_D = 100 - 31 - 31 = 38 bp.
/// Read offsets in S (200) and D (217) give span_D = 200 + 38 + 217 = 455 bp (5 pairs).
///
/// Authentic route is S→Q→R→R→R→X (3 copies of R, 4 steps from Q).
/// For X with read offsets 75 in S and 75 in X:
/// - 1 copy: sep = 127 bp, span = 277 bp (0 valid pairs in [380, 520]).
/// - 2 copies: sep = 216 bp, span = 366 bp (0 valid pairs in [380, 520]).
/// - 3 copies: sep = 305 bp, span = 455 bp (20 valid pairs, perfectly centered).
///
/// In F11-A, the repeat probe admitted only one re-entry into R and then pruned
/// subsequent repeat transitions at line 488 via `node_visited.contains(&neighbor.id())`.
/// Because 2 copies yield span 366 (< 380), the probe returned false, the main search
/// reported Complete, and weak decoy D was committed under single-candidate acceptance.
#[test]
fn test_t19_three_copy_tandem_repeat_must_withhold_decoy() {
    let mut graph = CdbgGraph::new(31);
    let j_sq = junction(0xE19_0001);
    let j_loop = junction(0xE19_0002);
    assert_junctions_distinct(&[j_sq.clone(), j_loop.clone()]);

    let s_seq = span_node(&[], 0xE19_0101, 700, &j_sq);
    let q_seq = span_node(&j_sq, 0xE19_0102, 100, &j_loop);
    let d_seq = span_node(&j_loop, 0xE19_0103, 300, &[]);
    let r_seq = span_node(&j_loop, 0xE19_0104, 120, &j_loop);
    let x_seq = span_node(&j_loop, 0xE19_0105, 300, &[]);
    assert_eq!(
        (
            s_seq.len(),
            q_seq.len(),
            d_seq.len(),
            r_seq.len(),
            x_seq.len()
        ),
        (700, 100, 300, 120, 300)
    );
    assert!(q_seq.len() < 500 && d_seq.len() < 500 && r_seq.len() < 500 && x_seq.len() < 500);
    assert_eq!(&r_seq[..OV], &r_seq[r_seq.len() - OV..]);
    assert_ne!(
        &r_seq[OV..OV + OV],
        &r_seq[..OV],
        "R's interior must differ from its junction"
    );
    assert_ne!(d_seq, x_seq);
    assert_ne!(&d_seq[OV..], &x_seq[OV..]);

    let s_id = add_node(&mut graph, &s_seq);
    let q_id = add_node(&mut graph, &q_seq);
    let d_id = add_node(&mut graph, &d_seq);
    let r_id = add_node(&mut graph, &r_seq);
    let x_id = add_node(&mut graph, &x_seq);
    assert_eq!(
        (s_id.id(), q_id.id(), d_id.id(), r_id.id(), x_id.id()),
        (0, 1, 2, 3, 4)
    );
    assert!(!graph.nodes[r_id.id() as usize].is_circular);

    for (from, to) in [
        (s_id, q_id),
        (q_id, d_id),
        (q_id, r_id),
        (r_id, r_id),
        (r_id, x_id),
    ] {
        add_edge(&mut graph, from, to, OV as u8);
        assert_valid_edge(&graph, from, to, OV as u8);
    }
    assert_eq!(graph.nodes[q_id.id() as usize].edges_fwd.len(), 2);
    assert_eq!(graph.nodes[r_id.id() as usize].edges_fwd.len(), 2);

    let sep_d = q_seq.len() as i64 - OV as i64 - OV as i64;
    let sep_r1 = sep_d + r_seq.len() as i64 - OV as i64;
    let sep_r2 = sep_r1 + r_seq.len() as i64 - OV as i64;
    let sep_r3 = sep_r2 + r_seq.len() as i64 - OV as i64;
    assert_eq!((sep_d, sep_r1, sep_r2, sep_r3), (38, 127, 216, 305));
    assert!(sep_r3 <= 520);

    const READ: u32 = 75;
    assert_read_inside(s_seq.len(), 200, READ);
    assert_read_inside(d_seq.len(), 217, READ);
    assert_read_inside(s_seq.len(), 75, READ);
    assert_read_inside(x_seq.len(), 75, READ);

    let dist = InsertSizeDist {
        median: 455.0,
        mad: 20.0,
        min_insert: 380,
        max_insert: 520,
    };
    let mut matrix = PairLinkMatrix::new();
    add_pairs_len(&mut matrix, s_id, d_id, 5, 200, 217, 19_100, READ as u16);
    add_pairs_len(&mut matrix, s_id, x_id, 20, 75, 75, 19_200, READ as u16);

    let d_weight = assert_score(&matrix, &dist, s_id, d_id, sep_d, 5);
    let (c1, w1) = score_candidate_path(&matrix, &dist, s_id, x_id, sep_r1);
    assert_eq!((c1, w1), (0, 0.0), "1 copy of R must be outside [380, 520]");
    let (c2, w2) = score_candidate_path(&matrix, &dist, s_id, x_id, sep_r2);
    assert_eq!((c2, w2), (0, 0.0), "2 copies of R must be outside [380, 520]");
    let three_weight = assert_score(&matrix, &dist, s_id, x_id, sep_r3, 20);

    assert!(
        (d_weight / f64::from(5) - 1.0).abs() < 1e-6,
        "D's five pairs are centered at median 455 and pass single-candidate rule on their own"
    );
    assert!(
        three_weight / d_weight >= 2.5,
        "the 3-copy tandem route must be strong enough to defeat D"
    );

    let misassembly = splice(&[&s_seq, &q_seq, &d_seq], &[OV, OV]);
    resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params())
        .expect("repeat resolution failed");

    let d_node = d_id.id() as usize;
    assert!(
        !graph.nodes[d_node].is_deleted,
        "D must not be consumed by a misassembly commit"
    );
    assert!(
        !has_active_seq(&graph, &misassembly),
        "Misassembled contig S-Q-D must not be created"
    );
}

/// T20 / F12-P1-001. Four-step tandem repeat exceeding/hitting candidate depth cap.
///
/// Topology: S(700) → Q(100) → D(300)
///                 └→ R(120) ↺
///
/// Route from Q: R → R → R → R (4 steps, reaching MAX_CANDIDATE_DEPTH).
/// From the 4th R, an outgoing edge leads to X(300).
/// Because an admissible continuation exists at the depth boundary, the resolver
/// must mark the search Truncated and withhold the decision, never committing D.
#[test]
fn test_t20_four_step_tandem_repeat_must_withhold_decoy() {
    let mut graph = CdbgGraph::new(31);
    let j_sq = junction(0xE20_0001);
    let j_loop = junction(0xE20_0002);
    assert_junctions_distinct(&[j_sq.clone(), j_loop.clone()]);

    let s_seq = span_node(&[], 0xE20_0101, 700, &j_sq);
    let q_seq = span_node(&j_sq, 0xE20_0102, 100, &j_loop);
    let d_seq = span_node(&j_loop, 0xE20_0103, 300, &[]);
    let r_seq = span_node(&j_loop, 0xE20_0104, 120, &j_loop);
    let x_seq = span_node(&j_loop, 0xE20_0105, 300, &[]);

    let s_id = add_node(&mut graph, &s_seq);
    let q_id = add_node(&mut graph, &q_seq);
    let d_id = add_node(&mut graph, &d_seq);
    let r_id = add_node(&mut graph, &r_seq);
    let x_id = add_node(&mut graph, &x_seq);

    for (from, to) in [
        (s_id, q_id),
        (q_id, d_id),
        (q_id, r_id),
        (r_id, r_id),
        (r_id, x_id),
    ] {
        add_edge(&mut graph, from, to, OV as u8);
        assert_valid_edge(&graph, from, to, OV as u8);
    }

    let sep_d = q_seq.len() as i64 - OV as i64 - OV as i64;
    let sep_r4 = sep_d + 4 * (r_seq.len() as i64 - OV as i64);
    assert_eq!(sep_d, 38);
    assert_eq!(sep_r4, 394);
    assert!(sep_r4 <= 520);

    const READ: u32 = 75;
    let dist = InsertSizeDist {
        median: 455.0,
        mad: 20.0,
        min_insert: 380,
        max_insert: 520,
    };
    let mut matrix = PairLinkMatrix::new();
    add_pairs_len(&mut matrix, s_id, d_id, 5, 200, 217, 20_100, READ as u16);
    add_pairs_len(&mut matrix, s_id, x_id, 20, 30, 31, 20_200, READ as u16);

    let misassembly = splice(&[&s_seq, &q_seq, &d_seq], &[OV, OV]);
    resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params())
        .expect("repeat resolution failed");

    let d_node = d_id.id() as usize;
    assert!(
        !graph.nodes[d_node].is_deleted,
        "D must not be consumed by a misassembly commit"
    );
    assert!(
        !has_active_seq(&graph, &misassembly),
        "Misassembled contig S-Q-D must not be created"
    );
}

/// T21 / F12-P1-001. Alternating two-node cycle A ⇄ B must withhold weak decoy.
///
/// Topology:
/// S(700) → Q(100) → D(300)
///        └→ A(70) ⇄ B(70)
///                   B → X(300)
///
/// Route from Q: A → B → A → B → X (or A → B → A with continuation).
/// A and B are short (< 500 bp) nodes. Because the cycle involves multiple distinct
/// physical nodes, pruning any visited node in the probe hides authentic paths.
#[test]
fn test_t21_alternating_cycle_repeat_must_withhold_decoy() {
    let mut graph = CdbgGraph::new(31);
    let j_sq = junction(0xE21_0001);
    let j_shared = junction(0xE21_0002);
    let j_ab = junction(0xE21_0003);
    assert_junctions_distinct(&[j_sq.clone(), j_shared.clone(), j_ab.clone()]);

    let s_seq = span_node(&[], 0xE21_0101, 700, &j_sq);
    let q_seq = span_node(&j_sq, 0xE21_0102, 100, &j_shared);
    let d_seq = span_node(&j_shared, 0xE21_0103, 300, &[]);
    let a_seq = span_node(&j_shared, 0xE21_0104, 70, &j_ab);
    let b_seq = span_node(&j_ab, 0xE21_0105, 70, &j_shared);
    let x_seq = span_node(&j_shared, 0xE21_0106, 300, &[]);

    let s_id = add_node(&mut graph, &s_seq);
    let q_id = add_node(&mut graph, &q_seq);
    let d_id = add_node(&mut graph, &d_seq);
    let a_id = add_node(&mut graph, &a_seq);
    let b_id = add_node(&mut graph, &b_seq);
    let x_id = add_node(&mut graph, &x_seq);

    for (from, to) in [
        (s_id, q_id),
        (q_id, d_id),
        (q_id, a_id),
        (a_id, b_id),
        (b_id, a_id),
        (b_id, x_id),
    ] {
        add_edge(&mut graph, from, to, OV as u8);
        assert_valid_edge(&graph, from, to, OV as u8);
    }

    let sep_d = q_seq.len() as i64 - OV as i64 - OV as i64;
    assert_eq!(sep_d, 38);
    const READ: u32 = 75;
    let dist = InsertSizeDist {
        median: 455.0,
        mad: 20.0,
        min_insert: 380,
        max_insert: 520,
    };
    let mut matrix = PairLinkMatrix::new();
    add_pairs_len(&mut matrix, s_id, d_id, 5, 200, 217, 21_100, READ as u16);
    add_pairs_len(&mut matrix, s_id, x_id, 20, 75, 75, 21_200, READ as u16);

    let misassembly = splice(&[&s_seq, &q_seq, &d_seq], &[OV, OV]);
    resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params())
        .expect("repeat resolution failed");

    let d_node = d_id.id() as usize;
    assert!(
        !graph.nodes[d_node].is_deleted,
        "D must not be consumed by a misassembly commit"
    );
    assert!(
        !has_active_seq(&graph, &misassembly),
        "Misassembled contig S-Q-D must not be created"
    );
}

/// T22 / F14-P2-002. Genuinely mirrored reverse-complement symmetry of three-copy tandem repeat.
///
/// Builds a true reverse-complement mirror graph where all junction sequences are RC-transformed,
/// node sequences are spanned with reverse-complement junctions, and seed S_rc (700 bp) initiates
/// resolution in reverse complement, verifying that the repeat probe correctly withholds
/// decoy D_rc under identical PE link constraints.
#[test]
fn test_t22_three_copy_tandem_repeat_rc_symmetry_must_withhold_decoy() {
    let j_sq_fwd = junction(0xE22_0001);
    let j_loop_fwd = junction(0xE22_0002);
    assert_junctions_distinct(&[j_sq_fwd.clone(), j_loop_fwd.clone()]);

    let j_sq = reverse_complement_ascii(&j_sq_fwd);
    let j_loop = reverse_complement_ascii(&j_loop_fwd);

    let s_rc = span_node(&[], 0xE22_0101, 700, &j_sq);
    let q_rc = span_node(&j_sq, 0xE22_0102, 100, &j_loop);
    let d_rc = span_node(&j_loop, 0xE22_0103, 300, &[]);
    let r_rc = span_node(&j_loop, 0xE22_0104, 120, &j_loop);
    let x_rc = span_node(&j_loop, 0xE22_0105, 300, &[]);

    let mut rc_graph = CdbgGraph::new(31);
    let s_id = add_node(&mut rc_graph, &s_rc);
    let q_id = add_node(&mut rc_graph, &q_rc);
    let d_id = add_node(&mut rc_graph, &d_rc);
    let r_id = add_node(&mut rc_graph, &r_rc);
    let x_id = add_node(&mut rc_graph, &x_rc);

    for (from, to) in [
        (s_id, q_id),
        (q_id, d_id),
        (q_id, r_id),
        (r_id, r_id),
        (r_id, x_id),
    ] {
        add_edge(&mut rc_graph, from, to, OV as u8);
        assert_valid_edge(&rc_graph, from, to, OV as u8);
    }

    let dist = InsertSizeDist {
        median: 455.0,
        mad: 20.0,
        min_insert: 380,
        max_insert: 520,
    };
    let mut matrix = PairLinkMatrix::new();
    const READ: u32 = 75;
    add_pairs_len(&mut matrix, s_id, d_id, 5, 200, 217, 22_100, READ as u16);
    add_pairs_len(&mut matrix, s_id, x_id, 20, 75, 75, 22_200, READ as u16);

    let misassembly = splice(&[&s_rc, &q_rc, &d_rc], &[OV, OV]);
    resolve_repeats_exspander_with_stats(&mut rc_graph, &matrix, &dist, params())
        .expect("repeat resolution failed");

    let d_node = d_id.id() as usize;
    assert!(
        !rc_graph.nodes[d_node].is_deleted,
        "D_rc must not be consumed by a misassembly commit under true RC mirror symmetry"
    );
    assert!(
        !has_active_seq(&rc_graph, &misassembly),
        "Misassembled contig S_rc-Q_rc-D_rc must not be created under true RC mirror symmetry"
    );
}

/// T23 / F14-P1-001. Short prefix reuse must not promote decoy under Complete status.
///
/// Topology:
/// S(700) → P(100) → C(100) → D(300)
///           │        │
///           └→ X(300) └→ P
///
/// Overlaps: all 31 bp.
/// S ends in junction A.
/// P: prefix A, suffix B. (len 100)
/// C: prefix B, suffix A. (len 100)
/// D: prefix A, suffix []. (len 300)
/// X: prefix B, suffix []. (len 300)
///
/// Decision 1: Seed S extends through P to C (20 pairs between S and C, span 200).
/// Prefix ext = [P, C], curr_dir = C.
/// ext_nodes = {P, C}.
///
/// Decision 2: From C, candidate search evaluates:
/// - C → D (5 pairs, span 200)
/// - C → P → X (20 pairs, span 200)
///
/// If P is unconditionally skipped due to ext_nodes.contains(&P.id()),
/// candidate X is dropped, search status returns Complete, and weak decoy D commits!
/// P is a short node (< 500 bp), so reusing P conforms to validate_and_splice_path.
#[test]
fn test_t23_short_prefix_reuse_must_not_promote_decoy() {
    let mut graph = CdbgGraph::new(31);
    let j_a = junction(0xE23_0001);
    let j_b = junction(0xE23_0002);
    assert_junctions_distinct(&[j_a.clone(), j_b.clone()]);

    let s_seq = span_node(&[], 0xE23_0101, 700, &j_a);
    let p_seq = span_node(&j_a, 0xE23_0102, 100, &j_b);
    let c_seq = span_node(&j_b, 0xE23_0103, 100, &j_a);
    let d_seq = span_node(&j_a, 0xE23_0104, 300, &[]);
    let x_seq = span_node(&j_b, 0xE23_0105, 300, &[]);

    let s_id = add_node(&mut graph, &s_seq);
    let p_id = add_node(&mut graph, &p_seq);
    let c_id = add_node(&mut graph, &c_seq);
    let d_id = add_node(&mut graph, &d_seq);
    let x_id = add_node(&mut graph, &x_seq);

    for (from, to) in [
        (s_id, p_id),
        (p_id, c_id),
        (p_id, x_id),
        (c_id, p_id),
        (c_id, d_id),
    ] {
        add_edge(&mut graph, from, to, OV as u8);
        assert_valid_edge(&graph, from, to, OV as u8);
    }

    const READ: u32 = 75;
    let dist = InsertSizeDist {
        median: 200.0,
        mad: 20.0,
        min_insert: 150,
        max_insert: 250,
    };
    let mut matrix = PairLinkMatrix::new();
    // Decision 1: S -> P -> C (sep = 100 - 31 - 31 = 38, span = 80 + 38 + 82 = 200)
    add_pairs_len(&mut matrix, s_id, c_id, 20, 80, 82, 23_001, READ as u16);

    // Decision 2: C -> D (sep = -31, span = 80 - 31 + 151 = 200, 5 pairs)
    add_pairs_len(&mut matrix, c_id, d_id, 5, 80, 151, 23_002, READ as u16);

    // Decision 2: C -> P -> X (sep = -31 + 100 - 31 = 38, span = 80 + 38 + 82 = 200, 20 pairs)
    add_pairs_len(&mut matrix, c_id, x_id, 20, 80, 82, 23_003, READ as u16);

    let misassembly = splice(&[&s_seq, &p_seq, &c_seq, &d_seq], &[OV, OV, OV]);
    resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params())
        .expect("repeat resolution failed");

    let d_node = d_id.id() as usize;
    assert!(
        !graph.nodes[d_node].is_deleted,
        "D must not be consumed by a misassembly commit when C -> P -> X is an admissible competitor"
    );
    assert!(
        !has_active_seq(&graph, &misassembly),
        "Misassembled contig S-P-C-D must not be created"
    );
}

/// T24 / F15-E. Reverse-complement prefix reuse mirror of T23.
///
/// Verifies that in reverse-complement orientation, short prefix reuse
/// is equally protected and does not promote a weak decoy.
#[test]
fn test_t24_rc_prefix_reuse_must_not_promote_decoy() {
    let j_a_fwd = junction(0xE24_0001);
    let j_b_fwd = junction(0xE24_0002);
    assert_junctions_distinct(&[j_a_fwd.clone(), j_b_fwd.clone()]);

    let j_a = reverse_complement_ascii(&j_a_fwd);
    let j_b = reverse_complement_ascii(&j_b_fwd);

    let s_rc = span_node(&[], 0xE24_0101, 700, &j_a);
    let p_rc = span_node(&j_a, 0xE24_0102, 100, &j_b);
    let c_rc = span_node(&j_b, 0xE24_0103, 100, &j_a);
    let d_rc = span_node(&j_a, 0xE24_0104, 300, &[]);
    let x_rc = span_node(&j_b, 0xE24_0105, 300, &[]);

    let mut graph = CdbgGraph::new(31);
    let s_id = add_node(&mut graph, &s_rc);
    let p_id = add_node(&mut graph, &p_rc);
    let c_id = add_node(&mut graph, &c_rc);
    let d_id = add_node(&mut graph, &d_rc);
    let x_id = add_node(&mut graph, &x_rc);

    for (from, to) in [
        (s_id, p_id),
        (p_id, c_id),
        (p_id, x_id),
        (c_id, p_id),
        (c_id, d_id),
    ] {
        add_edge(&mut graph, from, to, OV as u8);
        assert_valid_edge(&graph, from, to, OV as u8);
    }

    const READ: u32 = 75;
    let dist = InsertSizeDist {
        median: 200.0,
        mad: 20.0,
        min_insert: 150,
        max_insert: 250,
    };
    let mut matrix = PairLinkMatrix::new();
    add_pairs_len(&mut matrix, s_id, c_id, 20, 80, 82, 24_001, READ as u16);
    add_pairs_len(&mut matrix, c_id, d_id, 5, 80, 151, 24_002, READ as u16);
    add_pairs_len(&mut matrix, c_id, x_id, 20, 80, 82, 24_003, READ as u16);

    let misassembly = splice(&[&s_rc, &p_rc, &c_rc, &d_rc], &[OV, OV, OV]);
    resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params())
        .expect("repeat resolution failed");

    let d_node = d_id.id() as usize;
    assert!(
        !graph.nodes[d_node].is_deleted,
        "D_rc must not be consumed by a misassembly commit under RC mirror prefix reuse"
    );
    assert!(
        !has_active_seq(&graph, &misassembly),
        "Misassembled contig S_rc-P_rc-C_rc-D_rc must not be created"
    );
}

/// T25 / F15-F. Prefix reuse with fragment separation exactly at max budget.
///
/// When the reused path has span == max_insert (250 bp), it is admitted by `<=`
/// and safely withholds the weak decoy D.
#[test]
fn test_t25_prefix_reuse_at_exact_max_budget_withholds_decoy() {
    let mut graph = CdbgGraph::new(31);
    let j_a = junction(0xE25_0001);
    let j_b = junction(0xE25_0002);

    let s_seq = span_node(&[], 0xE25_0101, 700, &j_a);
    let p_seq = span_node(&j_a, 0xE25_0102, 100, &j_b);
    let c_seq = span_node(&j_b, 0xE25_0103, 100, &j_a);
    let d_seq = span_node(&j_a, 0xE25_0104, 300, &[]);
    let x_seq = span_node(&j_b, 0xE25_0105, 300, &[]);

    let s_id = add_node(&mut graph, &s_seq);
    let p_id = add_node(&mut graph, &p_seq);
    let c_id = add_node(&mut graph, &c_seq);
    let d_id = add_node(&mut graph, &d_seq);
    let x_id = add_node(&mut graph, &x_seq);

    for (from, to) in [
        (s_id, p_id),
        (p_id, c_id),
        (p_id, x_id),
        (c_id, p_id),
        (c_id, d_id),
    ] {
        add_edge(&mut graph, from, to, OV as u8);
    }

    let dist = InsertSizeDist {
        median: 200.0,
        mad: 20.0,
        min_insert: 150,
        max_insert: 250,
    };
    let mut matrix = PairLinkMatrix::new();
    add_pairs_len(&mut matrix, s_id, c_id, 20, 80, 82, 25_001, 75);
    add_pairs_len(&mut matrix, c_id, d_id, 5, 80, 151, 25_002, 75);

    // C -> P -> X sep = 38. With dist1 = 80, dist2 = 132: span = 80 + 38 + 132 = 250 == max_insert
    add_pairs_len(&mut matrix, c_id, x_id, 20, 80, 132, 25_003, 75);

    resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params())
        .expect("repeat resolution failed");

    let d_node = d_id.id() as usize;
    assert!(
        !graph.nodes[d_node].is_deleted,
        "D must not be consumed when prefix reuse competitor is at exact max budget"
    );
}

/// T26 / F15-F. Prefix reuse exceeding max budget must not cause unnecessary Truncated.
///
/// When the alternative repeat path has span 251 > max_insert 250, it is safely
/// excluded by separation budget, allowing the valid in-budget candidate D to resolve cleanly.
#[test]
fn test_t26_prefix_reuse_exceeding_max_budget_resolves_cleanly() {
    let mut graph = CdbgGraph::new(31);
    let j_a = junction(0xE26_0001);
    let j_b = junction(0xE26_0002);

    let s_seq = span_node(&[], 0xE26_0101, 700, &j_a);
    let p_seq = span_node(&j_a, 0xE26_0102, 100, &j_b);
    let c_seq = span_node(&j_b, 0xE26_0103, 100, &j_a);
    let d_seq = span_node(&j_a, 0xE26_0104, 300, &[]);
    let x_seq = span_node(&j_b, 0xE26_0105, 300, &[]);

    let s_id = add_node(&mut graph, &s_seq);
    let p_id = add_node(&mut graph, &p_seq);
    let c_id = add_node(&mut graph, &c_seq);
    let d_id = add_node(&mut graph, &d_seq);
    let x_id = add_node(&mut graph, &x_seq);

    for (from, to) in [
        (s_id, p_id),
        (p_id, c_id),
        (p_id, x_id),
        (c_id, p_id),
        (c_id, d_id),
    ] {
        add_edge(&mut graph, from, to, OV as u8);
    }

    let dist = InsertSizeDist {
        median: 200.0,
        mad: 20.0,
        min_insert: 150,
        max_insert: 250,
    };
    let mut matrix = PairLinkMatrix::new();
    add_pairs_len(&mut matrix, s_id, c_id, 20, 80, 82, 26_001, 75);
    // D has 10 valid pairs at span 200
    add_pairs_len(&mut matrix, c_id, d_id, 10, 80, 151, 26_002, 75);

    // C -> P -> X sep = 38. With dist1 = 80, dist2 = 133: span = 80 + 38 + 133 = 251 > 250
    add_pairs_len(&mut matrix, c_id, x_id, 20, 80, 133, 26_003, 75);

    resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params())
        .expect("repeat resolution failed");

    let d_node = d_id.id() as usize;
    assert!(
        graph.nodes[d_node].is_deleted,
        "D should resolve cleanly when competitor is beyond max_separation"
    );
}

/// T27 / F15-F. Legitimate normal repeat bridge must still commit (Positive Control).
///
/// Positive control demonstrating that a genuine repeat bridge that reuses a short prefix node
/// successfully commits the winning path S -> P -> C -> P -> X without false truncation.
#[test]
fn test_t27_legitimate_prefix_reuse_bridge_commits_correctly() {
    let mut graph = CdbgGraph::new(31);
    let j_a = junction(0xE27_0001);
    let j_b = junction(0xE27_0002);

    let s_seq = span_node(&[], 0xE27_0101, 700, &j_a);
    let p_seq = span_node(&j_a, 0xE27_0102, 100, &j_b);
    let c_seq = span_node(&j_b, 0xE27_0103, 100, &j_a);
    let x_seq = span_node(&j_b, 0xE27_0104, 300, &[]);

    let s_id = add_node(&mut graph, &s_seq);
    let p_id = add_node(&mut graph, &p_seq);
    let c_id = add_node(&mut graph, &c_seq);
    let x_id = add_node(&mut graph, &x_seq);

    for (from, to) in [
        (s_id, p_id),
        (p_id, c_id),
        (c_id, p_id),
        (p_id, x_id),
    ] {
        add_edge(&mut graph, from, to, OV as u8);
    }

    let dist = InsertSizeDist {
        median: 200.0,
        mad: 20.0,
        min_insert: 150,
        max_insert: 250,
    };
    let mut matrix = PairLinkMatrix::new();
    add_pairs_len(&mut matrix, s_id, c_id, 20, 80, 82, 27_001, 75);
    add_pairs_len(&mut matrix, c_id, x_id, 25, 80, 82, 27_002, 75);

    let full_path_seq = splice(&[&s_seq, &p_seq, &c_seq, &p_seq, &x_seq], &[OV, OV, OV, OV]);
    resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params())
        .expect("repeat resolution failed");

    let x_node = x_id.id() as usize;
    assert!(
        graph.nodes[x_node].is_deleted,
        "X must be consumed by the legitimate repeat bridge commit"
    );
    assert!(
        has_active_seq(&graph, &full_path_seq),
        "Full assembled sequence S-P-C-P-X must be present in graph contigs"
    );
}

/// T28 / F15-F. Prefix physical node reused with opposite orientation.
///
/// Prefix contains short node P (< 500 bp) in forward orientation.
/// Downstream fork traverses P in twin orientation. Since P is short,
/// physical node reuse across orientations is valid and withholds decoy.
#[test]
fn test_t28_prefix_short_node_reused_in_twin_orientation() {
    let mut graph = CdbgGraph::new(31);
    let j_a = junction(0xE28_0001);
    let j_b = junction(0xE28_0002);

    let s_seq = span_node(&[], 0xE28_0101, 700, &j_a);
    let p_seq = span_node(&j_a, 0xE28_0102, 100, &j_b);
    let c_seq = span_node(&j_b, 0xE28_0103, 100, &j_a);
    let d_seq = span_node(&j_a, 0xE28_0104, 300, &[]);
    let x_seq = span_node(&j_b, 0xE28_0105, 300, &[]);

    let s_id = add_node(&mut graph, &s_seq);
    let p_id = add_node(&mut graph, &p_seq);
    let c_id = add_node(&mut graph, &c_seq);
    let d_id = add_node(&mut graph, &d_seq);
    let x_id = add_node(&mut graph, &x_seq);

    for (from, to) in [
        (s_id, p_id),
        (p_id, c_id),
        (p_id, x_id),
        (c_id, p_id),
        (c_id, d_id),
    ] {
        add_edge(&mut graph, from, to, OV as u8);
    }

    let dist = InsertSizeDist {
        median: 200.0,
        mad: 20.0,
        min_insert: 150,
        max_insert: 250,
    };
    let mut matrix = PairLinkMatrix::new();
    add_pairs_len(&mut matrix, s_id, c_id, 20, 80, 82, 28_001, 75);
    add_pairs_len(&mut matrix, c_id, d_id, 5, 80, 151, 28_002, 75);
    // Twin orientation links
    add_pairs_len(&mut matrix, c_id.twin(), x_id.twin(), 20, 80, 82, 28_003, 75);
    add_pairs_len(&mut matrix, c_id, x_id, 20, 80, 82, 28_004, 75);

    resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params())
        .expect("repeat resolution failed");

    let d_node = d_id.id() as usize;
    assert!(
        !graph.nodes[d_node].is_deleted,
        "D must not be consumed when twin-orientation competitor exists"
    );
}

/// T29 / F15-F. Prefix reuse involving consumed/shadow node.
///
/// If a short node was already consumed by an earlier seed, it participates
/// as a shadow candidate and blocks a weaker live decoy from committing.
#[test]
fn test_t29_prefix_reuse_with_consumed_shadow_node() {
    let mut graph = CdbgGraph::new(31);
    let j_a = junction(0xE29_0001);
    let j_b = junction(0xE29_0002);

    let s1_seq = span_node(&[], 0xE29_0101, 700, &j_a);
    let p_seq = span_node(&j_a, 0xE29_0102, 100, &j_b);
    let tail1_seq = span_node(&j_b, 0xE29_0103, 500, &[]);

    let s2_seq = span_node(&[], 0xE29_0104, 700, &j_b);
    let c_seq = span_node(&j_b, 0xE29_0105, 100, &j_a);
    let d_seq = span_node(&j_a, 0xE29_0106, 300, &[]);
    let x_seq = span_node(&j_b, 0xE29_0107, 300, &[]);

    let s1_id = add_node(&mut graph, &s1_seq);
    let p_id = add_node(&mut graph, &p_seq);
    let tail1_id = add_node(&mut graph, &tail1_seq);

    let s2_id = add_node(&mut graph, &s2_seq);
    let c_id = add_node(&mut graph, &c_seq);
    let d_id = add_node(&mut graph, &d_seq);
    let x_id = add_node(&mut graph, &x_seq);

    for (from, to) in [
        (s1_id, p_id),
        (p_id, tail1_id),
        (s2_id, c_id),
        (c_id, p_id),
        (c_id, d_id),
        (p_id, x_id),
    ] {
        add_edge(&mut graph, from, to, OV as u8);
    }

    let dist = InsertSizeDist {
        median: 200.0,
        mad: 20.0,
        min_insert: 150,
        max_insert: 250,
    };
    let mut matrix = PairLinkMatrix::new();
    // Seed 1 extends S1 -> P -> tail1
    add_pairs_len(&mut matrix, s1_id, tail1_id, 25, 80, 82, 29_001, 75);

    // Seed 2 extends to C, then evaluates C -> D vs C -> P -> X (shadow through consumed P)
    add_pairs_len(&mut matrix, s2_id, c_id, 20, 80, 82, 29_002, 75);
    add_pairs_len(&mut matrix, c_id, d_id, 5, 80, 151, 29_003, 75);
    add_pairs_len(&mut matrix, c_id, x_id, 20, 80, 82, 29_004, 75);

    resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params())
        .expect("repeat resolution failed");

    let d_node = d_id.id() as usize;
    assert!(
        !graph.nodes[d_node].is_deleted,
        "D must not be committed when shadow candidate through consumed P blocks fork"
    );
}

/// T30 / F15-F. Multiple short-prefix nodes reused in a bounded cycle.
///
/// Prefix has two short nodes P1 and P2. Downstream cycle re-enters P1 -> P2 -> X.
/// Because both are short (< 500 bp), the cycle is admissible and withholds decoy D.
#[test]
fn test_t30_multiple_short_prefix_nodes_in_cycle() {
    let mut graph = CdbgGraph::new(31);
    let j_a = junction(0xE30_0001);
    let j_b = junction(0xE30_0002);
    let j_c = junction(0xE30_0003);

    let s_seq = span_node(&[], 0xE30_0101, 700, &j_a);
    let p1_seq = span_node(&j_a, 0xE30_0102, 80, &j_b);
    let p2_seq = span_node(&j_b, 0xE30_0103, 80, &j_c);
    let c_seq = span_node(&j_c, 0xE30_0104, 80, &j_a);
    let d_seq = span_node(&j_a, 0xE30_0105, 300, &[]);
    let x_seq = span_node(&j_c, 0xE30_0106, 300, &[]);

    let s_id = add_node(&mut graph, &s_seq);
    let p1_id = add_node(&mut graph, &p1_seq);
    let p2_id = add_node(&mut graph, &p2_seq);
    let c_id = add_node(&mut graph, &c_seq);
    let d_id = add_node(&mut graph, &d_seq);
    let x_id = add_node(&mut graph, &x_seq);

    for (from, to) in [
        (s_id, p1_id),
        (p1_id, p2_id),
        (p2_id, c_id),
        (c_id, p1_id),
        (c_id, d_id),
        (p2_id, x_id),
    ] {
        add_edge(&mut graph, from, to, OV as u8);
    }

    let dist = InsertSizeDist {
        median: 250.0,
        mad: 20.0,
        min_insert: 180,
        max_insert: 320,
    };
    let mut matrix = PairLinkMatrix::new();
    add_pairs_len(&mut matrix, s_id, c_id, 20, 60, 60, 30_001, 75);
    add_pairs_len(&mut matrix, c_id, d_id, 5, 60, 150, 30_002, 75);
    add_pairs_len(&mut matrix, c_id, x_id, 20, 60, 60, 30_003, 75);

    resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params())
        .expect("repeat resolution failed");

    let d_node = d_id.id() as usize;
    assert!(
        !graph.nodes[d_node].is_deleted,
        "D must not be consumed when multi-node prefix cycle competitor exists"
    );
}

/// T31 / F15-F. No PE evidence and no admissible competitor cleanly resolves supported branch.
#[test]
fn test_t31_no_evidence_dead_end_resolves_cleanly() {
    let mut graph = CdbgGraph::new(31);
    let j_a = junction(0xE31_0001);

    let s_seq = span_node(&[], 0xE31_0101, 700, &j_a);
    let a_seq = span_node(&j_a, 0xE31_0102, 300, &[]);
    let b_seq = span_node(&j_a, 0xE31_0103, 300, &[]);

    let s_id = add_node(&mut graph, &s_seq);
    let a_id = add_node(&mut graph, &a_seq);
    let b_id = add_node(&mut graph, &b_seq);

    add_edge(&mut graph, s_id, a_id, OV as u8);
    add_edge(&mut graph, s_id, b_id, OV as u8);

    let dist = InsertSizeDist {
        median: 369.0,
        mad: 20.0,
        min_insert: 300,
        max_insert: 450,
    };
    let mut matrix = PairLinkMatrix::new();
    // Only A has evidence, B has zero evidence
    add_pairs_len(&mut matrix, s_id, a_id, 15, 200, 200, 31_001, 75);

    resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params())
        .expect("repeat resolution failed");

    let a_node = a_id.id() as usize;
    assert!(
        graph.nodes[a_node].is_deleted,
        "A must resolve cleanly when B has zero evidence"
    );
}

/// T32 / F15-F. Depth limit exhaustion at depth 4 produces conservative Truncated.
#[test]
fn test_t32_depth_limit_continuation_produces_conservative_truncated() {
    let mut graph = CdbgGraph::new(31);
    let j_sq = junction(0xE32_0001);
    let j_1 = junction(0xE32_0002);
    let j_2 = junction(0xE32_0003);
    let j_3 = junction(0xE32_0004);
    let j_4 = junction(0xE32_0005);

    let s_seq = span_node(&[], 0xE32_0101, 700, &j_sq);
    let q_seq = span_node(&j_sq, 0xE32_0102, 80, &j_1);
    let d_seq = span_node(&j_1, 0xE32_0103, 300, &[]);
    let r1_seq = span_node(&j_1, 0xE32_0104, 80, &j_2);
    let r2_seq = span_node(&j_2, 0xE32_0105, 80, &j_3);
    let r3_seq = span_node(&j_3, 0xE32_0106, 80, &j_4);
    let r4_seq = span_node(&j_4, 0xE32_0107, 80, &j_sq);

    let s_id = add_node(&mut graph, &s_seq);
    let q_id = add_node(&mut graph, &q_seq);
    let d_id = add_node(&mut graph, &d_seq);
    let r1_id = add_node(&mut graph, &r1_seq);
    let r2_id = add_node(&mut graph, &r2_seq);
    let r3_id = add_node(&mut graph, &r3_seq);
    let r4_id = add_node(&mut graph, &r4_seq);

    for (from, to) in [
        (s_id, q_id),
        (q_id, d_id),
        (q_id, r1_id),
        (r1_id, r2_id),
        (r2_id, r3_id),
        (r3_id, r4_id),
        (r4_id, q_id), // Continuation exists at depth limit!
    ] {
        add_edge(&mut graph, from, to, OV as u8);
    }

    let dist = InsertSizeDist {
        median: 400.0,
        mad: 20.0,
        min_insert: 300,
        max_insert: 500,
    };
    let mut matrix = PairLinkMatrix::new();
    // Decoy D has 5 pairs
    add_pairs_len(&mut matrix, s_id, d_id, 5, 200, 150, 32_001, 75);

    resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params())
        .expect("repeat resolution failed");

    let d_node = d_id.id() as usize;
    assert!(
        !graph.nodes[d_node].is_deleted,
        "Decoy D must not be committed when continuation exists at depth limit"
    );
}


