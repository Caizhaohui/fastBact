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

use std::collections::BTreeSet;

use fastbact::common::bitseq::reverse_complement_ascii;
use fastbact::common::sketch::splitmix64;
use fastbact::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, UnitigNode};
use fastbact::polish::output::curate_contigs;
use fastbact::resolve::exspander::{
    resolve_edge_between, resolve_repeats_exspander_with_stats, score_candidate_path,
    EdgeResolveResult, ExspanderParams,
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
    resolve_repeats_exspander_with_stats(graph, matrix, dist, params());
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
    for offset in 0..n {
        matrix.add_evidence(
            from,
            to,
            ReadPairEvidence {
                pair_id: id_base + offset,
                dist_u1_exit: dist1,
                dist_u2_entry: dist2,
                r1_len: 150,
                r2_len: 150,
            },
        );
    }
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
    );
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
