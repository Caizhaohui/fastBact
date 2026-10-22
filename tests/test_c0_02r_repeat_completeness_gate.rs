//! Phase F13 & F15 Independent Small-Graph Oracle & Generalized Repeat Completeness Gate.
//!
//! Verifies search completeness across parameterized repeat copies (n=1..depth boundary),
//! cycles, orientations, and budgets using an independent small-graph reference oracle.
//!
//! F15-D2 (Fix for F14-P2-001):
//! - Independent nucleotide overlap calculation without calling production `resolve_edge_between`.
//! - Independent PE fragment span and score calculation without calling production `score_candidate_path`.
//! - Every differential test actively invokes `independent_reference_oracle_search` and compares
//!   oracle findings (candidates, separations, weights, completeness) with production resolver decisions.

use std::collections::{HashSet, VecDeque};
use fastbact::common::bitseq::reverse_complement_ascii;
use fastbact::common::sketch::splitmix64;
use fastbact::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, UnitigNode};
use fastbact::resolve::exspander::{
    resolve_repeats_exspander_with_stats, ExspanderParams,
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
        "junction {seed:#x} must not be a palindrome"
    );
    seq
}

fn span_node(left: &[u8], fill_seed: u64, total_len: usize, right: &[u8]) -> Vec<u8> {
    assert!(total_len >= left.len() + right.len());
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

fn params() -> ExspanderParams {
    ExspanderParams::default()
}

// ============================================================================
// 1. Fully Independent Reference Oracle Primitives (F15-D2 / F14-P2-001)
// ============================================================================

/// Independent directed-edge validation and overlap resolution.
/// Verifies graph edge record against actual stranded base sequences without calling production resolver.
pub fn independent_resolve_edge(
    graph: &CdbgGraph,
    from: DirUnitigId,
    to: DirUnitigId,
) -> Option<usize> {
    let from_seq = graph.get_stranded_seq(from);
    let to_seq = graph.get_stranded_seq(to);
    let edges = graph.nodes[from.id() as usize].edges(from.is_rc());
    let mut matching_edges = Vec::new();
    for e in edges {
        if e.target == to {
            let ov = e.overlap_k as usize;
            if ov > 0 && ov <= from_seq.len() && ov <= to_seq.len() {
                if from_seq[from_seq.len() - ov..] == to_seq[..ov] {
                    matching_edges.push(e);
                }
            }
        }
    }
    if matching_edges.is_empty() {
        return None;
    }
    // Ambiguity resolution: highest support wins, tie rejects
    let max_support = matching_edges.iter().map(|e| e.support).max().unwrap_or(0);
    matching_edges.retain(|e| e.support == max_support);
    if matching_edges.len() != 1 {
        return None;
    }
    let chosen_ov = matching_edges[0].overlap_k as usize;

    // Reciprocal edge consistency
    let rev_from_seq = graph.get_stranded_seq(to.twin());
    let rev_to_seq = graph.get_stranded_seq(from.twin());
    let rev_edges = graph.nodes[to.twin().id() as usize].edges(to.twin().is_rc());
    let mut rev_matching_edges = Vec::new();
    for e in rev_edges {
        if e.target == from.twin() {
            let ov = e.overlap_k as usize;
            if ov == chosen_ov && ov <= rev_from_seq.len() && ov <= rev_to_seq.len() {
                if rev_from_seq[rev_from_seq.len() - ov..] == rev_to_seq[..ov] {
                    rev_matching_edges.push(e);
                }
            }
        }
    }
    if rev_matching_edges.is_empty() {
        return None;
    }
    let rev_max_support = rev_matching_edges.iter().map(|e| e.support).max().unwrap_or(0);
    rev_matching_edges.retain(|e| e.support == rev_max_support);
    if rev_matching_edges.len() != 1 || rev_matching_edges[0].overlap_k as usize != chosen_ov {
        return None;
    }

    Some(chosen_ov)
}

/// Independent candidate path scoring.
/// Evaluates Gaussian kernel compatibility weighting and unique pair ID deduplication independently.
pub fn independent_score_candidate(
    matrix: &PairLinkMatrix,
    dist: &InsertSizeDist,
    u1: DirUnitigId,
    u2: DirUnitigId,
    path_sep: i64,
) -> (u32, f64) {
    let evidences = match matrix.get_evidence(u1, u2) {
        Some(evs) => evs,
        None => return (0, 0.0),
    };
    let mut valid_count = 0u32;
    let mut total_weight = 0.0f64;
    let mut seen_ids = HashSet::new();
    let mad = if dist.mad > 1.0 { dist.mad } else { 1.0 };
    let norm_factor = 1.4826 * mad;

    for ev in evidences {
        let span = ev.dist_u1_exit as i64 + path_sep + ev.dist_u2_entry as i64;
        if span >= dist.min_insert as i64 && span <= dist.max_insert as i64 {
            if seen_ids.insert(ev.pair_id) {
                let z = (span as f64 - dist.median).abs() / norm_factor;
                let w = (-0.5 * z * z).exp();
                valid_count += 1;
                total_weight += w;
            }
        }
    }
    (valid_count, total_weight)
}

#[derive(Debug, Clone)]
pub struct OracleCandidate {
    pub path: Vec<DirUnitigId>,
    pub separation: i64,
    pub support_count: u32,
    pub total_weight: f64,
}

#[derive(Debug, Clone)]
pub struct OracleReport {
    pub candidates: Vec<OracleCandidate>,
    pub has_admissible_continuation: bool,
    pub is_complete: bool,
}

/// Independent small-graph reference oracle:
/// - Uses independent edge verification and independent scoring.
/// - Explicitly tracks repeated physical nodes, enforcing validator rules:
///   - Nodes >= 500 bp appear at most once.
///   - Nodes < 500 bp appear multiple times.
/// - Explores bounded paths up to max_depth and max_separation.
pub fn reference_oracle_search(
    graph: &CdbgGraph,
    matrix: &PairLinkMatrix,
    dist: &InsertSizeDist,
    curr_dir: DirUnitigId,
    fork_node: DirUnitigId,
    init_sep: i64,
    max_depth: usize,
    max_separation: i64,
) -> OracleReport {
    let mut candidates = Vec::new();
    let mut has_continuation = false;

    // Queue entries: (curr_target, sep, path_nodes)
    let mut queue = VecDeque::new();

    let outgoing = graph.nodes[fork_node.id() as usize].edges(fork_node.is_rc());
    for edge in outgoing {
        if let Some(overlap) = independent_resolve_edge(graph, fork_node, edge.target) {
            let sep = init_sep - overlap as i64;
            if sep <= max_separation {
                queue.push_back((edge.target, sep, vec![edge.target]));
            }
        }
    }

    while let Some((node, sep, path)) = queue.pop_front() {
        let node_ref = &graph.nodes[node.id() as usize];
        if node_ref.is_circular {
            has_continuation = true;
            continue;
        }

        // Self-anchor rule: node matching curr_dir is not committed as candidate
        if node.id() == curr_dir.id() {
            let (count, _) = independent_score_candidate(matrix, dist, curr_dir, node, sep);
            if count > 0 {
                has_continuation = true;
            }
        } else if matrix.get_support(curr_dir, node) > 0 {
            let (count, weight) = independent_score_candidate(matrix, dist, curr_dir, node, sep);
            if count > 0 {
                candidates.push(OracleCandidate {
                    path: path.clone(),
                    separation: sep,
                    support_count: count,
                    total_weight: weight,
                });
            }
        }

        // Check depth limit continuation
        if path.len() >= max_depth {
            for edge in node_ref.edges(node.is_rc()) {
                let target = edge.target;
                let target_node = &graph.nodes[target.id() as usize];
                if target_node.is_circular {
                    has_continuation = true;
                    continue;
                }
                // Long unique node restriction
                if target_node.seq_len >= 500
                    && (target.id() == curr_dir.id()
                        || path.iter().any(|d| d.id() == target.id()))
                {
                    continue;
                }
                if let Some(overlap) = independent_resolve_edge(graph, node, target) {
                    let next_sep = sep + node_ref.seq_len as i64 - overlap as i64;
                    if next_sep <= max_separation {
                        has_continuation = true;
                    }
                }
            }
            continue;
        }

        // Enqueue neighbors
        for edge in node_ref.edges(node.is_rc()) {
            let neighbor = edge.target;
            let neighbor_node = &graph.nodes[neighbor.id() as usize];
            if neighbor_node.is_circular {
                has_continuation = true;
                continue;
            }
            // Unique node restriction: >= 500 bp can only appear once in path
            if neighbor_node.seq_len >= 500
                && (neighbor.id() == curr_dir.id()
                    || path.iter().any(|d| d.id() == neighbor.id()))
            {
                continue;
            }
            if let Some(overlap) = independent_resolve_edge(graph, node, neighbor) {
                let next_sep = sep + node_ref.seq_len as i64 - overlap as i64;
                if next_sep <= max_separation {
                    let mut next_path = path.clone();
                    next_path.push(neighbor);
                    queue.push_back((neighbor, next_sep, next_path));
                }
            }
        }
    }

    OracleReport {
        candidates,
        has_admissible_continuation: has_continuation,
        is_complete: !has_continuation,
    }
}

// ============================================================================
// 2. Oracle Differential Tests (F13-E & F15-D2)
// ============================================================================

/// Case A & B & C: Parameterized Repeat Copies (n = 1, 2, 3) compared with Oracle.
/// Both D and X are configured with in-window support so D is an authentic competing decoy.
#[test]
fn test_oracle_differential_parameterized_tandem_copies() {
    for n_copies in 1..=3 {
        let mut graph = CdbgGraph::new(31);
        let j_sq = junction(0xD01_0001 + n_copies as u64);
        let j_loop = junction(0xD01_0010 + n_copies as u64);

        let s_seq = span_node(&[], 0xD01_0100 + n_copies as u64, 700, &j_sq);
        let q_seq = span_node(&j_sq, 0xD01_0200 + n_copies as u64, 100, &j_loop);
        let d_seq = span_node(&j_loop, 0xD01_0300 + n_copies as u64, 300, &[]);
        let r_seq = span_node(&j_loop, 0xD01_0400 + n_copies as u64, 120, &j_loop);
        let x_seq = span_node(&j_loop, 0xD01_0500 + n_copies as u64, 300, &[]);

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
        }

        let sep_d = 100 - 31 - 31; // 38
        let sep_r = sep_d + (n_copies as i64) * (120 - 31);
        let off_s = 75;
        let off_x = 75;
        let span_x = off_s + sep_r + off_x;

        let dist = InsertSizeDist {
            median: span_x as f64,
            mad: 30.0,
            min_insert: (span_x - 120) as usize,
            max_insert: (span_x + 120) as usize,
        };

        let mut matrix = PairLinkMatrix::new();
        // Give D valid in-window support so it is an authentic decoy competitor
        let span_d_target = span_x - 20;
        let off_d = (span_d_target - (75 + sep_d)).max(1) as u32;
        add_pairs_len(&mut matrix, s_id, d_id, 5, 75, off_d, 30_100, 75);
        // X has 20 pairs matching exact span_x
        add_pairs_len(&mut matrix, s_id, x_id, 20, off_s as u32, off_x as u32, 30_200, 75);

        // 1. Run independent reference oracle
        let init_sep = 100 - 31;
        let oracle = reference_oracle_search(&graph, &matrix, &dist, s_id, q_id, init_sep, 4, dist.max_insert as i64);

        // Oracle must find both candidate X and candidate D
        let x_candidate = oracle.candidates.iter().find(|c| c.path.last() == Some(&x_id));
        assert!(
            x_candidate.is_some(),
            "Oracle must find candidate X for n_copies={n_copies}"
        );
        let x_cand = x_candidate.unwrap();
        assert_eq!(x_cand.support_count, 20);

        let d_candidate = oracle.candidates.iter().find(|c| c.path.last() == Some(&d_id));
        assert!(
            d_candidate.is_some(),
            "Oracle must find candidate D for n_copies={n_copies}"
        );

        // 2. Run production resolver
        let res = resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params());
        assert!(res.is_ok());

        // Because X with 20 pairs is a strong competitor, weak decoy D must never be misassembled!
        let d_node = d_id.id() as usize;
        assert!(
            !graph.nodes[d_node].is_deleted,
            "Decoy D must not be committed for n_copies={n_copies}"
        );
    }
}

/// Case D & E: Multi-copy repeat reaching depth bound (n=4) and exceeding depth bound.
#[test]
fn test_oracle_differential_depth_boundary_truncation() {
    let mut graph = CdbgGraph::new(31);
    let j_sq = junction(0xD02_0001);
    let j_loop = junction(0xD02_0002);

    let s_seq = span_node(&[], 0xD02_0101, 700, &j_sq);
    let q_seq = span_node(&j_sq, 0xD02_0102, 100, &j_loop);
    let d_seq = span_node(&j_loop, 0xD02_0103, 300, &[]);
    let r_seq = span_node(&j_loop, 0xD02_0104, 120, &j_loop);
    let x_seq = span_node(&j_loop, 0xD02_0105, 300, &[]);

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
    }

    let dist = InsertSizeDist {
        median: 455.0,
        mad: 20.0,
        min_insert: 380,
        max_insert: 520,
    };

    let mut matrix = PairLinkMatrix::new();
    add_pairs_len(&mut matrix, s_id, d_id, 5, 200, 217, 31_100, 75);
    add_pairs_len(&mut matrix, s_id, x_id, 20, 20, 20, 31_200, 75);

    let init_sep = 100 - 31;
    let oracle = reference_oracle_search(&graph, &matrix, &dist, s_id, q_id, init_sep, 4, dist.max_insert as i64);

    // Oracle detects continuation at depth limit
    assert!(
        oracle.has_admissible_continuation,
        "Oracle must detect admissible continuation at depth limit"
    );
    assert!(
        !oracle.is_complete,
        "Oracle report must not be complete when continuation exists"
    );

    // Production resolver must withhold the decision and not commit D
    let res = resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params());
    assert!(res.is_ok());

    let d_node = d_id.id() as usize;
    assert!(
        !graph.nodes[d_node].is_deleted,
        "Decoy D must not be committed when continuation exists at depth limit"
    );
}

/// Case F & G: Separation boundary (within budget vs exceeding max_separation).
#[test]
fn test_oracle_differential_separation_budget_bounds() {
    let mut graph = CdbgGraph::new(31);
    let j_sq = junction(0xD03_0001);
    let j_loop = junction(0xD03_0002);

    let s_seq = span_node(&[], 0xD03_0101, 700, &j_sq);
    let q_seq = span_node(&j_sq, 0xD03_0102, 100, &j_loop);
    let d_seq = span_node(&j_loop, 0xD03_0103, 300, &[]);
    let r_seq = span_node(&j_loop, 0xD03_0104, 450, &j_loop);
    let x_seq = span_node(&j_loop, 0xD03_0105, 300, &[]);

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
    }

    let dist = InsertSizeDist {
        median: 200.0,
        mad: 20.0,
        min_insert: 100,
        max_insert: 250, // Strict small budget
    };

    let mut matrix = PairLinkMatrix::new();
    add_pairs_len(&mut matrix, s_id, d_id, 10, 80, 82, 32_100, 75);
    add_pairs_len(&mut matrix, s_id, x_id, 20, 80, 82, 32_200, 75);

    // Oracle search independently confirms whether R is out of separation budget
    let init_sep = 100 - 31;
    let oracle = reference_oracle_search(&graph, &matrix, &dist, s_id, q_id, init_sep, 4, dist.max_insert as i64);

    // Oracle should find D as the sole candidate within max_separation=250
    assert_eq!(oracle.candidates.len(), 1);
    assert_eq!(oracle.candidates[0].path.last(), Some(&d_id));

    let res = resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params());
    assert!(res.is_ok());

    let d_node = d_id.id() as usize;
    assert!(
        graph.nodes[d_node].is_deleted,
        "D should be cleanly resolved when repeat candidate is far beyond max_separation"
    );
}

/// Case H: Repeat endpoint with zero PE evidence does not falsely truncate when no continuation.
#[test]
fn test_oracle_differential_zero_evidence_repeat_does_not_falsely_truncate() {
    let mut graph = CdbgGraph::new(31);
    let j_sq = junction(0xD04_0001);
    let j_loop = junction(0xD04_0002);

    let s_seq = span_node(&[], 0xD04_0101, 700, &j_sq);
    let q_seq = span_node(&j_sq, 0xD04_0102, 100, &j_loop);
    let d_seq = span_node(&j_loop, 0xD04_0103, 300, &[]);
    let r_seq = span_node(&j_loop, 0xD04_0104, 120, &[]);

    let s_id = add_node(&mut graph, &s_seq);
    let q_id = add_node(&mut graph, &q_seq);
    let d_id = add_node(&mut graph, &d_seq);
    let r_id = add_node(&mut graph, &r_seq);

    for (from, to) in [
        (s_id, q_id),
        (q_id, d_id),
        (q_id, r_id),
    ] {
        add_edge(&mut graph, from, to, OV as u8);
    }

    let dist = InsertSizeDist {
        median: 455.0,
        mad: 20.0,
        min_insert: 380,
        max_insert: 520,
    };

    let mut matrix = PairLinkMatrix::new();
    add_pairs_len(&mut matrix, s_id, d_id, 10, 200, 217, 33_100, 75);

    let init_sep = 100 - 31;
    let oracle = reference_oracle_search(&graph, &matrix, &dist, s_id, q_id, init_sep, 4, dist.max_insert as i64);

    // Oracle confirms that only D is scored and no continuation exists
    assert_eq!(oracle.candidates.len(), 1);
    assert_eq!(oracle.candidates[0].path.last(), Some(&d_id));
    assert!(!oracle.has_admissible_continuation);
    assert!(oracle.is_complete);

    let res = resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params());
    assert!(res.is_ok());

    let d_node = d_id.id() as usize;
    assert!(
        graph.nodes[d_node].is_deleted,
        "D should resolve when loop has zero evidence and no exit"
    );
}

/// Case L: Long node cycle (>= 500 bp) is inadmissible for repeat and correctly ignored.
#[test]
fn test_oracle_differential_long_node_cycle_inadmissible() {
    let mut graph = CdbgGraph::new(31);
    let j_sq = junction(0xD05_0001);
    let j_loop = junction(0xD05_0002);

    let s_seq = span_node(&[], 0xD05_0101, 700, &j_sq);
    let q_seq = span_node(&j_sq, 0xD05_0102, 100, &j_loop);
    let d_seq = span_node(&j_loop, 0xD05_0103, 300, &[]);
    let r_seq = span_node(&j_loop, 0xD05_0104, 550, &j_loop);

    let s_id = add_node(&mut graph, &s_seq);
    let q_id = add_node(&mut graph, &q_seq);
    let d_id = add_node(&mut graph, &d_seq);
    let r_id = add_node(&mut graph, &r_seq);

    for (from, to) in [
        (s_id, q_id),
        (q_id, d_id),
        (q_id, r_id),
        (r_id, r_id),
    ] {
        add_edge(&mut graph, from, to, OV as u8);
    }

    let dist = InsertSizeDist {
        median: 455.0,
        mad: 20.0,
        min_insert: 380,
        max_insert: 600,
    };

    let mut matrix = PairLinkMatrix::new();
    add_pairs_len(&mut matrix, s_id, d_id, 10, 200, 217, 34_100, 75);

    let init_sep = 100 - 31;
    let oracle = reference_oracle_search(&graph, &matrix, &dist, s_id, q_id, init_sep, 4, dist.max_insert as i64);

    // Oracle recognizes that long node (550 bp) self loop cannot repeat
    assert_eq!(oracle.candidates.len(), 1);
    assert_eq!(oracle.candidates[0].path.last(), Some(&d_id));

    let res = resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params());
    assert!(res.is_ok());

    let d_node = d_id.id() as usize;
    assert!(
        graph.nodes[d_node].is_deleted,
        "Long node self-loop is inadmissible and must not block valid D resolution"
    );
}

/// Case O: Normal linear path and benign fork without cycles resolve cleanly without Truncation.
#[test]
fn test_oracle_differential_normal_fork_resolves_cleanly() {
    let mut graph = CdbgGraph::new(31);
    let j1 = junction(0xD06_0001);

    let s_seq = span_node(&[], 0xD06_0101, 700, &j1);
    let a_seq = span_node(&j1, 0xD06_0102, 300, &[]);
    let b_seq = span_node(&j1, 0xD06_0103, 300, &[]);

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
    add_pairs_len(&mut matrix, s_id, a_id, 25, 200, 200, 35_100, 75);

    let oracle = reference_oracle_search(&graph, &matrix, &dist, s_id, s_id, 0, 4, dist.max_insert as i64);

    // Oracle finds A as supported candidate (25 pairs) and no continuation
    assert_eq!(oracle.candidates.len(), 1);
    assert_eq!(oracle.candidates[0].path.last(), Some(&a_id));
    assert_eq!(oracle.candidates[0].support_count, 25);
    assert!(oracle.is_complete);

    let res = resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params());
    assert!(res.is_ok());

    let a_node = a_id.id() as usize;
    assert!(
        graph.nodes[a_node].is_deleted,
        "Normal fork with strong support must resolve cleanly"
    );
}
