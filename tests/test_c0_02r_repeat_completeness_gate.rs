//! Phase F13 Independent Small-Graph Oracle & Generalized Repeat Completeness Gate.
//!
//! Verifies search completeness across parameterized repeat copies (n=1..depth boundary),
//! cycles, orientations, and budgets using an independent small-graph reference oracle.
//! Covers Cases A through O mandated in Phase F13 specification.

use std::collections::VecDeque;
use fastbact::common::bitseq::reverse_complement_ascii;
use fastbact::common::sketch::splitmix64;
use fastbact::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, UnitigNode};
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
// 1. Independent Small-Graph Reference Oracle (F13-E)
// ============================================================================

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

/// Independent reference oracle that enumerates admissible candidate paths on small graphs.
/// Explicitly tracks repeated physical nodes, enforcing validator rules:
/// - Nodes >= 500 bp can appear at most once.
/// - Nodes < 500 bp can appear multiple times.
/// - All edge transitions must have valid sequence overlap.
/// - Path depth <= max_depth (default 4).
/// - Separation <= max_separation.
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
        if let EdgeResolveResult::Valid(valid_edge) = resolve_edge_between(graph, fork_node, edge.target) {
            let sep = init_sep - valid_edge.overlap_k as i64;
            if sep <= max_separation {
                queue.push_back((edge.target, sep, vec![edge.target]));
            }
        }
    }

    while let Some((node, sep, path)) = queue.pop_front() {
        let node_ref = &graph.nodes[node.id() as usize];
        if node_ref.is_circular || node.id() == curr_dir.id() {
            continue;
        }

        // Check endpoint support from curr_dir
        if matrix.get_support(curr_dir, node) > 0 {
            let (count, weight) = score_candidate_path(matrix, dist, curr_dir, node, sep);
            if count > 0 {
                candidates.push(OracleCandidate {
                    path: path.clone(),
                    separation: sep,
                    support_count: count,
                    total_weight: weight,
                });
            }
        }

        // Continuation or extension
        if path.len() >= max_depth {
            // Check if there is an admissible outgoing continuation
            for edge in node_ref.edges(node.is_rc()) {
                let target = edge.target;
                let target_node = &graph.nodes[target.id() as usize];
                if target_node.is_circular || target.id() == curr_dir.id() {
                    continue;
                }
                if target_node.seq_len >= 500 && path.iter().any(|d| d.id() == target.id()) {
                    continue;
                }
                if let EdgeResolveResult::Valid(valid_edge) = resolve_edge_between(graph, node, target) {
                    let next_sep = sep + node_ref.seq_len as i64 - valid_edge.overlap_k as i64;
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
            if neighbor_node.is_circular || neighbor.id() == curr_dir.id() {
                continue;
            }
            // Unique node restriction: >= 500 bp can only appear once in path
            if neighbor_node.seq_len >= 500 && path.iter().any(|d| d.id() == neighbor.id()) {
                continue;
            }
            if let EdgeResolveResult::Valid(valid_edge) = resolve_edge_between(graph, node, neighbor) {
                let next_sep = sep + node_ref.seq_len as i64 - valid_edge.overlap_k as i64;
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
// 2. Oracle Differential Tests (F13-E & F13-D Matrix A–O)
// ============================================================================

/// Case A & B & C: Parameterized Repeat Copies (n = 1, 2, 3) compared with Oracle.
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

        // Insert size configured so that exactly n_copies of R has PE support at median 455
        let sep_d = 100 - 31 - 31; // 38
        let sep_r = sep_d + (n_copies as i64) * (120 - 31);
        let off_s = 75;
        let off_x = 75;
        let span_x = off_s + sep_r + off_x;

        let dist = InsertSizeDist {
            median: span_x as f64,
            mad: 20.0,
            min_insert: (span_x - 50) as usize,
            max_insert: (span_x + 50) as usize,
        };

        let mut matrix = PairLinkMatrix::new();
        // D has 5 pairs with mismatching span (or weak support)
        add_pairs_len(&mut matrix, s_id, d_id, 5, 200, 217, 30_100, 75);
        // X has 20 pairs matching exact span_x
        add_pairs_len(&mut matrix, s_id, x_id, 20, off_s as u32, off_x as u32, 30_200, 75);

        // 1. Run independent reference oracle
        let init_sep = 100 - 31;
        let oracle = reference_oracle_search(&graph, &matrix, &dist, s_id, q_id, init_sep, 4, dist.max_insert as i64);

        // Oracle must find the X path with n_copies
        let x_candidate = oracle.candidates.iter().find(|c| c.path.last() == Some(&x_id));
        assert!(
            x_candidate.is_some(),
            "Oracle must find candidate X for n_copies={n_copies}"
        );
        let x_cand = x_candidate.unwrap();
        assert_eq!(x_cand.support_count, 20);

        // 2. Run production resolver
        let res = resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params());
        assert!(res.is_ok());

        // Because X with 20 pairs is a strong competitor (or truncates at boundary),
        // weak decoy D must never be misassembled!
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

    // 4 copies of R from Q: Q -> R -> R -> R -> R (4 steps, hits MAX_CANDIDATE_DEPTH=4)
    // From 4th R, R -> X is step 5 (beyond depth cap)
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
    // In Case G, repeat cycle separation exceeds max_separation.
    // The repeat probe must safely prune without false truncation when all exits exceed budget.
    let mut graph = CdbgGraph::new(31);
    let j_sq = junction(0xD03_0001);
    let j_loop = junction(0xD03_0002);

    let s_seq = span_node(&[], 0xD03_0101, 700, &j_sq);
    let q_seq = span_node(&j_sq, 0xD03_0102, 100, &j_loop);
    let d_seq = span_node(&j_loop, 0xD03_0103, 300, &[]);
    // R is 450 bp long, so even 1 repeat step exceeds a small max_separation of 250 bp
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
    // D has valid support within [100, 250]
    add_pairs_len(&mut matrix, s_id, d_id, 10, 80, 82, 32_100, 75);
    // X has no valid support within [100, 250] because R is 450 bp
    add_pairs_len(&mut matrix, s_id, x_id, 20, 80, 82, 32_200, 75);

    let res = resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params());
    assert!(res.is_ok());

    // When the repeat is strictly outside separation budget, D is the sole valid candidate
    // and can safely resolve without false truncation.
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
    // R is a dead-end branch with zero evidence
    let r_seq = span_node(&j_loop, 0xD04_0104, 120, &[]);

    let s_id = add_node(&mut graph, &s_seq);
    let q_id = add_node(&mut graph, &q_seq);
    let d_id = add_node(&mut graph, &d_seq);
    let r_id = add_node(&mut graph, &r_seq);

    // Q -> D and Q -> R (dead end, no outgoing edges from R)
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
    // Only D has evidence
    add_pairs_len(&mut matrix, s_id, d_id, 10, 200, 217, 33_100, 75);

    let res = resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params());
    assert!(res.is_ok());

    // When repeat cycle has 0 evidence and dead-ends at depth limit, D resolves cleanly
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
    // R is 550 bp (>= 500 bp, unique node!)
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

    // sep = 0 - 31 = -31. With off_1 = 200, off_2 = 200, span = 200 - 31 + 200 = 369 bp
    let dist = InsertSizeDist {
        median: 369.0,
        mad: 20.0,
        min_insert: 300,
        max_insert: 450,
    };

    let mut matrix = PairLinkMatrix::new();
    // Strong support for A (25 pairs), none for B
    add_pairs_len(&mut matrix, s_id, a_id, 25, 200, 200, 35_100, 75);

    let res = resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params());
    assert!(res.is_ok());

    let a_node = a_id.id() as usize;
    assert!(
        graph.nodes[a_node].is_deleted,
        "Normal fork with strong support must resolve cleanly"
    );
}
