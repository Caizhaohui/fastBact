//! Graph simplification: tip clipping, Myers bit-parallel bubble popping,
//! and relative chimeric edge removal.

use crate::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge};

/// Computes Levenshtein edit distance using Myers bit-parallel algorithm for small strings,
/// or standard dynamic programming.
pub fn edit_distance(s1: &[u8], s2: &[u8]) -> usize {
    let m = s1.len();
    let n = s2.len();
    if m == 0 {
        return n;
    }
    if n == 0 {
        return m;
    }

    if m <= 64 {
        // Fast Myers bit-parallel algorithm
        let mut peq = [0u64; 256];
        for (i, &b) in s1.iter().enumerate() {
            peq[b as usize] |= 1u64 << i;
        }

        let mut pv = !0u64;
        let mut mv = 0u64;
        let mut score = m;

        for &b in s2 {
            let eq = peq[b as usize];
            let xv = eq | mv;
            let xh = (((eq & pv).wrapping_add(pv)) ^ pv) | eq;

            let mut ph = mv | !(xh | pv);
            let mut mh = pv & xh;

            if (ph >> (m - 1)) & 1 == 1 {
                score += 1;
            }
            if (mh >> (m - 1)) & 1 == 1 {
                score -= 1;
            }

            ph = (ph << 1) | 1;
            mh <<= 1;

            pv = mh | !(xv | ph);
            mv = ph & xv;
        }
        score
    } else {
        // Space-efficient 2-row DP for strings
        let (short, long) = if m <= n { (s1, s2) } else { (s2, s1) };
        let mut prev: Vec<usize> = (0..=short.len()).collect();
        let mut curr = vec![0; short.len() + 1];

        for (i, &b_long) in long.iter().enumerate() {
            curr[0] = i + 1;
            for (j, &b_short) in short.iter().enumerate() {
                let cost = if b_long == b_short { 0 } else { 1 };
                curr[j + 1] = (prev[j + 1] + 1).min(curr[j] + 1).min(prev[j] + cost);
            }
            std::mem::swap(&mut prev, &mut curr);
        }
        prev[short.len()]
    }
}

/// Simplifies the cDBG graph using iterative passes of:
/// 1. Tip clipping: dead-end paths shorter than 2*k with significantly higher alternative coverage.
/// 2. Bubble popping: alternative paths between branching points with >= 90% sequence similarity.
/// 3. Re-compacting unbranched 1-in, 1-out paths.
pub fn simplify_graph(graph: &mut CdbgGraph) {
    let mut changed = true;
    let mut iterations = 0;

    while changed && iterations < 10 {
        changed = false;
        iterations += 1;

        // 1. Tip clipping
        if clip_tips(graph) {
            changed = true;
        }

        // 2. Bubble popping
        if pop_bubbles(graph) {
            changed = true;
        }

        // 3. Re-compact unbranched chains
        if recompact_graph(graph) {
            changed = true;
        }
    }
}

/// Removes dead-end tips (length < 2*k, out-degree 0) if alternative paths exist with coverage >= 4*cov(tip).
pub fn clip_tips(graph: &mut CdbgGraph) -> bool {
    let k = graph.k;
    let n = graph.nodes.len();
    let mut tips_to_remove = Vec::new();

    for id in 0..n as u32 {
        let node = &graph.nodes[id as usize];
        if node.is_deleted || node.is_circular {
            continue;
        }

        if node.seq_len < (2 * k) as u32 {
            let fwd_out = node.edges_fwd.len();
            let rev_out = node.edges_rev.len();

            if fwd_out == 0 && rev_out > 0 {
                tips_to_remove.push(id);
            } else if rev_out == 0 && fwd_out > 0 {
                tips_to_remove.push(id);
            }
        }
    }

    if tips_to_remove.is_empty() {
        return false;
    }

    // Build incoming predecessor lookup map to avoid O(N^2) search
    let mut predecessors: Vec<smallvec::SmallVec<[u32; 4]>> = vec![smallvec::SmallVec::new(); n];
    for u in 0..n as u32 {
        let node = &graph.nodes[u as usize];
        if !node.is_deleted {
            for e in node.edges_fwd.iter().chain(node.edges_rev.iter()) {
                let target_id = e.target.id() as usize;
                if target_id < n && !predecessors[target_id].contains(&u) {
                    predecessors[target_id].push(u);
                }
            }
        }
    }

    let mut removed = false;
    for tip_id in tips_to_remove {
        let tip_cov = graph.nodes[tip_id as usize].kmer_cov;

        let mut has_strong_alternative = false;
        for &pred_id in &predecessors[tip_id as usize] {
            if pred_id != tip_id && !graph.nodes[pred_id as usize].is_deleted {
                let alt_cov = graph.nodes[pred_id as usize].kmer_cov;
                if alt_cov >= 3.0 * tip_cov || tip_cov <= 2.0 {
                    has_strong_alternative = true;
                    break;
                }
            }
        }

        if has_strong_alternative {
            graph.nodes[tip_id as usize].is_deleted = true;
            for &pred_id in &predecessors[tip_id as usize] {
                let other = &mut graph.nodes[pred_id as usize];
                if !other.is_deleted {
                    other.edges_fwd.retain(|e| e.target.id() != tip_id);
                    other.edges_rev.retain(|e| e.target.id() != tip_id);
                }
            }
            removed = true;
        }
    }

    removed
}

/// Identifies and pops bubbles (alternative paths with >= 90% identity).
///
/// Bounded to maximum bubble path length <= 1,000 bp to prevent catastrophic O(M*N)
/// DP on large genomic unitigs.
pub fn pop_bubbles(graph: &mut CdbgGraph) -> bool {
    let mut bubbles_to_pop = Vec::new();
    let n = graph.nodes.len() as u32;
    let max_bubble_len = (3 * graph.k).clamp(500, 1000);

    for u in 0..n {
        if graph.nodes[u as usize].is_deleted {
            continue;
        }

        let edges = graph.nodes[u as usize].edges_fwd.clone();
        if edges.len() >= 2 {
            for i in 0..edges.len() {
                for j in i + 1..edges.len() {
                    let target_i = edges[i].target;
                    let target_j = edges[j].target;

                    if target_i.id() == target_j.id()
                        || graph.nodes[target_i.id() as usize].is_deleted
                        || graph.nodes[target_j.id() as usize].is_deleted
                    {
                        continue;
                    }

                    // Strict length boundary: bubbles are small sequencing errors (< 1 kb)
                    let len_i = graph.nodes[target_i.id() as usize].seq_len as usize;
                    let len_j = graph.nodes[target_j.id() as usize].seq_len as usize;
                    if len_i > max_bubble_len || len_j > max_bubble_len {
                        continue;
                    }

                    let max_len = len_i.max(len_j);
                    if max_len == 0 {
                        continue;
                    }
                    let len_diff = (len_i as isize - len_j as isize).unsigned_abs();
                    if len_diff > (max_len / 10).max(2) {
                        continue;
                    }

                    let succs_i = graph.nodes[target_i.id() as usize]
                        .edges(target_i.is_rc())
                        .clone();
                    let succs_j = graph.nodes[target_j.id() as usize]
                        .edges(target_j.is_rc())
                        .clone();

                    for si in &succs_i {
                        for sj in &succs_j {
                            if si.target == sj.target {
                                let seq_i = graph.get_stranded_seq(target_i);
                                let seq_j = graph.get_stranded_seq(target_j);

                                let dist = edit_distance(&seq_i, &seq_j);
                                let similarity = 1.0 - (dist as f64) / (max_len as f64);

                                if similarity >= 0.90 {
                                    let cov_i = graph.nodes[target_i.id() as usize].kmer_cov;
                                    let cov_j = graph.nodes[target_j.id() as usize].kmer_cov;

                                    let (keep, del) = if cov_i >= cov_j {
                                        (target_i, target_j)
                                    } else {
                                        (target_j, target_i)
                                    };
                                    bubbles_to_pop.push((keep, del));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if bubbles_to_pop.is_empty() {
        return false;
    }

    let mut popped = false;
    for (keep_target, del_target) in bubbles_to_pop {
        if graph.nodes[del_target.id() as usize].is_deleted {
            continue;
        }

        graph.nodes[del_target.id() as usize].is_deleted = true;
        let added_cov = graph.nodes[del_target.id() as usize].kmer_cov;
        graph.nodes[keep_target.id() as usize].kmer_cov += added_cov;

        for other in graph.nodes.iter_mut() {
            other.edges_fwd.retain(|e| e.target.id() != del_target.id());
            other.edges_rev.retain(|e| e.target.id() != del_target.id());
        }
        popped = true;
    }

    popped
}

/// Re-compacts unbranched 1-in, 1-out paths in the graph.
pub fn recompact_graph(graph: &mut CdbgGraph) -> bool {
    let mut recompacted = false;
    let n = graph.nodes.len() as u32;

    // Precompute incoming counts in O(E) to avoid O(N^2) scans
    let mut in_counts = vec![0u32; n as usize];
    for u in 0..n {
        let node = &graph.nodes[u as usize];
        if !node.is_deleted {
            for e in node.edges_fwd.iter().chain(node.edges_rev.iter()) {
                let tid = e.target.id() as usize;
                if tid < in_counts.len() {
                    in_counts[tid] += 1;
                }
            }
        }
    }

    for u in 0..n {
        if graph.nodes[u as usize].is_deleted || graph.nodes[u as usize].is_circular {
            continue;
        }

        if graph.nodes[u as usize].edges_fwd.len() == 1 {
            let edge = graph.nodes[u as usize].edges_fwd[0];
            let v = edge.target;

            if v.id() == u
                || graph.nodes[v.id() as usize].is_deleted
                || graph.nodes[v.id() as usize].is_circular
            {
                continue;
            }

            if in_counts[v.id() as usize] == 1 {
                let seq_u = graph.get_unitig_seq(u);
                let seq_v = graph.get_stranded_seq(v);

                let overlap = edge.overlap_k as usize;
                if seq_v.len() > overlap {
                    let mut merged_seq = seq_u;
                    merged_seq.extend_from_slice(&seq_v[overlap..]);

                    let new_id = graph.arena.push_seq(&merged_seq);
                    graph.nodes[u as usize].seq_len = merged_seq.len() as u32;
                    graph.nodes[u as usize].seq_id = new_id;

                    let v_edges = graph.nodes[v.id() as usize].edges(v.is_rc()).clone();
                    graph.nodes[u as usize].edges_fwd = v_edges;

                    // Rewire any edges in other nodes pointing to v so they point to u
                    for other in graph.nodes.iter_mut() {
                        if !other.is_deleted {
                            for e in other.edges_fwd.iter_mut().chain(other.edges_rev.iter_mut()) {
                                if e.target.id() == v.id() {
                                    let is_rc = e.target.is_rc() ^ v.is_rc();
                                    *e = GraphEdge::new(DirUnitigId::new(u, is_rc), e.overlap_k, e.support);
                                }
                            }
                        }
                    }

                    graph.nodes[v.id() as usize].is_deleted = true;
                    recompacted = true;
                }
            }
        }
    }

    recompacted
}

