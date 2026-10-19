//! SPAdes-style exSPAnder paired-end repeat resolution engine.
//!
//! Traverses across repeats longer than k_max using paired-end read insert size constraints,
//! connecting unique flanking unitigs across repetitive regions.

use crate::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, UnitigNode};
use crate::resolve::mapper::{InsertSizeDist, PairLinkMatrix};
use hashbrown::HashSet;

/// Configuration parameters for exSPAnder repeat resolution.
#[derive(Clone, Copy, Debug)]
pub struct ExspanderParams {
    /// Minimum raw supporting read-pair links to extend path across repeat
    pub min_support: u32,
    /// Winner-take-all confidence ratio: score(best) >= decision_threshold * score(second)
    pub decision_threshold: f64,
}

impl Default for ExspanderParams {
    fn default() -> Self {
        Self {
            min_support: 5,
            decision_threshold: 2.5,
        }
    }
}

/// Statistics tracking junction operations and correctness decisions in exSPAnder.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct JunctionStats {
    /// Successfully spliced joins with exact k overlap
    pub valid_k_joins: usize,
    /// Successfully spliced joins with exact (k - 1) overlap
    pub valid_k_minus_1_joins: usize,
    /// Successfully spliced joins with other validated overlap
    pub valid_other_joins: usize,
    /// Joins rejected due to ambiguous edge selection (e.g. tied support)
    pub rejected_ambiguous_joins: usize,
    /// Joins rejected due to sequence mismatch or invalid overlap
    pub rejected_invalid_joins: usize,
}

/// Internal path step preserving directed unitig identity, chosen incoming edge,
/// validated overlap, and supporting evidence.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PathStep {
    /// Directed stranded unitig ID.
    pub dir: DirUnitigId,
    /// Edge chosen when entering this unitig from the preceding step (None for seed/start).
    pub incoming_edge: Option<GraphEdge>,
    /// Validated overlap length in bp between preceding unitig and this unitig.
    pub overlap: usize,
    /// Supporting evidence score (paired-end link support or transition count).
    pub support: u32,
}

/// Result of evaluating edge candidates between two directed unitigs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EdgeResolveResult {
    /// Exactly one unambiguous edge matches actual sequence suffix/prefix.
    Valid(GraphEdge),
    /// Multiple edges match sequence with tied support (insufficient evidence).
    Ambiguous,
    /// Edge(s) exist but sequence suffix/prefix does not match recorded overlap.
    InvalidSequenceOverlap,
    /// No edge exists in graph.
    NoEdge,
}

/// Resolves the directed edge from `curr_dir` to `target_dir` with explicit, auditable criteria:
/// 1. Finds all candidate edges from `curr_dir` to `target_dir`.
/// 2. Verifies whether each candidate edge's recorded overlap matches the actual stranded sequence suffix and prefix.
/// 3. If multiple edges match at sequence level (e.g. k and k-1 in homopolymers or tandem motifs):
///    - Compares edge transition support counts.
///    - If one edge has strictly higher support, selects it.
///    - If top edges have equal support, rejects as Ambiguous (no arbitrary choice based on list order).
pub fn resolve_edge_between(
    graph: &CdbgGraph,
    curr_dir: DirUnitigId,
    target_dir: DirUnitigId,
) -> EdgeResolveResult {
    let curr_node = &graph.nodes[curr_dir.id() as usize];
    let candidate_edges: Vec<&GraphEdge> = curr_node
        .edges(curr_dir.is_rc())
        .iter()
        .filter(|e| e.target == target_dir)
        .collect();

    if candidate_edges.is_empty() {
        return EdgeResolveResult::NoEdge;
    }

    let u_seq = graph.get_stranded_seq(curr_dir);
    let v_seq = graph.get_stranded_seq(target_dir);

    let valid_seq_edges: Vec<&GraphEdge> = candidate_edges
        .into_iter()
        .filter(|e| {
            let ov = e.overlap_k as usize;
            if ov == 0 || ov > u_seq.len() || ov > v_seq.len() {
                false
            } else {
                u_seq[u_seq.len() - ov..] == v_seq[..ov]
            }
        })
        .collect();

    if valid_seq_edges.is_empty() {
        return EdgeResolveResult::InvalidSequenceOverlap;
    }

    if valid_seq_edges.len() == 1 {
        return EdgeResolveResult::Valid(*valid_seq_edges[0]);
    }

    // Multiple edges match sequence: inspect read support
    let max_support = valid_seq_edges.iter().map(|e| e.support).max().unwrap_or(0);
    let best_edges: Vec<&GraphEdge> = valid_seq_edges
        .into_iter()
        .filter(|e| e.support == max_support)
        .collect();

    if best_edges.len() == 1 {
        EdgeResolveResult::Valid(*best_edges[0])
    } else {
        EdgeResolveResult::Ambiguous
    }
}

/// Extends path from `curr_dir` in one direction using paired-end read constraints
/// while recording the exact `PathStep` for every transition.
fn extend_one_direction(
    mut curr_dir: DirUnitigId,
    graph: &CdbgGraph,
    matrix: &PairLinkMatrix,
    dist: &InsertSizeDist,
    params: ExspanderParams,
    used_unique: &mut HashSet<u32>,
    stats: &mut JunctionStats,
    rejected_nodes: &mut HashSet<u32>,
) -> Vec<PathStep> {
    let mut ext = Vec::new();

    loop {
        let curr_node = &graph.nodes[curr_dir.id() as usize];
        let outgoing = curr_node.edges(curr_dir.is_rc()).clone();

        if outgoing.is_empty() {
            break;
        }

        // Collect distinct target nodes to prevent self-competition on duplicate overlap edges
        let mut unique_targets = Vec::new();
        for e in &outgoing {
            if !unique_targets.contains(&e.target) {
                unique_targets.push(e.target);
            }
        }

        let mut valid_continuations = Vec::new();
        for &target in &unique_targets {
            match resolve_edge_between(graph, curr_dir, target) {
                EdgeResolveResult::Valid(edge) => {
                    valid_continuations.push((target, edge));
                }
                EdgeResolveResult::Ambiguous => {
                    stats.rejected_ambiguous_joins += 1;
                    rejected_nodes.insert(target.id());
                }
                EdgeResolveResult::InvalidSequenceOverlap => {
                    stats.rejected_invalid_joins += 1;
                    rejected_nodes.insert(target.id());
                }
                EdgeResolveResult::NoEdge => {}
            }
        }

        if valid_continuations.is_empty() {
            break;
        }

        // Case A: 1 valid outgoing continuation
        if valid_continuations.len() == 1 {
            let (target, edge_to_target) = valid_continuations[0];
            let target_node = &graph.nodes[target.id() as usize];
            if target_node.is_deleted || target_node.is_circular {
                break;
            }

            let next_outgoing = target_node.edges(target.is_rc()).clone();
            let mut next_unique_targets = Vec::new();
            for e in &next_outgoing {
                if !next_unique_targets.contains(&e.target) {
                    next_unique_targets.push(e.target);
                }
            }

            let mut next_continuations = Vec::new();
            for &cand in &next_unique_targets {
                match resolve_edge_between(graph, target, cand) {
                    EdgeResolveResult::Valid(edge) => {
                        next_continuations.push((cand, edge));
                    }
                    EdgeResolveResult::Ambiguous => {
                        stats.rejected_ambiguous_joins += 1;
                        rejected_nodes.insert(cand.id());
                    }
                    EdgeResolveResult::InvalidSequenceOverlap => {
                        stats.rejected_invalid_joins += 1;
                        rejected_nodes.insert(cand.id());
                    }
                    EdgeResolveResult::NoEdge => {}
                }
            }

            if next_continuations.len() > 1 {
                // Target is a branching repeat node! Evaluate candidates exiting the repeat
                let mut best_cand: Option<(DirUnitigId, GraphEdge)> = None;
                let mut best_support = 0u32;
                let mut second_support = 0u32;

                for &(cand, edge) in &next_continuations {
                    let supp = matrix.get_support(curr_dir, cand);
                    if supp > best_support {
                        second_support = best_support;
                        best_support = supp;
                        best_cand = Some((cand, edge));
                    } else if supp > second_support {
                        second_support = supp;
                    }
                }

                if let Some((cand, edge_target_to_cand)) = best_cand {
                    let ratio = if second_support > 0 {
                        (best_support as f64) / (second_support as f64)
                    } else {
                        100.0
                    };

                    // A pair linking curr_dir to cand cannot span a repeat wider than the insert size.
                    let gap = target_node.seq_len as i64
                        - edge_to_target.overlap_k as i64
                        - edge_target_to_cand.overlap_k as i64;
                    if best_support >= params.min_support && ratio >= params.decision_threshold && gap <= dist.max_insert as i64 {
                        if !used_unique.contains(&cand.id()) && !graph.nodes[cand.id() as usize].is_deleted {
                            used_unique.insert(cand.id());
                            ext.push(PathStep {
                                dir: target,
                                incoming_edge: Some(edge_to_target),
                                overlap: edge_to_target.overlap_k as usize,
                                support: edge_to_target.support as u32,
                            });
                            ext.push(PathStep {
                                dir: cand,
                                incoming_edge: Some(edge_target_to_cand),
                                overlap: edge_target_to_cand.overlap_k as usize,
                                support: best_support,
                            });
                            curr_dir = cand;
                            continue;
                        }
                    }
                }
                // Exit repeat decision failed: protect target and candidate exit nodes
                rejected_nodes.insert(target.id());
                for &(cand, _) in &next_continuations {
                    rejected_nodes.insert(cand.id());
                }
                break;
            } else if next_continuations.len() == 1 {
                // Unbranched continuation
                if target_node.seq_len >= 500 {
                    if used_unique.contains(&target.id()) {
                        break;
                    }
                    used_unique.insert(target.id());
                }
                ext.push(PathStep {
                    dir: target,
                    incoming_edge: Some(edge_to_target),
                    overlap: edge_to_target.overlap_k as usize,
                    support: edge_to_target.support as u32,
                });
                curr_dir = target;
                continue;
            } else {
                // Dead end at target
                ext.push(PathStep {
                    dir: target,
                    incoming_edge: Some(edge_to_target),
                    overlap: edge_to_target.overlap_k as usize,
                    support: edge_to_target.support as u32,
                });
                break;
            }
        }

        // Case B: multiple outgoing edges directly from curr_dir
        let mut best_cand: Option<(DirUnitigId, GraphEdge)> = None;
        let mut best_support = 0u32;
        let mut second_support = 0u32;

        for &(cand, edge) in &valid_continuations {
            let supp = matrix.get_support(curr_dir, cand);
            if supp > best_support {
                second_support = best_support;
                best_support = supp;
                best_cand = Some((cand, edge));
            } else if supp > second_support {
                second_support = supp;
            }
        }

        if let Some((cand, edge_to_cand)) = best_cand {
            let ratio = if second_support > 0 {
                (best_support as f64) / (second_support as f64)
            } else {
                100.0
            };

            if best_support >= params.min_support && ratio >= params.decision_threshold {
                if !used_unique.contains(&cand.id()) && !graph.nodes[cand.id() as usize].is_deleted {
                    if graph.nodes[cand.id() as usize].seq_len >= 500 {
                        used_unique.insert(cand.id());
                    }
                    ext.push(PathStep {
                        dir: cand,
                        incoming_edge: Some(edge_to_cand),
                        overlap: edge_to_cand.overlap_k as usize,
                        support: best_support,
                    });
                    curr_dir = cand;
                    continue;
                }
            }
        }

        // Branching decision failed: protect candidate targets
        for &(cand, _) in &valid_continuations {
            rejected_nodes.insert(cand.id());
        }
        break;
    }

    ext
}

/// Verifies path validity, node survival, edge orientations, and actual sequence overlaps,
/// returning the fully spliced DNA sequence without mutating the graph.
///
/// Strict verification contract:
/// 1. First step: must conform to seed/path-start contract (incoming_edge = None, overlap = 0).
/// 2. For every step: node must be non-deleted, within bounds, and no unauthorized duplicate reuse.
/// 3. For every non-initial step (idx >= 1):
///    - incoming_edge must be Some(edge).
///    - edge.target must equal step.dir.
///    - edge.overlap_k as usize must equal step.overlap.
///    - Edge must truly exist in predecessor's outgoing edges for the correct strand.
///    - Edge support field must exactly match the authentic graph edge (no forged metadata).
///    - Suffix/prefix base-by-base real sequence must match.
pub fn validate_and_splice_path(
    graph: &CdbgGraph,
    full_path: &[PathStep],
) -> Result<Vec<u8>, &'static str> {
    if full_path.len() <= 1 {
        return Err("Path too short to splice");
    }

    // Seed/path-start contract for first step
    let first = &full_path[0];
    if first.incoming_edge.is_some() {
        return Err("First step must not have incoming edge");
    }
    if first.overlap != 0 {
        return Err("First step overlap must be 0");
    }

    // 1. Verify all constituent nodes are non-deleted, within bounds, and unique nodes are not duplicated
    let mut seen_unique_nodes = hashbrown::HashSet::with_capacity(full_path.len());
    for step in full_path {
        let node_id = step.dir.id() as usize;
        if node_id >= graph.nodes.len() {
            return Err("Node ID out of bounds");
        }
        let node = &graph.nodes[node_id];
        if node.is_deleted {
            return Err("Path contains deleted node");
        }
        if node.seq_len >= 500 && !seen_unique_nodes.insert(step.dir.id()) {
            return Err("Unauthorized unique node reuse in path");
        }
    }

    // 2. Validate edge existence, orientation, support, and base-by-base sequence overlaps
    let mut assembled_seq = Vec::new();
    for (idx, step) in full_path.iter().enumerate() {
        let curr_seq = graph.get_stranded_seq(step.dir);
        if idx == 0 {
            assembled_seq.extend_from_slice(&curr_seq);
        } else {
            let prev = &full_path[idx - 1];
            let prev_node = &graph.nodes[prev.dir.id() as usize];
            let outgoing = prev_node.edges(prev.dir.is_rc());

            // 1. incoming_edge exists
            let edge = step
                .incoming_edge
                .as_ref()
                .ok_or("Non-initial step missing incoming edge")?;

            // 3. incoming_edge.target equals current step.dir
            if edge.target != step.dir {
                return Err("Incoming edge target does not match step dir");
            }

            // 4. incoming_edge.overlap_k equals step.overlap
            if edge.overlap_k as usize != step.overlap {
                return Err("Incoming edge overlap_k does not match step overlap");
            }

            // 2. Edge truly exists in predecessor outgoing edges
            let candidate_edges: Vec<&GraphEdge> = outgoing
                .iter()
                .filter(|e| e.target == step.dir)
                .collect();
            if candidate_edges.is_empty() {
                return Err("Incoming edge does not exist in predecessor outgoing edges");
            }

            let matching_edge = candidate_edges
                .iter()
                .find(|e| e.overlap_k as usize == step.overlap)
                .ok_or("Incoming edge overlap_k does not exist in predecessor outgoing edges")?;

            // 5. Support field matches graph edge
            if edge.support != matching_edge.support {
                return Err("Incoming edge support does not match graph edge support");
            }

            match resolve_edge_between(graph, prev.dir, step.dir) {
                EdgeResolveResult::Valid(selected) if selected == *edge => {}
                EdgeResolveResult::Valid(_) => return Err("Incoming edge is not the selected edge"),
                EdgeResolveResult::Ambiguous => return Err("Ambiguous incoming edge"),
                EdgeResolveResult::InvalidSequenceOverlap => return Err("Sequence overlap mismatch at junction"),
                EdgeResolveResult::NoEdge => return Err("Incoming edge does not exist in predecessor outgoing edges"),
            }
            match resolve_edge_between(graph, step.dir.twin(), prev.dir.twin()) {
                EdgeResolveResult::Valid(reciprocal)
                    if reciprocal.overlap_k == edge.overlap_k
                        && reciprocal.support == edge.support => {}
                EdgeResolveResult::Valid(_)
                | EdgeResolveResult::Ambiguous
                | EdgeResolveResult::InvalidSequenceOverlap
                | EdgeResolveResult::NoEdge => return Err("Missing or inconsistent reciprocal edge"),
            }

            // 7. Real sequence base-by-base match
            let ov = step.overlap;
            let prev_seq = graph.get_stranded_seq(prev.dir);
            if ov == 0 || ov > prev_seq.len() || ov > curr_seq.len() {
                return Err("Invalid overlap length");
            }
            if prev_seq[prev_seq.len() - ov..] != curr_seq[..ov] {
                return Err("Sequence overlap mismatch at junction");
            }
            if assembled_seq.len() < ov
                || assembled_seq[assembled_seq.len() - ov..] != curr_seq[..ov]
            {
                return Err("Sequence overlap mismatch at junction");
            }

            assembled_seq.extend_from_slice(&curr_seq[ov..]);
        }
    }

    Ok(assembled_seq)
}

/// Resolves repeats in `graph` using paired-end link constraints from `PairLinkMatrix`,
/// returning comprehensive junction statistics.
pub fn resolve_repeats_exspander_with_stats(
    graph: &mut CdbgGraph,
    matrix: &PairLinkMatrix,
    dist: &InsertSizeDist,
    params: ExspanderParams,
) -> JunctionStats {
    let k = graph.k;
    let n = graph.nodes.len() as u32;
    let mut stats = JunctionStats::default();

    // 1. Identify long unique seed unitigs (length >= 500 bp, non-deleted, non-circular)
    let mut seeds: Vec<u32> = (0..n)
        .filter(|&id| {
            let node = &graph.nodes[id as usize];
            !node.is_deleted && !node.is_circular && node.seq_len >= 500
        })
        .collect();

    // Sort seeds by length descending
    seeds.sort_by_key(|&id| std::cmp::Reverse(graph.nodes[id as usize].seq_len));

    let mut used_in_path: HashSet<u32> = HashSet::new();
    let mut rejected_nodes: HashSet<u32> = HashSet::new();

    for seed_id in seeds {
        if used_in_path.contains(&seed_id) || graph.nodes[seed_id as usize].is_deleted {
            continue;
        }

        let seed_fwd = DirUnitigId::new(seed_id, false);

        // Exploration phase: use isolated local copy of used set
        let mut local_used = used_in_path.clone();
        local_used.insert(seed_id);

        // 1. Extend forward from 3' end
        let fwd_ext = extend_one_direction(seed_fwd, graph, matrix, dist, params, &mut local_used, &mut stats, &mut rejected_nodes);

        // 2. Extend backward from 5' end (forward from seed.twin())
        let bwd_ext = extend_one_direction(seed_fwd.twin(), graph, matrix, dist, params, &mut local_used, &mut stats, &mut rejected_nodes);

        // Construct full path with explicit PathSteps:
        // bwd_ext is reversed and twinned
        let mut full_path: Vec<PathStep> = Vec::new();
        let m = bwd_ext.len();
        if m > 0 {
            // First node in reversed path: bwd_ext[m-1].dir.twin() (start node, no incoming edge)
            full_path.push(PathStep {
                dir: bwd_ext[m - 1].dir.twin(),
                incoming_edge: None,
                overlap: 0,
                support: 0,
            });

            // Intermediate nodes from bwd_ext: look up authentic reciprocal edges in graph
            let mut bwd_failed = false;
            for i in 1..m {
                let step_idx = m - i;
                let step = &bwd_ext[step_idx];
                let source_dir = bwd_ext[step_idx].dir.twin();
                let target_dir = bwd_ext[step_idx - 1].dir.twin();

                let source_node = &graph.nodes[source_dir.id() as usize];
                let reciprocal = source_node.find_edge(source_dir.is_rc(), target_dir, step.overlap as u8);

                if let Some(&edge) = reciprocal {
                    full_path.push(PathStep {
                        dir: target_dir,
                        incoming_edge: Some(edge),
                        overlap: step.overlap,
                        support: edge.support as u32,
                    });
                } else {
                    tracing::debug!(
                        "[exSPAnder] Missing reciprocal edge from {:?} to {:?} (expected ov={}): rejecting backward path",
                        source_dir, target_dir, step.overlap
                    );
                    bwd_failed = true;
                    break;
                }
            }

            if !bwd_failed {
                // Transition from bwd_ext[0] to seed_fwd: look up reciprocal edge in graph
                let step0 = &bwd_ext[0];
                let source_dir = step0.dir.twin();
                let target_dir = seed_fwd;

                let source_node = &graph.nodes[source_dir.id() as usize];
                let reciprocal = source_node.find_edge(source_dir.is_rc(), target_dir, step0.overlap as u8);

                if let Some(&edge) = reciprocal {
                    full_path.push(PathStep {
                        dir: target_dir,
                        incoming_edge: Some(edge),
                        overlap: step0.overlap,
                        support: edge.support as u32,
                    });
                } else {
                    tracing::debug!(
                        "[exSPAnder] Missing reciprocal edge from {:?} to seed {:?} (expected ov={}): rejecting backward path",
                        source_dir, target_dir, step0.overlap
                    );
                    bwd_failed = true;
                }
            }

            if bwd_failed {
                // Conservatively reject candidate path when reciprocal edge is missing
                stats.rejected_invalid_joins += 1;
                rejected_nodes.insert(seed_id);
                for step in &bwd_ext {
                    rejected_nodes.insert(step.dir.id());
                }
                for step in &fwd_ext {
                    rejected_nodes.insert(step.dir.id());
                }
                continue;
            }
        } else {
            full_path.push(PathStep {
                dir: seed_fwd,
                incoming_edge: None,
                overlap: 0,
                support: 0,
            });
        }

        // Append forward extension steps
        full_path.extend(fwd_ext);

        // Splicing & Atomic Commit Phase
        if full_path.len() > 1 {
            match validate_and_splice_path(graph, &full_path) {
                Ok(assembled_seq) => {
                    // Record successful join metrics
                    for step in &full_path[1..] {
                        if step.overlap == k {
                            stats.valid_k_joins += 1;
                        } else if step.overlap == k.saturating_sub(1) {
                            stats.valid_k_minus_1_joins += 1;
                        } else {
                            stats.valid_other_joins += 1;
                        }
                    }

                    // Commit: mark non-circular constituent nodes as deleted
                    for step in &full_path {
                        let nid = step.dir.id() as usize;
                        if !graph.nodes[nid].is_circular {
                            graph.nodes[nid].is_deleted = true;
                        }
                    }

                    // Push assembled contig
                    let new_seq_id = graph.arena.push_seq(&assembled_seq);
                    let new_node_id = graph.nodes.len() as u32;
                    let new_len = assembled_seq.len() as u32;
                    let node = UnitigNode::new(new_node_id, new_seq_id, new_len, 20.0);
                    graph.nodes.push(node);

                    // Commit used nodes
                    used_in_path.extend(local_used);
                }
                Err(err) => {
                    stats.rejected_invalid_joins += 1;
                    tracing::debug!("[exSPAnder] Path rejected for seed {}: {}", seed_id, err);
                    for &nid in &local_used {
                        rejected_nodes.insert(nid);
                    }
                    for step in &full_path {
                        rejected_nodes.insert(step.dir.id());
                    }
                    // Contract 1: Do NOT commit: used_in_path, node deletions, and graph stay untouched
                }
            }
        }
    }

    // Post-exSPAnder Cleanup Phase:
    // Defined behavioral contracts:
    // Contract 1: Candidate paths rejected during exploration or validation make ZERO modifications to graph topology.
    // Contract 2: Short unitigs (< 500 bp) involved in rejected candidate paths are protected via `rejected_nodes`
    //             and MUST NOT be deleted by cleanup, preserving graph connectivity for downstream stages.
    // Contract 3: Unrelated short unitigs (< 500 bp, non-circular, not protected by `rejected_nodes` and never incorporated
    //             into any contig) are safely cleaned up to eliminate orphan sequencing artifacts.
    let mut connected = HashSet::new();
    for node in graph.nodes.iter().filter(|node| !node.is_deleted) {
        for edge in node.edges_fwd.iter().chain(&node.edges_rev) {
            if graph.nodes.get(edge.target.id() as usize).is_some_and(|target| !target.is_deleted) {
                connected.insert(node.id);
                connected.insert(edge.target.id());
            }
        }
    }
    for node in graph.nodes.iter_mut() {
        if !node.is_circular && node.seq_len < 500 && !rejected_nodes.contains(&node.id) && !connected.contains(&node.id) {
            node.is_deleted = true;
        }
    }

    stats
}

/// Resolves repeats in `graph` using paired-end link constraints from `PairLinkMatrix`.
pub fn resolve_repeats_exspander(
    graph: &mut CdbgGraph,
    matrix: &PairLinkMatrix,
    dist: &InsertSizeDist,
    params: ExspanderParams,
) {
    let stats = resolve_repeats_exspander_with_stats(graph, matrix, dist, params);
    tracing::info!(
        "[exSPAnder] Junction stats: valid k joins = {}, valid k-1 joins = {}, ambiguous rejected = {}, invalid rejected = {}",
        stats.valid_k_joins,
        stats.valid_k_minus_1_joins,
        stats.rejected_ambiguous_joins,
        stats.rejected_invalid_joins
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::bitseq::{reverse_complement_ascii, ReadBank};
    use crate::common::sketch::splitmix64;
    use crate::graph::build::build_cdbg;
    use crate::resolve::mapper::map_paired_reads;

    #[test]
    fn test_verification_gate_exspander_350bp_repeat_resolution() {
        // Verification Gate for Phase 5:
        // Simulate a 30 kb chromosome containing two copies of a 350 bp IS-like repeat
        // sequenced with 2 x 150 bp reads at 450 bp insert size:
        // Verify exspander accurately traverses across both copies of the 350 bp repeat
        // without misassembly, producing a single 30 kb contig.

        let k = 31;
        let alphabet = [b'A', b'C', b'G', b'T'];
        let gen_seq = |seed: u64, len: usize| -> Vec<u8> {
            let mut buf = Vec::with_capacity(len);
            for i in 0..len {
                let r = splitmix64(seed + i as u64);
                buf.push(alphabet[(r % 4) as usize]);
            }
            buf
        };

        // 350 bp IS-like repeat
        let rep350 = gen_seq(0x1111_3500, 350);
        let u1 = gen_seq(0xAAAA_1000, 9400);
        let u2 = gen_seq(0xBBBB_2000, 9700);
        let u3 = gen_seq(0xCCCC_3000, 10200);

        // 30 kb chromosome: U1 -> Rep -> U2 -> Rep -> U3
        let mut chrom = Vec::new();
        chrom.extend_from_slice(&u1);
        chrom.extend_from_slice(&rep350);
        chrom.extend_from_slice(&u2);
        chrom.extend_from_slice(&rep350);
        chrom.extend_from_slice(&u3);

        let total_chrom_len = chrom.len(); // 9400 + 350 + 9700 + 350 + 10200 = 30,000 bp
        assert_eq!(total_chrom_len, 30_000);

        // Simulate 40x depth of 2 x 150bp paired-end reads with insert size ~450 bp
        let mut bank = ReadBank::new();
        let num_pairs = (total_chrom_len * 40) / 300;

        let insert_size = 450usize;
        let read_len = 150usize;

        let mut rng: u64 = 0x5555_6666_7777_8888;
        for i in 0..num_pairs {
            rng = splitmix64(rng + i as u64);
            if total_chrom_len <= insert_size {
                break;
            }
            let start = (rng as usize) % (total_chrom_len - insert_size);

            let r1 = &chrom[start..start + read_len];
            let r2_frag = &chrom[start + insert_size - read_len..start + insert_size];
            let r2 = reverse_complement_ascii(r2_frag);

            bank.pe_forward.push_seq(r1);
            bank.pe_reverse.push_seq(&r2);
        }

        // 1. Build cDBG graph at k=31
        let mut graph = build_cdbg::<u64>(&bank, k, Some(3));

        // 2. Map paired reads to construct link matrix
        let (matrix, dist) = map_paired_reads::<u64>(&graph, &bank, k);

        // 3. Run exSPAnder
        let params = ExspanderParams {
            min_support: 5,
            decision_threshold: 2.5,
        };
        resolve_repeats_exspander(&mut graph, &matrix, &dist, params);

        // Verify: The 30 kb chromosome assembles into a single continuous contig
        let long_contigs: Vec<&UnitigNode> = graph
            .nodes
            .iter()
            .filter(|n| !n.is_deleted && n.seq_len >= 29_000)
            .collect();

        assert_eq!(
            long_contigs.len(),
            1,
            "exSPAnder must resolve both 350 bp repeats to produce a single ~30 kb contig! Found {} contigs >= 29kb",
            long_contigs.len()
        );
        assert!(
            long_contigs[0].seq_len >= 29_500,
            "Contig length should be ~30,000 bp, found {}",
            long_contigs[0].seq_len
        );
    }
}
