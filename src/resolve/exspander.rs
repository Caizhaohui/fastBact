//! SPAdes-style exSPAnder paired-end repeat resolution engine.
//!
//! Traverses across repeats longer than k_max using paired-end read insert size constraints,
//! connecting unique flanking unitigs across repetitive regions.

use crate::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, ResolverLifecycleState, UnitigNode};
use crate::resolve::mapper::{InsertSizeDist, PairLinkMatrix};
use hashbrown::HashSet;
use smallvec::SmallVec;
use std::collections::VecDeque;

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

/// Evaluates candidate-specific paired-end support for a path between `u1` and `u2` with known `path_separation`.
///
/// Span(p, P) = dist_u1_exit + path_separation + dist_u2_entry
/// A read pair is valid if Span in [min_insert, max_insert].
/// Normalized score Z = |Span - median| / (1.4826 * MAD).
/// Gaussian kernel weight w = exp(-0.5 * Z^2).
///
/// Returns (valid_count, total_weight). Each `pair_id` contributes at most once;
/// the first in-range row for that id supplies the weight.
pub fn score_candidate_path(
    matrix: &PairLinkMatrix,
    dist: &InsertSizeDist,
    u1: DirUnitigId,
    u2: DirUnitigId,
    path_separation: i64,
) -> (u32, f64) {
    let evidences = match matrix.get_evidence(u1, u2) {
        Some(evs) => evs,
        None => return (0, 0.0),
    };

    let mut valid_count = 0u32;
    let mut total_weight = 0.0f64;
    let mut seen_pair_ids = HashSet::with_capacity(evidences.len());

    let mad = if dist.mad > 1.0 { dist.mad } else { 1.0 };
    let norm_factor = 1.4826 * mad;

    for ev in evidences {
        let span = ev.dist_u1_exit as i64 + path_separation + ev.dist_u2_entry as i64;
        if span >= dist.min_insert as i64 && span <= dist.max_insert as i64 {
            if !seen_pair_ids.insert(ev.pair_id) {
                continue;
            }
            let z = (span as f64 - dist.median).abs() / norm_factor;
            let w = (-0.5 * z * z).exp();
            valid_count += 1;
            total_weight += w;
        }
    }

    (valid_count, total_weight)
}

const MAX_CANDIDATE_DEPTH: usize = 4;
const MAX_CANDIDATE_PATHS: usize = 32;
const MAX_DFS_STATES: usize = 256;
const MAX_EXTENSION_STEPS: usize = 64;
const MAX_QUEUE_CAPACITY: usize = 1024;

#[derive(Clone, Debug)]
struct CandidatePath {
    steps: SmallVec<[(DirUnitigId, GraphEdge); 4]>,
    endpoint: DirUnitigId,
    #[allow(dead_code)]
    separation: i64,
    support_count: u32,
    total_weight: f64,
    /// Route enters a node spliced by an earlier seed, or a unique node this
    /// seed already uses. Never committed; its evidence still competes.
    through_consumed: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SearchStatus {
    Complete,
    Truncated,
}

/// Admission class of one directed hop. Every search walker consumes this and
/// nothing else, so the hop rule cannot drift between them.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum HopClass {
    Live {
        edge: GraphEdge,
        sep: i64,
        shadow: bool,
    },
    Probe {
        edge: GraphEdge,
        sep: i64,
        shadow: bool,
    },
    SafePrune,
    Inadmissible,
    Unresolved,
}

/// Minimum `dist_u1_exit + dist_u2_entry` over evidence anchored at `curr_dir`.
///
/// Any later endpoint has `span >= sep + beta`. A hop with `sep + beta > max_insert`
/// cannot hide an in-window competitor. No anchored evidence yields 0, which
/// keeps the historical `sep <= max_insert` test.
fn separation_beta(matrix: &PairLinkMatrix, curr_dir: DirUnitigId) -> i64 {
    matrix
        .links
        .iter()
        .filter(|((u1, _), _)| *u1 == curr_dir)
        .flat_map(|(_, evs)| {
            evs.iter()
                .map(|ev| ev.dist_u1_exit as i64 + ev.dist_u2_entry as i64)
        })
        .min()
        .unwrap_or(0)
}

fn distinct_outgoing(graph: &CdbgGraph, node: DirUnitigId) -> Vec<DirUnitigId> {
    let mut targets = Vec::new();
    for edge in graph.nodes[node.id() as usize].edges(node.is_rc()) {
        if !targets.contains(&edge.target) {
            targets.push(edge.target);
        }
    }
    targets
}

/// Stateless hop admission. `step_depth == 0` means `u` is the fork and `sep_u`
/// is already `init_sep`; later hops add `len(u) - overlap`.
fn classify_hop(
    graph: &CdbgGraph,
    u: DirUnitigId,
    v: DirUnitigId,
    sep_u: i64,
    visited: &[u32],
    step_depth: usize,
    shadow_u: bool,
    curr_dir: DirUnitigId,
    ext_nodes: &HashSet<u32>,
    used_unique: &HashSet<u32>,
    consumed_nodes: &HashSet<u32>,
    max_separation: i64,
    beta: i64,
) -> HopClass {
    // Self-anchor scoring is applied by the walkers. `curr_dir` stays in the
    // contract so every call site passes the same anchor the bound was built for.
    let _ = curr_dir;
    let EdgeResolveResult::Valid(edge) = resolve_edge_between(graph, u, v) else {
        return HopClass::Inadmissible;
    };

    let sep_v = (if step_depth == 0 {
        sep_u
    } else {
        sep_u + graph.nodes[u.id() as usize].seq_len as i64
    }) - edge.overlap_k as i64;

    let v_node = &graph.nodes[v.id() as usize];
    // Raw insert budget. Equality remains admissible.
    if sep_v > max_separation {
        return HopClass::SafePrune;
    }
    // Evidence floor on short re-entries only. A repeat hop with
    // `sep + β > max_insert` cannot hide an in-window endpoint, so the probe
    // must not walk it (this keeps T26/T27 from false truncation). A first
    // visit stays inside the raw budget: the depth cap still has to see that
    // chain, and `β` is passed as 0 by the depth-cap check itself.
    let repeat_reentry = v_node.seq_len < 500 && visited.contains(&v.id()) && step_depth >= 1;
    if repeat_reentry && sep_v + beta > max_separation {
        return HopClass::SafePrune;
    }

    if v_node.is_deleted && !consumed_nodes.contains(&v.id()) {
        return HopClass::Inadmissible;
    }

    if v_node.is_circular {
        return HopClass::Unresolved;
    }

    if v_node.seq_len >= 500 && (ext_nodes.contains(&v.id()) || visited.contains(&v.id())) {
        return HopClass::SafePrune;
    }

    let shadow_v = shadow_u
        || (v_node.is_deleted && consumed_nodes.contains(&v.id()))
        || (v_node.seq_len >= 500 && used_unique.contains(&v.id()));

    if v_node.seq_len < 500 && visited.contains(&v.id()) && step_depth >= 1 {
        return HopClass::Probe {
            edge,
            sep: sep_v,
            shadow: shadow_v,
        };
    }

    HopClass::Live {
        edge,
        sep: sep_v,
        shadow: shadow_v,
    }
}

/// Bounded multi-hop candidate path search from `fork_node` to endpoints supported by read pairs with `curr_dir`.
fn find_candidate_paths(
    graph: &CdbgGraph,
    matrix: &PairLinkMatrix,
    dist: &InsertSizeDist,
    curr_dir: DirUnitigId,
    fork_node: DirUnitigId,
    init_sep: i64,
    used_unique: &HashSet<u32>,
    consumed_nodes: &HashSet<u32>,
    ext_nodes: &HashSet<u32>,
) -> (Vec<CandidatePath>, SearchStatus) {
    let max_separation = dist.max_insert as i64;
    let beta = separation_beta(matrix, curr_dir);
    let mut candidate_paths = Vec::new();
    let mut search_status = SearchStatus::Complete;
    let mut states_visited = 0usize;

    let mut queue = VecDeque::new();
    for target in distinct_outgoing(graph, fork_node) {
        let mut visited: SmallVec<[u32; 8]> = SmallVec::new();
        visited.push(curr_dir.id());
        visited.push(fork_node.id());
        match classify_hop(
            graph,
            fork_node,
            target,
            init_sep,
            &visited,
            0,
            false,
            curr_dir,
            ext_nodes,
            used_unique,
            consumed_nodes,
            max_separation,
            beta,
        ) {
            HopClass::Live { edge, sep, shadow } => {
                if queue.len() >= MAX_QUEUE_CAPACITY {
                    search_status = SearchStatus::Truncated;
                    break;
                }
                let mut steps: SmallVec<[(DirUnitigId, GraphEdge); 4]> = SmallVec::new();
                steps.push((target, edge));
                queue.push_back((target, sep, steps, visited, shadow));
            }
            HopClass::Unresolved => {
                search_status = SearchStatus::Truncated;
                break;
            }
            HopClass::Probe { .. } | HopClass::SafePrune | HopClass::Inadmissible => {}
        }
    }

    'bfs: while let Some((u, sep_u, steps, mut visited, through_consumed)) = queue.pop_front() {
        states_visited += 1;
        if states_visited > MAX_DFS_STATES {
            search_status = SearchStatus::Truncated;
            break;
        }

        visited.push(u.id());

        // A return to the anchor is evidence about this anchor, not an extension
        // endpoint. In-window support makes the decision unresolved; either way
        // the walk continues so a later exit can still be scored.
        if u.id() == curr_dir.id() {
            let (support_count, _) = score_candidate_path(matrix, dist, curr_dir, u, sep_u);
            if support_count > 0 {
                search_status = SearchStatus::Truncated;
            }
        } else {
            let (support_count, total_weight) =
                score_candidate_path(matrix, dist, curr_dir, u, sep_u);
            if support_count > 0 {
                candidate_paths.push(CandidatePath {
                    steps: steps.clone(),
                    endpoint: u,
                    separation: sep_u,
                    support_count,
                    total_weight,
                    through_consumed,
                });
                if candidate_paths.len() >= MAX_CANDIDATE_PATHS {
                    search_status = SearchStatus::Truncated;
                    break;
                }
            }
        }

        if steps.len() < MAX_CANDIDATE_DEPTH {
            for neighbor in distinct_outgoing(graph, u) {
                match classify_hop(
                    graph,
                    u,
                    neighbor,
                    sep_u,
                    &visited,
                    steps.len(),
                    through_consumed,
                    curr_dir,
                    ext_nodes,
                    used_unique,
                    consumed_nodes,
                    max_separation,
                    beta,
                ) {
                    HopClass::Probe { .. } => {
                        if pruned_repeat_hides_supported_endpoint(
                            graph,
                            matrix,
                            dist,
                            curr_dir,
                            u,
                            neighbor,
                            sep_u,
                            steps.len(),
                            &visited,
                            ext_nodes,
                            used_unique,
                            max_separation,
                            consumed_nodes,
                            beta,
                            through_consumed,
                        ) {
                            search_status = SearchStatus::Truncated;
                            break 'bfs;
                        }
                    }
                    HopClass::Live { edge, sep, shadow } => {
                        if queue.len() >= MAX_QUEUE_CAPACITY {
                            search_status = SearchStatus::Truncated;
                            break 'bfs;
                        }
                        let mut next_steps = steps.clone();
                        next_steps.push((neighbor, edge));
                        queue.push_back((neighbor, sep, next_steps, visited.clone(), shadow));
                    }
                    HopClass::Unresolved => {
                        search_status = SearchStatus::Truncated;
                        break 'bfs;
                    }
                    HopClass::SafePrune | HopClass::Inadmissible => {}
                }
            }
        } else if depth_limit_hides_continuation(
            graph,
            u,
            sep_u,
            &visited,
            steps.len(),
            through_consumed,
            curr_dir,
            ext_nodes,
            used_unique,
            consumed_nodes,
            max_separation,
            // The cap asks whether any raw-budget hop is still unexplored.
            // The evidence floor applies inside the repeat probe, which is
            // what decides whether that hop can actually score.
            0,
        ) {
            search_status = SearchStatus::Truncated;
            break;
        }
    }

    (candidate_paths, search_status)
}

/// True when unrolling pruned short-repeat hops would score an endpoint the
/// acyclic search cannot see, or would run into an unresolved hop.
///
/// The walk uses `classify_hop`, so a short node may recur and a node of at
/// least 500 bp may appear at most once. Landing inside the insert window with
/// paired-end support, including support for the anchor itself, abstains.
fn pruned_repeat_hides_supported_endpoint(
    graph: &CdbgGraph,
    matrix: &PairLinkMatrix,
    dist: &InsertSizeDist,
    curr_dir: DirUnitigId,
    from: DirUnitigId,
    next_target: DirUnitigId,
    sep_u: i64,
    steps_before: usize,
    visited: &[u32],
    ext_nodes: &HashSet<u32>,
    used_unique: &HashSet<u32>,
    max_separation: i64,
    consumed_nodes: &HashSet<u32>,
    beta: i64,
    shadow_u: bool,
) -> bool {
    let (entry_sep, entry_shadow) = match classify_hop(
        graph,
        from,
        next_target,
        sep_u,
        visited,
        steps_before,
        shadow_u,
        curr_dir,
        ext_nodes,
        used_unique,
        consumed_nodes,
        max_separation,
        beta,
    ) {
        HopClass::Live { sep, shadow, .. } | HopClass::Probe { sep, shadow, .. } => (sep, shadow),
        HopClass::Unresolved => return true,
        HopClass::SafePrune | HopClass::Inadmissible => return false,
    };

    let mut queue = VecDeque::new();
    let entry_visited: SmallVec<[u32; 8]> = SmallVec::from_slice(visited);
    queue.push_back((
        next_target,
        entry_sep,
        steps_before + 1,
        entry_visited,
        entry_shadow,
    ));
    let mut states = 0usize;

    while let Some((node, sep, steps, mut node_visited, shadow)) = queue.pop_front() {
        states += 1;
        if states > MAX_DFS_STATES || queue.len() >= MAX_QUEUE_CAPACITY {
            return true;
        }

        node_visited.push(node.id());

        let (count, _) = score_candidate_path(matrix, dist, curr_dir, node, sep);
        if count > 0 {
            return true;
        }

        if steps >= MAX_CANDIDATE_DEPTH {
            if depth_limit_hides_continuation(
                graph,
                node,
                sep,
                &node_visited,
                steps,
                shadow,
                curr_dir,
                ext_nodes,
                used_unique,
                consumed_nodes,
                max_separation,
                beta,
            ) {
                return true;
            }
            continue;
        }

        for neighbor in distinct_outgoing(graph, node) {
            match classify_hop(
                graph,
                node,
                neighbor,
                sep,
                &node_visited,
                steps,
                shadow,
                curr_dir,
                ext_nodes,
                used_unique,
                consumed_nodes,
                max_separation,
                beta,
            ) {
                HopClass::Live {
                    sep: sep_next,
                    shadow: shadow_next,
                    ..
                }
                | HopClass::Probe {
                    sep: sep_next,
                    shadow: shadow_next,
                    ..
                } => {
                    if queue.len() >= MAX_QUEUE_CAPACITY {
                        return true;
                    }
                    queue.push_back((
                        neighbor,
                        sep_next,
                        steps + 1,
                        node_visited.clone(),
                        shadow_next,
                    ));
                }
                HopClass::Unresolved => return true,
                HopClass::SafePrune | HopClass::Inadmissible => {}
            }
        }
    }
    false
}

/// True when `u` still has a hop the depth cap did not enqueue.
///
/// Live, probe, and unresolved hops are unfinished business. Safe prunes and
/// inadmissible hops are not.
fn depth_limit_hides_continuation(
    graph: &CdbgGraph,
    u: DirUnitigId,
    sep_u: i64,
    visited: &[u32],
    step_depth: usize,
    shadow_u: bool,
    curr_dir: DirUnitigId,
    ext_nodes: &HashSet<u32>,
    used_unique: &HashSet<u32>,
    consumed_nodes: &HashSet<u32>,
    max_separation: i64,
    beta: i64,
) -> bool {
    for next_target in distinct_outgoing(graph, u) {
        match classify_hop(
            graph,
            u,
            next_target,
            sep_u,
            visited,
            step_depth,
            shadow_u,
            curr_dir,
            ext_nodes,
            used_unique,
            consumed_nodes,
            max_separation,
            beta,
        ) {
            HopClass::Live { .. } | HopClass::Probe { .. } | HopClass::Unresolved => return true,
            HopClass::SafePrune | HopClass::Inadmissible => {}
        }
    }
    false
}

/// True when a consumed route that leaves `best`'s exit has enough evidence to
/// withhold the decision. A shadow that continues through `best` is the same
/// route and does not compete with it.
fn shadow_blocks_decision(best: &CandidatePath, shadows: &[CandidatePath], threshold: f64) -> bool {
    shadows.iter().any(|shadow| {
        let is_continuation = shadow.steps.len() >= best.steps.len()
            && shadow.steps[..best.steps.len()]
                .iter()
                .map(|step| step.0)
                .eq(best.steps.iter().map(|step| step.0));
        !is_continuation && best.total_weight < threshold * shadow.total_weight
    })
}

/// Keep outgoing neighbors that the step cap did not merge.
///
/// After the prefix contig deletes `curr_dir`, a short neighbor whose only active
/// edge ran back to that node would otherwise be removed as an orphan.
fn protect_unmerged_frontier(
    graph: &CdbgGraph,
    curr_dir: DirUnitigId,
    ext: &[PathStep],
    cleanup_protected_nodes: &mut HashSet<u32>,
) {
    let curr_node = &graph.nodes[curr_dir.id() as usize];
    for edge in curr_node.edges(curr_dir.is_rc()) {
        let frontier_id = edge.target.id();
        if ext.iter().any(|step| step.dir.id() == frontier_id) {
            continue;
        }
        cleanup_protected_nodes.insert(frontier_id);
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
    cleanup_protected_nodes: &mut HashSet<u32>,
    consumed_nodes: &HashSet<u32>,
) -> Vec<PathStep> {
    let mut ext: Vec<PathStep> = Vec::new();

    loop {
        if ext.len() >= MAX_EXTENSION_STEPS {
            // The node in hand is the last merged step. Its outgoing neighbors were
            // not consumed; keep them out of orphan cleanup after the prefix commit.
            protect_unmerged_frontier(graph, curr_dir, &ext, cleanup_protected_nodes);
            break;
        }

        // Pairs that close on this anchor are not an extension endpoint. An
        // in-window hit at separation 0 means splicing forward would commit the
        // anchor as if that evidence were not there.
        let (self_support, _) = score_candidate_path(matrix, dist, curr_dir, curr_dir, 0);
        if self_support > 0 {
            if !ext.is_empty() {
                protect_unmerged_frontier(graph, curr_dir, &ext, cleanup_protected_nodes);
            }
            break;
        }

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
                    cleanup_protected_nodes.insert(target.id());
                }
                EdgeResolveResult::InvalidSequenceOverlap => {
                    stats.rejected_invalid_joins += 1;
                    cleanup_protected_nodes.insert(target.id());
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
                        cleanup_protected_nodes.insert(cand.id());
                    }
                    EdgeResolveResult::InvalidSequenceOverlap => {
                        stats.rejected_invalid_joins += 1;
                        cleanup_protected_nodes.insert(cand.id());
                    }
                    EdgeResolveResult::NoEdge => {}
                }
            }

            if next_continuations.len() > 1 {
                // Target is a branching repeat node! Evaluate complete candidate paths exiting the repeat
                let ext_nodes: HashSet<u32> = ext.iter().map(|s: &PathStep| s.dir.id()).collect();
                let init_sep = target_node.seq_len as i64 - edge_to_target.overlap_k as i64;
                let (candidate_paths, search_status) = find_candidate_paths(
                    graph,
                    matrix,
                    dist,
                    curr_dir,
                    target,
                    init_sep,
                    used_unique,
                    consumed_nodes,
                    &ext_nodes,
                );

                if let SearchStatus::Truncated = search_status {
                    cleanup_protected_nodes.insert(target.id());
                    for &(cand, _) in &next_continuations {
                        cleanup_protected_nodes.insert(cand.id());
                    }
                    break;
                }

                let (shadow_paths, mut candidate_paths): (Vec<_>, Vec<_>) = candidate_paths
                    .into_iter()
                    .partition(|candidate| candidate.through_consumed);

                if candidate_paths.is_empty() {
                    cleanup_protected_nodes.insert(target.id());
                    for &(cand, _) in &next_continuations {
                        cleanup_protected_nodes.insert(cand.id());
                    }
                    break;
                }

                candidate_paths.sort_by(|a, b| {
                    b.total_weight
                        .partial_cmp(&a.total_weight)
                        .unwrap_or(std::cmp::Ordering::Equal)
                        .then_with(|| b.support_count.cmp(&a.support_count))
                });

                let best = &candidate_paths[0];
                let second_count = if candidate_paths.len() > 1 {
                    candidate_paths[1].support_count
                } else {
                    0
                };
                let second_weight = if candidate_paths.len() > 1 {
                    candidate_paths[1].total_weight
                } else {
                    0.0
                };

                let rivals_agree = if best.support_count < params.min_support {
                    false
                } else if second_count == 0 || second_weight < 1e-9 {
                    // Single candidate: require average weight >= 0.05 (Z <= 2.45) to reject tail outliers
                    (best.total_weight / best.support_count as f64) >= 0.05
                } else {
                    (best.total_weight / second_weight) >= params.decision_threshold
                };
                let is_decision_valid = rivals_agree
                    && !shadow_blocks_decision(best, &shadow_paths, params.decision_threshold);

                let all_steps_valid = best.steps.iter().all(|(step_dir, _)| {
                    let step_node = &graph.nodes[step_dir.id() as usize];
                    !step_node.is_deleted
                        && !(step_node.seq_len >= 500 && used_unique.contains(&step_dir.id()))
                });

                if is_decision_valid && all_steps_valid {
                    ext.push(PathStep {
                        dir: target,
                        incoming_edge: Some(edge_to_target),
                        overlap: edge_to_target.overlap_k as usize,
                        support: edge_to_target.support as u32,
                    });
                    if target_node.seq_len >= 500 {
                        used_unique.insert(target.id());
                    }

                    for (step_dir, step_edge) in &best.steps {
                        let step_node = &graph.nodes[step_dir.id() as usize];
                        if step_node.seq_len >= 500 {
                            used_unique.insert(step_dir.id());
                        }
                        let support = if *step_dir == best.endpoint {
                            best.support_count
                        } else {
                            step_edge.support as u32
                        };
                        ext.push(PathStep {
                            dir: *step_dir,
                            incoming_edge: Some(*step_edge),
                            overlap: step_edge.overlap_k as usize,
                            support,
                        });
                    }

                    curr_dir = best.endpoint;
                    continue;
                }

                // Exit repeat decision failed: protect target and candidate exit nodes
                cleanup_protected_nodes.insert(target.id());
                for &(cand, _) in &next_continuations {
                    cleanup_protected_nodes.insert(cand.id());
                }
                break;
            } else if next_continuations.len() == 1 {
                // Unbranched continuation: protect against infinite repeat cycle loops (F1-P2-003)
                if target.id() == curr_dir.id() || ext.iter().any(|s| s.dir.id() == target.id()) {
                    break;
                }
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
        let ext_nodes: HashSet<u32> = ext.iter().map(|s: &PathStep| s.dir.id()).collect();
        let (candidate_paths, search_status) = find_candidate_paths(
            graph,
            matrix,
            dist,
            curr_dir,
            curr_dir,
            0,
            used_unique,
            consumed_nodes,
            &ext_nodes,
        );

        if let SearchStatus::Truncated = search_status {
            for &(cand, _) in &valid_continuations {
                cleanup_protected_nodes.insert(cand.id());
            }
            break;
        }

        let (shadow_paths, mut candidate_paths): (Vec<_>, Vec<_>) = candidate_paths
            .into_iter()
            .partition(|candidate| candidate.through_consumed);

        if candidate_paths.is_empty() {
            for &(cand, _) in &valid_continuations {
                cleanup_protected_nodes.insert(cand.id());
            }
            break;
        }

        candidate_paths.sort_by(|a, b| {
            b.total_weight
                .partial_cmp(&a.total_weight)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| b.support_count.cmp(&a.support_count))
        });

        let best = &candidate_paths[0];
        let second_count = if candidate_paths.len() > 1 {
            candidate_paths[1].support_count
        } else {
            0
        };
        let second_weight = if candidate_paths.len() > 1 {
            candidate_paths[1].total_weight
        } else {
            0.0
        };

        let rivals_agree = if best.support_count < params.min_support {
            false
        } else if second_count == 0 || second_weight < 1e-9 {
            // Single candidate: require average weight >= 0.05 (Z <= 2.45) to reject tail outliers
            (best.total_weight / best.support_count as f64) >= 0.05
        } else {
            (best.total_weight / second_weight) >= params.decision_threshold
        };
        let is_decision_valid =
            rivals_agree && !shadow_blocks_decision(best, &shadow_paths, params.decision_threshold);

        let all_steps_valid = best.steps.iter().all(|(step_dir, _)| {
            let step_node = &graph.nodes[step_dir.id() as usize];
            !step_node.is_deleted
                && !(step_node.seq_len >= 500 && used_unique.contains(&step_dir.id()))
        });

        if is_decision_valid && all_steps_valid {
            for (step_dir, step_edge) in &best.steps {
                let step_node = &graph.nodes[step_dir.id() as usize];
                if step_node.seq_len >= 500 {
                    used_unique.insert(step_dir.id());
                }
                let support = if *step_dir == best.endpoint {
                    best.support_count
                } else {
                    step_edge.support as u32
                };
                ext.push(PathStep {
                    dir: *step_dir,
                    incoming_edge: Some(*step_edge),
                    overlap: step_edge.overlap_k as usize,
                    support,
                });
            }

            curr_dir = best.endpoint;
            continue;
        }

        // Branching decision failed: protect candidate targets
        for &(cand, _) in &valid_continuations {
            cleanup_protected_nodes.insert(cand.id());
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
            let candidate_edges: Vec<&GraphEdge> =
                outgoing.iter().filter(|e| e.target == step.dir).collect();
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
                EdgeResolveResult::Valid(_) => {
                    return Err("Incoming edge is not the selected edge")
                }
                EdgeResolveResult::Ambiguous => return Err("Ambiguous incoming edge"),
                EdgeResolveResult::InvalidSequenceOverlap => {
                    return Err("Sequence overlap mismatch at junction")
                }
                EdgeResolveResult::NoEdge => {
                    return Err("Incoming edge does not exist in predecessor outgoing edges")
                }
            }
            match resolve_edge_between(graph, step.dir.twin(), prev.dir.twin()) {
                EdgeResolveResult::Valid(reciprocal)
                    if reciprocal.overlap_k == edge.overlap_k
                        && reciprocal.support == edge.support => {}
                EdgeResolveResult::Valid(_)
                | EdgeResolveResult::Ambiguous
                | EdgeResolveResult::InvalidSequenceOverlap
                | EdgeResolveResult::NoEdge => {
                    return Err("Missing or inconsistent reciprocal edge")
                }
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

/// Failure of the checked single-use repeat-resolution contract.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExspanderError {
    AlreadyResolved,
    ResolutionInProgress,
    ResolutionPreviouslyFailed,
}

impl std::fmt::Display for ExspanderError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::AlreadyResolved => {
                write!(
                    f,
                    "repeat resolution has already been executed on this graph"
                )
            }
            Self::ResolutionInProgress => {
                write!(
                    f,
                    "repeat resolution is currently in progress on this graph"
                )
            }
            Self::ResolutionPreviouslyFailed => {
                write!(f, "repeat resolution previously failed on this graph")
            }
        }
    }
}

impl std::error::Error for ExspanderError {}

/// Marks `Failed` when resolution unwinds before a normal completion.
///
/// Partial graph edits are left in place. Panic and OOM are not a transaction.
struct ResolverLifecycleGuard<'a> {
    graph: &'a mut CdbgGraph,
    terminal: ResolverLifecycleState,
}

impl<'a> ResolverLifecycleGuard<'a> {
    fn enter(graph: &'a mut CdbgGraph) -> Self {
        graph.set_resolver_state(ResolverLifecycleState::Running);
        Self {
            graph,
            terminal: ResolverLifecycleState::Failed,
        }
    }

    fn graph_mut(&mut self) -> &mut CdbgGraph {
        self.graph
    }

    fn finish_ok(&mut self) {
        self.terminal = ResolverLifecycleState::Completed;
        self.graph
            .set_resolver_state(ResolverLifecycleState::Completed);
    }
}

impl Drop for ResolverLifecycleGuard<'_> {
    fn drop(&mut self) {
        if self.graph.resolver_state() == ResolverLifecycleState::Running {
            self.graph.set_resolver_state(self.terminal);
        }
    }
}

/// Resolves repeats in `graph` using paired-end link constraints from `PairLinkMatrix`,
/// returning comprehensive junction statistics.
///
/// The first call on a [`ResolverLifecycleState::NeverRun`] graph runs Stage 5 and,
/// on success, records [`ResolverLifecycleState::Completed`]. Any later call returns
/// [`ExspanderError`] before exploration, mutation, or Stage 5 cleanup.
pub fn resolve_repeats_exspander_with_stats(
    graph: &mut CdbgGraph,
    matrix: &PairLinkMatrix,
    dist: &InsertSizeDist,
    params: ExspanderParams,
) -> Result<JunctionStats, ExspanderError> {
    match graph.resolver_state() {
        ResolverLifecycleState::NeverRun => {}
        ResolverLifecycleState::Completed => return Err(ExspanderError::AlreadyResolved),
        ResolverLifecycleState::Running => return Err(ExspanderError::ResolutionInProgress),
        ResolverLifecycleState::Failed => return Err(ExspanderError::ResolutionPreviouslyFailed),
    }

    let mut guard = ResolverLifecycleGuard::enter(graph);
    let stats = execute_repeat_resolution(guard.graph_mut(), matrix, dist, params);
    guard.finish_ok();
    Ok(stats)
}

fn execute_repeat_resolution(
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
    let mut cleanup_protected_nodes: HashSet<u32> = HashSet::new();
    let mut consumed_nodes: HashSet<u32> = HashSet::new();

    for seed_id in seeds {
        if used_in_path.contains(&seed_id) || graph.nodes[seed_id as usize].is_deleted {
            continue;
        }

        let seed_fwd = DirUnitigId::new(seed_id, false);

        // Exploration phase: use isolated local copy of used set
        let mut local_used = used_in_path.clone();
        local_used.insert(seed_id);

        // 1. Extend forward from 3' end
        let fwd_ext = extend_one_direction(
            seed_fwd,
            graph,
            matrix,
            dist,
            params,
            &mut local_used,
            &mut stats,
            &mut cleanup_protected_nodes,
            &consumed_nodes,
        );

        // 2. Extend backward from 5' end (forward from seed.twin())
        let bwd_ext = extend_one_direction(
            seed_fwd.twin(),
            graph,
            matrix,
            dist,
            params,
            &mut local_used,
            &mut stats,
            &mut cleanup_protected_nodes,
            &consumed_nodes,
        );

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
                let reciprocal =
                    source_node.find_edge(source_dir.is_rc(), target_dir, step.overlap as u8);

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
                let reciprocal =
                    source_node.find_edge(source_dir.is_rc(), target_dir, step0.overlap as u8);

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
                cleanup_protected_nodes.insert(seed_id);
                for step in &bwd_ext {
                    cleanup_protected_nodes.insert(step.dir.id());
                }
                for step in &fwd_ext {
                    cleanup_protected_nodes.insert(step.dir.id());
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

                    // Protect every active neighbor that is not itself consumed.
                    // Incoming side branches, intermediate joins, outgoing frontiers,
                    // and reverse-orientation neighbors lose their edge into the path
                    // once constituent nodes are deleted. Register them before that
                    // deletion so Stage 5 orphan cleanup cannot drop their sequence.
                    // Validation failure does not reach here, so external neighbors
                    // stay out of `cleanup_protected_nodes` when the commit is rolled back.
                    // These ids are the constituent nodes deleted below. Later seeds
                    // may traverse them as shadows, but must not splice them again.
                    consumed_nodes.extend(full_path.iter().map(|step| step.dir.id()));
                    for step in &full_path {
                        let neighbor_ids: Vec<u32> = {
                            let node = &graph.nodes[step.dir.id() as usize];
                            node.edges_fwd
                                .iter()
                                .chain(node.edges_rev.iter())
                                .map(|edge| edge.target.id())
                                .collect()
                        };
                        for neighbor_id in neighbor_ids {
                            if consumed_nodes.contains(&neighbor_id) {
                                continue;
                            }
                            if graph.nodes[neighbor_id as usize].is_deleted {
                                continue;
                            }
                            cleanup_protected_nodes.insert(neighbor_id);
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
                        cleanup_protected_nodes.insert(nid);
                    }
                    for step in &full_path {
                        cleanup_protected_nodes.insert(step.dir.id());
                    }
                    // Contract 1: Do NOT commit: used_in_path, node deletions, and graph stay untouched
                }
            }
        }
    }

    // Post-exSPAnder Cleanup Phase:
    // Defined behavioral contracts:
    // Contract 1: Candidate paths rejected during exploration or validation make ZERO modifications to graph topology.
    // Contract 2: Short unitigs (< 500 bp) involved in rejected candidate paths are protected via `cleanup_protected_nodes`
    //             and MUST NOT be deleted by cleanup, preserving graph connectivity for downstream stages.
    // Contract 3: Unrelated short unitigs (< 500 bp, non-circular, not protected by `cleanup_protected_nodes` and never incorporated
    //             into any contig) are safely cleaned up to eliminate orphan sequencing artifacts.
    let mut connected = HashSet::new();
    for node in graph.nodes.iter().filter(|node| !node.is_deleted) {
        for edge in node.edges_fwd.iter().chain(&node.edges_rev) {
            if graph
                .nodes
                .get(edge.target.id() as usize)
                .is_some_and(|target| !target.is_deleted)
            {
                connected.insert(node.id);
                connected.insert(edge.target.id());
            }
        }
    }
    for node in graph.nodes.iter_mut() {
        if !node.is_circular
            && node.seq_len < 500
            && !cleanup_protected_nodes.contains(&node.id)
            && !connected.contains(&node.id)
        {
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
) -> Result<(), ExspanderError> {
    let stats = resolve_repeats_exspander_with_stats(graph, matrix, dist, params)?;
    tracing::info!(
        "[exSPAnder] Junction stats: valid k joins = {}, valid k-1 joins = {}, ambiguous rejected = {}, invalid rejected = {}",
        stats.valid_k_joins,
        stats.valid_k_minus_1_joins,
        stats.rejected_ambiguous_joins,
        stats.rejected_invalid_joins
    );
    Ok(())
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
        resolve_repeats_exspander(&mut graph, &matrix, &dist, params)
            .expect("repeat resolution failed");

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
        assert_eq!(graph.resolver_state(), ResolverLifecycleState::Completed);
    }

    #[test]
    fn non_never_run_states_reject_before_mutation() {
        let params = ExspanderParams::default();
        let matrix = PairLinkMatrix::new();
        let dist = InsertSizeDist::default();

        for (state, expected) in [
            (
                ResolverLifecycleState::Running,
                ExspanderError::ResolutionInProgress,
            ),
            (
                ResolverLifecycleState::Failed,
                ExspanderError::ResolutionPreviouslyFailed,
            ),
            (
                ResolverLifecycleState::Completed,
                ExspanderError::AlreadyResolved,
            ),
        ] {
            let mut graph = CdbgGraph::new(31);
            let seq_id = graph.arena.push_seq(b"ACGTACGTACGT");
            graph.nodes.push(UnitigNode::new(0, seq_id, 12, 10.0));
            graph.set_resolver_state(state);
            let before = format!("{graph:?}");
            let err = resolve_repeats_exspander_with_stats(&mut graph, &matrix, &dist, params);
            assert_eq!(err, Err(expected));
            assert_eq!(graph.resolver_state(), state);
            assert_eq!(format!("{graph:?}"), before);
        }
    }
}
