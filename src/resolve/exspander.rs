//! SPAdes-style exSPAnder paired-end repeat resolution engine.
//!
//! Traverses across repeats longer than k_max using paired-end read insert size constraints,
//! connecting unique flanking unitigs across repetitive regions.

use crate::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, ResolverLifecycleState, UnitigNode};
use crate::resolve::mapper::{InsertSizeDist, PairLinkMatrix, ReadPairEvidence};
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
    /// Fork decisions and seed walks withheld by the F23 provenance guard:
    /// the walk's last long-unique anchor was beyond the insert horizon at a
    /// branching repeat (horizon rejection, CS-3), no attributable anchor
    /// covered a fork in a contaminated zone (PAG), a branch had no
    /// admissible pairs under PAG (§3.10.2), or the seed itself was not
    /// anchor-eligible (SI-2 seed rejection).
    pub withheld_no_unique_anchor: usize,
    /// Breakdown: SI-2 unanchored seed rejected before exploration
    pub withheld_si2_seed: usize,
    /// Breakdown: fork distance exceeds max_insert
    pub withheld_horizon: usize,
    /// Breakdown: contaminated fork with no attributable anchor
    pub withheld_no_attributable_anchor: usize,
    /// Breakdown: contaminated branch with no admissible pairs under PAG
    pub withheld_pag_empty: usize,
}

// ==========================================================================
// F23 provenance architecture — Opus Hybrid Option A+B
// (DCC + BCI + PB, enforced at runtime by ARG and PAG)
//
// Reference: F23_FINAL_PROVENANCE_ARCHITECTURE.md (Claude Opus 5.5 ruling):
//   DCC  — Directed Context Count per node orientation (§3.1, Q1/Q2)
//   BCI  — Bounded Contamination Index via bounded reverse search (§Q6, LB-1)
//   PB   — Provenance Boundary scalar on composite unitigs (§Q4, MR-1/MR-2)
//   ARG  — Anchor Reset Gate (§3.7, AR-1..AR-4)
//   PAG  — PE Admissibility Gate: in a shared/contaminated zone, evidence
//          keyed at the walk head has zero deciding weight and candidates
//          are scored from the unique anchor instead (§3.6)
// ==========================================================================

/// Provenance metadata for a single unitig node, kept in a side table keyed by
/// node id (the graph node itself is untouched).
///
/// All classification fields are orientation-specific (DSC-3): a walk that
/// enters a node in the reverse-complement direction reads the `_rev` fields,
/// computed independently on the reverse-complement strand graph.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct NodeProvenanceMeta {
    /// DCC: number of distinct external incoming lineages (graph context
    /// sources) that reach this node's forward orientation. Capped at 2.
    /// 0 means no external context is established (a cyclic source SCC with
    /// no 0-in-degree vertex): attribution is unknown, so the node is not
    /// contaminating but is never anchor-eligible (DCC strict, AR-2).
    pub external_context_count_fwd: u8,
    /// DCC: same, for the reverse-complement orientation.
    pub external_context_count_rev: u8,
    /// BCI: a multi-context confluence is reachable from the forward
    /// orientation within `max_insert` bp (pooled-locus zone lies ahead).
    pub bci_fwd: bool,
    /// BCI: same, for the reverse-complement orientation.
    pub bci_rev: bool,
    /// BCI downstream: this node is within `max_insert` bp downstream of a
    /// multi-context confluence in the walk direction.
    pub bci_downstream_fwd: bool,
    /// BCI downstream: same, for reverse-complement orientation.
    pub bci_downstream_rev: bool,
    /// Maximum 5' entry-shadow reach in forward orientation (coordinate b).
    /// Evidence-source read starts p must satisfy p > b + H (i.e. dist_u1_exit < seq_len - (b + H)).
    pub entry_shadow_reach_fwd: u32,
    /// Maximum 5' entry-shadow reach in reverse-complement orientation.
    pub entry_shadow_reach_rev: u32,
    /// Exit isolation (R-I): the 3' exit of this node in forward orientation
    /// is physically separated from any upstream multi-context confluence (> max_insert).
    pub exit_isolated_fwd: bool,
    /// Exit isolation (R-I): same, for reverse-complement orientation.
    pub exit_isolated_rev: bool,
    /// Structural in-degree in forward orientation at analysis snapshot time.
    pub structural_in_degree_fwd: u16,
    /// Structural in-degree in reverse-complement orientation at analysis snapshot time.
    pub structural_in_degree_rev: u16,
    /// Structural out-degree in forward orientation at analysis snapshot time.
    pub structural_out_degree_fwd: u16,
    /// Structural out-degree in reverse-complement orientation at analysis snapshot time.
    pub structural_out_degree_rev: u16,
    /// PB: for composite unitigs, the spliced 5'-offset at which the sequence
    /// shifts from a unique-locus prefix to a pooled-locus suffix
    /// (`None` = uniformly classified). Merge contract MR-1. A read at
    /// 5'-offset `p` has `dist_u1_exit = seq_len - p`, so pairs from the
    /// pooled suffix — those with `dist_u1_exit <= seq_len - boundary` —
    /// are inadmissible as unique-locus evidence (MR-2).
    pub provenance_boundary: Option<u32>,
}

impl NodeProvenanceMeta {
    /// DCC context count in the given walk direction.
    #[inline]
    pub fn external_context_count(&self, is_rc: bool) -> u8 {
        if is_rc {
            self.external_context_count_rev
        } else {
            self.external_context_count_fwd
        }
    }

    /// BCI flag in the given walk direction.
    #[inline]
    pub fn bci(&self, is_rc: bool) -> bool {
        if is_rc {
            self.bci_rev
        } else {
            self.bci_fwd
        }
    }

    /// BCI downstream flag in the given walk direction.
    #[inline]
    pub fn bci_downstream(&self, is_rc: bool) -> bool {
        if is_rc {
            self.bci_downstream_rev
        } else {
            self.bci_downstream_fwd
        }
    }

    /// Exit isolation flag in the given walk direction.
    #[inline]
    pub fn exit_isolated(&self, is_rc: bool) -> bool {
        if is_rc {
            self.exit_isolated_rev
        } else {
            self.exit_isolated_fwd
        }
    }

    /// Maximum 5' entry-shadow reach in the given walk direction (coordinate b).
    #[inline]
    pub fn entry_shadow_reach(&self, is_rc: bool) -> u32 {
        if is_rc {
            self.entry_shadow_reach_rev
        } else {
            self.entry_shadow_reach_fwd
        }
    }

    /// Structural in-degree at analysis snapshot time in the given walk direction.
    #[inline]
    pub fn structural_in_degree(&self, is_rc: bool) -> u16 {
        if is_rc {
            self.structural_in_degree_rev
        } else {
            self.structural_in_degree_fwd
        }
    }

    /// Structural out-degree at analysis snapshot time in the given walk direction.
    #[inline]
    pub fn structural_out_degree(&self, is_rc: bool) -> u16 {
        if is_rc {
            self.structural_out_degree_rev
        } else {
            self.structural_out_degree_fwd
        }
    }

    /// True when this orientation of the node is a pooled (multi-context) locus.
    #[inline]
    pub fn is_multi_context(&self, is_rc: bool) -> bool {
        self.external_context_count(is_rc) >= 2
    }

    /// True when traversing this orientation taints the walk: shared locus or
    /// bounded contamination ahead (§3.4 step 1).
    #[inline]
    pub fn contaminating(&self, is_rc: bool) -> bool {
        self.is_multi_context(is_rc) || self.bci(is_rc)
    }
}

/// Side table of [`NodeProvenanceMeta`], one entry per graph node (including
/// composite unitigs pushed during resolution).
#[derive(Clone, Debug, Default)]
pub struct ProvenanceTable {
    metas: Vec<NodeProvenanceMeta>,
}

impl ProvenanceTable {
    /// Static provenance analysis over the whole graph: DCC and BCI, computed
    /// once per resolution in O(V + E) on the strand-directed walk graph.
    ///
    /// The walk graph has one vertex per (node, strand) orientation; every
    /// directed edge `u -> v` appears together with its reverse-complement
    /// twin `v.twin() -> u.twin()`, so forward and reverse walks read
    /// symmetric, independently computed fields (DSC-1..DSC-3).
    ///
    /// DCC: context sources are the graph's external entries — strand
    /// vertices with no incoming external edge. Self-loops are internal to
    /// the node's equivalence class and never introduce a context (AR-4).
    /// Source lineages propagate over the SCC condensation: simple bubbles
    /// diverging from a single ancestor reconverge with one context, while
    /// nodes with two or more independent external incoming lineages get
    /// `external_context_count == 2` (capped).
    ///
    /// BCI: bounded multi-source reverse search from every multi-context
    /// vertex, with `max_insert` budget and edge cost `seq_len - overlap`
    /// (matching the walk's overlap-adjusted anchor distance). Any node that
    /// leads into a multi-context confluence within the horizon is flagged.
    /// A long (>= 500 bp) single-context predecessor absorbs the contamination
    /// instead of propagating it: the region behind a clean unique anchor
    /// stays attributable, so the anchor itself is never flagged (Q8
    /// over-truncation avoidance; F5/F6 contracts).
    pub fn analyze(graph: &CdbgGraph, max_insert: usize) -> Self {
        let n = graph.nodes.len();
        let nv = n * 2;

        // Strand-directed CSR adjacency (skipping deleted nodes and AR-4
        // self-loops).
        let mut out_deg = vec![0u32; nv];
        let mut in_deg = vec![0u32; nv];
        let mut edges: Vec<(u32, u32, u8)> = Vec::new();
        for id in 0..n {
            if graph.nodes[id].is_deleted {
                continue;
            }
            for rc in [false, true] {
                let from = DirUnitigId::new(id as u32, rc);
                for e in graph.nodes[id].edges(rc) {
                    if e.target == from || e.target.id() as usize >= n {
                        continue;
                    }
                    if graph.nodes[e.target.id() as usize].is_deleted {
                        continue;
                    }
                    edges.push((from.0, e.target.0, e.overlap_k));
                    out_deg[from.0 as usize] += 1;
                    in_deg[e.target.0 as usize] += 1;
                }
            }
        }

        let build_csr = |deg: &[u32]| -> (Vec<u32>, Vec<(u32, u8)>) {
            let mut off = Vec::with_capacity(nv + 1);
            let mut acc = 0u32;
            off.push(0);
            for &d in deg {
                acc += d;
                off.push(acc);
            }
            (off, vec![(0u32, 0u8); acc as usize])
        };
        let (out_off, mut out_adj) = build_csr(&out_deg);
        let (in_off, mut in_adj) = build_csr(&in_deg);
        {
            let mut cursor = out_off[..nv].to_vec();
            for &(from, to, ov) in &edges {
                out_adj[cursor[from as usize] as usize] = (to, ov);
                cursor[from as usize] += 1;
            }
            let mut cursor = in_off[..nv].to_vec();
            for &(from, to, ov) in &edges {
                in_adj[cursor[to as usize] as usize] = (from, ov);
                cursor[to as usize] += 1;
            }
        }
        drop(edges);

        // Iterative Tarjan SCC over the walk graph. Tarjan emits SCC ids in
        // reverse topological order of the condensation, so a smaller id is
        // downstream: descending id order is source-first (topological).
        let mut scc = vec![u32::MAX; nv];
        let mut index = vec![u32::MAX; nv];
        let mut low = vec![0u32; nv];
        let mut on_stack = vec![false; nv];
        let mut vertex_stack: Vec<u32> = Vec::new();
        let mut call_stack: Vec<(u32, u32)> = Vec::new();
        let mut next_index = 0u32;
        let mut scc_count = 0u32;
        for start in 0..nv {
            if index[start] != u32::MAX {
                continue;
            }
            index[start] = next_index;
            low[start] = next_index;
            next_index += 1;
            vertex_stack.push(start as u32);
            on_stack[start] = true;
            call_stack.push((start as u32, out_off[start]));
            while let Some(&(v, next_child)) = call_stack.last() {
                if next_child < out_off[v as usize + 1] {
                    call_stack.last_mut().unwrap().1 += 1;
                    let (w, _) = out_adj[next_child as usize];
                    let wu = w as usize;
                    if index[wu] == u32::MAX {
                        index[wu] = next_index;
                        low[wu] = next_index;
                        next_index += 1;
                        vertex_stack.push(w);
                        on_stack[wu] = true;
                        call_stack.push((w, out_off[wu]));
                    } else if on_stack[wu] {
                        low[v as usize] = low[v as usize].min(index[wu]);
                    }
                } else {
                    call_stack.pop();
                    if let Some(&(parent, _)) = call_stack.last() {
                        low[parent as usize] = low[parent as usize].min(low[v as usize]);
                    }
                    if low[v as usize] == index[v as usize] {
                        loop {
                            let w = vertex_stack
                                .pop()
                                .expect("Tarjan vertex stack cannot be empty");
                            on_stack[w as usize] = false;
                            scc[w as usize] = scc_count;
                            if w == v {
                                break;
                            }
                        }
                        scc_count += 1;
                    }
                }
            }
        }

        // Astra DCC (§1 Precise distinction rule & §2 Shared-context trigger):
        // An SCC is a multi-context repeat confluence if and only if 2 or more
        // distinct upstream locus classes converge into it.
        // Certified same-locus bubble alternatives and error tips share the same
        // upstream origin and do NOT introduce shared context.
        let mut scc_nodes: Vec<SmallVec<[u32; 2]>> = vec![SmallVec::new(); scc_count as usize];
        for v in 0..nv {
            scc_nodes[scc[v] as usize].push(v as u32);
        }

        // Invariant F24-I08 & Audit Verification:
        // Condensation graph is a DAG where directed edges between distinct SCCs (u -> v)
        // must satisfy scc[u] > scc[v], guaranteeing that descending SCC iteration
        // (scc_count - 1 down to 0) strictly visits predecessors before successors.
        #[cfg(debug_assertions)]
        for from in 0..nv {
            let su = scc[from] as usize;
            for ei in out_off[from]..out_off[from + 1] {
                let (to, _) = out_adj[ei as usize];
                let sv = scc[to as usize] as usize;
                if su != sv {
                    debug_assert!(
                        su > sv,
                        "Condensation DAG topological ordering invariant violated: edge {} -> {} with scc {} -> {}",
                        from, to, su, sv
                    );
                }
            }
        }

        // Astra DCC (§1 Precise distinction rule & §2 Shared-context trigger, Sol Round 8):
        // Single-pass topological origin propagation in strictly O(V + E) time.
        // Replaces quadratic predecessor search: an unbranched predecessor chain
        // inherits its root origin in O(1) per node.
        let mut origin = vec![0u32; nv];
        let mut origin_dist = vec![0u32; nv];
        for s in (0..scc_count).rev() {
            let su = s as usize;
            for &v in &scc_nodes[su] {
                let vu = v as usize;
                let v_len = graph.nodes[vu >> 1].seq_len;
                let in_start = in_off[vu];
                let in_end = in_off[vu + 1];
                let deg = in_end - in_start;
                if v_len >= 500 || deg != 1 {
                    origin[vu] = v;
                    origin_dist[vu] = 0;
                    continue;
                }
                let (p, ov) = in_adj[in_start as usize];
                let pu = p as usize;
                let p_len = graph.nodes[pu >> 1].seq_len;
                let p_in_deg = in_off[pu + 1] - in_off[pu];
                if p_len >= 500 || p_in_deg != 1 {
                    origin[vu] = p;
                    origin_dist[vu] = (v_len as usize).saturating_sub(ov as usize) as u32;
                } else {
                    let step_dist = (v_len as usize).saturating_sub(ov as usize) as u32;
                    let tot_dist = origin_dist[pu].saturating_add(step_dist);
                    if tot_dist > 2000 {
                        origin[vu] = v;
                        origin_dist[vu] = 0;
                    } else {
                        origin[vu] = origin[pu];
                        origin_dist[vu] = tot_dist;
                    }
                }
            }
        }

        // Provenance Attribution Lattice (Astra §3.1 & Ruling):
        // 0 = Unknown (Empty/ungrounded)
        // 1 = Single(lineage_id)
        // 2 = Shared (Many / absorbing)
        // Evaluated across condensation DAG in topological order (0..scc_count).rev().
        // For an in-degree d junction, processing takes O(d) time with O(1) state per edge.
        #[derive(Clone, Copy, PartialEq, Eq)]
        enum LineageAcc {
            Empty,
            One(u32),
            Many,
        }

        let mut scc_ctx = vec![0u8; scc_count as usize];
        let mut scc_lineage = vec![u32::MAX; scc_count as usize];

        for s in (0..scc_count).rev() {
            let su = s as usize;
            let mut acc = LineageAcc::Empty;
            let mut has_external_in = false;

            for &v in &scc_nodes[su] {
                let vu = v as usize;
                for ei in in_off[vu]..in_off[vu + 1] {
                    let (p, _) = in_adj[ei as usize];
                    let sp = scc[p as usize] as usize;
                    if sp != su {
                        has_external_in = true;
                        let p_ctx = scc_ctx[sp];
                        if p_ctx >= 2 {
                            acc = LineageAcc::Many;
                            break;
                        } else if p_ctx == 1 {
                            let orig = origin[p as usize];
                            match acc {
                                LineageAcc::Empty => acc = LineageAcc::One(orig),
                                LineageAcc::One(prev) => {
                                    if prev != orig {
                                        acc = LineageAcc::Many;
                                        break;
                                    }
                                }
                                LineageAcc::Many => break,
                            }
                        } else {
                            // Predecessor is Unknown (0). Unknown cannot certify Single.
                            acc = LineageAcc::Many;
                            break;
                        }
                    }
                }
                if acc == LineageAcc::Many {
                    break;
                }
            }

            if !has_external_in {
                // Root component: if any node in SCC has length >= 500 bp, it defines its own unique Single lineage.
                // If all nodes are short without external entry, it has no witness -> Unknown (0).
                let long_node = scc_nodes[su].iter().find(|&&v| graph.nodes[(v as usize) >> 1].seq_len >= 500);
                if let Some(&root_v) = long_node {
                    acc = LineageAcc::One(root_v);
                } else {
                    acc = LineageAcc::Empty;
                }
            }

            match acc {
                LineageAcc::Empty => {
                    scc_ctx[su] = 0;
                    scc_lineage[su] = u32::MAX - 1;
                }
                LineageAcc::One(orig) => {
                    scc_ctx[su] = 1;
                    scc_lineage[su] = orig;
                }
                LineageAcc::Many => {
                    scc_ctx[su] = 2;
                    scc_lineage[su] = u32::MAX;
                }
            }

            for &v in &scc_nodes[su] {
                origin[v as usize] = scc_lineage[su];
            }
        }

        let mut ctx_count = vec![0u8; nv];
        for v in 0..nv {
            ctx_count[v] = scc_ctx[scc[v] as usize];
        }

        // BCI: linear reverse reach propagation from all multi-context vertices.
        // Sol Round 8/9: Topological DAG dynamic programming over reversed condensation graph.
        // Processes SCCs in reverse topological order (0..scc_count).
        // Cyclic components (scc_nodes > 1) are conservatively collapsed to minimum reach;
        // only inter-SCC edges across condensation DAG are relaxed, strictly O(V + E).
        let budget = max_insert as i64;
        let mut bci = vec![false; nv];
        if budget > 0 {
            let mut best = vec![u64::MAX; nv];
            for v in 0..nv {
                if ctx_count[v] >= 2 {
                    best[v] = 0;
                    bci[v] = true;
                }
            }

            for s in 0..scc_count {
                let su = s as usize;
                let is_cyclic = scc_nodes[su].len() > 1 || {
                    let v = scc_nodes[su][0] as usize;
                    out_adj[out_off[v] as usize..out_off[v + 1] as usize]
                        .iter()
                        .any(|&(to, _)| to == scc_nodes[su][0])
                };
                let comp_d = if is_cyclic {
                    let mut min_d = u64::MAX;
                    for &v in &scc_nodes[su] {
                        min_d = min_d.min(best[v as usize]);
                    }
                    if min_d <= budget as u64 {
                        for &v in &scc_nodes[su] {
                            best[v as usize] = min_d;
                            bci[v as usize] = true;
                        }
                    }
                    min_d
                } else {
                    best[scc_nodes[su][0] as usize]
                };

                for &v in &scc_nodes[su] {
                    let vu = v as usize;
                    let d = comp_d;
                    if d == u64::MAX {
                        continue;
                    }
                    let v_len = graph.nodes[vu >> 1].seq_len as i64;
                    for ei in in_off[vu]..in_off[vu + 1] {
                        let (p, ov) = in_adj[ei as usize];
                        let pu = p as usize;
                        if scc[pu] as usize == su {
                            // Intra-SCC edge: handled by conservative cyclic component collapse above.
                            continue;
                        }
                        let nd = d + (v_len - ov as i64).max(0) as u64;
                        if nd as i64 > budget {
                            continue;
                        }
                        // Anchor absorption (Q8): a long single-context
                        // predecessor is a valid unique anchor; the pooled zone
                        // ahead of it does not flag it or anything behind it.
                        // Unknown (0) does NOT absorb; only Single (1) absorbs.
                        if ctx_count[pu] == 1 && graph.nodes[pu >> 1].seq_len >= 500 {
                            continue;
                        }
                        if nd < best[pu] {
                            best[pu] = nd;
                            bci[pu] = true;
                        }
                    }
                }
            }
        }

        // Linear forward shadow reach propagation (O(V + E), Astra §1/§3.2-§3.3):
        // Sol Round 5/8/9: Topological DAG dynamic programming over condensation graph.
        // Guarantees shortest distance (maximum conservative shadow reach) in strictly linear time.
        // Cyclic components (scc_nodes > 1) are conservatively collapsed to maximum reach;
        // only inter-SCC edges across condensation DAG are relaxed, strictly O(V + E).
        let mut bci_downstream = vec![false; nv];
        let mut exit_isolated = vec![false; nv];
        let mut entry_shadow_reach = vec![0u32; nv];
        if budget > 0 {
            let mut best_fwd = vec![u64::MAX; nv];
            for v in 0..nv {
                if ctx_count[v] >= 2 {
                    best_fwd[v] = 0;
                    bci_downstream[v] = true;
                    entry_shadow_reach[v] = budget as u32;
                }
            }

            // Iterate SCC condensation in topological order (scc_count - 1 down to 0).
            // Each vertex is processed exactly once; each directed inter-SCC edge is relaxed exactly once.
            for s in (0..scc_count).rev() {
                let su = s as usize;
                let is_cyclic = scc_nodes[su].len() > 1 || {
                    let v = scc_nodes[su][0] as usize;
                    out_adj[out_off[v] as usize..out_off[v + 1] as usize]
                        .iter()
                        .any(|&(to, _)| to == scc_nodes[su][0])
                };
                let comp_d = if is_cyclic {
                    let mut min_d = u64::MAX;
                    for &v in &scc_nodes[su] {
                        min_d = min_d.min(best_fwd[v as usize]);
                    }
                    if min_d <= budget as u64 {
                        let reach = (budget - min_d as i64).max(0) as u32;
                        for &v in &scc_nodes[su] {
                            best_fwd[v as usize] = min_d;
                            bci_downstream[v as usize] = true;
                            entry_shadow_reach[v as usize] = entry_shadow_reach[v as usize].max(reach);
                        }
                    }
                    min_d
                } else {
                    best_fwd[scc_nodes[su][0] as usize]
                };

                for &v in &scc_nodes[su] {
                    let vu = v as usize;
                    let d = comp_d;
                    if d == u64::MAX {
                        continue;
                    }
                    let reach = (budget - d as i64).max(0) as u32;
                    entry_shadow_reach[vu] = entry_shadow_reach[vu].max(reach);
                    bci_downstream[vu] = true;

                    let v_len = graph.nodes[vu >> 1].seq_len as u64;
                    if ctx_count[vu] == 1 && v_len > (reach as u64 + max_insert as u64) {
                        // Forward anchor absorption (R-I, Astra §3.2-§3.3):
                        // 3' exit of v is certified physically isolated (> max_insert beyond reach b).
                        exit_isolated[vu] = true;
                        // Do not propagate entry shadow downstream of v's isolated exit.
                        continue;
                    }

                    for ei in out_off[vu]..out_off[vu + 1] {
                        let (w, ov) = out_adj[ei as usize];
                        let wu = w as usize;
                        if scc[wu] as usize == su {
                            // Intra-SCC edge: handled by conservative cyclic component collapse above.
                            continue;
                        }
                        let nd = d + (v_len as i64 - ov as i64).max(0) as u64;
                        if nd as i64 > budget {
                            continue;
                        }
                        if nd < best_fwd[wu] {
                            best_fwd[wu] = nd;
                            bci_downstream[wu] = true;
                        }
                    }
                }
            }
        }

        let mut metas = Vec::with_capacity(n);
        for id in 0..n {
            metas.push(NodeProvenanceMeta {
                external_context_count_fwd: ctx_count[id * 2],
                external_context_count_rev: ctx_count[id * 2 + 1],
                bci_fwd: bci[id * 2],
                bci_rev: bci[id * 2 + 1],
                bci_downstream_fwd: bci_downstream[id * 2],
                bci_downstream_rev: bci_downstream[id * 2 + 1],
                entry_shadow_reach_fwd: entry_shadow_reach[id * 2],
                entry_shadow_reach_rev: entry_shadow_reach[id * 2 + 1],
                exit_isolated_fwd: exit_isolated[id * 2],
                exit_isolated_rev: exit_isolated[id * 2 + 1],
                structural_in_degree_fwd: in_deg[id * 2].min(u16::MAX as u32) as u16,
                structural_in_degree_rev: in_deg[id * 2 + 1].min(u16::MAX as u32) as u16,
                structural_out_degree_fwd: out_deg[id * 2].min(u16::MAX as u32) as u16,
                structural_out_degree_rev: out_deg[id * 2 + 1].min(u16::MAX as u32) as u16,
                provenance_boundary: None,
            });
        }
        let mut count_0 = 0usize;
        let mut count_1 = 0usize;
        let mut count_2plus = 0usize;
        let mut bci_count = 0usize;
        for v in 0..nv {
            if ctx_count[v] == 0 {
                count_0 += 1;
            } else if ctx_count[v] == 1 {
                count_1 += 1;
            } else {
                count_2plus += 1;
            }
            if bci[v] {
                bci_count += 1;
            }
        }
        tracing::info!(
            "[ProvenanceTable] nv={}, scc_count={}, count_0={}, count_1={}, count_2plus={}, bci_nodes={}",
            nv, scc_count, count_0, count_1, count_2plus, bci_count
        );

        Self { metas }
    }

    #[inline]
    pub fn meta(&self, id: u32) -> &NodeProvenanceMeta {
        &self.metas[id as usize]
    }

    /// Provenance cleanliness of a directed node in its walk direction:
    /// exactly one established external lineage (DCC strict: count == 1,
    /// originated either from a 0-in-degree root or an entry source SCC in
    /// the condensation graph) and no bounded contamination ahead (BCI),
    /// with entry shadow either absent or physically absorbed at an isolated exit (R-I).
    /// Count >= 2 indicates a multi-context pooled locus.
    #[inline]
    pub fn anchor_clean(&self, d: DirUnitigId) -> bool {
        let m = self.meta(d.id());
        m.external_context_count(d.is_rc()) == 1
            && !m.bci(d.is_rc())
            && (!m.bci_downstream(d.is_rc()) || m.exit_isolated(d.is_rc()))
    }

    /// Appends metadata for a newly pushed composite unitig, keeping the
    /// table in lockstep with `graph.nodes`.
    pub fn push_composite(&mut self, meta: NodeProvenanceMeta) {
        self.metas.push(meta);
    }
}

/// ARG (Anchor Reset Gate, AR-1..AR-4): a traversed node may become the
/// walk's unique anchor only if it is long (>= 500 bp), single-context with
/// an established lineage in the walk direction (DCC count exactly 1 —
/// count 0 is unknown attribution, e.g. cyclic source SCCs, and must not
/// establish uniqueness), uncontaminated ahead (BCI), and the walk itself is
/// still clean — contamination is monotonically increasing and never resets
/// (AR-3).
fn anchor_reset_eligible(
    graph: &CdbgGraph,
    prov: &ProvenanceTable,
    node: DirUnitigId,
    contaminated: bool,
) -> bool {
    graph.nodes[node.id() as usize].seq_len >= 500 && prov.anchor_clean(node) && !contaminated
}

/// SI-1/SI-2 (Seed Integrity gate): a walk may originate from `seed` only if it
/// is a long unique locus (>= 500 bp), carries exactly one established incoming
/// context (external_context_count == 1), and is either free of upstream confluence
/// shadows (!bci_downstream) or has its exit physically isolated from all upstream
/// confluences (exit_isolated, R-I).
#[inline]
fn seed_anchor_eligible(
    graph: &CdbgGraph,
    prov: &ProvenanceTable,
    node: DirUnitigId,
) -> bool {
    let m = prov.meta(node.id());
    graph.nodes[node.id() as usize].seq_len >= 500
        && m.external_context_count(node.is_rc()) == 1
        && (!m.bci_downstream(node.is_rc()) || m.exit_isolated(node.is_rc()))
}

/// PAG withhold with no attributable anchor covering the fork: everything the
/// walk spliced past its last valid anchor is locus-ambiguous and must not be
/// committed. Truncates `ext` back to the anchor (`keep == 0` means the
/// anchor is the seed itself) and undoes this call's `used_unique`
/// registrations for the dropped steps so later seeds can still use them.
fn discard_unattributable_walk(
    ext: &mut Vec<PathStep>,
    keep: usize,
    used_unique: &mut HashSet<u32>,
    inserted_now: &mut Vec<u32>,
    cleanup_protected_nodes: &mut HashSet<u32>,
) {
    for step in ext.drain(keep..) {
        cleanup_protected_nodes.insert(step.dir.id());
    }
    let retained: SmallVec<[u32; 8]> = ext.iter().map(|s| s.dir.id()).collect();
    let mut i = 0;
    while i < inserted_now.len() {
        if retained.contains(&inserted_now[i]) {
            i += 1;
        } else {
            let id = inserted_now.swap_remove(i);
            used_unique.remove(&id);
        }
    }
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
    match matrix.get_evidence(u1, u2) {
        Some(evidences) => score_evidence_rows(evidences, dist, path_separation, None, 0, 0),
        None => (0, 0.0),
    }
}

/// Kernel scoring over evidence rows, optionally restricted by a composite
/// unitig's provenance boundary (PB / MR-2) and 5' entry-shadow reach (Sol R1-F2, R2-F1).
/// `pooled_suffix_exit` is the PB boundary in exit coordinates (`u1_len - boundary`).
/// `entry_shadow_reach` is the maximum 5' entry reach b from upstream confluences.
/// The read's start coordinate p = u1_len - dist_u1_exit must lie strictly beyond the
/// shadow interval [0, entry_shadow_reach]: p > entry_shadow_reach.
fn score_evidence_rows(
    evidences: &[ReadPairEvidence],
    dist: &InsertSizeDist,
    path_separation: i64,
    pooled_suffix_exit: Option<u32>,
    u1_len: u32,
    entry_shadow_reach: u32,
) -> (u32, f64) {
    let mut valid_count = 0u32;
    let mut total_weight = 0.0f64;
    let mut seen_pair_ids = HashSet::with_capacity(evidences.len());

    let mad = if dist.mad > 1.0 { dist.mad } else { 1.0 };
    let norm_factor = 1.4826 * mad;

    for ev in evidences {
        if let Some(suffix_exit) = pooled_suffix_exit {
            if (ev.dist_u1_exit as u64) <= suffix_exit as u64 {
                continue;
            }
        }
        if entry_shadow_reach > 0 {
            // Coordinate-faithful start p = u1_len - dist_u1_exit (Sol R2-F1).
            // Under Contract §3.2 & Astra ruling: p must strictly satisfy p > b + H,
            // where b = entry_shadow_reach and H = dist.max_insert.
            let p = u1_len.saturating_sub(ev.dist_u1_exit);
            if (p as u64) <= (entry_shadow_reach as u64 + dist.max_insert as u64) {
                continue;
            }
        }
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

/// PB-aware scoring used by the resolver: converts the evidence-source
/// composite's provenance boundary (a 5'-offset) into the exit-coordinate
/// cutoff `u1_len - boundary` before kernel scoring (PAG attribution, §3.6
/// rule 3). Also enforces multi-context zero-weight (Sol R1-F3) and
/// entry-shadow reach exclusion (Sol R1-F2, R2-F1).
pub fn score_candidate_path_prov(
    matrix: &PairLinkMatrix,
    graph: &CdbgGraph,
    prov: &ProvenanceTable,
    dist: &InsertSizeDist,
    u1: DirUnitigId,
    u2: DirUnitigId,
    path_separation: i64,
) -> (u32, f64) {
    let m = prov.meta(u1.id());
    // Astra Ruling: Reject Unknown (0) and Shared (>=2) in deciding-source paths.
    // Only Single (1) may decide forks.
    if m.external_context_count(u1.is_rc()) != 1 {
        return (0, 0.0);
    }
    let u1_len = graph.nodes[u1.id() as usize].seq_len;
    let boundary = m
        .provenance_boundary
        .map(|b| u1_len.saturating_sub(b));
    let reach = m.entry_shadow_reach(u1.is_rc());
    match matrix.get_evidence(u1, u2) {
        Some(evidences) => score_evidence_rows(evidences, dist, path_separation, boundary, u1_len, reach),
        None => (0, 0.0),
    }
}

const MAX_CANDIDATE_DEPTH: usize = 4;
const MAX_CANDIDATE_PATHS: usize = 32;
const MAX_DFS_STATES: usize = 256;
const MAX_EXTENSION_STEPS: usize = 64;
const MAX_QUEUE_CAPACITY: usize = 1024;

#[derive(Clone, Debug)]
struct CandidatePath {
    steps: SmallVec<[(DirUnitigId, GraphEdge); 4]>,
    /// Anchor-to-node separation for each step (parallel to `steps`), used
    /// by F6b differential scoring over path-local evidence.
    step_seps: SmallVec<[i64; 4]>,
    endpoint: DirUnitigId,
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

/// Bounded multi-hop candidate path search from `fork_node` to endpoints
/// supported by read pairs.
///
/// Candidate scoring is keyed at `evidence_anchor` with separation offset
/// `sep_offset` (PAG): for an uncontaminated walk the anchor is the walk head
/// `curr_dir` with offset 0 — the pre-F23 semantics; for a walk in a
/// shared/contaminated zone the anchor is the walk's unique anchor and the
/// offset is the overlap-adjusted anchor-to-head distance, so only pairs
/// attributable to the unique locus can decide the fork.
#[allow(clippy::too_many_arguments)]
fn find_candidate_paths(
    graph: &CdbgGraph,
    matrix: &PairLinkMatrix,
    dist: &InsertSizeDist,
    prov: &ProvenanceTable,
    evidence_anchor: DirUnitigId,
    sep_offset: i64,
    fork_node: DirUnitigId,
    init_sep: i64,
    used_unique: &HashSet<u32>,
    consumed_nodes: &HashSet<u32>,
    ext_nodes: &HashSet<u32>,
) -> (Vec<CandidatePath>, SearchStatus) {
    let max_separation = dist.max_insert as i64;
    let beta = separation_beta(matrix, evidence_anchor);
    let mut candidate_paths = Vec::new();
    let mut search_status = SearchStatus::Complete;
    let mut states_visited = 0usize;

    let mut queue = VecDeque::new();
    for target in distinct_outgoing(graph, fork_node) {
        let mut visited: SmallVec<[u32; 8]> = SmallVec::new();
        visited.push(evidence_anchor.id());
        visited.push(fork_node.id());
        match classify_hop(
            graph,
            fork_node,
            target,
            init_sep,
            &visited,
            0,
            false,
            evidence_anchor,
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
                let mut step_seps: SmallVec<[i64; 4]> = SmallVec::new();
                steps.push((target, edge));
                step_seps.push(sep);
                queue.push_back((target, sep, steps, step_seps, visited, shadow));
            }
            HopClass::Unresolved => {
                search_status = SearchStatus::Truncated;
                break;
            }
            HopClass::Probe { .. } | HopClass::SafePrune | HopClass::Inadmissible => {}
        }
    }

    'bfs: while let Some((u, sep_u, steps, step_seps, mut visited, through_consumed)) =
        queue.pop_front()
    {
        states_visited += 1;
        if states_visited > MAX_DFS_STATES {
            search_status = SearchStatus::Truncated;
            break;
        }

        visited.push(u.id());

        // A return to the anchor is evidence about this anchor, not an
        // extension endpoint. In-window support makes the decision unresolved;
        // either way the walk continues so a later exit can still be scored.
        if u.id() == evidence_anchor.id() {
            let (support_count, _) = score_candidate_path_prov(
                matrix,
                graph,
                prov,
                dist,
                evidence_anchor,
                u,
                sep_offset + sep_u,
            );
            if support_count > 0 {
                search_status = SearchStatus::Truncated;
            }
        } else {
            let (support_count, total_weight) = score_candidate_path_prov(
                matrix,
                graph,
                prov,
                dist,
                evidence_anchor,
                u,
                sep_offset + sep_u,
            );
            if support_count > 0 {
                candidate_paths.push(CandidatePath {
                    steps: steps.clone(),
                    step_seps: step_seps.clone(),
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
                    evidence_anchor,
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
                            prov,
                            evidence_anchor,
                            sep_offset,
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
                        let mut next_seps = step_seps.clone();
                        next_seps.push(sep);
                        queue.push_back((
                            neighbor,
                            sep,
                            next_steps,
                            next_seps,
                            visited.clone(),
                            shadow,
                        ));
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
            evidence_anchor,
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
#[allow(clippy::too_many_arguments)]
fn pruned_repeat_hides_supported_endpoint(
    graph: &CdbgGraph,
    matrix: &PairLinkMatrix,
    dist: &InsertSizeDist,
    prov: &ProvenanceTable,
    evidence_anchor: DirUnitigId,
    sep_offset: i64,
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
        evidence_anchor,
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

        let (count, _) = score_candidate_path_prov(
            matrix,
            graph,
            prov,
            dist,
            evidence_anchor,
            node,
            sep_offset + sep,
        );
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
                evidence_anchor,
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
                evidence_anchor,
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

/// Per-pair kernel weights of the in-window evidence rows keyed at
/// `evidence_anchor` about `node` at `path_separation` — the same window,
/// kernel, and PB restriction (MR-2) as [`score_evidence_rows`]: when the
/// anchor is a composite carrying a provenance boundary, rows keyed inside
/// its pooled-locus suffix (`dist_u1_exit <= pooled_suffix_exit`, with
/// `pooled_suffix_exit` already in exit coordinates as `u1_len - boundary`)
/// are inadmissible as attributable evidence and are dropped before kernel
/// weighting.
fn in_window_pair_weights(
    matrix: &PairLinkMatrix,
    dist: &InsertSizeDist,
    evidence_anchor: DirUnitigId,
    node: DirUnitigId,
    path_separation: i64,
    pooled_suffix_exit: Option<u32>,
    anchor_len: u32,
    entry_shadow_reach: u32,
) -> SmallVec<[(u32, f64); 8]> {
    let mut out: SmallVec<[(u32, f64); 8]> = SmallVec::new();
    let Some(evidences) = matrix.get_evidence(evidence_anchor, node) else {
        return out;
    };
    let mad = if dist.mad > 1.0 { dist.mad } else { 1.0 };
    let norm_factor = 1.4826 * mad;
    for ev in evidences {
        if let Some(suffix_exit) = pooled_suffix_exit {
            if (ev.dist_u1_exit as u64) <= suffix_exit as u64 {
                continue;
            }
        }
        if entry_shadow_reach > 0 {
            // Coordinate-faithful start p = anchor_len - dist_u1_exit (Sol R2-F1).
            // Under Contract §3.2 & Astra ruling: p must strictly satisfy p > b + H,
            // where b = entry_shadow_reach and H = dist.max_insert.
            let p = anchor_len.saturating_sub(ev.dist_u1_exit);
            if (p as u64) <= (entry_shadow_reach as u64 + dist.max_insert as u64) {
                continue;
            }
        }
        let span = ev.dist_u1_exit as i64 + path_separation + ev.dist_u2_entry as i64;
        if span < dist.min_insert as i64 || span > dist.max_insert as i64 {
            continue;
        }
        let z = (span as f64 - dist.median).abs() / norm_factor;
        let w = (-0.5 * z * z).exp();
        match out.iter_mut().find(|(pid, _)| *pid == ev.pair_id) {
            Some(slot) => {
                if w > slot.1 {
                    slot.1 = w;
                }
            }
            None => out.push((ev.pair_id, w)),
        }
    }
    out
}

/// Admissible pair weights over one candidate route: every node the route
/// enters — each step at its own anchor separation, plus the endpoint —
/// counted once per pair at its best in-window weight (F6b). A direct pair
/// link into an intermediate node (e.g. the S -> B1 allele bridge across a
/// same-length SNP bubble) is thereby part of the route's own support.
/// Rows are PB-restricted at `evidence_anchor` (MR-2): when the anchor is a
/// composite carrying a provenance boundary, rows keyed inside its
/// pooled-locus suffix are inadmissible attributable evidence.
/// Sol R1-F3: An evidence anchor that is multi-context has zero deciding weight.
/// Sol R1-F2, R2-F1: Entry-shadow reach excludes rows originating within [0, reach].
fn candidate_pair_weights(
    matrix: &PairLinkMatrix,
    graph: &CdbgGraph,
    prov: &ProvenanceTable,
    dist: &InsertSizeDist,
    evidence_anchor: DirUnitigId,
    sep_offset: i64,
    candidate: &CandidatePath,
) -> SmallVec<[(u32, f64); 8]> {
    let m = prov.meta(evidence_anchor.id());
    // Astra Ruling: Reject Unknown (0) and Shared (>=2) in deciding-source paths.
    // Only Single (1) may decide forks.
    if m.external_context_count(evidence_anchor.is_rc()) != 1 {
        return SmallVec::new();
    }
    let anc_len = graph.nodes[evidence_anchor.id() as usize].seq_len;
    let pooled_suffix_exit = m
        .provenance_boundary
        .map(|b| anc_len.saturating_sub(b));
    let reach = m.entry_shadow_reach(evidence_anchor.is_rc());
    let mut acc: SmallVec<[(u32, f64); 8]> = SmallVec::new();
    let mut absorb = |rows: SmallVec<[(u32, f64); 8]>| {
        for (pid, w) in rows {
            match acc.iter_mut().find(|(p, _)| *p == pid) {
                Some(slot) => {
                    if w > slot.1 {
                        slot.1 = w;
                    }
                }
                None => acc.push((pid, w)),
            }
        }
    };
    absorb(in_window_pair_weights(
        matrix,
        dist,
        evidence_anchor,
        candidate.endpoint,
        sep_offset + candidate.separation,
        pooled_suffix_exit,
        anc_len,
        reach,
    ));
    for (&(node, _), &sep) in candidate.steps.iter().zip(candidate.step_seps.iter()) {
        absorb(in_window_pair_weights(
            matrix,
            dist,
            evidence_anchor,
            node,
            sep_offset + sep,
            pooled_suffix_exit,
            anc_len,
            reach,
        ));
    }
    acc
}

/// F6b differential rivals test. Two rival routes through a same-length
/// bubble reconverge on a shared endpoint: pairs that bridge the whole
/// bubble and land in the shared node are in-window for both arms at the
/// SAME kernel weight and carry no evidence about which arm the walk should
/// take, whatever their count. The rivals ratio therefore compares each
/// candidate's *differential* support over every node of its own route: a
/// pair's differential contribution is its weight surplus over the same
/// pair's best admissible weight in the rival route (zero when the rival
/// explains it at least as well), so a direct pair link into one arm (e.g.
/// S -> B1 vs S -> B2) decides the fork cleanly, while routes of different
/// length still separate through the span shift of their shared pairs.
/// When the candidates share no admissible pairs this reduces to the
/// pre-F23 weight ratio.
fn differential_rivals_agree(
    matrix: &PairLinkMatrix,
    graph: &CdbgGraph,
    prov: &ProvenanceTable,
    dist: &InsertSizeDist,
    evidence_anchor: DirUnitigId,
    sep_offset: i64,
    best: &CandidatePath,
    second: &CandidatePath,
    min_support: u32,
    decision_threshold: f64,
) -> bool {
    let best_pairs = candidate_pair_weights(
        matrix, graph, prov, dist, evidence_anchor, sep_offset, best,
    );
    let second_pairs = candidate_pair_weights(
        matrix, graph, prov, dist, evidence_anchor, sep_offset, second,
    );
    let rival_weight = |pairs: &SmallVec<[(u32, f64); 8]>, pid: u32| {
        pairs
            .iter()
            .find(|&&(p, _)| p == pid)
            .map(|&(_, w)| w)
            .unwrap_or(0.0)
    };

    let mut best_exclusive = (0u32, 0.0f64);
    for &(pid, w) in &best_pairs {
        let rival_w = rival_weight(&second_pairs, pid);
        if w > rival_w {
            best_exclusive.0 += 1;
            best_exclusive.1 += w - rival_w;
        }
    }
    let mut second_exclusive = (0u32, 0.0f64);
    for &(pid, w) in &second_pairs {
        let best_w = rival_weight(&best_pairs, pid);
        if w > best_w {
            second_exclusive.0 += 1;
            second_exclusive.1 += w - best_w;
        }
    }

    if second_exclusive.1 < 1e-9 {
        // The runner-up is explained at least as well as the winner by every
        // pair it has: only the winner's differential pairs may decide, and
        // they must stand on their own — count >= min_support and the
        // single-candidate average-weight quality bar (Z <= 2.45) against
        // tail outliers.
        best_exclusive.0 >= min_support.max(1) && best_exclusive.1 / best_exclusive.0 as f64 >= 0.05
    } else {
        best_exclusive.1 / second_exclusive.1 >= decision_threshold
    }
}

/// Two candidate paths are compatible (the same route) if one is a prefix
/// of the other — they agree on every step up to the shorter path's length.
/// Compatible paths do not compete with each other at a fork; only paths
/// that diverge (choose different outgoing branches or different downstream
/// exits) are rival alternatives.
fn are_divergent_rivals(a: &CandidatePath, b: &CandidatePath) -> bool {
    let min_len = a.steps.len().min(b.steps.len());
    !a.steps[..min_len]
        .iter()
        .map(|s| s.0)
        .eq(b.steps[..min_len].iter().map(|s| s.0))
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
///
/// F23 provenance walk state (ARG + PAG):
/// - `unique_anchor` is the walk's last ARG-eligible anchor; the seed itself
///   must be ARG-eligible or the walk is abandoned before it starts (SI-2:
///   an anchor-less walk must never commit, not even linear paths from a
///   shared seed such as pooled L in A,B -> L -> F -> T).
/// - `contaminated` is set when any traversed node is a shared locus
///   (multi-context) or BCI-flagged, and is monotonically increasing (AR-3).
/// - At a fork in a contaminated zone, evidence keyed at the walk head has
///   zero deciding weight: candidates are scored from `unique_anchor` when it
///   covers the fork within `max_insert`, otherwise the fork is withheld and
///   the unattributable portion of the walk is rolled back.
fn extend_one_direction(
    mut curr_dir: DirUnitigId,
    graph: &CdbgGraph,
    matrix: &PairLinkMatrix,
    dist: &InsertSizeDist,
    params: ExspanderParams,
    prov: &ProvenanceTable,
    used_unique: &mut HashSet<u32>,
    stats: &mut JunctionStats,
    cleanup_protected_nodes: &mut HashSet<u32>,
    consumed_nodes: &HashSet<u32>,
) -> Vec<PathStep> {
    let mut ext: Vec<PathStep> = Vec::new();

    // F23 horizon guard (narrowed Option A): track the walk's most recent long
    // (>= 500 bp) node and the overlap-adjusted distance from its 3' end. When
    // the walk stands beyond that anchor's insert horizon at a branching
    // repeat, anchored evidence is arithmetically impossible and the only
    // remaining links are repeat-internal ones with no locus specificity
    // (runtime decision DEC_000022, R1-001 contig_16 chimeric join): the fork
    // is withheld instead of decided. Within the horizon, forks keep the
    // pre-F23 curr_dir-anchored semantics unchanged (accepted T26/T27/T30/T34
    // contracts). Entry nodes are always seeds (>= 500 bp).
    //
    // SI-1/SI-2 (seed initialization): a seed anchors the walk only if it is
    // ARG-eligible in the walk direction — single-context with an
    // established lineage (DCC count exactly 1) and BCI-clean. A shared long
    // seed (e.g. pooled L in A,B -> L -> F -> T) must not commit anything,
    // not even unbranched linear splices: an anchor-less walk could commit
    // linear paths from the shared seed while relying on every later fork to
    // withhold. The walk is therefore abandoned before it starts (SI-2). A
    // seed with no outgoing edge in this direction had nothing to withhold
    // (pre-SI-2 walks terminated immediately there), so the rejection is
    // counted only when an unanchored walk could otherwise have committed.
    // SI-1/SI-2 / F24 R-S:
    // A seed anchors the walk (unique_anchor = Some(curr_dir)) only if it is
    // ARG-eligible in the walk direction (single-context, established lineage,
    // and BCI-clean).
    // An unanchored seed (shared or shadowed) cannot decide forks (unique_anchor = None),
    // but may traverse verified unbranched 1-out/1-in structural corridors (R-S, Astra §1/§3.5).
    let (mut unique_anchor, mut contaminated) = if seed_anchor_eligible(graph, prov, curr_dir) {
        (Some(curr_dir), false)
    } else {
        (None, true)
    };
    // Number of ext steps (inclusive) that end at the current anchor;
    // 0 means the anchor is the seed itself.
    let mut anchor_steps = 0usize;
    let mut anchor_gap = 0i64;
    // Node ids this call newly registered in `used_unique`, so a PAG rollback
    // can undo exactly its own registrations.
    let mut inserted_now: Vec<u32> = Vec::new();

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

            let target_in_deg = graph.nodes[target.id() as usize]
                .edges(target.twin().is_rc())
                .len();
            if unique_anchor.is_none() {
                let curr_m = prov.meta(curr_dir.id());
                let target_m = prov.meta(target.id());
                if curr_m.structural_out_degree(curr_dir.is_rc()) != 1
                    || target_m.structural_in_degree(target.is_rc()) != 1
                    || outgoing.len() != 1
                    || target_in_deg != 1
                {
                    // Structural corridor traversal (R-S, Astra §1/§3.5, Sol R7 finding):
                    // Unanchored walk requires strict 1-out/1-in corridor both structurally in the
                    // immutable analysis snapshot and dynamically in the live graph.
                    cleanup_protected_nodes.insert(target.id());
                    break;
                }
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
                // Target is a branching repeat node! Evaluate complete candidate paths exiting the repeat.
                let fork_dist =
                    anchor_gap + target_node.seq_len as i64 - edge_to_target.overlap_k as i64;

                // PAG (contaminated zone first): evidence keyed at the walk
                // head is pooled and has zero deciding weight; only pairs
                // attributable to the walk's unique anchor may decide.
                let (evidence_anchor, sep_offset) = if contaminated {
                    match unique_anchor {
                        Some(anchor) if anchor_gap <= dist.max_insert as i64 => {
                            (anchor, anchor_gap)
                        }
                        _ => {
                            // No attributable anchor covers this fork: withhold
                            // and roll the walk back to its anchor (§3.10.1/4).
                            stats.withheld_no_unique_anchor += 1;
                            stats.withheld_no_attributable_anchor += 1;
                            discard_unattributable_walk(
                                &mut ext,
                                anchor_steps,
                                used_unique,
                                &mut inserted_now,
                                cleanup_protected_nodes,
                            );
                            cleanup_protected_nodes.insert(target.id());
                            for &(cand, _) in &next_continuations {
                                cleanup_protected_nodes.insert(cand.id());
                            }
                            break;
                        }
                    }
                } else {
                    // F23 horizon guard: an uncontaminated fork beyond the
                    // last unique anchor's insert horizon cannot be decided by
                    // locus-specific evidence — withhold (overlap-adjusted
                    // anchor-to-fork distance, Codex-reviewed formulas). The
                    // committed prefix stands (CS-3).
                    if unique_anchor != Some(curr_dir) && fork_dist > dist.max_insert as i64 {
                        stats.withheld_no_unique_anchor += 1;
                        stats.withheld_horizon += 1;
                        cleanup_protected_nodes.insert(target.id());
                        for &(cand, _) in &next_continuations {
                            cleanup_protected_nodes.insert(cand.id());
                        }
                        break;
                    }
                    (curr_dir, 0)
                };

                let ext_nodes: HashSet<u32> = ext.iter().map(|s: &PathStep| s.dir.id()).collect();
                let init_sep = target_node.seq_len as i64 - edge_to_target.overlap_k as i64;
                let (candidate_paths, search_status) = find_candidate_paths(
                    graph,
                    matrix,
                    dist,
                    prov,
                    evidence_anchor,
                    sep_offset,
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

                let (mut shadow_paths, mut candidate_paths): (Vec<_>, Vec<_>) = candidate_paths
                    .into_iter()
                    .partition(|candidate| candidate.through_consumed);

                if candidate_paths.is_empty() {
                    if contaminated {
                        // §3.10.2: branch with no admissible pairs under PAG.
                        stats.withheld_no_unique_anchor += 1;
                        stats.withheld_pag_empty += 1;
                    }
                    cleanup_protected_nodes.insert(target.id());
                    for &(cand, _) in &next_continuations {
                        cleanup_protected_nodes.insert(cand.id());
                    }
                    break;
                }

                // Route-wide support: the DFS admits a route on endpoint
                // evidence alone, but the fork must rank and compare whole
                // routes (F6b) — every node the route enters contributes its
                // PB-restricted in-window pairs. Shadows get the same
                // recomputation so `shadow_blocks_decision` compares like
                // for like.
                for cand in &mut candidate_paths {
                    let pairs = candidate_pair_weights(
                        matrix, graph, prov, dist, evidence_anchor, sep_offset, cand,
                    );
                    cand.support_count = pairs.len() as u32;
                    cand.total_weight = pairs.iter().map(|&(_, w)| w).sum();
                }
                for shadow in &mut shadow_paths {
                    let pairs = candidate_pair_weights(
                        matrix, graph, prov, dist, evidence_anchor, sep_offset, shadow,
                    );
                    shadow.support_count = pairs.len() as u32;
                    shadow.total_weight = pairs.iter().map(|&(_, w)| w).sum();
                }

                candidate_paths.sort_by(|a, b| {
                    b.total_weight
                        .partial_cmp(&a.total_weight)
                        .unwrap_or(std::cmp::Ordering::Equal)
                        .then_with(|| b.support_count.cmp(&a.support_count))
                });

                let best = &candidate_paths[0];

                // The winner must beat every DIVERGENT rival route, not
                // only the runner-up. A prefix-compatible candidate agrees
                // with `best` on every step up to its own length (the same
                // route, just a shorter exit), so it is not an alternative
                // at this fork and does not compete; only divergent routes
                // — a different outgoing branch or a different downstream
                // exit — can veto the decision.
                let divergent_rivals: Vec<&CandidatePath> = candidate_paths[1..]
                    .iter()
                    .filter(|rival| are_divergent_rivals(best, rival))
                    .collect();

                let rivals_agree = if best.support_count < params.min_support {
                    false
                } else if divergent_rivals.is_empty() {
                    // Every other candidate is prefix-compatible with `best`:
                    // only the single-candidate quality bar applies (average
                    // weight >= 0.05, Z <= 2.45, rejects tail outliers).
                    (best.total_weight / best.support_count as f64) >= 0.05
                } else if prov.meta(evidence_anchor.id()).bci_downstream(evidence_anchor.is_rc())
                    && divergent_rivals.iter().any(|r| r.support_count >= params.min_support)
                {
                    // F24 Invariant F24-I04 & N01: An evidence anchor carrying an upstream
                    // confluence shadow (bci_downstream) cannot decide between divergent rivals
                    // when multiple exits carry substantial evidence: competing evidence indicates
                    // pooled/homologous repeat reads rather than authentic locus-specific support.
                    stats.withheld_no_unique_anchor += 1;
                    stats.withheld_no_attributable_anchor += 1;
                    false
                } else {
                    // F6b: pairs admissible for both rivals (e.g. rows that
                    // bridge a whole same-length bubble into the shared
                    // endpoint) carry no differential evidence, so the ratio
                    // compares exclusive support only. A rival with no
                    // exclusive pairs degenerates to the single-candidate
                    // quality bar inside `differential_rivals_agree`.
                    divergent_rivals.iter().all(|rival| {
                        differential_rivals_agree(
                            matrix,
                            graph,
                            prov,
                            dist,
                            evidence_anchor,
                            sep_offset,
                            best,
                            rival,
                            params.min_support,
                            params.decision_threshold,
                        )
                    })
                };
                let is_decision_valid = rivals_agree
                    && !shadow_blocks_decision(best, &shadow_paths, params.decision_threshold);

                let all_steps_valid = best.steps.iter().all(|(step_dir, _)| {
                    let step_node = &graph.nodes[step_dir.id() as usize];
                    !step_node.is_deleted
                        && !(step_node.seq_len >= 500 && used_unique.contains(&step_dir.id()))
                });

                if is_decision_valid && all_steps_valid {
                    // Commit the fork: push the fork node, then the winning
                    // candidate steps, updating walk contamination
                    // (monotonic, AR-3) and the ARG anchor.
                    if prov.meta(target.id()).contaminating(target.is_rc()) {
                        contaminated = true;
                    }
                    if target_node.seq_len >= 500 && used_unique.insert(target.id()) {
                        inserted_now.push(target.id());
                    }
                    ext.push(PathStep {
                        dir: target,
                        incoming_edge: Some(edge_to_target),
                        overlap: edge_to_target.overlap_k as usize,
                        support: edge_to_target.support as u32,
                    });
                    anchor_gap += target_node.seq_len as i64 - edge_to_target.overlap_k as i64;
                    if anchor_reset_eligible(graph, prov, target, contaminated) {
                        unique_anchor = Some(target);
                        anchor_steps = ext.len();
                        anchor_gap = 0;
                    }

                    for (step_dir, step_edge) in &best.steps {
                        let step_node = &graph.nodes[step_dir.id() as usize];
                        if prov.meta(step_dir.id()).contaminating(step_dir.is_rc()) {
                            contaminated = true;
                        }
                        if step_node.seq_len >= 500 && used_unique.insert(step_dir.id()) {
                            inserted_now.push(step_dir.id());
                        }
                        anchor_gap += step_node.seq_len as i64 - step_edge.overlap_k as i64;
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
                        if anchor_reset_eligible(graph, prov, *step_dir, contaminated) {
                            unique_anchor = Some(*step_dir);
                            anchor_steps = ext.len();
                            anchor_gap = 0;
                        }
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
                    if used_unique.insert(target.id()) {
                        inserted_now.push(target.id());
                    }
                }
                if prov.meta(target.id()).contaminating(target.is_rc()) {
                    contaminated = true;
                }
                ext.push(PathStep {
                    dir: target,
                    incoming_edge: Some(edge_to_target),
                    overlap: edge_to_target.overlap_k as usize,
                    support: edge_to_target.support as u32,
                });
                anchor_gap += target_node.seq_len as i64 - edge_to_target.overlap_k as i64;
                if anchor_reset_eligible(graph, prov, target, contaminated) {
                    unique_anchor = Some(target);
                    anchor_steps = ext.len();
                    anchor_gap = 0;
                } else if unique_anchor.is_none() {
                    // Structural corridor traversal (R-S, Astra §1/§3.5):
                    // Advance anchor_steps so verified 1-out/1-in structural steps are preserved.
                    anchor_steps = ext.len();
                }
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

        // Case B: multiple outgoing edges directly from curr_dir (the walk
        // stands on the fork). PAG (contaminated zone) first, then the F23
        // horizon guard (Case B formula: fork distance == anchor_gap).
        let (evidence_anchor, sep_offset) = if contaminated {
            match unique_anchor {
                Some(anchor) if anchor_gap <= dist.max_insert as i64 => (anchor, anchor_gap),
                _ => {
                    // No attributable anchor covers this fork: withhold and
                    // roll the walk back to its anchor (§3.10.1/4).
                    stats.withheld_no_unique_anchor += 1;
                    stats.withheld_no_attributable_anchor += 1;
                    discard_unattributable_walk(
                        &mut ext,
                        anchor_steps,
                        used_unique,
                        &mut inserted_now,
                        cleanup_protected_nodes,
                    );
                    for &(cand, _) in &valid_continuations {
                        cleanup_protected_nodes.insert(cand.id());
                    }
                    break;
                }
            }
        } else {
            if unique_anchor != Some(curr_dir) && anchor_gap > dist.max_insert as i64 {
                stats.withheld_no_unique_anchor += 1;
                stats.withheld_horizon += 1;
                for &(cand, _) in &valid_continuations {
                    cleanup_protected_nodes.insert(cand.id());
                }
                break;
            }
            (curr_dir, 0)
        };
        let ext_nodes: HashSet<u32> = ext.iter().map(|s: &PathStep| s.dir.id()).collect();
        let (candidate_paths, search_status) = find_candidate_paths(
            graph,
            matrix,
            dist,
            prov,
            evidence_anchor,
            sep_offset,
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

        let (mut shadow_paths, mut candidate_paths): (Vec<_>, Vec<_>) = candidate_paths
            .into_iter()
            .partition(|candidate| candidate.through_consumed);

        if candidate_paths.is_empty() {
            if contaminated {
                // §3.10.2: branch with no admissible pairs under PAG.
                stats.withheld_no_unique_anchor += 1;
                stats.withheld_pag_empty += 1;
            }
            for &(cand, _) in &valid_continuations {
                cleanup_protected_nodes.insert(cand.id());
            }
            break;
        }

        // Route-wide support: the DFS admits a route on endpoint evidence
        // alone, but the fork must rank and compare whole routes (F6b) —
        // every node the route enters contributes its PB-restricted
        // in-window pairs. Shadows get the same recomputation so
        // `shadow_blocks_decision` compares like for like.
        for cand in &mut candidate_paths {
            let pairs =
                candidate_pair_weights(matrix, graph, prov, dist, evidence_anchor, sep_offset, cand);
            cand.support_count = pairs.len() as u32;
            cand.total_weight = pairs.iter().map(|&(_, w)| w).sum();
        }
        for shadow in &mut shadow_paths {
            let pairs = candidate_pair_weights(
                matrix, graph, prov, dist, evidence_anchor, sep_offset, shadow,
            );
            shadow.support_count = pairs.len() as u32;
            shadow.total_weight = pairs.iter().map(|&(_, w)| w).sum();
        }

        candidate_paths.sort_by(|a, b| {
            b.total_weight
                .partial_cmp(&a.total_weight)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| b.support_count.cmp(&a.support_count))
        });

        let best = &candidate_paths[0];

        // The winner must beat every DIVERGENT rival route, not only the
        // runner-up. A prefix-compatible candidate agrees with `best` on
        // every step up to its own length (the same route, just a shorter
        // exit), so it is not an alternative at this fork and does not
        // compete; only divergent routes — a different outgoing branch or a
        // different downstream exit — can veto the decision.
        let divergent_rivals: Vec<&CandidatePath> = candidate_paths[1..]
            .iter()
            .filter(|rival| are_divergent_rivals(best, rival))
            .collect();

        let rivals_agree = if best.support_count < params.min_support {
            false
        } else if divergent_rivals.is_empty() {
            // Every other candidate is prefix-compatible with `best`: only
            // the single-candidate quality bar applies (average weight
            // >= 0.05, Z <= 2.45, rejects tail outliers).
            (best.total_weight / best.support_count as f64) >= 0.05
        } else if prov.meta(evidence_anchor.id()).bci_downstream(evidence_anchor.is_rc())
            && divergent_rivals.iter().any(|r| r.support_count >= params.min_support)
        {
            // F24 Invariant F24-I04 & N01: An evidence anchor carrying an upstream
            // confluence shadow (bci_downstream) cannot decide between divergent rivals
            // when multiple exits carry substantial evidence: competing evidence indicates
            // pooled/homologous repeat reads rather than authentic locus-specific support.
            stats.withheld_no_unique_anchor += 1;
            stats.withheld_no_attributable_anchor += 1;
            false
        } else {
            // F6b: pairs admissible for both rivals (e.g. rows that bridge a
            // whole same-length bubble into the shared endpoint) carry no
            // differential evidence, so the ratio compares exclusive
            // support only. A rival with no exclusive pairs degenerates to
            // the single-candidate quality bar inside
            // `differential_rivals_agree`.
            divergent_rivals.iter().all(|rival| {
                differential_rivals_agree(
                    matrix,
                    graph,
                    prov,
                    dist,
                    evidence_anchor,
                    sep_offset,
                    best,
                    rival,
                    params.min_support,
                    params.decision_threshold,
                )
            })
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
                if prov.meta(step_dir.id()).contaminating(step_dir.is_rc()) {
                    contaminated = true;
                }
                if step_node.seq_len >= 500 && used_unique.insert(step_dir.id()) {
                    inserted_now.push(step_dir.id());
                }
                anchor_gap += step_node.seq_len as i64 - step_edge.overlap_k as i64;
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
                if anchor_reset_eligible(graph, prov, *step_dir, contaminated) {
                    unique_anchor = Some(*step_dir);
                    anchor_steps = ext.len();
                    anchor_gap = 0;
                }
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

/// PB (provenance boundary): builds the provenance meta of the composite
/// unitig about to be pushed for a committed `full_path`.
///
/// `provenance_boundary` is the spliced 5'-offset at which the committed
/// path enters its first pooled-locus node (multi-context or BCI-flagged in
/// the walked orientation). A pair read at 5'-offset `p` sits at exit
/// distance `dist_u1_exit = seq_len - p`, so it is attributable to the
/// unique locus only while `p < boundary`, i.e. while
/// `dist_u1_exit > seq_len - boundary`; rows at or below that cutoff belong
/// to the shared zone and are inadmissible as unique-locus evidence (MR-2).
/// Constituent context counts merge by max (MR-1), read in each
/// constituent's walked orientation (DSC-3): an rc-entered constituent
/// contributes its `_rev` scalars to the composite's `_fwd` side. A
/// constituent that already carries its own composite boundary preserves
/// it in composite spliced coordinates: in reverse orientation
/// `[0, seq_len - b)` is pooled, so the composite boundary begins
/// immediately at the constituent's entry offset (`step_start_offset`);
/// in forward orientation `[0, b)` is unique, so the composite boundary
/// is `step_start_offset + b`. A pristine constituent is marked at its
/// entry offset when its walked orientation is contaminating. `None` means
/// the whole composite is uniformly classified (no pooled-locus
/// constituent).
fn composite_provenance(
    graph: &CdbgGraph,
    prov: &ProvenanceTable,
    full_path: &[PathStep],
) -> NodeProvenanceMeta {
    let mut meta = NodeProvenanceMeta::default();
    let mut assembled_len: u32 = 0;
    for (idx, step) in full_path.iter().enumerate() {
        let node = &graph.nodes[step.dir.id() as usize];
        let m = prov.meta(step.dir.id());
        let (fwd_cnt, rev_cnt, bci_fwd, bci_rev, bci_downstream_fwd, bci_downstream_rev, exit_iso_fwd, exit_iso_rev) = if step.dir.is_rc() {
            (
                m.external_context_count_rev,
                m.external_context_count_fwd,
                m.bci_rev,
                m.bci_fwd,
                m.bci_downstream_rev,
                m.bci_downstream_fwd,
                m.exit_isolated_rev,
                m.exit_isolated_fwd,
            )
        } else {
            (
                m.external_context_count_fwd,
                m.external_context_count_rev,
                m.bci_fwd,
                m.bci_rev,
                m.bci_downstream_fwd,
                m.bci_downstream_rev,
                m.exit_isolated_fwd,
                m.exit_isolated_rev,
            )
        };
        meta.external_context_count_fwd = meta.external_context_count_fwd.max(fwd_cnt);
        meta.external_context_count_rev = meta.external_context_count_rev.max(rev_cnt);
        meta.bci_fwd |= bci_fwd;
        meta.bci_rev |= bci_rev;
        if idx == 0 {
            meta.bci_downstream_fwd = bci_downstream_fwd;
            meta.exit_isolated_rev = exit_iso_rev;
            meta.entry_shadow_reach_fwd = m.entry_shadow_reach(step.dir.is_rc());
        }
        if idx + 1 == full_path.len() {
            meta.bci_downstream_rev = bci_downstream_rev;
            meta.exit_isolated_fwd = exit_iso_fwd;
            meta.entry_shadow_reach_rev = m.entry_shadow_reach(step.dir.is_rc());
        }

        let step_start_offset = if idx == 0 {
            0
        } else {
            assembled_len.saturating_sub(step.overlap as u32)
        };

        if meta.provenance_boundary.is_none() {
            if let Some(b) = m.provenance_boundary {
                if step.dir.is_rc() {
                    // In reverse orientation, [0, seq_len - b) is pooled, so the
                    // constituent enters pooled territory immediately at step_start_offset.
                    meta.provenance_boundary = Some(step_start_offset);
                } else {
                    // In forward orientation, [0, b) is unique, so the pooled boundary
                    // occurs at step_start_offset + b.
                    meta.provenance_boundary = Some(step_start_offset.saturating_add(b));
                }
            } else if m.contaminating(step.dir.is_rc()) {
                meta.provenance_boundary = Some(step_start_offset);
            }
        }

        if idx == 0 {
            assembled_len = node.seq_len;
        } else {
            assembled_len = step_start_offset + node.seq_len;
        }
    }
    if let Some(first_step) = full_path.first() {
        meta.structural_in_degree_fwd = prov
            .meta(first_step.dir.id())
            .structural_in_degree(first_step.dir.is_rc());
        meta.structural_out_degree_rev = prov
            .meta(first_step.dir.id())
            .structural_out_degree(first_step.dir.twin().is_rc());
    }
    if let Some(last_step) = full_path.last() {
        meta.structural_out_degree_fwd = prov
            .meta(last_step.dir.id())
            .structural_out_degree(last_step.dir.is_rc());
        meta.structural_in_degree_rev = prov
            .meta(last_step.dir.id())
            .structural_in_degree(last_step.dir.twin().is_rc());
    }
    meta
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

    // F23 provenance precompute (DCC + BCI): static graph analysis in
    // O(V + E), read-only over the graph, executed once per resolution.
    let mut prov = ProvenanceTable::analyze(graph, dist.max_insert);

    // 1. Identify long unique seed unitigs (length >= 500 bp, non-deleted, non-circular)
    let mut seeds: Vec<u32> = (0..n)
        .filter(|&id| {
            let node = &graph.nodes[id as usize];
            !node.is_deleted && !node.is_circular && node.seq_len >= 500
        })
        .collect();

    // Sort seeds by length descending
    seeds.sort_by_key(|&id| std::cmp::Reverse(graph.nodes[id as usize].seq_len));

    let mut fail_ctx = 0usize;
    let mut fail_bci_ds = 0usize;
    let mut fail_bci = 0usize;
    for &seed_id in &seeds {
        let m = prov.meta(seed_id);
        if m.external_context_count_fwd != 1 {
            fail_ctx += 1;
        }
        if m.bci_downstream_fwd {
            fail_bci_ds += 1;
        }
        if m.bci_fwd {
            fail_bci += 1;
        }
    }
    tracing::info!(
        "[SeedProfile] total_seeds={}, fail_ctx={}, fail_bci_ds={}, fail_bci={}",
        seeds.len(),
        fail_ctx,
        fail_bci_ds,
        fail_bci
    );

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
            &prov,
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
            &prov,
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

                    // PB (MR-1): the composite unitig records where its
                    // sequence shifts from the unique-locus prefix into the
                    // pooled-locus suffix, keeping the provenance side table
                    // in lockstep with the graph.
                    let composite_meta = composite_provenance(graph, &prov, &full_path);
                    prov.push_composite(composite_meta);

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
    tracing::info!(
        "[exSPAnder] F23 Horizon/PAG stats: withheld no unique anchor = {} (si2_seed = {}, horizon = {}, no_attr_anchor = {}, pag_empty = {})",
        stats.withheld_no_unique_anchor,
        stats.withheld_si2_seed,
        stats.withheld_horizon,
        stats.withheld_no_attributable_anchor,
        stats.withheld_pag_empty
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

    fn push_plain_node(graph: &mut CdbgGraph, len: usize) -> u32 {
        let id = graph.nodes.len() as u32;
        let seq: Vec<u8> = (0..len).map(|i| b"ACGT"[(i + id as usize) % 4]).collect();
        let seq_id = graph.arena.push_seq(&seq);
        graph
            .nodes
            .push(UnitigNode::new(id, seq_id, len as u32, 25.0));
        id
    }

    fn link(graph: &mut CdbgGraph, from: DirUnitigId, to: DirUnitigId, ov: u8) {
        graph.nodes[from.id() as usize].add_edge(from.is_rc(), GraphEdge::new(to, ov, 10));
        graph.nodes[to.id() as usize]
            .add_edge(to.twin().is_rc(), GraphEdge::new(from.twin(), ov, 10));
    }

    /// PB coordinate orientation (MR-2): `dist_u1_exit` is the distance from
    /// the read to u1's 3' exit, so a read at 5'-offset `p` has
    /// `dist_u1_exit = u1_len - p` and falls in the unique prefix
    /// `[0, boundary)` iff `dist_u1_exit > u1_len - boundary`. With
    /// u1_len = 470 and boundary = 300, a prefix read at offset 50
    /// (dist_u1_exit = 420) must stay admissible while a suffix read at
    /// offset 350 (dist_u1_exit = 120) must be excluded.
    #[test]
    fn pb_boundary_filters_pooled_suffix_in_exit_coordinates() {
        let mut graph = CdbgGraph::new(31);
        let u1_id = push_plain_node(&mut graph, 470);
        let u2_id = push_plain_node(&mut graph, 400);
        let u1 = DirUnitigId::new(u1_id, false);
        let u2 = DirUnitigId::new(u2_id, false);

        let mut prov = ProvenanceTable::default();
        prov.push_composite(NodeProvenanceMeta {
            provenance_boundary: Some(300),
            external_context_count_fwd: 1,
            external_context_count_rev: 1,
            ..NodeProvenanceMeta::default()
        });
        prov.push_composite(NodeProvenanceMeta {
            external_context_count_fwd: 1,
            external_context_count_rev: 1,
            ..NodeProvenanceMeta::default()
        });

        let mut matrix = PairLinkMatrix::new();
        // Both spans (620 and 320) are inside the default [250, 700] window,
        // so the only filter that can separate the rows is the boundary.
        matrix.add_evidence(
            u1,
            u2,
            ReadPairEvidence {
                pair_id: 0,
                dist_u1_exit: 420,
                dist_u2_entry: 100,
                r1_len: 100,
                r2_len: 100,
            },
        );
        matrix.add_evidence(
            u1,
            u2,
            ReadPairEvidence {
                pair_id: 1,
                dist_u1_exit: 120,
                dist_u2_entry: 100,
                r1_len: 100,
                r2_len: 100,
            },
        );

        let dist = InsertSizeDist::default();
        let kernel = |span: i64| -> f64 {
            let z = (span as f64 - dist.median).abs() / (1.4826 * dist.mad);
            (-0.5 * z * z).exp()
        };

        let (count, weight) = score_candidate_path_prov(&matrix, &graph, &prov, &dist, u1, u2, 100);
        assert_eq!(count, 1, "only the unique-prefix row may count");
        assert!(
            (weight - kernel(620)).abs() < 1e-12,
            "weight must be the prefix row's kernel weight, got {weight}"
        );

        // Without a boundary both rows are admissible (control orientation).
        let mut unbounded = ProvenanceTable::default();
        unbounded.push_composite(NodeProvenanceMeta {
            external_context_count_fwd: 1,
            external_context_count_rev: 1,
            ..NodeProvenanceMeta::default()
        });
        unbounded.push_composite(NodeProvenanceMeta {
            external_context_count_fwd: 1,
            external_context_count_rev: 1,
            ..NodeProvenanceMeta::default()
        });
        let (count, weight) =
            score_candidate_path_prov(&matrix, &graph, &unbounded, &dist, u1, u2, 100);
        assert_eq!(count, 2);
        assert!(
            (weight - (kernel(620) + kernel(320))).abs() < 1e-12,
            "unbounded control must score both rows, got {weight}"
        );
    }

    /// DCC source SCC origin (AR-2): a source node of a linear chain has
    /// exactly one established lineage and stays anchor-eligible; similarly,
    /// an isolated circular chromosome or plasmid cycle (a source SCC in the
    /// condensation graph) originates its own authentic genomic lineage and
    /// is anchor-clean.
    #[test]
    fn dcc_cyclic_source_scc_originates_established_lineage() {
        let mut chain = CdbgGraph::new(31);
        let s = push_plain_node(&mut chain, 600);
        let c = push_plain_node(&mut chain, 300);
        link(
            &mut chain,
            DirUnitigId::new(s, false),
            DirUnitigId::new(c, false),
            30,
        );
        let prov = ProvenanceTable::analyze(&chain, 700);
        assert_eq!(
            prov.meta(s).external_context_count(false),
            1,
            "linear source carries its own established lineage"
        );
        assert!(
            prov.anchor_clean(DirUnitigId::new(s, false)),
            "single-lineage source must remain anchor-eligible"
        );

        let mut cycle = CdbgGraph::new(31);
        let a = push_plain_node(&mut cycle, 600);
        let b = push_plain_node(&mut cycle, 600);
        link(
            &mut cycle,
            DirUnitigId::new(a, false),
            DirUnitigId::new(b, false),
            30,
        );
        link(
            &mut cycle,
            DirUnitigId::new(b, false),
            DirUnitigId::new(a, false),
            30,
        );
        let prov = ProvenanceTable::analyze(&cycle, 700);
        assert_eq!(
            prov.meta(a).external_context_count(false),
            1,
            "cyclic source SCC originates its own authentic lineage"
        );
        assert!(
            !prov.meta(a).is_multi_context(false),
            "single source SCC is not multi-context"
        );
        assert!(
            prov.anchor_clean(DirUnitigId::new(a, false)),
            "cyclic source SCC must be anchor-clean"
        );
    }
}
