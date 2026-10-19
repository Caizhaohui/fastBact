# fastBact — Design Proposal: Lightweight Fragment Span Verification for Repeat Resolution (P1)

**Document ID**: `DESIGN-20261019-P1-FRAGMENT-SPAN`  
**Target Subsystem**: `src/resolve/mapper.rs`, `src/resolve/exspander.rs`  
**Author**: Antigravity CLI (Main Implementer & Coordinator)  
**Status**: DRAFT / PROPOSED (Phase E Deliverable; non-invasive, no production refactor in Phase E)  
**Context**: Verification of Codex Finding P1 (`AUDIT-20261019-001`, verified empirically in `tests/test_c0_02r_p1_counterexample.rs`)  

---

## 1. Executive Summary

In commit `f241d9d`, exSPAnder repeat resolution added a necessary scalar constraint:
```rust
let gap = target_node.seq_len as i64 - edge_to_target.overlap_k as i64 - edge_target_to_cand.overlap_k as i64;
if gap <= dist.max_insert as i64 { ... }
```
As revealed by Codex audit `AUDIT-20261019-001` and empirically verified in `test_p1_shortcut_vulnerability_counterexample`:
- `gap` measures only the repeat node separation between flanking unitigs, **not** the physical span of the paired-end fragments.
- `PairLinkMatrix` only aggregates link counts `(u1 -> u2)` without retaining read coordinates.
- When an erroneous shortcut (e.g. $A \to R \to B$, span 391 bp) competes with a true longer path ($A \to R \to C \to B$, span 600 bp), authentic fragments originating from the longer path produce `(A, B)` link counts. Because $R$ has `gap = 38 <= 700`, exSPAnder falsely chooses the shortcut, deleting $C$ from the assembly.

This document proposes a **minimal, zero-copy, lightweight fragment coordinate evidence architecture** that enables exSPAnder to verify candidate-specific path spans without storing raw reads, quality scores, or full read alignments.

---

## 2. Technical Architecture & Answers to Core Requirements

### 2.1 Lightweight Representation of Read-Pair Coordinates

A paired-end fragment $F$ spans from the 5' end of read $R_1$ on unitig $u_1$ to the 5' end of reverse read $R_2$ on unitig $u_2$.
Along any candidate path $P = (u_1, \dots, u_2)$, the total physical fragment span is:
$$\text{Span}(P) = d_1 + \text{PathSeparation}(P) + d_2$$
where:
- $d_1 = \text{len}(u_1) - s_1$: Distance from read 1 start to the 3' exit junction of $u_1$.
- $d_2 = e_2$: Distance from the 5' entry junction of $u_2$ to read 2 end.
- $\text{PathSeparation}(P) = \sum_{v \in P_{\text{internal}}} \text{len}(v) - \sum_{e \in P} \text{overlap}(e)$: Net sequence length separating the 3' junction of $u_1$ and 5' junction of $u_2$.

Because `index_terminal_kmers` indexes only the terminal 1000 bp of each unitig, $d_1 \le 1000$ and $d_2 \le 1000$.
Therefore, each supporting read pair requires only **4 bytes** (two `u16` fields):

```rust
/// Compact evidence of a paired-end read fragment linking two unitigs.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct FragmentTerminalDistance {
    /// Distance in bp from R1 start to u1 3' exit: len(u1) - s1 (<= 1000 bp)
    pub dist_from_u1_exit: u16,
    /// Distance in bp from u2 5' entry to R2 end: e2 (<= 1000 bp)
    pub dist_from_u2_entry: u16,
}
```

In `PairLinkMatrix`:
```rust
#[derive(Clone, Debug, Default)]
pub struct PairLinkMatrix {
    /// Directed unitig pair links with exact terminal distance evidence
    pub links: FxHashMap<(DirUnitigId, DirUnitigId), SmallVec<[FragmentTerminalDistance; 8]>>,
}
```

### 2.2 Memory Overhead Analysis

- **Current Implementation**:
  `links: FxHashMap<(DirUnitigId, DirUnitigId), u32>`  
  $\approx 24\text{ bytes}$ per unique unitig pair entry.
- **Proposed Implementation**:
  `links: FxHashMap<(DirUnitigId, DirUnitigId), SmallVec<[FragmentTerminalDistance; 8]>>`
  - Inline storage: `SmallVec<[FragmentTerminalDistance; 8]>` occupies 32 bytes inline.
  - Up to 8 supporting pairs are stored inline without heap allocation.
  - Pairs with $>8$ links allocate a small heap buffer (4 bytes per additional pair, capped at e.g. 64 pairs).
- **Scale in Bacterial Genomes**:
  - Typical bacterial assembly: $200 \sim 1,500$ unitig pairs with non-zero link support.
  - $1,500 \times 48\text{ bytes} \approx 72\text{ KB}$.
  - Compared to the hundreds of megabytes used by the k-mer index and de Bruijn graph, the memory overhead of storing terminal distances is **completely negligible (< 0.1% of total peak RSS)**.

### 2.3 Computational Overhead Analysis

- **During `map_paired_reads`**:
  - $s_1$ and $e_2$ are already computed in the P2 estimator loop.
  - Computing $d_1 = \text{len}(u_1) - s_1$ and $d_2 = e_2$ requires two integer subtractions.
  - Pushing 4 bytes into `SmallVec`: $O(1)$.
  - **Estimated time overhead**: $< 1.5\%$ increase in `map_paired_reads` runtime.
- **During `exSPAnder` Path Extension**:
  - For each candidate outgoing edge $cand$ from repeat $R$:
    Compute $\text{PathSeparation}(P) = \text{len}(R) - ov_1 - ov_2$.
  - Evaluate support:
    ```rust
    let valid_support = matrix.get_span_support(curr_dir, cand, path_separation, dist);
    ```
  - For a typical support depth of $10 \sim 30$ pairs, evaluating $\text{span} = d_1 + S(P) + d_2 \in [\text{min\_insert}, \text{max\_insert}]$ takes $< 50\text{ ns}$ per candidate.
  - **Estimated time overhead**: Undetectable (< 0.05% of assembly time).

### 2.4 Interface Compatibility

The API of `PairLinkMatrix` remains fully backwards-compatible:
```rust
impl PairLinkMatrix {
    /// Returns total raw observed pairs (preserves legacy API).
    pub fn get_support(&self, u1: DirUnitigId, u2: DirUnitigId) -> u32 {
        self.links.get(&(u1, u2)).map_or(0, |v| v.len() as u32)
    }

    /// Returns the number of pairs whose implied fragment span across
    /// the specific candidate path falls within [dist.min_insert, dist.max_insert].
    pub fn get_span_support(
        &self,
        u1: DirUnitigId,
        u2: DirUnitigId,
        path_separation: i64,
        dist: &InsertSizeDist,
    ) -> u32 {
        let Some(evidence_list) = self.links.get(&(u1, u2)) else { return 0 };
        evidence_list
            .iter()
            .filter(|ev| {
                let span = ev.dist_from_u1_exit as i64 + path_separation + ev.dist_from_u2_entry as i64;
                span >= dist.min_insert as i64 && span <= dist.max_insert as i64
            })
            .count() as u32
    }
}
```

### 2.5 Protection Against Duplicate Counting

- **Fragment Independence**: Each read pair $i$ in `ReadBank` is iterated once. If both reads uniquely anchor to terminal k-mers on $(u_1, u_2)$, exactly one `FragmentTerminalDistance` entry is appended to `(u1, u2)` and its reciprocal `(u2.twin(), u1.twin())`.
- **Hypothesis Isolation**: When comparing candidate extension paths $P_1 = (curr \to R \to B)$ and $P_2 = (curr \to R \to C \to B)$, the support for each hypothesis is computed independently against the relevant pair endpoints:
  - For $P_1$: evaluated against `(curr, B)` links with separation $S(P_1) = \text{len}(R) - 2k$.
  - For $P_2$: evaluated against `(curr, C)` links with separation $S(curr \to R \to C) = \text{len}(R) - 2k$, and subsequently `(R, B)` or `(curr, B)` with separation $S(P_2) = \text{len}(R) + \text{len}(C) - 3k$.
- Because each physical read pair exists as a single event, no pair can be counted multiple times toward the same path hypothesis.

### 2.6 Repeat-Copy Ambiguity Evasion

- **K-mer Uniqueness**: In `index_terminal_kmers`, any k-mer appearing $>1$ time in the graph is discarded (`kmer_pos.retain(|h, _| *kmer_counts.get(h) == 1)`).
- **Flanking Specificity**: Reads are never anchored to repetitive internal sequence. They are only anchored in unique terminal regions flanking the repeat.
- **Span Specificity**: If two identical repeat copies have different lengths or connect different genomic loci, read pairs spanning locus 1 will have physical distances incompatible with locus 2. The explicit test `span in [min_insert, max_insert]` rejects chimeric or cross-copy support.

### 2.7 Fallback Rules when Distance Evidence is Missing

- **Principle: Safe Non-Resolution**.
- If a repeat has candidate outgoing edges, but none of the candidates achieve `valid_span_support >= params.min_support` with confidence ratio $\ge 2.5$:
  - **Action**: Do **NOT** resolve the repeat.
  - The branch point in the cDBG is left intact.
  - In downstream graph cleaning or output contig generation, the unitigs remain separate contigs rather than being falsely spliced into a misassembled chimeric contig.
- **Prohibition**: Never fall back to unconstrained raw link counts (`get_support()`), as that re-introduces the P1 shortcut vulnerability.

### 2.8 Required Deterministic Test Suites for Production Readiness

Before merging this feature into any release branch, the following deterministic tests must be implemented:
1. `test_p1_shortcut_rejection`: The exact counterexample in `tests/test_c0_02r_p1_counterexample.rs`, asserting that $A \to R \to B$ is rejected and the true path $A \to R \to C \to B$ is assembled.
2. `test_span_boundary_precision`: Synthetic fragments constructed at `min_insert - 1` (rejected), `min_insert` (accepted), `max_insert` (accepted), `max_insert + 1` (rejected).
3. `test_variable_repeat_lengths`: Repeats of varying lengths ($50, 150, 300, 500\text{ bp}$) correctly resolved when within library span.
4. `test_twin_direction_reciprocity`: Verifying that backward extension from `seed.twin()` evaluates identical valid span support as forward extension.
5. `test_quast_regression_gate`: Full 6-dataset Slurm benchmark demonstrating zero increase in mismatches or indels.

---

## 3. Implementation Roadmap

- **Phase E** (Current): Proposal drafted and preserved; empirical counterexample committed on feature branch.
- **Phase F** (Future Work upon User Approval):
  1. Create isolated feature branch `feature/c0-02r-p1-span-support`.
  2. Implement `FragmentTerminalDistance` and `get_span_support()` in `mapper.rs` and `exspander.rs`.
  3. Validate against deterministic test suites on Slurm compute node `vnode1`.
  4. Submit full 3-way A/B/C benchmark job on Slurm cluster.
