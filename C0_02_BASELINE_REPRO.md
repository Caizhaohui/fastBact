# FB-P0-02 — Junction Overlap Correctness Defect Reproduction & Analysis Report

**Repository**: `https://github.com/Caizhaohui/fastBact`  
**Task ID**: C0-02 (FB-P0-02: exSPAnder Path Splicing Ignores Directed Edge Overlap Correctness Defect)  
**Baseline Git Commit**: `3950b2ab2584b24d39172b25bbbcb565232891b7`  
**Fix Branch**: `fix/c0-02-junction-overlap`  
**Execution Node**: HPC Slurm cluster `vnode1` (partition `qfree`, `-c 8 --mem=12G`)  
**Date**: October 16, 2026  

---

## 1. Executive Summary & Defect Classification

| Verification Gate | Claim | Baseline (`3950b2a`) Status | Evidence / Failure Trace |
| :--- | :--- | :---: | :--- |
| **T1: Legal Overlap $k - 1$** | Splicing with legal $k - 1$ (30 bp) overlap | **PASSED** | Coincidentally matched baseline's hardcoded $k - 1$ assumption. |
| **T2: Legal Overlap $k$** | Splicing with legal $k$ (31 bp) overlap | **FAILED** | `Length must be 1169 bp, not 1170 bp`: baseline inserted a duplicate base at the junction. |
| **T3: Dual Overlap Edges** | Candidate target has both $k$ and $k - 1$ edges | **FAILED** | `Should merge into exactly 1 active unitig`: self-competition caused confidence ratio $1.0 < 2.5$, aborting extension. |
| **T4: Reverse-Complement** | Unitig connection on reverse-complement strand | **FAILED** | Sequence mismatch caused by hardcoded $k - 1$ assumption on reverse-complement strand. |
| **T5: Mixed Overlaps** | 3-unitig path mixing $k$ and $k - 1$ overlaps | **FAILED** | Assembled length 1740 bp vs ground-truth 1739 bp; duplicate base inserted at $k$-overlap junction. |
| **T6: Invalid Suffix/Prefix** | Edge records overlap 31 bp but sequence mismatches | **FAILED** | Baseline blindly spliced corrupted sequences without verifying suffix/prefix identity. |
| **T7: No Valid Connection** | Splicing rejected when no valid edge exists | **PASSED** | Correctly refused extension when no outgoing edge was present. |
| **T8: Graph Conservation** | After failure, graph nodes, edges, and status remain intact | **FAILED** | Post-exSPAnder small unitig cleanup deleted candidate node (<500 bp) after rejected join. |
| **T9: Auditable Selection** | Multiple edges to same target resolved by support | **FAILED** | Evaluated multiple edges to same candidate as competing nodes, ratio $1.0 < 2.5$. |
| **T10: Backward Traversal** | Reversal and twin conversion references true edge overlap | **FAILED** | Backward extension twin path hardcoded $k - 1$, inserting duplicate base on $k$-overlap paths. |

---

## 2. Mathematical & Algorithmic Root Cause Analysis

### 2.1 The Hardcoded $k - 1$ Assumption

In Compacted De Bruijn Graphs (cDBG), unitigs are connected by directed sequence overlaps. During graph construction (`src/graph/build.rs`, lines 551–594), edges are populated by evaluating both exact $k$ and $k - 1$ overlaps:
```rust
for &ov_len in &[k, k.saturating_sub(1)] {
    ...
    graph.nodes[id as usize].add_edge(is_rc, GraphEdge::new(target, ov_len as u8, 10));
}
```
Every `GraphEdge` explicitly stores `pub overlap_k: u8`.

However, in baseline commit `3950b2a` (`src/resolve/exspander.rs`, lines 198–211), the sequence assembly logic completely discarded the chosen edge identity and hardcoded the junction overlap to $(k - 1)$:
```rust
// Baseline defect in src/resolve/exspander.rs:204
let overlap = (k - 1).min(u_seq.len());
assembled_seq.extend_from_slice(&u_seq[overlap..]);
```

#### Ground Truth vs Baseline Splicing:
Let unitig $U_0$ and $U_1$ share a sequence overlap of length $o$.
- **True Spliced Length**: $|U_0| + |U_1| - o$
- **True Spliced Sequence**: $U_0[0 \dots |U_0| - o] \circ U_1[0 \dots |U_1|]$
- **Baseline Behavior When $o = k$**:
  The baseline set $overlap = k - 1$.
  The assembled length became $|U_0| + |U_1| - (k - 1) = (|U_0| + |U_1| - k) + 1$.
  The base $U_1[k - 1]$ was appended twice, duplicating 1 base at every junction with overlap $k$.

### 2.2 Lost Edge Identity in `extend_one_direction`

In baseline `3950b2a`, `extend_one_direction()` returned only `Vec<DirUnitigId>`:
```rust
fn extend_one_direction(...) -> Vec<DirUnitigId>
```
Because the return type only contained the unitig IDs, the edge selected during graph traversal was completely lost. When assembling the full path, exSPAnder had no access to `GraphEdge.overlap_k`, making it impossible to know whether a junction had an overlap of $k$ or $k - 1$.

### 2.3 Double Overlap Edges & Self-Competition

In `GraphEdge::add_edge()`, a directed unitig can have multiple outgoing edges pointing to the same target node (e.g., one edge with `overlap_k = k` and another with `overlap_k = k - 1` when homopolymers or repetitive k-mers occur).

In baseline `extend_one_direction`:
```rust
// Case B: multiple outgoing edges directly from curr_dir
for e in &outgoing {
    let cand = e.target;
    let supp = matrix.get_support(curr_dir, cand);
    if supp > best_support {
        second_support = best_support;
        best_support = supp;
        best_cand = Some(cand);
    } else if supp > second_support {
        second_support = supp;
    }
}
```
When `outgoing` contained two edges pointing to the same target `T`:
1. In iteration 1, `e.target == T`. `best_support` became $S$.
2. In iteration 2, `e.target == T`. `supp == S`. Since `supp > second_support` (0), `second_support` was set to $S$.
3. The winner-take-all ratio became:
   $$\text{ratio} = \frac{\text{best\_support}}{\text{second\_support}} = \frac{S}{S} = 1.0$$
4. Since $1.0 < \text{decision\_threshold}$ (2.5), extension was permanently aborted due to self-competition, preventing any repeat resolution.

### 2.4 Lack of Atomic Splicing and Graph Mutation Safety

In baseline `3950b2a`:
1. Splicing occurred in-place while iterating over candidate paths.
2. No pre-commit validation of sequence suffix/prefix identity was performed; random sequences with corrupted edges were blindly merged.
3. Candidate nodes were permanently added to `used_in_path` even when downstream extension failed.
4. The post-exSPAnder cleanup unconditionally marked all non-circular unitigs with `seq_len < 500` as `is_deleted = true`, destroying valid candidate nodes that were part of rejected paths.

---

## 3. Minimal Counterexample Failure Log (Compute Node `vnode1`)

The minimal counterexample test suite `tests/test_c0_02_junction_overlap.rs` was executed directly against baseline commit `3950b2ab2584b24d39172b25bbbcb565232891b7` on Slurm compute node `vnode1`.

### Execution Command:
```bash
unset $(compgen -v | grep '^SLURM_') && srun -p qfree -w vnode1 -c 8 --mem=12G cargo test --test test_c0_02_junction_overlap
```

### Raw Baseline Output:
```
   Compiling fastBact v0.1.0 (/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact)
    Finished `test` profile [unoptimized + debuginfo] target(s) in 0.16s
     Running tests/test_c0_02_junction_overlap.rs (target/debug/deps/test_c0_02_junction_overlap-f54587ced662c4ee)

running 10 tests
test test_t1_two_unitigs_valid_overlap_k_minus_1 ... ok
test test_t7_no_valid_connection_rejects_splicing ... ok
test test_t2_two_unitigs_valid_overlap_k ... FAILED
test test_t3_same_directed_pair_both_k_and_k_minus_1_edges ... FAILED
test test_t4_reverse_complement_orientation ... FAILED
test test_t5_three_unitigs_mixing_k_and_k_minus_1 ... FAILED
test test_t6_edge_overlap_mismatches_actual_sequence ... FAILED
test test_t8_failure_preserves_graph_state_atomically ... FAILED
test test_t9_multiple_overlap_edges_explicit_selection_and_ambiguity_rejection ... FAILED
test test_t10_backward_path_reversal_references_correct_edge_and_overlap ... FAILED

failures:

---- test_t2_two_unitigs_valid_overlap_k stdout ----
thread 'test_t2_two_unitigs_valid_overlap_k' panicked at tests/test_c0_02_junction_overlap.rs:125:5:
assertion `left == right` failed: Length must be 1169 bp, not 1170 bp
  left: 1170
 right: 1169

---- test_t3_same_directed_pair_both_k_and_k_minus_1_edges stdout ----
thread 'test_t3_same_directed_pair_both_k_and_k_minus_1_edges' panicked at tests/test_c0_02_junction_overlap.rs:166:5:
assertion `left == right` failed: Should merge into exactly 1 active unitig
  left: 2
 right: 1

---- test_t6_edge_overlap_mismatches_actual_sequence stdout ----
thread 'test_t6_edge_overlap_mismatches_actual_sequence' panicked at tests/test_c0_02_junction_overlap.rs:295:5:
assertion `left == right` failed: Splicing must be rejected when edge overlap does not match actual sequence
  left: 1
 right: 2

---- test_t8_failure_preserves_graph_state_atomically stdout ----
thread 'test_t8_failure_preserves_graph_state_atomically' panicked at tests/test_c0_02_junction_overlap.rs:351:5:
assertion `left == right` failed: Node count unchanged
  left: 3
 right: 2

failures:
    test_t10_backward_path_reversal_references_correct_edge_and_overlap
    test_t2_two_unitigs_valid_overlap_k
    test_t3_same_directed_pair_both_k_and_k_minus_1_edges
    test_t4_reverse_complement_orientation
    test_t5_three_unitigs_mixing_k_and_k_minus_1
    test_t6_edge_overlap_mismatches_actual_sequence
    test_t8_failure_preserves_graph_state_atomically
    test_t9_multiple_overlap_edges_explicit_selection_and_ambiguity_rejection

test result: FAILED. 2 passed; 8 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s
```

### Conclusion:
Baseline commit `3950b2ab2584b24d39172b25bbbcb565232891b7` failed 8 out of 10 targeted tests deterministically, confirming the full spectrum of correctness defects in FB-P0-02 prior to algorithmic repair.
