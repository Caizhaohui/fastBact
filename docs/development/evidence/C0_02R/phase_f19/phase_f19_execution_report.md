# fastBact — Phase F19 Execution & Evidence Report: Codex-Led P2 Verification Repair

**Phase**: F19  
**Audit ID Reference**: AUDIT-20261022-013 (Baseline), AUDIT-20261022-014 (Audit)  
**Target Branch**: `feature/c0-02r-phase-f19-codex`  
**Target Worktree**: `fastBact-phase-f19-codex-worktree`  
**Baseline Commit**: `39f11ca90a5de13a5b3d5572d29c191c40d9c03d`  
**Test Repair Commit**: `a9930471a73612159a83694d3d661f060042efc8`  
**Orchestrator**: Antigravity  
**Primary Implementer & Auditor**: Codex CLI (`gpt-6.1-sol`, medium reasoning effort, read-only sandbox)  
**C0-02R Status**: **NOT ACCEPTED**  

---

## 1. Executive Summary

Phase F19 addressed the 5 secondary test and oracle findings identified in the Phase F18 Independent Codex Audit (`AUDIT-20261022-013`):
- `F18-P2-001`: Oracle reciprocal selection ambiguity & tie rejection
- `F18-P2-002`: Oracle initialization long-fork reuse
- `F18-P2-003`: Test T30 decision reachability & identical-topology sensitivity control
- `F18-P2-004`: Test T37 $\beta$ boundary sensitivity ($\text{sep} + \beta == \text{max\_insert}$ vs $\text{max\_insert} + 1$)
- `F18-P2-005`: Test T39 queue capacity limit saturation ($Q_{\text{MAX}} = 1024$)

In accordance with Phase F19 governance:
- **Zero production code changes** were made in `src/` (`git diff 39f11ca..a993047 -- src/` is strictly 0 lines).
- All compilation, testing, and regression runs were executed on Slurm `vnode1` (`qfree` partition).
- An independent audit was conducted by Codex CLI (`gpt-6.1-sol`, medium reasoning effort, read-only sandbox), recorded in `AUDIT-20261022-014`.
- **4 of 5 findings are fully RESOLVED and CLOSED**.
- `F18-P2-005` is partially resolved (queue saturation reached at $N = 1025$), but remains open at P2 due to state limit ($S_{\text{MAX}} = 256$) masking queue capacity guard deletion in black-box testing.

---

## 2. Implementation Summary

### F19-A (F18-P2-001): Oracle Reciprocal Selection & Tie Rejection
- Rewrote `independent_resolve_edge` in `tests/test_c0_02r_repeat_completeness_gate.rs`:
  - Removed premature filtering by `ov == chosen_ov` before finding reverse maximum support.
  - Required all reverse edges targeting `from.twin()` to be evaluated.
  - Enforced exact reciprocal support equality: `rev_max_support == chosen_support`.
  - Enforced strict unique winner on reverse strand: `rev_candidates_with_max_support == 1` (rejecting ties).
- Added `test_oracle_reciprocal_edge_controls` with 5 controls (symmetric positive, asymmetric support, tied forward support, competing reverse overlap with higher support, tied reverse support).

### F19-B (F18-P2-002): Oracle Initialization Long-Fork Rejection
- Updated `reference_oracle_search` in `tests/test_c0_02r_repeat_completeness_gate.rs`:
  - Rejected initial hop if `target.id() == fork_node.id() && target_node.seq_len >= 500` or `target.id() == curr_dir.id()`.
  - Strictly enforced long unique node invariant across `curr_dir.id()`, `fork_node.id()`, and candidate path.
- Added `test_oracle_differential_long_fork_self_loop_rejected` covering 600 bp long fork self-loop (rejected) vs 100 bp repeat loop (admitted).

### F19-C (F18-P2-003): Test T30 Decision Reachability & Sensitivity Control
- Updated `test_t30_multiple_short_prefix_nodes_in_cycle` in `tests/test_c0_02r_state_semantics_gate.rs`:
  - Repaired junction sequence mismatch: $X$ starts with `j_b` matching $P_1$ exit `j_b`.
  - Adjusted node $C$ length (129 bp) and distance calibration ($d_1 = 110, d_2 = 110 \implies \text{span} = 238 \in [180, 320]$) so Decision 1 commits cleanly and reaches Decision 2 at $P_2$.
  - At $P_2$, competitor $X$ via $P_2 \to C \to P_1 \to X$ ($\text{span} = 246$, 20 pairs) actively withholds decoy $D$ ($\text{span} = 249$, 5 pairs).
  - Built sensitivity control using exact cloned 6-edge graph (`ctrl_graph = graph.clone()`), omitting only $P_2 \to X$ evidence, asserting $D$ commits.

### F19-D (F18-P2-004): Test T37 $\beta$ Boundary Sensitivity
- Updated `test_t37_beta_distance_boundary` in `tests/test_c0_02r_state_semantics_gate.rs`:
  - Fixture: Seed $S$ (700 bp), short repeat $P$ ($L_P$ bp, self-loop $P \to P$), decoy $D$ (300 bp).
  - Seed-anchored evidence: $D$ ($\text{span} = 209$, support 5), $P$ (support 20, $d_1 + d_2 = 162 \implies \beta(S) = 162$).
  - Evaluates exact boundary $\text{sep}_v + \beta == \text{max\_insert}$:
    - At $L_P = 150$: $\text{sep}_v + \beta = 250 == \text{max\_insert} \implies \text{Probe} \implies \text{Truncated} \implies D$ withheld.
    - At $L_P = 151$: $\text{sep}_v + \beta = 251 > \text{max\_insert} \implies \text{SafePrune} \implies \text{Complete} \implies D$ commits.

### F19-E (F18-P2-005): Test T39 Queue Capacity Limit Saturation
- Updated `test_t39_queue_capacity_limit_truncates` in `tests/test_c0_02r_state_semantics_gate.rs`:
  - Tested $N = 1025$ branches ($> Q_{\text{MAX}} = 1024$), reaching the queue saturation condition.
  - Tested $N = 50$ branches ($\le 1024$), confirming $D$ commits when queue is not saturated.

---

## 3. Slurm Execution Evidence

All test suites were executed on dedicated Slurm compute node `vnode1` (`qfree` partition).

### Job 3001447: Gate Verification (test_c0_02r_repeat_completeness_gate + test_c0_02r_state_semantics_gate)
```text
=== Running test_c0_02r_repeat_completeness_gate ===
running 8 tests
test test_oracle_differential_depth_boundary_truncation ... ok
test test_oracle_differential_normal_fork_resolves_cleanly ... ok
test test_oracle_differential_long_node_cycle_inadmissible ... ok
test test_oracle_differential_zero_evidence_repeat_does_not_falsely_truncate ... ok
test test_oracle_differential_separation_budget_bounds ... ok
test test_oracle_differential_long_fork_self_loop_rejected ... ok
test test_oracle_reciprocal_edge_controls ... ok
test test_oracle_differential_parameterized_tandem_copies ... ok

test result: ok. 8 passed; 0 failed; 0 ignored; finished in 0.00s

=== Running test_c0_02r_state_semantics_gate ===
running 42 tests
test test_t14_tandem_repeat_cycle_must_not_allow_decoy ... ok
...
test test_t30_multiple_short_prefix_nodes_in_cycle ... ok
test test_t37_beta_distance_boundary ... ok
test test_t39_queue_capacity_limit_truncates ... ok
test result: ok. 42 passed; 0 failed; 0 ignored; finished in 0.10s
```

### Job 3001450: Full Release Regression (cargo test --release across all 16 test targets)
```text
running unittests src/lib.rs: 34 passed
running integration test suites:
- test_c0_02r_gate: 1 passed
- test_c0_02r_unitig_graph: 29 passed
- test_c0_02r_pe_matrix: 10 passed
- test_c0_02r_kmer_dispatch: 6 passed
- test_c0_02r_depth_continuation: 8 passed
- test_c0_02r_circular_topology: 5 passed
- test_c0_02r_decision_weights: 5 passed
- test_c0_02r_cycle_termination: 6 passed
- test_c0_02r_pe_indexing: 9 passed
- test_c0_02r_read_boundary: 7 passed
- test_c0_02r_edge_reciprocity: 6 passed
- test_c0_02r_pe_evidence: 6 passed
- test_c0_02r_polish_residual: 1 passed
- test_c0_02r_repeat_completeness_gate: 8 passed
- test_c0_02r_state_semantics_gate: 42 passed
Total: 183 passed, 0 failed, 0 ignored.
```

---

## 4. Independent Codex Audit (AUDIT-20261022-014) Summary

| Finding | Pre-F19 Status | F19 Disposition | Auditor Assessment |
|---|---|---|---|
| `F18-P2-001` | Open P2 | **CLOSED** | Reciprocal selection, reverse overlap inspection, support equality, and all 5 controls verified sound. |
| `F18-P2-002` | Open P2 | **CLOSED** | Long-fork self-loop rejected at initialization; short repeat loop unrolling preserved. |
| `F18-P2-003` | Open P2 | **CLOSED** | T30 reaches Decision 2, isolates competitor evidence, and sensitivity control uses identical topology. |
| `F18-P2-004` | Open P2 | **CLOSED** | T37 sharp boundary verified ($L_P=150 \implies \text{Probe} \implies \text{Withheld}$, $L_P=151 \implies \text{SafePrune} \implies \text{Committed}$). |
| `F18-P2-005` | Open P2 | **PARTIALLY RESOLVED (Open P2)** | Queue saturation ($N=1025$) reached, but state limit ($S_{\text{MAX}}=256$) produces identical truncation if queue guard is removed in black-box test. |

---

## 5. Findings Ledger

- **Confirmed P0**: 0
- **Confirmed P1**: 0
- **Remaining P2**: 1 (`F18-P2-005` / `F19-P2-001`)
- **Stage C0-02R Acceptance**: **NOT ACCEPTED**
