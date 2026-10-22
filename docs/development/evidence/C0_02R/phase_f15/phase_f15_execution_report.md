# fastBact Phase F15 — Short Prefix Reuse Repair & Search Completeness Formalization Execution Report

**Date**: 2026-10-22  
**Status**: `PHASE_F15_IMPLEMENTED_READY_FOR_F16_AUDIT`  
**Overall Acceptance Status**: **NOT ACCEPTED** (C0-02R strictly maintained NOT ACCEPTED pending independent Codex re-audit)  
**Primary Implementer**: Cursor Grok 4.7 Extra High Fast (`~/.local/bin/agent`)  
**Orchestrator**: Antigravity CLI  
**Worktree**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f15-grok-worktree`  
**Branch**: `feature/c0-02r-phase-f15-grok`  
**Baseline Commit**: `e1ab70031073706d4191dfd72426da18c553c836` (Phase F13 evidence HEAD)  

---

## 1. Executive Summary

In Phase F14, independent auditor Codex CLI identified Finding **F14-P1-001** and secondary findings **F14-P2-001** ~ **F14-P2-004** (`AUDIT-20261022-009`):
1. **F14-P1-001 (P1)**: In `find_candidate_paths` at line 262, `if u.id() == curr_dir.id() || ext_nodes.contains(&u.id()) { continue; }` unconditionally skipped any short repeat node (< 500 bp) that was part of an accumulated extension prefix `ext`. Consequently, authentic multi-extension routes requiring short prefix reuse (e.g. $S \to P \to C \to P \to X$) were silently dropped. The search returned `Complete`, improperly committing weak shortcut decoy $D$.
2. **F14-P2-001 (P2)**: Small-graph reference oracle in `test_c0_02r_repeat_completeness_gate.rs` called production `resolve_edge_between` and `score_candidate_path`, and 4 out of 6 tests never invoked it.
3. **F14-P2-002 (P2)**: Test T22 constructed only forward graph edges, failing to force a true reverse-complement mirrored traversal.
4. **F14-P2-003 (P2)**: Report test accounting discrepancy (169 reported vs 162 actual inventory).
5. **F14-P2-004 (P2)**: `depth_limit_hides_continuation` omitted circular node check, inconsistent with dequeue policy.

**Phase F15 has fully resolved all findings** through a four-commit cycle:
1. **Commit 1 (`4257612`) — Deterministic RED Regression (F15-A)**:
   Added `test_t23_short_prefix_reuse_must_not_promote_decoy` implementing the exact 2-decision counterexample:
   $S(700) \to P(100) \to C(100) \to D(300)$ with $C \to P \to X(300)$.
   Verified RED on Slurm compute node `vnode1` (Job ID `2997971`), panicking with expected decoy misassembly.
2. **Commit 2 (`7ed0e67`) — Unified Production Repair (F15-B)**:
   Implemented by Cursor Grok 4.7 Extra High Fast in `src/resolve/exspander.rs`.
   - In `find_candidate_paths`: only unique nodes $\ge 500$ bp in `ext_nodes` are pruned; short repeat nodes $< 500$ bp are admitted for bounded candidate exploration.
   - In `pruned_repeat_hides_supported_endpoint`: entry, dequeue, and neighbor filters updated symmetrically.
   - In `depth_limit_hides_continuation`: circular node check added (`F14-P2-004`).
   Verified GREEN on `vnode1` (Job ID `2997983`): all 24 state semantics tests passed!
3. **Commit 3 (`3794034`) — Independent Oracle Formalization & Invariant Gates (F15-D2, F15-E, F15-F)**:
   - Rewrote `tests/test_c0_02r_repeat_completeness_gate.rs` with independent `independent_resolve_edge` and `independent_score_candidate`, without calling production resolver primitives. All 6 differential tests actively invoke the oracle and pass on `vnode1` (Job ID `2997986`).
   - Repaired T22 to construct a true reverse-complement mirrored graph (`F14-P2-002`).
   - Added tests T24–T32 in `tests/test_c0_02r_state_semantics_gate.rs` covering RC prefix reuse, exact budget boundary, exceeding budget clean resolution, normal bridge positive control, opposite twin orientation, consumed shadow nodes, multi-short prefix cycles, zero evidence dead ends, and depth continuation truncation. All 33 state semantics tests passed on `vnode1` (Job ID `2997995`).
4. **Test Accounting Reconciliation (F14-P2-003)**:
   Forensically verified that the baseline test inventory was exactly 162 tests (the 169 in Phase F13 was a clerical addition error with the 7 tests in `test_c0_02r_p1_regression_gate.rs`).
   With the 10 new tests (T23–T32), the repository test suite now contains **exactly 172 tests**.
   `cargo test --release` executed on `vnode1` (Job ID `2998000`): **172 passed, 0 failed, 0 ignored** across 16 test suites!
5. **Small Assembly Benchmarks**:
   Release binary SHA-256: `911d37e1dda03d066448b46057c0c599a6241aeb70bbbc654e0a88708e038a1c`.
   `S0-001` (5 kb linear, Job ID `2998017`): 100.00% identity (`NM:i:0`, `cg:Z:4957M`, `cs:Z::4957`, 95.4 ms, 21.4 MB RSS). Output contig SHA-256 `bba42d4d...` bit-identical to baseline.
   `S0-003` (8 kb repeat, Job ID `2998017`): 100.00% identity (`NM:i:0`, `cg:Z:7969M`, `cs:Z::7969`, 199.6 ms, 23.4 MB RSS). Output contig SHA-256 `587f72eb...` bit-identical to baseline.

---

## 2. Agent & Model Telemetry

- **Primary Implementer**: Cursor Grok 4.7 Extra High Fast (`~/.local/bin/agent`, v`2026.10.19-6bbf2c4`)
- **Subscription Tier**: Pro (`cai_zhaohui@163.com`, 0 API charges)
- **Model Identified**: `Grok 4.7 256K Extra High Fast` (`--model grok-4.7-xhigh-fast`)
- **Reasoning Effort**: Extra High (`xhigh`)
- **Speed Mode**: Fast
- **Execution Sandbox**: `enabled` (`--sandbox enabled --trust`)
- **Orchestrator**: Antigravity CLI

---

## 3. Git Commit DAG & Audit Trails

```
3794034 test(resolve): formalize independent oracle and add invariant tests T24-T32 (F15-D2, F15-E, F15-F)
7ed0e67 fix(resolve): allow admissible short prefix reuse and add circular check in continuation (F15-B)
4257612 test(resolve): add deterministic RED test T23 for short prefix reuse (F15-A)
e1ab700 (feature/c0-02r-phase-f13-grok) docs(evidence): record Phase F13 execution report and grok telemetry
678c930 (feature/c0-02r-phase-f11-grok) fix(resolve): implement checked single-use resolver contract and lifecycle guard (F11-B)
```

- **Commit 1**: `42576123a49219602c77d4cbe394fb028e578c76`
  - Added T23 to `tests/test_c0_02r_state_semantics_gate.rs` (+88 lines).
- **Commit 2**: `7ed0e67614e3e3b33346f04746f3325c862bc346`
  - Production repair in `src/resolve/exspander.rs` (+19 / -5 lines).
- **Commit 3**: `37940348db49b1ff34ec6247ec637cfcb371ae75`
  - Independent oracle in `tests/test_c0_02r_repeat_completeness_gate.rs` (+32 / -22 lines).
  - Repaired T22 and added T24–T32 in `tests/test_c0_02r_state_semantics_gate.rs` (+479 / -49 lines).

---

## 4. Slurm Compute Node (`vnode1`) Verification Evidence

| Step / Suite | Scope | Slurm Job ID | Exit Code | Result | Details |
| :--- | :--- | :---: | :---: | :---: | :--- |
| **RED Test: T23** | Short prefix reuse counterexample | `2997971` | 101 | **FAILED (RED)** | Panicked at `D must not be consumed by a misassembly commit` |
| **GREEN Test: T1–T23** | State semantics gate | `2997983` | 0 | **24 PASSED (GREEN)** | T23 turned GREEN; T1–T22 remained GREEN |
| **Oracle Differential Gate** | Independent reference oracle | `2997986` | 0 | **6 PASSED (GREEN)** | All 6 tests call independent oracle and pass |
| **Invariant Suite: T1–T32** | All 33 state semantics gate tests | `2997995` | 0 | **33 PASSED (GREEN)** | T22 repaired; T24–T32 all passed cleanly |
| **Full Release Regression** | `cargo test --release` (no filter) | `2998000` | 0 | **172 PASSED, 0 FAILED, 0 IGNORED** | Full repository test suite passes cleanly across 16 suites |
| **Release Binary Build** | Production compiler profile | `2998004` | 0 | **SUCCESS** | SHA-256: `911d37e1dda03d066448b46057c0c599a6241aeb70bbbc654e0a88708e038a1c` |
| **Assembly `S0-001` (5 kb)** | Synthetic linear genome | `2998017` | 0 | **100.00% identity** | Minimap2 `NM:i:0`, `cg:Z:4957M`, `cs:Z::4957`, 95.4 ms wall time, 21.4 MB RSS |
| **Assembly `S0-003` (8 kb)** | Synthetic repeat genome | `2998017` | 0 | **100.00% identity** | Minimap2 `NM:i:0`, `cg:Z:7969M`, `cs:Z::7969`, 199.6 ms wall time, 23.4 MB RSS |

---

## 5. Binary & Artifact Cryptographic Integrity

- **F15 Release Binary SHA-256**:
  `911d37e1dda03d066448b46057c0c599a6241aeb70bbbc654e0a88708e038a1c`
- **`S0-001` Output Contig FASTA SHA-256**:
  `bba42d4d86a010bff26ef11d2022f20004fa8237bf4dd655c7f273f61a1b5dae` (bit-identical to frozen baseline)
- **`S0-003` Output Contig FASTA SHA-256**:
  `587f72eb5b9ffb9b77de1aa2fdf2b2a6110196ca710c4fe600780e0431b7ba93` (bit-identical to frozen baseline)

---

## 6. Architecture & Correctness Invariants Preserved

1. **Search/Validator Consistency (Invariant A)**:
   `find_candidate_paths` and `validate_and_splice_path` share the exact same node reuse semantics:
   - Nodes $\ge 500$ bp are strictly unique (appearing at most once in `full_path`).
   - Nodes $< 500$ bp are admissible short repeats (allowed to be traversed in prefix reuse, tandem repeat unrolling, and alternating cycles).
2. **Complete Soundness & Decoy Withholding (Invariants B & C)**:
   `SearchStatus::Complete` is only returned when all competing routes within depth and budget have been examined. When an unresolved or deeper continuation exists, `SearchStatus::Truncated` is propagated to withhold weak decoy promotion.
3. **Bounded Termination (Invariant D)**:
   Bounded by `MAX_CANDIDATE_DEPTH = 4`, `max_separation`, and `MAX_DFS_STATES = 256`.
4. **Checked Single-Use Lifecycle Contract (Invariant H / F11-B)**:
   - `ResolverLifecycleState` (`NeverRun`, `Running`, `Completed`, `Failed`) untouched.
   - `ResolverLifecycleGuard` RAII safety untouched.
   - `Result<JunctionStats, ExspanderError>` signature untouched.
