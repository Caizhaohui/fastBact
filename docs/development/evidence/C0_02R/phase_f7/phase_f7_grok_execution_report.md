# fastBact Phase F7 — Cursor Grok 4.7 Extra High Fast Execution Report

**Date**: 2026-10-21  
**Status**: `PHASE_F7_COMPLETED` (C0-02R maintained `NOT ACCEPTED` pending independent re-audit)  
**Implementer**: Cursor Grok 4.7 Extra High Fast (`~/.local/bin/agent`)  
**Orchestrator**: Antigravity CLI  
**Worktree**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f7-grok-worktree`  
**Branch**: `feature/c0-02r-phase-f7-grok`  
**Base Commit (F5B)**: `0e513a3d4c232e501322f2ae3da0df3e7784a8e4`  
**RED Tests Commit**: `4ed4f8c6b90757a41982b6b27928236d6e7fca4b`  
**Production Repair Commit**: `37a7c0baeb251e6211beafc60a2608304916a0f6`  

---

## 1. Agent & Model Telemetry

- **CLI Tool**: Cursor Agent (`~/.local/bin/agent`, v`2026.10.19-6bbf2c4`)
- **Subscription Tier**: Pro (`cai_zhaohui@163.com`, `apiKeySource: login`, 0 API key charges)
- **Model Identified**: `Grok 4.7 256K Extra High Fast` (`--model grok-4.7-xhigh-fast`)
- **Reasoning Effort**: Extra High (`xhigh`)
- **Speed Mode**: Fast
- **Execution Sandbox**: `enabled` (`--sandbox enabled --trust`)
- **Session ID**: `1c044a72-e978-401f-bf89-2b0d4062c878`
- **Cumulative Token Usage**:
  - Step 1 (RED tests): Input 194,086 tokens, Output 97,109 tokens, Cache Read 2,005,504 tokens
  - Step 2 (Production repairs): Input 38,936 tokens, Output 16,566 tokens, Cache Read 4,516,608 tokens

---

## 2. Findings Resolved & Engineering Solutions

### F5-P1-001 (P1 Correctness): Incoming Unmerged Side-Branch Sequence Loss
- **Defect**: When path $P$ was committed and constituent nodes deleted, `protect_unmerged_frontier` only inspected outgoing edges from the last extension node at the cap. Active incoming side branches (e.g. $X \to R$, intermediate joins $X \to R_i$, reverse-orientation side branches) lost their edges into $P$. Being $< 500$ bp and absent from `rejected_nodes`, Stage 5 cleanup deleted $X$ as an orphan, permanently dropping authentic genomic sequence.
- **RED Evidence**: In `tests/test_c0_02r_incoming_side_branch_gate.rs`, Cases A, B, C, D, and E all failed on unpatched code (`4ed4f8c`). For example, Case A panicked with: `F5-P1-001 case A RED: S→R spliced=true, X alive=false, unique 69 bp in graph=false, in raw curation=false... Unpatched cleanup deletes the incoming 100 bp branch after R is consumed.`
- **Implementation**: In `src/resolve/exspander.rs` (`resolve_repeats_exspander_with_stats`), immediately after `validate_and_splice_path` succeeds and before constituent nodes are marked deleted:
  1. Built `consumed_nodes: HashSet<u32>` of all node IDs in `full_path`.
  2. Executed a pre-commit neighbor snapshot across every constituent node in both orientations (`edges_fwd` and `edges_rev`).
  3. Registered every unconsumed, active neighbor (`!consumed_nodes.contains(&neighbor_id) && !graph.nodes[neighbor_id].is_deleted`) into `rejected_nodes`.
  4. On validation failure (`Err`), no external neighbors are registered into `rejected_nodes`, preserving atomic rollback and avoiding debris pollution.
- **GREEN Evidence**: Cases A, B, C, D, E, F, and G in `test_c0_02r_incoming_side_branch_gate` all turned **GREEN** on compute node `vnode1`. All side branches and sequences survive intact.

### F5-P2-002 (P2 Robustness): Search Distance Budget Overlap Arithmetic
- **Defect**: `find_candidate_paths` (line 282) and `depth_limit_hides_continuation` (line 327) rejected paths where `sep_u_exit = sep_u + u_node.seq_len > max_separation` before subtracting `edge.overlap_k`. Transitions where $sep_u + len(u) > max\_separation$ but $sep_u + len(u) - overlap \le max\_separation$ were prematurely rejected.
- **RED Evidence**: `test_f7b_overlap_adjusted_transition_is_explored` and `test_f7b_depth_limit_overlap_adjusted_continuation_stays_hidden` panicked on unpatched code (`4ed4f8c`).
- **Implementation**: In `find_candidate_paths` and `depth_limit_hides_continuation`:
  - Removed premature rejection on `sep_u_exit > max_separation`.
  - For each sequence-valid edge, calculated `sep_next = sep_u + u_node.seq_len as i64 - edge.overlap_k as i64;`
  - Checked `if sep_next <= max_separation` before enqueueing or signaling truncation.
- **GREEN Evidence**: Both tests turned **GREEN** on compute node `vnode1`.

### F5-P2-003 (P2 Test Rigor): Strengthened 66-Step Boundary Gate
- **Defect**: $N=66$ test in `tests/test_c0_02r_frontier_preservation_gate.rs` used a weak assertion `short_survivors >= 1`.
- **Implementation**:
  - Strictly asserted `short_survivors == 2`.
  - Verified exact identity and sequence retention of both $R_{65}$ and $R_{66}$.
  - Asserted exact 100.00% reference coverage union across emitted raw contigs under `min_len = 0`.
- **GREEN Evidence**: Verified **GREEN** on compute node `vnode1`.

---

## 3. Slurm `vnode1` Verification Telemetry

### Full Release Regression Suite
- **Command**: `cargo test --release` executed on `vnode1`
- **Result**: **138 passed, 0 failed, 0 ignored** across all 14 test suites
  - `src/lib.rs`: 32 tests
  - `integration_test.rs`: 1 test
  - `test_c0_01_polish_orientation.rs`: 29 tests
  - `test_c0_02_junction_overlap.rs`: 17 tests
  - `test_c0_02r_cleanup_gate.rs`: 6 tests
  - `test_c0_02r_directed_edge_identity.rs`: 8 tests
  - `test_c0_02r_edge_gate.rs`: 5 tests
  - `test_c0_02r_frontier_preservation_gate.rs`: 5 tests
  - `test_c0_02r_graph_atomicity.rs`: 6 tests
  - `test_c0_02r_incoming_side_branch_gate.rs`: 9 tests (Cases A–G, F7-B budget tests)
  - `test_c0_02r_p1_regression_gate.rs`: 7 tests
  - `test_c0_02r_p2_estimator.rs`: 6 tests
  - `test_c0_02r_pe_evidence.rs`: 6 tests
  - `test_c0_02r_polish_residual.rs`: 1 test

### Small-Scale Synthetic Assembly Benchmarks
- **Binary SHA-256**: `5dd39fa1d7146eaace9d812bcf197d7284bb18146525201e80b77f107ebc7bfd`
- **`S0-001` (5 kb linear)**:
  - Input Ref SHA-256: `562317e0951d25b360c8710a98d06e542dbb2f277eca9b3a4a0a5a01bc859144`
  - Output FASTA SHA-256: `bba42d4d86a010bff26ef11d2022f20004fa8237bf4dd655c7f273f61a1b5dae`
  - Minimap2 Alignment: `NM:i:0`, `cg:Z:4957M`, `cs:Z::4957` (**100.00% identity**, 0 mismatches, 0 indels)
  - Telemetry: Wall time 91.64 ms, Peak RSS 23.9 MB, Slurm Job ID `2992270`
- **`S0-003` (8 kb with repeat)**:
  - Input Ref SHA-256: `b6b97494bc9044bce194f53c821a1dd7b464a4c2ccc2d8f0881763e70f8389a1`
  - Output FASTA SHA-256: `587f72eb5b9ffb9b77de1aa2fdf2b2a6110196ca710c4fe600780e0431b7ba93`
  - Minimap2 Alignment: `NM:i:0`, `cg:Z:7969M`, `cs:Z::7969` (**100.00% identity**, 0 mismatches, 0 indels)
  - Telemetry: Wall time 181.99 ms, Peak RSS 21.9 MB, Slurm Job ID `2992271`
