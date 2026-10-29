# ACTIVE TASK: fastBact F24 — Ready for Final Acceptance

**F24 status (2026-10-29):** Phase F24 Provenance-Safe Contiguity Recovery has completed sequential implementation across all 12 architectural invariants (F24-I01 through F24-I12), independent review with explicit `ACCEPT` by Codex GPT-6.1 Sol Medium (Round 11, task-4105), 100% deterministic verification across all 8 frozen test suites (Slurm Job `3065181`, 106/106 passed, 0 failures), and real-data A/B QUAST benchmark on R1-001 and R1-003 vs frozen baseline F23 (Slurm Job `3065210`). Contiguity recovered significantly (+78.0% NGA50 on R1-001, +164.1% NGA50 on R1-003) with zero new misassemblies and zero new local misassemblies. Status is **`F24_READY_FOR_FINAL_ACCEPTANCE`**. Local merge into `main`, remote push, and tag creation remain pending user authorization as `FINAL_AUTHORIZATION_AUTHORITY`. Full documentation: [`F24_FINAL_TECHNICAL_ACCEPTANCE_REPORT.md`](./F24_FINAL_TECHNICAL_ACCEPTANCE_REPORT.md).

## Phase F24 Acceptance Record

- **Overall Status**: **`F24_READY_FOR_FINAL_ACCEPTANCE`**
- **Active Worktree**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f24-contiguity-worktree`
- **Active Branch**: `feature/f24-provenance-safe-contiguity`
- **Candidate Release Binary**: `target/release/fastBact` SHA-256 `8348395b9de19f7a3d21bf1179737d3b10427ee782c6fa6461f9c0b3608f1820` (Job `3065182`)
- **Baseline Release Binary**: `01_fastBact/target/release/fastBact` SHA-256 `97adb3c2ae386d453ba3f50fa61920836fd6c37d83924c0e0bc04d6763caaeb7` (Commit `68b19f06d22537cf199e8fe124ea6b62bb0bb50a`)
- **Deterministic Verification**: 106/106 tests passed (Job `3065181`, exit 0, 1.76s on `vnode1`)
- **Independent Review**: Codex GPT-6.1 Sol Medium Round 11: **`ACCEPT`**
- **Real-Data A/B Benchmark (QUAST v5.2.0, Job `3065210` on `vnode1`)**:
  - `R1-001` (*E. coli* K-12): N50 28.1 kb $\to$ **48.6 kb** (+72.8%), NGA50 27.0 kb $\to$ **48.1 kb** (+78.0%), Misassemblies: strictly **5** (zero new misassemblies, AUD-01 chimera remains eliminated), Local misassemblies: strictly **0**, Valid joins: 544 $\to$ **840** (+296 joins).
  - `R1-003` (*E. coli* BW25113): N50 17.9 kb $\to$ **46.2 kb** (+158.0%), NGA50 17.5 kb $\to$ **46.2 kb** (+164.1%), Misassemblies: strictly **0**, Local misassemblies: strictly **1**, Valid joins: 63 $\to$ **347** (+450% joins).
- **Execution Safety**: Contained strictly in feature worktree. Zero merge to main, zero remote push, zero git tags.

## Frozen F23 acceptance record

> **Project status (2026-10-29; general rules are Multi-AI Collaboration V3.2 in `~/.config/ai-collaboration/POLICY.md`, which is not restated here):** Phase F23 Provenance-Aware Repeat Resolution has achieved complete technical completion and final scientific acceptance. Following comprehensive acceptance evidence closure (`FASTBACT_F23_ACCEPTANCE_EVIDENCE_CLOSURE.md`) and independent review PASS by Codex GPT-6.1 Sol Medium (R1–R6), the user as `FINAL_AUTHORIZATION_AUTHORITY` formally granted final merge and acceptance authorization (`USER_FINAL_MERGE_AUTHORIZATION_GRANTED` on 2026-10-29). The candidate implementation from `feature/c0-02r-phase-f23-unique-anchor` (Commit `0b567c6`) has been integrated into `main`. The scientific contract `C0-02R` is formally **`ACCEPTED`**, and Phase F23 is transitioned to **`F23_COMPLETED`**.

**Overall Acceptance Status**: **ACCEPTED** (C0-02R formally ACCEPTED; Phase F23 COMPLETED; integrated into local `main` under user authorization on 2026-10-29)
**Final Integration SHA**: `0b567c62b1f680f8c24d9fe4fc95c905f28b9666` (feature branch HEAD)
**Production Source**: `src/resolve/exspander.rs` SHA-256 `fd59e8bdc231dbb83468c711f3d87834c51638d75941f587eb0aa1c89298790c`
**Release Binary**: `target/release/fastBact` SHA-256 `2f1b6229c241479e9e3480adc6d60ab009d96f7ca8f7ff78307412a74b1c9a11` (Slurm Job `3052676`)
**Deterministic Quality Gates**: 118/118 passed / 0 failed (Slurm Job `3052557`, exit `0:0`, `ALL_GATES_PASSED`)
**Full Release Suite**: 194/194 passed / 0 failed (Slurm Job `3049336`, exit `0:0`)
**Real-Data A/B Benchmark (Slurm Job `3052717`)**:
- `R1-003`: 0 extensive misassemblies (vs 2 in Baseline A), 0.00 indels/100kb (vs 0.04), 0.59 mismatches/100kb (vs 1.51).
- `R1-001`: Historical AUD-01 2.1 Mb `contig_16` chimera completely eliminated. The 5th misassembly (`contig_48`) attributed to QUAST ambiguity filtering threshold artifact over a pre-existing baseline precursor.
**Accepted Limitations**: NGA50 contiguity trade-off (27.0 kb on R1-001, 17.5 kb on R1-003) due to SI-2 seed withholding when repeat lengths exceed paired-end library insert size ($L_{repeat} > I_{max} \approx 700$ bp).
**Execution Safety**: Local main integrated. Zero remote push, zero GitHub release, zero git tags created.
**F23 Implementation**: Feature worktree: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f23-unique-anchor-worktree`, branch `feature/c0-02r-phase-f23-unique-anchor`. Candidate binary `target/release/fastBact` SHA-256 `2f1b6229c241479e9e3480adc6d60ab009d96f7ca8f7ff78307412a74b1c9a11`.
**F23 Test Status**:
- Targeted Gate Suite (Job `3052557`): 118 passed / 0 failed, exit `0:0`, `ALL_GATES_PASSED`.
- End-to-end Integration Benchmark (Job `3049334`): 1 passed in 1.60s (<3.5s limit), exit `0:0`.
- Full Workspace Release Suite (Job `3049336`): 194 passed / 0 failed across all targets, exit `0:0`.
**Slurm Benchmark Jobs**: Job `3052557` (F23 Gate Matrix, exit `0:0`), Job `3049334` (Integration, exit `0:0`), Job `3049336` (Workspace Release, exit `0:0`), Job `3052717` (Real Data A/B QUAST Benchmark, exit `0:0`) on `vnode1` (`qfree`).
**F22-D Diagnostic Build**: Binary SHA-256 `efc97b451549de1ca732e5411a8747926ee190989e9a3a62b5a0844c9a3e917c`, `cargo build --release --features resolver-trace`, rustc 1.89.0
**Active Worktree (F22-D)**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f22-diagnostics-worktree`
**Active Branch (F22-D)**: `feature/c0-02r-phase-f22-diagnostics` (base `2ae4acc` = F21-B evidence HEAD)
**Key Evidence Directories**:
- Phase F21-B: [`docs/development/evidence/C0_02R/phase_f21b/`](./evidence/C0_02R/phase_f21b/) — report `FORENSIC-20261023-002`, independent audit [`codex_independent_trace_audit.md`](./evidence/C0_02R/phase_f21b/codex_independent_trace_audit.md) (`AUDIT-20261023-017`, **FAIL / QUALIFIED FINDINGS**)
- Phase F22-D: [`docs/development/evidence/C0_02R/phase_f22/`](./evidence/C0_02R/phase_f22/) — trace design, decision traces, pair evidence, equivalence report
- Phase F21-A (prior): [`docs/development/evidence/C0_02R/phase_f21/`](./evidence/C0_02R/phase_f21/) (`AUDIT-20261023-016`)
**Controller / Handoff**: Follow current controller and handoff information in [`HANDOFF.md`](./HANDOFF.md). The main worktree remains on `fix/c0-02r-final-gate`; F23 work is isolated in the worktree and branch listed above. Do not infer a merge or cherry-pick between them.
**Last Updated**: 2026-10-23

---

## 1. Current Phase Summary (F21-B → F22-D)

### 1.1 Phase F21-B — Repeat-Join Decision Provenance (FORENSIC-20261023-002)
1. Static/derived evidence reconstruction of the chimeric join decision across the 7-node repeat chain (`u1926 → u3275 → u1475 → u1476 → u4349 → u4465 → u3919`, raw 2,092 bp, overlap-adjusted span 1,076 bp).
2. Independent Codex audit `AUDIT-20261023-017`: **FAIL / QUALIFIED FINDINGS** — mechanism plausible but not certifiable from derived evidence alone: `u1476` locus specificity contradicted (annotated `TWO_COPY_HOMOLOGOUS_REPEAT` at both 1.57 Mb and 3.66 Mb), candidate-specific decision weights not reproducible from the pair ledger, polishing stage arithmetic and unique-loss accounting inconsistencies. Disposition: `PHASE_F21B_DIAGNOSTIC_EVIDENCE_REQUIRED` → authoritative runtime evidence required.
3. Gate 0 re-verified `BENCHMARK_CONFIGURATION_VERIFIED`: A/B identical CLI, k-mer schedule `[31, 63, 127]`, 500 bp assembler cutoff, QUAST `--min-contig 200`; F20 divergence is algorithmic (P2 insert-size estimator: B `374 ± 53` vs A `476 ± 54`).

### 1.2 Phase F22-D — Runtime Repeat-Join Decision Trace (Diagnostic Non-Interference Verified)
1. Added compile-time-gated `resolver-trace` feature (default off) recording real resolver fork decisions (`decision_trace.tsv`) and read-pair provenance (`pair_evidence.tsv`); `tests/test_c0_02r_resolver_trace.rs` covers trace equivalence, dedup, rank/score consistency, rejection reasons.
2. **Diagnostic non-interference verified bit-for-bit**: instrumented release build reproduces Phase F20 Candidate B `contigs.fa` and `stage_05_exspander.fa` SHA-256 exactly on R1-001 and R1-003 ([`trace_equivalence_report.md`](./evidence/C0_02R/phase_f22/trace_equivalence_report.md), `EQUIV-20261023-F22D-TRACE`).
3. **Root cause recovered from runtime decision `DEC_000022`** (R1-001, seed `u6867`): extending through the unbranched short-repeat chain advanced the anchor step-by-step until the search anchor became the collapsed shared repeat node `u2102` (`unitig_1476`, 255 bp) instead of the upstream long unique anchor (`u241`/`unitig_169`, 12.6 kb). `score_pair_evidence` then counted ectopic read pairs originating from the 3.66 Mb homologous copy of `u2102`: 16 reads supported the wrong exit `u2732` (CAND_01) vs 6 for the authentic exit `u4675` (CAND_02); ratio `2.819 > 1.2` wrongly committed the 2.1 Mb chimeric join.
4. Counterfactual strategy matrix ([`counterfactual_strategy_matrix.md`](./evidence/C0_02R/phase_f22/counterfactual_strategy_matrix.md), `STRATEGY-20261023-F22D-MATRIX`): **Option A (Unique Anchor Retention) RECOMMENDED primary; Option C (bilateral unique flank check) secondary guard; Options B/D/E REJECTED**.

### 1.3 F21-A Forensic Facts (retained)
- R1-001 extra misassembly is an authentic chimeric relocation (`contig_16` = Stage 5 `unitig_1631` bit-identical), 3,666,175 bp $\leftrightarrow$ 1,570,641 bp (-2,095,535 bp).
- 84.6% (132/156) of B's final error events concentrate in the 101–350 bp terminal polishing boundary zone (Stage 5 pre-polish: B fewer mismatches than A, 24 vs 29; indels equal). Stage 6 polishing boundary repair remains OUT of F23 scope (future F24 candidate).
- 1,561 short unitigs (235 kb) fall below the 500 bp output cutoff; Stage 5 genome fraction is higher in B (+0.38 pp R1-001, +0.41 pp R1-003).

---

## 2. Prior Phase Milestone Registry

| Phase | Scope | Baseline / Target | Status | Key Milestone / Commit |
| :--- | :--- | :--- | :---: | :--- |
| **Phase F19** | Codex-Led P2 Verification Repair | `39f11ca` | **AUDITED (4/5 CLOSED, 1 P2 OPEN)** | Repaired oracle and sensitivity tests. Audit `AUDIT-20261022-014`. Commit `0606be4`. 183 tests pass on `vnode1`. |
| **Phase F20** | Codex-Led Selective Real-Data A/B Benchmark (R1-001, R1-003) | `0606be4` vs `b9727c3` | **COMPLETED (REGRESSION DETECTED)** | Slurm job `3001561` on `vnode1`. R1-003 improved; R1-001 regressed (+1 misassembly, indels 0.48 $\to$ 1.03). Audit `AUDIT-20261022-015`. C0-02R NOT ACCEPTED. |
| **Phase F21-A** | Targeted Forensic Investigation of Real-Data Regression | `4b956a4` | **FORENSIC COMPLETED** | Proved R1-001 extra misassembly is real chimeric join (`contig_16`). Proved base errors are 84.6% boundary polishing artifacts. Audit `AUDIT-20261023-016`. C0-02R NOT ACCEPTED. |
| **Phase F21-B** | Repeat-Join Decision Provenance & Mechanism Verification | `07e0e5f` | **AUDITED (FAIL / QUALIFIED)** | Reconstructed static decision evidence across 7-node repeat chain. Audit `AUDIT-20261023-017`: mechanism plausible but derived evidence not certifiable. Commit `2ae4acc`. C0-02R NOT ACCEPTED. |
| **Phase F22-D** | Runtime Repeat-Join Decision Trace (resolver-trace diagnostics) | `2ae4acc` | **DIAGNOSTIC COMPLETE (non-interference verified)** | Runtime decision `DEC_000022` recovered; anchor-loss + ectopic-evidence hijack root cause established. Diagnostics commit `9c543a0` on `feature/c0-02r-phase-f22-diagnostics` (feature-gated, not merged to production). Independent audit pending (scheduled post-F23 benchmark). C0-02R NOT ACCEPTED. |
| **Phase F23** | Provenance-Aware Repeat Resolution | `2ae4acc` vs `68b19f0` | **ACCEPTED (COMPLETED)** | Option A+B unique anchor retention; Sol Round 6 PASS; Job `3052717` eliminated AUD-01 chimera; user authorized merge to main on 2026-10-29. |
| **Phase F24** | Provenance-Safe Contiguity Recovery | `68b19f0` (baseline) | **READY_FOR_FINAL_ACCEPTANCE** | 12 architectural invariants implemented; Sol Round 11 ACCEPT; Slurm `3065181` 106/106 PASS; Slurm `3065210` R1-001 NGA50 +78.0%, R1-003 NGA50 +164.1%, 0 new misassemblies. |
