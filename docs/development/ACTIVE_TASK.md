# ACTIVE TASK: fastBact C0-02R Phase F23 — EVIDENCE CLOSURE COMPLETE; READY FOR USER ACCEPTANCE

> **Project status (2026-10-28; general rules are Multi-AI Collaboration V3.2 in `~/.config/ai-collaboration/POLICY.md`, which is not restated here):** Phase F23 Provenance-Aware Repeat Resolution has completed its final acceptance evidence closure (`FASTBACT_F23_ACCEPTANCE_EVIDENCE_CLOSURE.md`). Forensic analysis of R1-001 `contig_48` proved that the chimeric junction is a shared pre-existing baseline precursor; QUAST reported it as a 5th misassembly solely because Candidate F23 extended 504 bp further, reducing primary coverage below the 99.0% threshold (`--ambiguity-score 0.99`) and revealing the junction. Assembly theoretical limits have been qualified under explicit non-spanning assumptions without unproven universal error claims. Codex GPT-6.1 Sol Medium completed targeted review R6, issuing an explicit **`PASS`** on Anchor Absorption safety under ARG/PAG. The complete provenance ledger reconciles all `exspander.rs` SHA-256 hashes from unpatched HEAD to Candidate source (`fd59e8bdc231dbb83468c711f3d87834c51638d75941f587eb0aa1c89298790c`) and release binary (`2f1b6229c241479e9e3480adc6d60ab009d96f7ca8f7ff78307412a74b1c9a11`). Production code in main repository remains untouched. C0-02R strictly remains **`NOT ACCEPTED`** awaiting user final review and authorization as `FINAL_AUTHORIZATION_AUTHORITY`.

**Overall Acceptance Status**: **NOT ACCEPTED** (C0-02R remains strictly NOT ACCEPTED; F23 evidence closure complete; READY_FOR_USER_ACCEPTANCE; awaiting user final authorization for merge, push, or release)
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
