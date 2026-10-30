# ACTIVE TASK: fastBact F25 — Assembly Correctness, Residual Misassembly Forensics, and Cross-Dataset Robustness

**F25 Status (2026-10-30):** Phase F25 has achieved technical completion, independent audit pass by Codex Sol, full deterministic verification (213/213 release tests), and four-dataset benchmark validation.
- **Disposition:** **`F25_READY_FOR_USER_ACCEPTANCE`**
- **Feature Branch:** `feature/f25-assembly-correctness` (Worktree: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f25-correctness-worktree`)
- **Base Commit:** `92be30a38b5779c1da73f27f05c31718ef5dc059`
- **Release Binary Checksum:** `target/release/fastBact` SHA-256 `7c70e3572996ec777b62b8be0080cd8a0a1528ab01f03d0d4b9f3f56640a8112`
- **Autonomous Controller:** Antigravity CLI / Gemini 3.8 Flash Medium
- **Independent Acceptance Auditor:** Codex CLI / `gpt-6.1-sol` (Medium Reasoning, Audit `AUD-F25-03`)
- **Key Scientific Achievements:**
  1. **F25-P1 (Forensic Trace & Classification)**: The five residual extensive misassemblies in R1-001 were conclusively proven to be `EVALUATION_ARTIFACT`:
     - Two events (`contig_53` Left/Right) are active biological phase variations in the e14 prophage (Pin recombinase DNA inversion), confirmed by raw read-pair support for both isoforms in `SRR28349173` (29 and 37 reads for inverted vs 18 and 14 for reference).
     - Three events (`contig_63`, `contig_200`, `contig_222`) are multi-copy repeat entries (IS5 and IS1) that terminate safely at the repeat boundary without false joins; QUAST arbitrarily aligns the terminal stubs to distant genomic copies.
     - Zero of the five events represent chimeric joins or assembly defects.
  2. **Independent Audit Pass**: Codex Sol verified read support directly on compute node `vnode2`, endorsed all five classifications, and verified all corrections in [`F25_ROOT_CAUSE_FORENSICS.md`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f25-correctness-worktree/F25_ROOT_CAUSE_FORENSICS.md).
  3. **F25-P2 (Safety Invariants)**: Conservative evidence isolation ($p > b + H$) and loop re-entrancy protections remain intact.
  4. **F25-P3 (Cross-Dataset Robustness & Benchmark)**: Evaluated on Slurm node `vnode1` (`qfree`, Job `3072539`):
     - R1-001 & R1-003: 100% mathematical determinism ($\Delta = 0, \sigma^2 = 0.000000$) against accepted F24 baseline.
     - R1-001-Val (`SRR28349174`): N50/NGA50 **54,554 bp** (+21.0% vs R1-001), 98.233% genome fraction, **0 extensive misassemblies**, **0 local misassemblies**.
     - R1-003-Val (`SRR28349165`, 270x coverage, 4.27M read pairs): **0 extensive misassemblies**, NGA50 39,314 bp in 296.4s wall time.
- **Verification Gates**: Slurm Job `3072530` passed 213/213 release tests (`cargo test --release`).
- **Core Deliverables Produced**:
  - [`F25_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f25-correctness-worktree/F25_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md)
  - [`F25_ROOT_CAUSE_FORENSICS.md`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f25-correctness-worktree/F25_ROOT_CAUSE_FORENSICS.md)
  - [`F25_IMPLEMENTATION_AND_REVIEW_REPORT.md`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f25-correctness-worktree/F25_IMPLEMENTATION_AND_REVIEW_REPORT.md)
  - [`F25_BENCHMARK_AND_REGRESSION_REPORT.md`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f25-correctness-worktree/F25_BENCHMARK_AND_REGRESSION_REPORT.md)
  - [`F25_FINAL_SCIENTIFIC_ACCEPTANCE_REPORT.md`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f25-correctness-worktree/F25_FINAL_SCIENTIFIC_ACCEPTANCE_REPORT.md)
- **Next Action**: Awaiting user decision as `FINAL_AUTHORIZATION_AUTHORITY` for final merge to `main`.

---

# ACTIVE TASK: fastBact F24 — Final Technical Acceptance (F24_ACCEPTED_WITH_P1-4_WAIVER)

**F24 Status (2026-10-30):** Phase F24 Provenance-Safe Contiguity Recovery has achieved formal technical and scientific acceptance under User Authorization: **`F24_ACCEPTED_WITH_P1-4_WAIVER`**.
- **Disposition**: `F24_ACCEPTED_WITH_P1-4_WAIVER` (Formal user decision on 2026-10-30).
- **P1-4 Floor Retained & Waiver Rationale**: The original static floor target of 46,800 bp for R1-003 NGA50 is formally retained in the historical contract record. The waiver of this floor in favor of the observed 44,702 bp is approved based on:
  1. Significant contiguity recovery: +155.3% NGA50 over the frozen F23 baseline (17,509 bp).
  2. Substantial base-level accuracy gains: mismatches dropped to 0.13 per 100 kbp (-78.0% below baseline 0.59) and indels strictly 0.00.
  3. Zero new extensive or local misassemblies reported relative to F23.
  4. Methodological correctness: conservative evidence isolation takes strict precedence over speculative repeat joining.
  5. Contiguity qualification: 44,702 bp is characterized strictly as the *empirically observed safe contiguity boundary under the implemented conservative evidence isolation criteria*, and is explicitly **not** described as a mathematically proven maximum.
- **P1-1 (R1-003 Base Errors)**: Fully resolved. Polisher grid-sampling artifact in `consensus.rs` eliminated via Pass 2 full sliding-window candidate frequency counting (`seq.windows(k)`).
- **P1-2 & P1-3**: Fully resolved via strict $p > b + H$ evidence filtering and monotonic 3-state lattice (`CleanUnique` $\to$ `Contaminated` $\to$ `Absorbed`).
- **Residual F23-Era Misassemblies Preserved**:
  - `R1-001` (*E. coli* K-12): Strictly 5 extensive misassemblies (identical to F23 baseline; AUD-01 2.1 Mb chimera remains eliminated; 5th misassembly on `contig_48` is a QUAST min-alignment/ambiguity threshold artifact over a pre-existing baseline precursor).
  - `R1-003` (*S. aureus* MW2): Strictly 0 extensive misassemblies, 1 residual local misassembly (identical to baseline).
- **Integrated Commit SHA**: `cfbd00503b3031853bb41441160717d1e5bc3e69` (tag: `f24-accepted-2026-10-30`).
- **Post-Integration Verification**: Slurm Job `3072360` on `vnode1` (`qfree`): 213/213 release tests passed, CLI smoke test OK, binary SHA-256 `53998c592cfe2155fca21a57fa1f39b36f5fc5fb1aeb11db970e760a4742b8fd`.
- **Benchmark Reproducibility**: Slurm Job `3066737`: 3 paired computational replicates (12 assemblies + 12 QUAST evaluations), empirical variance strictly $\sigma^2 = 0.000000$.

## Phase F24 Acceptance Record

- **Overall Status**: **`F24_ACCEPTED_WITH_P1-4_WAIVER`**
- **User Authorization**: Formal Technical Acceptance approved with conditions on 2026-10-30.
- **Integrated Main Commit SHA**: `cfbd00503b3031853bb41441160717d1e5bc3e69`
- **Immutable Freeze Tag**: `f24-accepted-2026-10-30`
- **Post-Integration Release Binary SHA-256**: `53998c592cfe2155fca21a57fa1f39b36f5fc5fb1aeb11db970e760a4742b8fd` (Slurm Job `3072360`)
- **Pre-Merge Candidate Binary SHA-256**: `1f1e5d8770c7657f880a6387cf221a2b530c02993b36ccb89d92f5391723886d` (Slurm Job `3066736`)
- **Baseline Release Binary SHA-256**: `97adb3c2ae386d453ba3f50fa61920836fd6c37d83924c0e0bc04d6763caaeb7` (Commit `68b19f06d22537cf199e8fe124ea6b62bb0bb50a`)
- **Post-Integration Tests**: 213/213 release tests passed (Slurm Job `3072360`), exit `0:0`, CLI smoke test OK.
- **Real-Data A/B 3-Replicate Benchmark (QUAST v5.2.0, Job `3066737` on `vnode1`)**:
  - `R1-001` (*E. coli* K-12): NGA50 27,031 $\to$ **45,092 bp** (+66.8%), Misassemblies: strictly **5** (zero new vs F23), Local misassemblies: strictly **0**, Mismatches: 2.24 $\to$ **0.85** (-62.1%), Indels: 1.09 $\to$ **0.53** (-51.4%), Valid k joins: 542 $\to$ **790**.
  - `R1-003` (*S. aureus* MW2): NGA50 17,509 $\to$ **44,702 bp** (+155.3%), Misassemblies: strictly **0** (zero new), Local misassemblies: strictly **1** (identical to baseline), Mismatches: 0.59 $\to$ **0.13** (-78.0%), Indels: 0.00 $\to$ **0.00**, Valid k joins: 63 $\to$ **341**.
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
