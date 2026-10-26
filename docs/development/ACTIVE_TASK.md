# ACTIVE TASK: fastBact C0-02R Phase F23 (Unique-Anchor Horizon Guard) — PENDING INDEPENDENT AUDIT

> **Current F23 worktree status:** F23 Unique-Anchor Horizon Guard is implemented at source HEAD `7bdc6f962f93d54cda192d47fb410e046fd5a627` and awaits an independent algorithm audit. C0-02R remains **NOT ACCEPTED**. The historical F21-B investigation below is retained for context; it is not the active phase. Prepare the independent review in a fresh OpenCode session with `zhipuai-coding-plan/glm-5.3`, read-only permissions, and `EFFORT_UNVERIFIED` until a reasoning variant is verified. Do not start the review as part of this governance synchronization.

**Current Controller / GLM Route**: Codex GPT-6 Luna Low is controller; OpenCode with `zhipuai-coding-plan` is the sole formal GLM route. See the root [`AGENTS.md`](../../AGENTS.md) and [`HANDOFF.md`](HANDOFF.md). Antigravity takes controller duties only if Codex quota is exhausted; Cursor is never controller.

### Historical F21-B Audit Record (retained)

**Historical Status**: `PHASE_F21B_DIAGNOSTIC_EVIDENCE_REQUIRED`
**Overall Acceptance Status**: **Strictly NOT ACCEPTED** (C0-02R maintained strictly NOT ACCEPTED)
**Independent Audit**: **COMPLETED (`AUDIT-20261023-017`, Verdict: FAIL / QUALIFIED FINDINGS)**
**Auditor**: Codex CLI (`gpt-6.1-sol`, medium reasoning, read-only)
**Key Diagnostic Conclusions**:
- `NEW_JOIN_ERROR_CONFIRMED`: Real chimeric relocation on R1-001 `contig_16` (-2,095,535 bp across 3.66 Mb $\leftrightarrow$ 1.57 Mb).
- `REPEAT_CHAIN_VERIFIED`: 7-node repeat chain (`u1926` to `u3919`), net physical span 1,076 bp; 0 spanning physical fragments between flanks `u1927` and `u169`.
- `R1_003_BASELINE_JOIN_CONFIRMED`: In R1-003, Baseline A already committed this exact chimeric relocation on `contig_16` (-2,094,638 bp); Candidate B reproduced it while resolving an independent 590 kb misassembly.
- `POLISHING_BOUNDARY_REGRESSION_CONFIRMED`: 84.6% of base errors (132/156) occur in 101–350 bp terminal boundary zone; pre-polish Stage 5 Candidate B had fewer errors than Baseline A (24 vs 29). Net Stage 5 $\to$ 6 increases: +85 mismatches and +43 indels for B vs +69 mismatches and +18 indels for A.
- `OUTPUT_FILTER_SEQUENCE_LOSS_RECONCILED`: 1,561 short unitigs < 500 bp; 200–499 bp unitigs contribute 25,574 bp unique reference loss, accounting for the QUAST reference coverage drop of 25,715 bp (0.554%).
- `DECISION_TRACE_QUALIFIED`: Structural geometry and 0-spanning pairs confirmed; offline decision-trace models require dedicated runtime diagnostic tracing instrumentation for full active-walker state verification.
**Production Code Status**: **100% BIT-IDENTICAL & IMMUTABLE** (Tree SHA `30967bf4940844aeb6b744732bcad6d2a72daab4`, 0 lines changed in `src/`)
**Slurm Benchmark Jobs**: Job `3001561` (F20 benchmark) & Job `3007620` (F21 Stage 5 QUAST) on `vnode1` (`qfree`)
**Active Worktree**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21b-trace-worktree`
**Historical Active Branch at F21-B**: `feature/c0-02r-phase-f21b-trace`
**Base Evidence Commit**: `07e0e5fd927cc3c8be96b22927e33b8ae6b9e138` (Phase F21-A forensic HEAD)
**Deliverables Directory**: [`docs/development/evidence/C0_02R/phase_f21b/`](./evidence/C0_02R/phase_f21b/)
- Effective Parameter Audit: [`effective_parameter_audit.tsv`](./evidence/C0_02R/phase_f21b/effective_parameter_audit.tsv)
- Repeat Locus Graph: [`repeat_locus_graph.tsv`](./evidence/C0_02R/phase_f21b/repeat_locus_graph.tsv)
- Directed Path Geometry: [`directed_path_geometry.tsv`](./evidence/C0_02R/phase_f21b/directed_path_geometry.tsv)
- Read-Pair Provenance: [`read_pair_provenance.tsv`](./evidence/C0_02R/phase_f21b/read_pair_provenance.tsv)
- Candidate Decision Trace: [`candidate_decision_trace.tsv`](./evidence/C0_02R/phase_f21b/candidate_decision_trace.tsv)
- Alternative Path Matrix: [`alternative_path_matrix.tsv`](./evidence/C0_02R/phase_f21b/alternative_path_matrix.tsv)
- R1-001 vs R1-003 Repeat Comparison: [`r1_001_r1_003_repeat_comparison.md`](./evidence/C0_02R/phase_f21b/r1_001_r1_003_repeat_comparison.md)
- Polish Boundary Evidence: [`polish_boundary_evidence.tsv`](./evidence/C0_02R/phase_f21b/polish_boundary_evidence.tsv)
- Unique Coverage Loss: [`unique_coverage_loss.tsv`](./evidence/C0_02R/phase_f21b/unique_coverage_loss.tsv)
- Mechanism Evidence Matrix: [`mechanism_evidence_matrix.md`](./evidence/C0_02R/phase_f21b/mechanism_evidence_matrix.md)
- Independent Trace Audit: [`codex_independent_trace_audit.md`](./evidence/C0_02R/phase_f21b/codex_independent_trace_audit.md) (`AUDIT-20261023-017`)
- Full Synthesis Report: [`phase_f21b_report.md`](./evidence/C0_02R/phase_f21b/phase_f21b_report.md)
**Coordinator / Orchestrator**: Antigravity CLI
**External AI Resource**: Codex CLI (`gpt-6.1-sol`, medium reasoning effort, read-only sandbox)
**Last Updated**: 2026-10-23

---

## Historical F21-B Investigation (retained)

### 1. Forensic & Mechanism Investigation Summary

### 1.1 Core Findings
1. **R1-001 Extensive Misassembly ($4 \to 5$)**:
   - Confirmed as an authentic chimeric relocation on `contig_16` (71,975 bp).
   - Joins reference 3,666,175 bp to 1,570,641 bp (-2,095,535 bp inconsistency).
   - Bridges 7 short intermediate nodes (`u1926`, `u3275`, `u1475`, `u1476`, `u4349`, `u4465`, `u3919`, total 2,092 bp).
   - Exists bit-identically in Stage 5 `unitig_1631` (pre-polishing).
2. **Base-Level Mismatch and Indel Elevation**:
   - In Stage 5 (pre-polishing), Candidate B had **fewer mismatches (24 vs 29)** and **identical indels (4 vs 4)** compared to Baseline A.
   - **84.6% of errors (132/156)** are concentrated in the 101–350 bp terminal boundary zone.
   - Net increases from Stage 5 to Stage 6: +85 mismatches and +43 indels for B vs +69 mismatches and +18 indels for A.
   - In deep internal sequence ($> 350\text{ bp}$), Candidate B has fewer errors than Baseline A (20 vs 22).
3. **Genome Fraction Drop**:
   - 1,561 short unitigs (totalling 235 kb) in Candidate B fell below the 500 bp default contig output cutoff.
   - Reconciled non-redundant reference loss in 200–499 bp unitigs is 25,574 bp (0.551%), accounting for the QUAST reference coverage drop of 25,715 bp (0.554%).
4. **R1-003 vs R1-001 Divergence**:
   - Both isolates share the 1.56 Mb $\leftrightarrow$ 3.66 Mb repeat family.
   - In R1-003, Baseline A already committed this exact chimeric relocation on `contig_16` (-2,094,638 bp); Candidate B reproduced it on `contig_19` while eliminating an independent 590 kb misassembly on `contig_4`.
   - In R1-001, Baseline A stopped short due to search depth limits, leaving the repeat unbridged.

---

## 2. Prior Phase Milestone Registry

| Phase | Scope | Baseline / Target | Status | Key Milestone / Commit |
| :--- | :--- | :--- | :---: | :--- |
| **Phase F19** | Codex-Led P2 Verification Repair | `39f11ca` | **AUDITED (4/5 CLOSED, 1 P2 OPEN)** | Repaired oracle and sensitivity tests. Audit `AUDIT-20261022-014`. Commit `0606be4`. 183 tests pass on `vnode1`. |
| **Phase F20** | Codex-Led Selective Real-Data A/B Benchmark (R1-001, R1-003) | `0606be4` vs `b9727c3` | **COMPLETED (REGRESSION DETECTED)** | Slurm job `3001561` on `vnode1`. R1-003 improved; R1-001 regressed (+1 misassembly, indels 0.48 $\to$ 1.03). Audit `AUDIT-20261022-015`. C0-02R NOT ACCEPTED. |
| **Phase F21-A** | Targeted Forensic Investigation of Real-Data Regression | `4b956a4` | **FORENSIC COMPLETED** | Proved R1-001 extra misassembly is real chimeric join (`contig_16`). Proved base errors are 84.6% boundary polishing artifacts. Audit `AUDIT-20261023-016`. C0-02R NOT ACCEPTED. Commit `07e0e5f`. |
| **Phase F21-B** | Repeat-Join Decision Provenance & Mechanism Verification | `07e0e5f` | **AUDITED (FAIL / QUALIFIED)** | Provenance ledger & graph geometry verified; Baseline A pre-existence confirmed in R1-003; coverage loss reconciled to 25.5 kb. Audit `AUDIT-20261023-017`. Status: `PHASE_F21B_DIAGNOSTIC_EVIDENCE_REQUIRED`. C0-02R NOT ACCEPTED. |
