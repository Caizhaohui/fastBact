# ACTIVE TASK: fastBact C0-02R Phase F21-A (Targeted Forensic Investigation of Real-Data Regression)

**Status**: PHASE_F21_FORENSIC_COMPLETE_PENDING_DECISION  
**Overall Acceptance Status**: **NOT ACCEPTED** (C0-02R maintained strictly NOT ACCEPTED pending design decision on chimeric repeat filter)  
**Independent Forensic Audit**: **COMPLETED (AUDIT-20261023-016, Verdict: AUDITED WITH CORRECTIONS; 100% measurement replication)**  
**Forensic Mechanism Classifications**:
- `NEW_JOIN_ERROR_CONFIRMED`: Real chimeric relocation on R1-001 `contig_16` (-2,095,535 bp across 3.66 Mb $\leftrightarrow$ 1.57 Mb).
- `POLISHING_BOUNDARY_REGRESSION_SUPPORTED`: 84.6% of errors in 101-350 bp terminal zone; Stage 5 pre-polish had fewer errors in B than A.
- `OUTPUT_FILTER_SEQUENCE_LOSS_SUPPORTED`: 1,561 short unitigs (235 kb) filtered out by 500 bp cutoff; Stage 5 genome fraction higher in B.
**Production Code Status**: **100% BIT-IDENTICAL & IMMUTABLE** (Tree SHA `30967bf4940844aeb6b744732bcad6d2a72daab4`, 0 lines changed in `src/`)  
**Slurm Benchmark Jobs**: Job `3001561` (F20 benchmark) & Job `3007620` (F21 Stage 5 QUAST) on `vnode1` (`qfree`)  
**Active Worktree**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21-forensics-worktree`  
**Active Branch**: `feature/c0-02r-phase-f21-forensics`  
**Base Evidence Commit**: `4b956a46935248779e17b8e3453f4bcb188f8e1f`  
**Deliverables Directory**: [`docs/development/evidence/C0_02R/phase_f21/`](./evidence/C0_02R/phase_f21/)  
- Forensic Manifest: [`forensic_manifest.tsv`](./evidence/C0_02R/phase_f21/forensic_manifest.tsv)  
- Misassembly Event Matrix: [`misassembly_event_matrix.tsv`](./evidence/C0_02R/phase_f21/misassembly_event_matrix.tsv)  
- Misassembly Locus Analysis: [`misassembly_locus_analysis.md`](./evidence/C0_02R/phase_f21/misassembly_locus_analysis.md)  
- Stage 5 vs 6 Comparison Table: [`stage5_stage6_comparison.tsv`](./evidence/C0_02R/phase_f21/stage5_stage6_comparison.tsv)  
- Base Error Loci: [`base_error_loci.tsv`](./evidence/C0_02R/phase_f21/base_error_loci.tsv)  
- Boundary Distribution: [`boundary_error_distribution.tsv`](./evidence/C0_02R/phase_f21/boundary_error_distribution.tsv)  
- Unitig Fate Summary: [`unitig_fate_summary.tsv`](./evidence/C0_02R/phase_f21/unitig_fate_summary.tsv)  
- R1-003 Contrast Analysis: [`r1_003_contrast.md`](./evidence/C0_02R/phase_f21/r1_003_contrast.md)  
- Causal Hypothesis Matrix: [`causal_hypothesis_matrix.md`](./evidence/C0_02R/phase_f21/causal_hypothesis_matrix.md)  
- Independent Forensic Audit: [`codex_independent_forensic_audit.md`](./evidence/C0_02R/phase_f21/codex_independent_forensic_audit.md) (`AUDIT-20261023-016`)  
- Full Forensic Report: [`phase_f21_report.md`](./evidence/C0_02R/phase_f21/phase_f21_report.md)  
**Coordinator / Orchestrator**: Antigravity CLI  
**External AI Resource**: Codex CLI (`gpt-6.1-sol`, medium reasoning effort, read-only sandbox)  
**Last Updated**: 2026-10-23  

---

## 1. Forensic Investigation Summary

### 1.1 Core Findings
1. **R1-001 Extensive Misassembly ($4 \to 5$)**:
   - Confirmed as an authentic chimeric relocation on `contig_16` (71,975 bp).
   - Joins reference 3,666,175 bp to 1,570,641 bp (-2,095,535 bp inconsistency).
   - Bridges 7 short intermediate nodes (`u1926`, `u3275`, `u1475`, `u1476`, `u4349`, `u4465`, `u3919`, total 2,092 bp).
   - Exists bit-identically in Stage 5 `unitig_1631` (pre-polishing).
2. **Base-Level Mismatch and Indel Elevation**:
   - In Stage 5 (pre-polishing), Candidate B had **fewer mismatches (24 vs 29)** and **identical indels (4 vs 4)** compared to Baseline A.
   - **84.6% of errors (132/156)** are concentrated in the 101–350 bp terminal boundary zone.
   - In deep internal sequence ($> 350\text{ bp}$), Candidate B has fewer errors than Baseline A (20 vs 22).
3. **Genome Fraction Drop**:
   - 1,561 short unitigs (totalling 235 kb) in Candidate B fell below the 500 bp default contig output cutoff.
   - At Stage 5 (where all contigs are evaluated), Candidate B had higher genome fraction (+0.38 pp on R1-001, +0.41 pp on R1-003).
4. **R1-003 vs R1-001 Divergence**:
   - Both isolates share the 1.56 Mb $\leftrightarrow$ 3.66 Mb repeat family.
   - R1-003 improved because Candidate B eliminated a second false join (+590 kb relocation on `contig_4`), while R1-001 regressed because Baseline A had avoided the 2.1 Mb repeat whereas Candidate B bridged it.

---

## 2. Prior Phase Milestone Registry

| Phase | Scope | Baseline / Target | Status | Key Milestone / Commit |
| :--- | :--- | :--- | :---: | :--- |
| **Phase F19** | Codex-Led P2 Verification Repair | `39f11ca` | **AUDITED (4/5 CLOSED, 1 P2 OPEN)** | Repaired oracle and sensitivity tests. Audit `AUDIT-20261022-014`. Commit `0606be4`. 183 tests pass on `vnode1`. |
| **Phase F20** | Codex-Led Selective Real-Data A/B Benchmark (R1-001, R1-003) | `0606be4` vs `b9727c3` | **COMPLETED (REGRESSION DETECTED)** | Slurm job `3001561` on `vnode1`. R1-003 improved; R1-001 regressed (+1 misassembly, indels 0.48 $\to$ 1.03). Audit `AUDIT-20261022-015`. C0-02R NOT ACCEPTED. |
| **Phase F21-A** | Targeted Forensic Investigation of Real-Data Regression | `4b956a4` | **FORENSIC COMPLETED** | Proved R1-001 extra misassembly is real chimeric join (`contig_16`). Proved base errors are 84.6% boundary polishing artifacts. Audit `AUDIT-20261023-016`. C0-02R NOT ACCEPTED. |
