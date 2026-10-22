# ACTIVE TASK: fastBact C0-02R Phase F20 (Codex-Led Selective Real-Data A/B Benchmark)

**Status**: BENCHMARK_REGRESSION_REVIEW_REQUIRED  
**Overall Acceptance Status**: **NOT ACCEPTED** (C0-02R maintained strictly NOT ACCEPTED; benchmark regression detected on R1-001)  
**F20 Independent Audit**: **COMPLETED (AUDIT-20261022-015, Verdict: FAIL; 7 PASS, 7 REGRESSION)**  
**Benchmark Scope**: Selective real-data A/B benchmark (4 cells: R1-001 A, R1-001 B, R1-003 B, R1-003 A)  
**Slurm Benchmark Job**: `3001561` on `vnode1` (`qfree`, 8 CPUs, 12 GB RAM, runtime `00:13:13`, exit `0:0`)  
**Production Code Status**: **100% BIT-IDENTICAL & IMMUTABLE** (Tree SHA `30967bf4940844aeb6b744732bcad6d2a72daab4`, 0 lines changed in `src/`)  
**Baseline Version A**: `b9727c3ab372f7b995e439d92b72c7a4bfd051cf` (Binary SHA-256: `2b184db0f108cb508ca7e591240eddb551bea8ebad90ac5eb8ee360bac78ea38`)  
**Candidate Version B**: `0606be423d187771f5925a83da2cc377c40c2a3d` (Binary SHA-256: `4fbafeeeffe94abeef42efaa5710f75dbd0af8b234edd2bffba43cc4f50198f0`)  
**Main Branch**: `fix/c0-02r-final-gate`  
**Feature Branch**: `feature/c0-02r-phase-f19-codex` (worktree `fastBact-phase-f19-codex-worktree`)  
**Deliverables Directory**: [`docs/development/evidence/C0_02R/phase_f20/`](./evidence/C0_02R/phase_f20/)  
- Benchmark Manifest: [`benchmark_manifest.tsv`](./evidence/C0_02R/phase_f20/benchmark_manifest.tsv)  
- Dataset Checksums: [`dataset_checksums.tsv`](./evidence/C0_02R/phase_f20/dataset_checksums.tsv)  
- Binary Checksums: [`binary_checksums.tsv`](./evidence/C0_02R/phase_f20/binary_checksums.tsv)  
- Benchmark Script: [`benchmark_commands.sh`](./evidence/C0_02R/phase_f20/benchmark_commands.sh)  
- Slurm Ledger: [`slurm_job_ledger.tsv`](./evidence/C0_02R/phase_f20/slurm_job_ledger.tsv)  
- Slurm Accounting: [`slurm_accounting.tsv`](./evidence/C0_02R/phase_f20/slurm_accounting.tsv)  
- Assembly Metrics (56 fields): [`assembly_metrics.tsv`](./evidence/C0_02R/phase_f20/assembly_metrics.tsv)  
- Acceptance Evaluation: [`acceptance.tsv`](./evidence/C0_02R/phase_f20/acceptance.tsv)  
- Planner Report: [`benchmark_planner_output.md`](./evidence/C0_02R/phase_f20/benchmark_planner_output.md)  
- Assembly Comparison Report: [`assembly_comparison.md`](./evidence/C0_02R/phase_f20/assembly_comparison.md)  
- Scientific Interpretation: [`codex_benchmark_analysis.md`](./evidence/C0_02R/phase_f20/codex_benchmark_analysis.md)  
- Independent Audit Report: [`codex_independent_benchmark_audit.md`](./evidence/C0_02R/phase_f20/codex_independent_benchmark_audit.md) (`AUDIT-20261022-015`)  
**Coordinator / Orchestrator**: Antigravity CLI  
**External AI Resource**: Codex CLI (`gpt-6.1-sol`) for planning, data summary, interpretation, and audit  
**Last Updated**: 2026-10-22  

---

## 1. Quantitative Benchmark Results Summary

### Baseline A vs Candidate B Comparison Matrix

| Dataset | Metric | Baseline A (`b9727c3`) | Candidate B (`0606be4`) | $\Delta$ (B - A) | Change (%) | Gate Status |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **R1-001** (*E. coli* K-12 MG1655) | N50 (bp) | 58,318 | 57,462 | -856 | -1.47% | **REGRESSION** |
| | NGA50 (bp) | 58,318 | 57,127 | -1,191 | -2.04% | **REGRESSION** |
| | Genome Fraction (%) | 98.364 | 98.281 | -0.083 pp | -0.08% | **REGRESSION** |
| | Misassemblies | 4 | 5 | +1 | +25.0% | **REGRESSION** |
| | Local Misassemblies | 0 | 0 | 0 | 0.0% | **PASS** |
| | Mismatches/100kbp | 2.14 | 2.38 | +0.24 | +11.2% | **REGRESSION** |
| | Indels/100kbp | 0.48 | 1.03 | +0.55 | +114.6% | **REGRESSION** |
| | Valid k / k-1 Joins | 987 / 18 | 947 / 10 | -40 / -8 | - | - |
| | Active Unitigs Residual | 315 | 1,728 | +1,413 | 5.49× | - |
| **R1-003** (*E. coli* BW25113) | N50 (bp) | 67,452 | 70,868 | +3,416 | +5.06% | **PASS** |
| | NGA50 (bp) | 63,648 | 70,010 | +6,362 | +10.00% | **PASS** |
| | Genome Fraction (%) | 98.330 | 98.233 | -0.097 pp | -0.10% | **REGRESSION** |
| | Misassemblies | 2 | 1 | -1 | -50.0% | **PASS (IMPROVED)** |
| | Local Misassemblies | 1 | 1 | 0 | 0.0% | **PASS** |
| | Mismatches/100kbp | 1.51 | 1.08 | -0.43 | -28.5% | **PASS (IMPROVED)** |
| | Indels/100kbp | 0.04 | 0.00 | -0.04 | -100.0% | **PASS (IMPROVED)** |
| | Valid k / k-1 Joins | 564 / 7 | 497 / 2 | -67 / -5 | - | - |
| | Active Unitigs Residual | 300 | 704 | +404 | 2.35× | - |

---

## 2. Key Scientific Findings & Interpretation

1. **R1-003 Marked Improvement**:
   - Extensive misassemblies cut in half (2 $\to$ 1, -50%).
   - Mismatch rate decreased by 28.5% (1.51 $\to$ 1.08 per 100 kbp).
   - Indel rate dropped to 0.00 per 100 kbp (100% elimination).
   - Contiguity significantly improved: N50 +5.06%, NGA50 +10.00%.
   - Confirms that the unified search admissibility contract ($\beta(\text{curr})$ lower bound and Complete/Truncated semantics) suppresses false-positive chimeric repeat bridging on this genome.
2. **R1-001 Quality Regression**:
   - One extra extensive misassembly introduced (4 $\to$ 5, +25%).
   - Indels increased from 0.48 to 1.03 per 100 kbp (+114.6%).
   - Mismatches rose slightly from 2.14 to 2.38 per 100 kbp (+11.2%).
   - Active unitigs retained post-exspander increased 5.49× (315 $\to$ 1,728).
   - Hypothesis: Conservative withholding on ambiguous repeat structures prevented false joins, leaving fragmented unitigs that fell below the 500 bp output cutoff. The shift in contig terminal boundaries altered read recruitment during consensus polishing, shifting base-level calls against the reference genome.
3. **Independent Audit Verdict (`AUDIT-20261022-015`)**:
   - Verdict: **FAIL (REGRESSION DETECTED)**.
   - Status: `BENCHMARK_REGRESSION_REVIEW_REQUIRED`.
   - C0-02R strictly maintained **NOT ACCEPTED**.
   - Prohibits expanding to the 18-dataset suite until the extra misassembly on R1-001 is forensically characterized.

---

## 3. Prior Phase Milestone Registry

| Phase | Scope | Baseline / Target | Status | Key Milestone / Commit |
| :--- | :--- | :--- | :---: | :--- |
| **Phase F18** | Independent Codex Audit of F17-I Unified Search Admissibility | `c37ef6e` | **CONDITIONAL PASS** | Audit `AUDIT-20261022-013`: 0 P0, 0 P1, 5 P2 test findings. Commit `39f11ca`. |
| **Phase F19** | Codex-Led P2 Verification Repair | `39f11ca` | **AUDITED (4/5 CLOSED, 1 P2 OPEN)** | Repaired oracle and sensitivity tests (`a993047`). Audit `AUDIT-20261022-014`. Commit `0606be4`. 183 tests pass on `vnode1`. |
| **Phase F20** | Codex-Led Selective Real-Data A/B Benchmark (R1-001, R1-003) | `0606be4` vs `b9727c3` | **COMPLETED (REGRESSION DETECTED)** | Slurm job `3001561` on `vnode1`. R1-003 improved; R1-001 regressed (+1 misassembly, indels 0.48 $\to$ 1.03). Audit `AUDIT-20261022-015`. C0-02R NOT ACCEPTED. |
