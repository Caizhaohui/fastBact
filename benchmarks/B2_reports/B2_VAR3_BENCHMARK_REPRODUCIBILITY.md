# B2-C1b-3B Var 3: Benchmark Provenance & Reproducibility Audit Record

**Date:** 2026-10-16  
**Execution Node:** `vnode1` (Linux Slurm HPC cluster, partition `qfree`, 8 dedicated cpus)  
**Base Repository Commit:** `6425e7e141871db5f283923ed9990b6487f5b026`  
**Software Stack:**
- **Python:** `3.12.3`
- **Rustc:** `rustc 1.89.0 (29483883e 2025-08-04)`
- **Cargo:** `cargo 1.89.0 (c24e10642 2025-06-23)`
- **QUAST:** `QUAST v5.2.0`
- **fastBact Binary:** `./target/release/fastBact`

---

## 1. Execution Architecture & Zero-Tolerance Assertions

All benchmark runs in this audit were executed through `./scripts/hpc_run.sh` exclusively on Slurm compute node `vnode1`. Zero computation was performed on the login node.

### Quality Control Gates & Hard Assertions:
1. **NameError and Stagger Regression Fixed:**
   `parse_fastbact_stdout()` was verified with dedicated unit test suite ([`tests/test_benchmark_parser.py`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/tests/test_benchmark_parser.py)). Zero parser crashes or unhandled keys.
2. **Strict QUAST Zero-Exit-Code & Required Metrics Enforcement:**
   Every QUAST run was executed with `--min-contig 200` and evaluated strictly by `parse_quast_strict()`. Any non-zero exit code or missing metric immediately aborts execution.
3. **No Hardcoded Values in Reports:**
   All tables, timing metrics, N50 values, and seed success ratios in `B2_VAR3_FINAL_GATE.md` are dynamically generated directly from the underlying structured run dictionaries and validated TSVs.

---

## 2. Command Lines Executed for Verification

### A. Source Code Gates Verification:
```bash
./scripts/hpc_run.sh "cargo run --release --bin verify_var3_gates"
```
Verified:
- Rejection of illegal threshold ratios (`-1.0`, `0.0`, `0.5`, `0.99`, `NaN`, `Inf`, `-Inf`).
- Exact behavioral consistency between explicit `--genome-size` and HLL estimator paths.
- Bit-identical determinism of `--subsample-seed 0` vs implicit default.
- Bit-identical repeatability of identical non-zero seeds (seed 42).

### B. Python Parser Unit Tests:
```bash
./scripts/hpc_run.sh "python3 tests/test_benchmark_parser.py"
```
Verified 4 unit tests covering active subsampling, bypassed subsampling, safety guard intercepts, and empty log edge cases.

### C. Full Policy & Seed Sensitivity Benchmark:
```bash
./scripts/hpc_run.sh "python3 benchmarks/scripts/run_var3_ab_benchmark.py"
```

---

## 3. Provenance Artifact Manifest

| Generated Artifact | File Path | Records Count | SHA-256 Checksum |
|:---|:---|:---:|:---:|
| **Verified Policy TSV** | [`benchmarks/B2_reports/B2_VAR3_POLICY_AB_VERIFIED.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/benchmarks/B2_reports/B2_VAR3_POLICY_AB_VERIFIED.tsv) | 34 rows | `ef5842f04eb9e75ce7ab036aa159918dda60a88ece1283f8ff9e4cc8251cdedc` |
| **Verified Seed TSV** | [`benchmarks/B2_reports/B2_VAR3_SEED_SENSITIVITY_VERIFIED.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/benchmarks/B2_reports/B2_VAR3_SEED_SENSITIVITY_VERIFIED.tsv) | 31 rows | `4b27a965c0b76faaee5e5d8e6c4d1238ceaf00df61d845a79dfa668047f758a4` |
| **Real Isolate Quality TSV** | [`benchmarks/B2_reports/B2_VAR3_REAL_ISOLATE_QUALITY.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/benchmarks/B2_reports/B2_VAR3_REAL_ISOLATE_QUALITY.tsv) | 10 rows | `b62bdedf3af4bfc57f91239f7d581822e01e1125943eca57faab52d841bb8b88` |
| **Final Gate Verdict MD** | [`benchmarks/B2_reports/B2_VAR3_FINAL_GATE.md`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/benchmarks/B2_reports/B2_VAR3_FINAL_GATE.md) | Comprehensive | Dynamic |
