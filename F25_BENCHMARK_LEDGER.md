# fastBact Phase F25 — Benchmark Ledger

**Phase**: F25 — Assembly Correctness, Residual Misassembly Forensics, and Cross-Dataset Robustness  
**Controller**: Antigravity CLI / Gemini 3.8 Flash Medium  
**Policy**: Multi-AI Collaboration V3.2  

---

## Benchmark Run Registry

| Run ID | Timestamp | Slurm Job ID | Partition / Node | Dataset | Target Software / Binary SHA-256 | QUAST NGA50 (bp) | Genome Fraction (%) | Extensive Misassemblies | Local Misassemblies | Mismatches / 100kb | Indels / 100kb | Determinism $\sigma^2$ | Status |
| :--- | :--- | :--- | :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **BM-F24-FINAL** | 2026-10-29 | 3066737 | `qfree` / `vnode1` | R1-001 (rep1-3) | `1f1e5d87...` | 45,092 | 98.115% | 5 | 0 | 0.85 | 0.53 | 0.000000 | **F24_ACCEPTED_BASELINE** |
| **BM-F24-FINAL** | 2026-10-29 | 3066737 | `qfree` / `vnode1` | R1-003 (rep1-3) | `1f1e5d87...` | 44,702 | 98.060% | 0 | 1 | 0.13 | 0.00 | 0.000000 | **F24_ACCEPTED_BASELINE** |
| **VAL-F24-MAIN** | 2026-10-30 | 3072360 | `qfree` / `vnode1` | Release Test Suite (213 tests) | `53998c59...` | N/A | N/A | N/A | N/A | N/A | N/A | N/A | **POST_INTEGRATION_PASS** |
| **VAL-F25-SUITE**| 2026-10-30 | 3072530 | `qfree` / `vnode1` | Release Test Suite (213 tests) | `7c70e357...` | N/A | N/A | N/A | N/A | N/A | N/A | N/A | **PASS_213_TESTS** |
| **BM-F25-R1-001** | 2026-10-30 | 3072539 | `qfree` / `vnode1` | R1-001 (`SRR28349173`) | `7c70e357...` | 45,092 | 98.115% | 5* | 0 | 0.85 | 0.53 | 0.000000 | **IDENTICAL_NON_REGRESSION** |
| **BM-F25-R1-003** | 2026-10-30 | 3072539 | `qfree` / `vnode1` | R1-003 (`SRR28349164`) | `7c70e357...` | 44,702 | 98.060% | 0 | 1 | 0.13 | 0.00 | 0.000000 | **IDENTICAL_NON_REGRESSION** |
| **BM-F25-R1-001-VAL** | 2026-10-30 | 3072539 | `qfree` / `vnode1` | R1-001-Val (`SRR28349174`) | `7c70e357...` | 54,554 | 98.233% | 0 | 0 | 0.81 | 0.24 | N/A | **GENERALIZATION_PASS** |
| **BM-F25-R1-003-VAL** | 2026-10-30 | 3072539 | `qfree` / `vnode1` | R1-003-Val (`SRR28349165`) | `7c70e357...` | 39,314 | 98.097% | 0 | 1 | 0.26 | 0.00 | N/A | **HIGH_DEPTH_PASS** |

*\*Audited evaluation artifacts (biological e14 phase inversion + IS5/IS1 repeat entry stubs).*
