# fastBact Phase F26 — Benchmark Ledger

**Phase:** F26 — Cross-Species Assembly Benchmarking, Robustness Hardening & Competitive Validation  
**Controller:** Antigravity CLI / Gemini 3.8 Flash Medium  
**Policy:** Multi-AI Collaboration V3.2  

---

## Benchmark Run Registry

| Run ID | Timestamp | Slurm Job ID | Partition / Node | Dataset | Tool & Version | N50 (bp) | NGA50 (bp) | Genome Fraction (%) | Extensive Misasm | Local Misasm | Mismatches / 100kb | Indels / 100kb | Wall Time (s) | Peak RSS (MB) | Status |
|:---|:---|:---|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|
| **BM-F25-ORACLE-01** | 2026-10-30 | 3072539 | `qfree` / `vnode1` | R1-001 (*E. coli* MG1655) | fastBact F25 | 48,115 | 45,092 | 98.115 | 5* | 0 | 0.85 | 0.53 | 170.0 | 1,638 | **F25_FROZEN_BASELINE** |
| **BM-F25-ORACLE-02** | 2026-10-30 | 3072539 | `qfree` / `vnode1` | R1-003 (*E. coli* BW25113) | fastBact F25 | 44,896 | 44,702 | 98.060 | 0 | 1 | 0.13 | 0.00 | 245.1 | 1,817 | **F25_FROZEN_BASELINE** |
| **BM-F25-ORACLE-03** | 2026-10-30 | 3072539 | `qfree` / `vnode1` | R1-001-Val (*E. coli* MG1655) | fastBact F25 | 54,554 | 54,554 | 98.233 | 0 | 0 | 0.81 | 0.24 | 195.7 | 1,791 | **F25_FROZEN_BASELINE** |
| **BM-F25-ORACLE-04** | 2026-10-30 | 3072539 | `qfree` / `vnode1` | R1-003-Val (*E. coli* BW25113) | fastBact F25 | 40,664 | 39,314 | 98.097 | 0 | 1 | 0.26 | 0.00 | 296.4 | 2,046 | **F25_FROZEN_BASELINE** |

*\*Audited evaluation artifacts (biological e14 phase variation + IS5/IS1 repeat entry stubs).*
