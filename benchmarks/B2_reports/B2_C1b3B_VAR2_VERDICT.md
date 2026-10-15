# FB-BENCH-B2-C1b-3B: Var 2 High-k CMS Genome Size Propagation Verdict

- **Date:** 2026-10-15
- **Branch:** `fix/b2-c1b3b-var2`
- **Base Commit:** `8cb6fc7` (on `fix/b2-c1b3b-var1`)
- **Production Baseline:** `dac6a57417b390cf7385ff07adf569e550cc8ca9` (`master` / `v0.1.0-baseline: 4259fe3`)
- **Remote Repository:** `https://github.com/Caizhaohui/fastBact`
- **Execution Platform:** Slurm HPC Cluster Node `vnode1` (Partition `qfree`)

---

## 1. Executive Summary

Milestone **`FB-BENCH-B2-C1b-3B Var 2` (High-k CMS Genome Size Propagation)** has been implemented, validated, and benchmarked across all 6 target datasets under strict single-variable isolation:

1. **Unified Genome Size Propagation (`VERIFIED`):** In `src/main.rs`, the calibrated `effective_g` (derived identically from user `--genome-size` or HLL estimator `estimate.genome_size`) is now explicitly passed from the initial $k_1$ stage into each subsequent high-$k$ upgrade stage via `upgrade_multik_ext::<K>(&mut graph, &bank, next_k, Some(effective_g))`.
2. **Elimination of Silent Heuristic Assumptions (`VERIFIED`):** `build_cdbg_with_contigs_ext` in `src/graph/build.rs` no longer silently assumes that the sum of intermediate unitig lengths equals the true bacterial genome size. When an explicit genome size is passed, it is used directly; if omitted, it falls back to a documented, warned `(bank.total_bases() / 60)` estimate.
3. **Preservation of Core Graph Dynamics (`VERIFIED`):** The algorithmic logic of local multi-$k$ graph progression—including branching node identification, early stopping (`NO_BRANCHING_REGIONS`), active boundary read filtering (`filter_active_reads`), prior linear unitig inheritance, and circular plasmid preservation—remains 100% unaltered.
4. **Zero Algorithmic Regression & Bitwise Equivalence across All 6 Datasets (`VERIFIED`):** Across all 6 benchmarked datasets (`S1-COV-150X`, `S1-COV-003`, `SYN-2.5M`, `S2-REAL-001`, `R1-001`, `R1-003`), every single dataset produced **100% bit-for-bit identical assembly FASTA checksums (`SHA-256`)** between Var 1 and Var 2.
5. **Memory and Capacity Sizing Stabilization (`VERIFIED`):**
   - On `S1-COV-003` (300X), Var 1's legacy contig sum heuristic ($490,779\text{ bp}$) accidentally over-allocated a 16.0 MB sketch at $k=63$, while Var 2 maintained a consistent 8.0 MB sketch matching $k=31$ and $k=127$, reducing Peak RSS from 276.2 MB to 259.5 MB.
   - On full genomes (`S2-REAL-001`, `R1-001`, `R1-003`), sketches were consistently and reliably sized at 256.0 MB ($128\text{M}$ counters) throughout all $k$ stages.

---

## 2. A/B Matrix Assembly Results

Summary of empirical results from [`benchmarks/B2_reports/B2_C1b3B_VAR2_ASSEMBLY_AB.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/benchmarks/B2_reports/B2_C1b3B_VAR2_ASSEMBLY_AB.tsv):

| Dataset ID | Variant Mode | N50 (bp) | NGA50 (bp) | Contigs ($\ge 500\text{ bp}$) | Genome Fraction | Misassemblies (Ext / Local) | Peak RSS | Wall Time | Assembly SHA-256 | Bit-for-Bit vs Var 1 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **S1-COV-150X** | Var 1 (Legacy) | 199,983 | 199,983 | 1 | 99.992% | 0 / 0 | 113.9 MB | 6.12 s | `f2081f26d00...` | Baseline |
| **S1-COV-150X** | Var 2 (Propagated) | 199,983 | 199,983 | 1 | 99.992% | 0 / 0 | 111.3 MB | 5.60 s | `f2081f26d00...` | **100% BIT-IDENTICAL** |
| **S1-COV-003** | Var 1 (Legacy) | 199,999 | 199,999 | 1 | 99.999% | 0 / 0 | 276.2 MB | 9.15 s | `e3272828d7a...` | Baseline |
| **S1-COV-003** | Var 2 (Propagated) | 199,999 | 199,999 | 1 | 99.999% | 0 / 0 | 259.5 MB | 9.16 s | `e3272828d7a...` | **100% BIT-IDENTICAL** |
| **SYN-2.5M** | Var 1 (Legacy) | 2,499,977 | 2,499,977 | 1 | 99.999% | 0 / 0 | 582.4 MB | 41.35 s | `764969be5e3...` | Baseline |
| **SYN-2.5M** | Var 2 (Propagated) | 2,499,977 | 2,499,977 | 1 | 99.999% | 0 / 0 | 623.5 MB | 40.48 s | `764969be5e3...` | **100% BIT-IDENTICAL** |
| **S2-REAL-001** | Var 1 (Legacy) | 285,687 | 285,687 | 28 | 99.182% | 1 / 0 | 1,191.9 MB | 82.60 s | `1447ac45e77...` | Baseline |
| **S2-REAL-001** | Var 2 (Propagated) | 285,687 | 285,687 | 28 | 99.182% | 1 / 0 | 1,150.8 MB | 91.00 s | `1447ac45e77...` | **100% BIT-IDENTICAL** |
| **R1-001** | Var 1 (Legacy) | 58,335 | 58,335 | 165 | 98.361% | 4 / 0 | 1,593.3 MB | 161.42 s | `c97ddb5cd09...` | Baseline |
| **R1-001** | Var 2 (Propagated) | 58,335 | 58,335 | 165 | 98.361% | 4 / 0 | 1,427.4 MB | 164.70 s | `c97ddb5cd09...` | **100% BIT-IDENTICAL** |
| **R1-003** | Var 1 (Legacy) | 67,459 | 63,655 | 137 | 98.329% | 2 / 1 | 1,849.8 MB | 238.30 s | `efd20c3b99f...` | Baseline |
| **R1-003** | Var 2 (Propagated) | 67,459 | 63,655 | 137 | 98.329% | 2 / 1 | 1,793.1 MB | 233.83 s | `efd20c3b99f...` | **100% BIT-IDENTICAL** |

---

## 3. Count-Min Sketch Configuration & Stage-by-Stage Telemetry

Empirical stage telemetry extracted from [`benchmarks/B2_reports/B2_C1b3B_VAR2_CMS_CONFIG.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/benchmarks/B2_reports/B2_C1b3B_VAR2_CMS_CONFIG.tsv):

| Dataset | Stage | $k$ | Sizing Source (Var 2) | Size Used | CMS Width | CMS Table Size | Unitigs / Edges | Upgrade Status |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **S1-COV-150X** | `k1_cdbg` | 31 | `ESTIMATOR_HLL` | 206,425 bp | 4,194,304 | 8.0 MB | 1,905 / 3,688 | Initial Build |
| | `multik_63` | 63 | `EXPLICIT_PASSED` | 206,425 bp | 4,194,304 | 8.0 MB | 17 / 10 | Executed |
| | `multik_127` | 127 | `EXPLICIT_PASSED` | 206,425 bp | 4,194,304 | 8.0 MB | 1 / 0 | Executed (Single Contig) |
| **S1-COV-003** | `k1_cdbg` | 31 | `ESTIMATOR_HLL` | 207,337 bp | 4,194,304 | 8.0 MB | 12,396 / 24,868 | Initial Build |
| | `multik_63` | 63 | `EXPLICIT_PASSED` | 207,337 bp | 4,194,304 | **8.0 MB** *(vs 16.0 MB legacy)* | 211 / 186 | Executed |
| | `multik_127` | 127 | `EXPLICIT_PASSED` | 207,337 bp | 4,194,304 | 8.0 MB | 2 / 0 | Executed |
| **SYN-2.5M** | `k1_cdbg` | 31 | `ESTIMATOR_HLL` | 2,696,316 bp | 67,108,864 | 128.0 MB | 6,725 / 12,790 | Initial Build |
| | `multik_63` | 63 | `EXPLICIT_PASSED` | 2,696,316 bp | 67,108,864 | 128.0 MB | 36 / 20 | Executed |
| | `multik_127` | 127 | `SKIPPED` | 0 bp | 0 | 0.0 MB | 36 / 20 | **Early Stop (`NO_BRANCHING`)** |
| **S2-REAL-001** | `k1_cdbg` | 31 | `ESTIMATOR_HLL` | 4,848,540 bp | 134,217,728 | 256.0 MB | 27,018 / 66,220 | Initial Build |
| | `multik_63` | 63 | `EXPLICIT_PASSED` | 4,848,540 bp | 134,217,728 | 256.0 MB | 870 / 906 | Executed |
| | `multik_127` | 127 | `EXPLICIT_PASSED` | 4,848,540 bp | 134,217,728 | 256.0 MB | 62 / 93 | Executed |
| **R1-001** | `k1_cdbg` | 31 | `ESTIMATOR_HLL` | 4,812,237 bp | 134,217,728 | 256.0 MB | 125,933 / 323,694 | Initial Build |
| | `multik_63` | 63 | `EXPLICIT_PASSED` | 4,812,237 bp | 134,217,728 | 256.0 MB | 16,064 / 26,868 | Executed |
| | `multik_127` | 127 | `EXPLICIT_PASSED` | 4,812,237 bp | 134,217,728 | 256.0 MB | 5,597 / 6,084 | Executed |
| **R1-003** | `k1_cdbg` | 31 | `ESTIMATOR_HLL` | 4,638,223 bp | 134,217,728 | 256.0 MB | 187,768 / 384,718 | Initial Build |
| | `multik_63` | 63 | `EXPLICIT_PASSED` | 4,638,223 bp | 134,217,728 | 256.0 MB | 12,690 / 19,838 | Executed |
| | `multik_127` | 127 | `EXPLICIT_PASSED` | 4,638,223 bp | 134,217,728 | 256.0 MB | 2,133 / 2,844 | Executed |

---

## 4. Key Findings & Architectural Observations

1. **Elimination of Uncontrolled Sketch Fluctuation:**
   - Under Var 1, in high-coverage samples like `S1-COV-003` (300X un-subsampled), the intermediate graph prior to curation retained un-curated repetitive error unitigs totaling $490,779\text{ bp}$. This caused the legacy heuristic to size the $k=63$ CMS at 16.0 MB instead of the true 8.0 MB required for a 200 kb bacterial sequence.
   - Under Var 2, passing the vetted `effective_g` ($207,337\text{ bp}$) eliminates this artifact, keeping the sketch stably pinned at 8.0 MB.
2. **Early Stopping Preservation:**
   - In `SYN-2.5M`, after $k=63$, all branching junctions were fully resolved. `upgrade_multik_ext` correctly detected `has_branching == false` and cleanly returned `early_stop_reason: "NO_BRANCHING_REGIONS"`, avoiding unnecessary sketch allocation and execution for $k=127$.
3. **Absence of Any Unintended Behavioral Mutation:**
   - Because Count-Min Sketch saturating counters $\ge 4$ retain exact topological classification across all tested genomes, fixing the capacity formula preserves 100% bit-for-bit equivalence with the certified baseline while eliminating the risk of under-sizing on sparse/fragmented intermediate graphs.

---

## 5. Verification Checklist & Gate Compliance

| Requirement ID | Verification Gate Item | Status | Verified Evidence |
| :--- | :--- | :--- | :--- |
| **G1** | High-$k$ CMS receives explicit, unified genome size configuration | **PASSED** | Telemetry logs and TSV confirm `EXPLICIT_PASSED` matching initial $k_1$ build. |
| **G2** | No modifications to multi-$k$ graph algorithms or simplification rules | **PASSED** | Source diff confirms identical simplification passes and early-stop logic. |
| **G3** | No newly introduced structural errors or misassemblies | **PASSED** | QUAST confirmed identical extensive/local misassemblies across all 6 datasets. |
| **G4** | No regression in Genome Fraction or NGA50 | **PASSED** | QUAST confirms identical GF (e.g. 99.182% on S2, 98.361% on R1) and NGA50. |
| **G5** | CMS memory consumption conforms strictly to power-of-two formula | **PASSED** | Table sizes verified: 8.0 MB (200 kb), 128.0 MB (2.5 Mb), 256.0 MB (4.6-5.4 Mb). |
| **G6** | All builds, executions, and benchmarks executed on Slurm compute node | **PASSED** | Executed synchronously on Slurm compute node `vnode1` (Jobs 2961042, 2961043, 2961046). |

---

## 6. Stop Rule Adherence

- [x] Independent fix branch `fix/b2-c1b3b-var2` created directly from verified Var 1 commit (`8cb6fc7`).
- [x] Strict single-variable isolation: only genome size propagation was modified.
- [x] All 6 benchmark datasets evaluated under Var 1 and Var 2 via Slurm compute node `vnode1`.
- [x] Outputs generated: `B2_C1b3B_VAR2_CMS_CONFIG.tsv`, `B2_C1b3B_VAR2_ASSEMBLY_AB.tsv`, and `B2_C1b3B_VAR2_VERDICT.md`.
- [x] No automatic merge to `main`, no force push, no releases created.
- [x] Execution halted immediately; awaiting user review before proceeding.
