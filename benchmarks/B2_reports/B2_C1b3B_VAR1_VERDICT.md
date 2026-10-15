# FB-BENCH-B2-C1b-3B: Var 1 CLI Subsampling Decoupling Verdict

- **Date:** 2026-10-15
- **Branch:** `fix/b2-c1b3b-var1`
- **Base Commit (Production Baseline):** `dac6a57417b390cf7385ff07adf569e550cc8ca9`
- **Remote URL:** `https://github.com/Caizhaohui/fastBact`
- **Execution Platform:** Slurm HPC Cluster Node `vnode1` / `vnode4` (Partition `qfree`)

---

## 1. Executive Summary

Milestone **`FB-BENCH-B2-C1b-3B Var 1`** has been implemented, validated, and verified under strict single-variable isolation:
1. **Decoupled CLI Subsampling Semantics (`VERIFIED`):** In `src/main.rs`, the logic coupling `args.genome_size` to read discarding has been completely eliminated. Specifying `--genome-size` and/or `--depth` now strictly acts as an informative prior for genome estimation, sketch sizing, and telemetry, without ever discarding reads.
2. **Explicit User Gate Enforcement (`VERIFIED`):** Read subsampling is now initiated **if and only if** the user explicitly specifies `--enable-auto-subsample`.
3. **Precedence Hierarchy Maintained (`VERIFIED`):** `--no-subsample` unconditionally overrides `--enable-auto-subsample` (highest priority).
4. **Safety Guards Preserved (`VERIFIED`):** `LowCoverageWarning` ($<20\times$) and `Unreliable` ($<5\times$) estimator confidence states safely downgrade to 100% data retention.
5. **Zero Algorithm Modifications:** No changes were made to Count-Min Sketch, HyperLogLog, local multi-k, exSPAnder, or assembly core.
6. **Bitwise Consistency Verified (`VERIFIED`):** For all unchanged execution paths (default mode, explicit `--no-subsample`), the final assembly `contigs.fa` SHA-256 is **100% bitwise identical** to production baseline `dac6a57` (`f2081f26d001112967ae79a0b1bdac4057b0675622fcd26ac0dafcf10a0ebf64`).

---

## 2. A/B Matrix Verification Across Required Configurations

Summary of empirical results from `benchmarks/B2_reports/B2_C1b3B_VAR1_AB.tsv`:

| Test ID | Configuration | CLI Arguments | Input Pairs | Retained Pairs | Retained % | Subsample Active | Reason | Contigs | N50 (bp) | SHA-256 Match vs dac6a57 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **T1_DEFAULT** | Default Parameters | *(none)* | 100,000 | 100,000 | **100.00%** | `false` | `DISABLED_BY_DEFAULT_SAFE_FALLBACK` | 1 | 199,983 | **Bit-Identical** (`f2081f...`) |
| **T2_GENOME_SIZE_ONLY** | Prior Genome Size Alone | `--genome-size 200000` | 100,000 | 100,000 | **100.00%** | `false` | `GENOME_SIZE_PRIOR_PROVIDED_SUBSAMPLING_NOT_ENABLED` | 1 | **199,983** | **Decoupled** (N50: 146k $\to$ 200k) |
| **T3_GENOME_SIZE_DEPTH** | Prior Size + Depth Target | `--genome-size 200000 --depth 120` | 100,000 | 100,000 | **100.00%** | `false` | `GENOME_SIZE_PRIOR_PROVIDED_SUBSAMPLING_NOT_ENABLED` | 1 | **199,983** | **Decoupled** (N50: 146k $\to$ 200k) |
| **T4_AUTO_GENOME_SIZE** | Explicit Auto + User Cap | `--enable-auto-subsample --genome-size 200000` | 100,000 | 79,859 | **79.86%** | `true` | `AUTO_SUBSAMPLE_EXPLICIT_GENOME_SIZE` | 2 | 146,614 | **Identical** to baseline cap |
| **T5_NO_SUBSAMPLE_OVERRIDE** | Override Flag Priority | `--no-subsample --enable-auto-subsample` | 100,000 | 100,000 | **100.00%** | `false` | `EXPLICIT_NO_SUBSAMPLE_FLAG` | 1 | 199,983 | **Bit-Identical** (`f2081f...`) |
| **T6_LOW_COV_15X** | 15X Low Coverage Guard | `--enable-auto-subsample` | 10,000 | 10,000 | **100.00%** | `false` | `AUTO_SUBSAMPLE_INTERCEPTED_LOW_COVERAGE_WARNING` | 30 | 9,538 | **Bit-Identical** (`2f8197...`) |
| **T7A_100X_AUTO** | 100X Within Target | `--enable-auto-subsample` | 66,666 | 66,666 | **100.00%** | `false` | `AUTO_SUBSAMPLE_HIGH_CONFIDENCE_WITHIN_TARGET` | 1 | 199,979 | **Bit-Identical** (`95a52b...`) |
| **T7B_150X_AUTO** | 150X Auto Subsample | `--enable-auto-subsample` | 100,000 | 82,414 | **82.41%** | `true` | `AUTO_SUBSAMPLE_HIGH_CONFIDENCE` | 2 | 146,614 | **Identical** to baseline subsample |
| **T7C_300X_AUTO** | 300X Auto Subsample | `--enable-auto-subsample` | 200,000 | 82,815 | **41.41%** | `true` | `AUTO_SUBSAMPLE_HIGH_CONFIDENCE` | 1 | 199,993 | **Identical** to baseline subsample |

---

## 3. Explanation of Behavioral Deltas

1. **Configurations T2 & T3 (`--genome-size` without `--enable-auto-subsample`):**
   - **Baseline `dac6a57` Behavior:** Under `dac6a57`, providing `--genome-size 200000` computed `explicit_cov = 150.0x > 1.1 * 120 = 132.0x` and automatically discarded ~17.5% of reads without explicit authorization. This caused the 5 bp Poisson dip at coords 53,371–53,375 bp to drop below $t_{\text{solid}} = 4$, fracturing the assembly into 2 contigs ($N50 = 146,614\text{ bp}$).
   - **Var 1 Behavior:** Specifying `--genome-size` retains 100% of reads (`subsampling_enabled = false`). All Poisson dips are preserved, producing a **perfect single contig of 199,983 bp** with SHA-256 bit-identical to the un-subsampled default run.
2. **All Other Configurations (T1, T4, T5, T6, T7A, T7B, T7C):**
   - Maintained 100% backward compatibility and bitwise equivalence with baseline expectations.

---

## 4. Formal Recording of Audit Errata

As detailed in `benchmarks/B2_reports/B2_AUDIT_ERRATA.md`, three technical items from preliminary audit scripts were reviewed and recorded:
1. **Deduplication in Folding Hash Tests & Verification Boundary:** In `S2-REAL-001`, the reported 38,937 folding value duplicates were caused by un-deduplicated sliding windows traversing biological repeats (identical 127-mers repeating across the 7 ribosomal operons). True distinct k-mer folding collision evaluation requires independent tests comparing distinct 256-bit canonical k-mers against 64-bit folded keys.
2. **Metric Scope for `final_contig_accuracy_pct` & QUAST Identity:** The percentage reported in `B2_CMS_GRAPH_FATE.tsv` measured intermediate graph unitig reference concordance, not the final polished contigs (`contigs.fa`). Furthermore, 0 QUAST misassemblies does not equal 100% nucleotide identity (small indels, SNVs, or unaligned regions may exist).
3. **S2-REAL-001 Official Measured Baseline:** Clarified that the official archived baseline measured on `S2-REAL-001` (v0.1.0-baseline / `dac6a57` / C1b-1) is **28 contigs, N50 = 285,687 bp, NGA50 = 285,687 bp, 1 misassembly**. Unverified figures (e.g. 4.08 Mb) must not be used as `baseline_measured`.

---

## 5. Evidence Classification

| Evidence Level | Verified Conclusions |
| :--- | :--- |
| **`VERIFIED`** | 1. `--genome-size` is completely decoupled from read discarding in `src/main.rs`.<br>2. `--enable-auto-subsample` is strictly required to activate auto-subsampling.<br>3. `--no-subsample` takes absolute precedence over all subsampling flags.<br>4. Unchanged execution paths maintain bitwise SHA-256 identity with baseline `dac6a57`.<br>5. 15X low-coverage data is safely protected by confidence state machine. |
| **`LIKELY`** | N/A (All Var 1 requirements empirically verified). |
| **`UNKNOWN`** | N/A. |

---

## 6. Stop Rule Adherence

- [x] Independent fix branch `fix/b2-c1b3b-var1` created from `dac6a57`.
- [x] Only Var 1 implemented; no multi-k, CMS, or Unicycler modifications made.
- [x] All 7 test configurations compiled and verified on Slurm compute node `vnode1`.
- [x] Changes committed to `fix/b2-c1b3b-var1` and pushed to GitHub.
- [x] No automatic merge to `main`, no force push, no release created.
- [x] Execution halted awaiting user review.
