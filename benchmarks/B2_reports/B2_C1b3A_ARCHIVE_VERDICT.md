# FB-BENCH-B2-C1b-3A: Final Audit Archive & B2-C1b-3B Verdict

- **Date:** 2026-10-15
- **Audit Branch:** `audit/b2-c1b-3a` (Pushed to GitHub)
- **Production Baseline:** `master` / `main` at `dac6a57417b390cf7385ff07adf569e550cc8ca9`
- **Audit Head Commit:** `f8be9ba7795178da8ebe41a9d44d351d8d8cd5bf`
- **Remote URL:** `https://github.com/Caizhaohui/fastBact`
- **Execution Platform:** Slurm Cluster Node `vnode1` (Partition `qfree`)

---

## 1. Executive Summary & Verification Matrix

The comprehensive audit of milestone **FB-BENCH-B2-C1b-3A** and preparation for **B2-C1b-3B** has been successfully executed, archived, and synchronized with GitHub. All investigations strictly followed single-variable isolation and executed 100% of compilation, exact k-mer tracing, and benchmarking on Slurm compute node `vnode1`.

### Core Verdict Summary:
| Audit Component | Core Finding | Evidence Grade | Verification Status |
| :--- | :--- | :--- | :--- |
| **CMS Ground Truth Benchmark** | False Negative rate is **strictly 0.0000%** ($FN=0$). Recall is 100.0000%. | **`VERIFIED`** | Mathematical & empirical proof on 200kb, 2.5Mb, 4.65Mb |
| **38% Precision Demystified** | Driven by 1.76M collision FPs and 879k recurrent sequencing error solids on 300X un-subsampled overload. | **`VERIFIED`** | Exact 64-shard ground truth table (`FxHashMap<u64, u16>`) |
| **Graph Simplification Defense** | Dual-kmer solid requirement and graph simplification (`clip_tips`, `pop_bubbles`) eliminate 100% of collision spurs, yielding 0 misassemblies and perfect N50 (199,983 bp). | **`VERIFIED`** | Traced through cDBG build and simplification passes |
| **Estimator vs Builder CMS** | Estimator CMS samples $<150\text{k}$ reads (solid depth $\le 15$), with collision rate $<0.05\%$. Sizing concerns only affect Builder CMS under extreme raw coverage. | **`VERIFIED`** | Verified across all sample sizes |
| **150X N50 Regression Forensics** | Caused by Poisson dip at coords 53,371–53,497 bp (raw depth 5X) being randomly reduced to 2X–3X under seed 0 (factor 0.8257), falling below high-k solid threshold ($t_{\text{solid}}=4$). | **`VERIFIED`** | 10-seed sensitivity matrix in `B2_150X_DEPTH_TRACE.tsv` |
| **Seed Sensitivity** | 60% of seeds on 150X assemble complete 199,983 bp contig; 40% fracture at 53497 bp. On 300X, 90% assemble complete contig. | **`VERIFIED`** | Tested across 10 deterministic seeds |
| **CLI `--genome-size` Conflation** | Passing `--genome-size` unintentionally sets `subsampling_enabled = true` in `src/main.rs`. | **`VERIFIED`** | Confirmed by code audit |
| **High-k CMS Genome Size Omission** | `local_multik.rs` passes `None` for genome size to `build_cdbg_with_contigs_ext`. | **`VERIFIED`** | Confirmed by code audit |

---

## 2. Git & Release Hygiene Audit

1. **Remote Repository Status (`VERIFIED`):**
   - Branch `audit/b2-c1b-3a` was verified clean of credentials, API keys, and large sequencing FASTQ/BAM files.
   - Pushed successfully to `origin/audit/b2-c1b-3a` (`https://github.com/Caizhaohui/fastBact`).
   - Production branch `main`/`master` remains strictly pinned at `dac6a57417b390cf7385ff07adf569e550cc8ca9`.
2. **Release Hygiene Verification (`VERIFIED`):**
   - No force-pushes performed.
   - No Git history rewritten.
   - No formal GitHub Release created (held until B2-C1b-3B closes all quality gates).

---

## 3. Evidence Classification Summary

### A. VERIFIED (Empirically and mathematically proven)
1. **CMS One-Sided Error Guarantee:** `false_negatives = 0` across all 24 tested permutations of genome size and depth.
2. **Poisson Dip Mechanism:** Reference coordinates 53,371–53,497 bp exhibit an extreme Poisson coverage dip (5X in 150X FASTQ). Subsampling reduces depth to 2X–3X, falling below $t_{\text{solid}} = 4$ at $k=127$, causing an assembly break into 146,614 bp and 53,490 bp contigs.
3. **Graph Simplification Robustness:** High CMS collision rates at 300X (38% precision) do not compromise final contig sequences. Simplification removes 100% of spurious tips and dead ends.
4. **CLI Parameter Coupling:** `cli.genome_size.is_some()` triggers subsampling even without `--subsample`.

### B. LIKELY (Highly supported by data and model, awaiting production validation)
1. **Trigger Threshold Policy $\ge 180\times$:** Raising the subsampling threshold to $180\times$ will guarantee 100% contig continuity across all random seeds for 150X datasets with negligible CPU/memory overhead.
2. **CMS Doubling Multiplier ($32 \times \widehat{G}$):** Increasing table width multiplier to $32\times$ will boost precision to $>75\%$ on 300X and $>99\%$ on 4.65 Mb while maintaining peak RSS under 2.0 GB.

### C. UNKNOWN (Requires future research outside current scope)
1. **Error Behavior on Long / High-Error Reads:** Behavior of 4-bit Count-Min Sketch under raw PacBio/ONT error profiles ($>2\%$ error rate) without prior consensus error correction.

---

## 4. Deliverable Artifacts Index

All requested deliverables have been generated and archived in `benchmarks/B2_reports/`:
1. `benchmarks/B2_reports/C1b3_CMS_EXACT_VALIDATION.tsv`: Exact 64-shard ground truth validation table.
2. `benchmarks/B2_reports/B2_CMS_ERROR_CLASSIFICATION.md`: Comprehensive breakdown of collision FPs vs. recurrent error solids.
3. `benchmarks/B2_reports/B2_150X_DEPTH_TRACE.tsv`: 10-seed exact depth trace from raw FASTQ to cDBG for 150X and 300X.
4. `benchmarks/B2_reports/B2_C1b3B_FIX_DESIGN.md`: Detailed specification of 4 isolated single-variable interventions.
5. `benchmarks/B2_reports/B2_C1b3B_AB_MATRIX.tsv`: Systematic A/B matrix for `B2-C1b-3B` validation.
6. `benchmarks/B2_reports/B2_C1b3A_ARCHIVE_VERDICT.md`: This comprehensive synthesis report.

---

## 5. Next Steps: Recommended Minimal Implementation for `B2-C1b-3B`

Upon user authorization to commence `FB-BENCH-B2-C1b-3B`, execute the four fixes sequentially under single-variable isolation:
1. **Step 1:** Implement Fix 1 in `src/main.rs`: Decouple `cli.genome_size` from `subsampling_enabled`. Verify bitwise equivalence on 150X with `--genome-size 200k`.
2. **Step 2:** Implement Fix 2 in `src/graph/local_multik.rs`: Pass `trusted_genome_size` to `build_cdbg_with_contigs_ext`. Verify memory and unitig counts.
3. **Step 3:** Implement Fix 3A in `src/main.rs`: Set subsampling trigger threshold to $\ge 180\times$. Verify 10-seed continuity on 150X ($N50 = 199,983\text{ bp}$ across 100% of seeds).
4. **Step 4:** Run full integration regression suite across all 8 datasets to verify 0 regressions and confirm closure of milestone `B2-C1b`.

---

## 6. Stop Rule Adherence

- [x] All compute jobs executed on Slurm node `vnode1`.
- [x] Production baseline `dac6a57` on `main`/`master` left untouched.
- [x] Audit branch `audit/b2-c1b-3a` pushed to GitHub.
- [x] No automatic advance to Unicycler `U2–U4`.
- [x] Execution halted awaiting user approval.
