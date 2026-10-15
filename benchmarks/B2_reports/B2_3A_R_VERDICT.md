# FB-BENCH-B2-C1b-3A-R: Final Audit Integrity Verdict & Milestone B2-C1b-3B Gateway

- **Date:** 2026-10-15
- **Audit Branch:** `audit/b2-c1b-3a`
- **Production Baseline (`master` / `main`):** `dac6a57417b390cf7385ff07adf569e550cc8ca9`
- **Remote Repository:** `https://github.com/Caizhaohui/fastBact`
- **Execution Platform:** Slurm HPC Cluster Node `vnode1` (Partition `qfree`)

---

## 1. Audit Integrity & Forensic Resolution Summary

The integrity audit **`FB-BENCH-B2-C1b-3A-R`** has rigorously resolved all discrepancies identified in the preliminary reports:
1. **Hardcoded Labels Eliminated (`VERIFIED`):** The initial copy-paste string `"NO (Complete single contig)"` in the 300X sensitivity analysis has been replaced with dynamic evaluation based on actual contig arrays, exact reference alignment, and N50 metrics.
2. **Poisson 5-Kmer Window Demystified (`VERIFIED`):** Auditing the full junction (53,371..53,497 bp) demonstrated that single-position depth is insufficient to predict graph continuity. In Seed 1 and Seed 100, while coordinate 53,371 bp had depth $\ge 4$, intermediate overlapping 127-mers in the 5 bp window dropped to depth $3\text{X} < t_{\text{solid}} = 4$, severing the de Bruijn path.
3. **Prior Contig Backbone Inheritance Proven (`VERIFIED`):** In Seeds 42 and 2026, the junction was fully resolved at $k=63$ (where depth is $>15\times$). Although interior reads were subsequently filtered by `filter_active_reads`, the sequence was preserved into $k=127$ through prior contig backbone inheritance (+10 count in CMS, +1000 edge support).
4. **`upgrade_multik` Early-Return Behavior Documented (`VERIFIED`):** `upgrade_multik` includes an explicit guard (`if !has_branching { return; }`). On unbranched graphs at $k=63$, high-k upgrades are safely bypassed, retaining unfragmented contigs.
5. **Pre-CMS Key Folding Collisions Ruled Out (`VERIFIED`):** Testing `KmerInt::to_u64()` on all reference 127-mers across 200 kb, 2.5 Mb, 4.65 Mb, and on 1,000,000 synthetic random 127-mers yielded **strictly 0 folding collisions**, consistent with the theoretical birthday paradox expectation ($P < 6.8 \times 10^{-7}$).
6. **Interface Standardization Completed (`VERIFIED`):** All fix designs and A/B matrices now strictly reference the verified production interfaces: `args.enable_auto_subsample`, `args.genome_size`, `args.depth`, `upgrade_multik()`, and `build_cdbg_with_contigs_ext()`. S2-REAL-001 expected metrics are properly bounded at its 5kb `rrn` ribosomal repeats (N50 ~4.08 Mb, NGA50 ~2.90 Mb, contigs 2–3).

---

## 2. Complete Evidence Classification

| Category | Finding / Metric | Status | Evidence Level |
| :--- | :--- | :--- | :--- |
| **CMS Recall** | False negative rate across all datasets and scales | $0.0000\%$ ($FN = 0$) | **`VERIFIED`** |
| **CMS 38% Precision** | Caused by 1.76M collision FPs and 879k recurrent sequencing error solids under 300X overload | 38.0460% exact precision | **`VERIFIED`** |
| **CMS Graph Fate** | Dual-solid requirement and simplification prune 100% of collision spurs; final contigs have $>99.9\%$ accuracy | Contig accuracy verified | **`VERIFIED`** |
| **150X Fracture** | 5 bp Poisson dip (53,371..53,375 bp) dropping to 2X–3X ($< t_{\text{solid}}=4$) under factor 0.8257 | 10-seed sensitivity trace | **`VERIFIED`** |
| **Folding Collisions** | Pre-CMS `to_u64()` folding collisions on bacterial genomes | 0 collisions ($0.0000\%$) | **`VERIFIED`** |
| **CLI Parameter Coupling**| `args.genome_size` unintentionally triggers subsampling when coverage $> 1.1 \times \text{target}$ | Code audit verified | **`VERIFIED`** |
| **Trigger Policy 3A** | Raising trigger to $\text{cov} > 1.5 \times \text{target}$ ($180\times$) avoids 150X dips without performance penalty | Supported by depth trace | **`LIKELY`** |
| **High-Error Long Reads** | 4-bit CMS behavior on raw Nanopore/PacBio reads ($>2\%$ error rate) | Needs further study | **`UNKNOWN`** |

---

## 3. Recommended Minimal Implementation for Milestone `B2-C1b-3B`

Upon user authorization to commence implementation, execute the following isolated steps sequentially:

1. **Phase 1 (Variable 1 - CLI Decoupling):**
   In `src/main.rs`: Ensure `args.genome_size` only overrides sketch sizing (`effective_g`), and read subsampling requires explicit `--enable-auto-subsample`.
   *Gate:* Run 150X with `--genome-size 200000` (without `--enable-auto-subsample`) and verify 100% data retention and bitwise-identical assembly ($N50 = 199,983\text{ bp}$).
2. **Phase 2 (Variable 2 - High-k Genome Size Propagation):**
   In `src/graph/local_multik.rs`: Provide `upgrade_multik_ext` accepting `trusted_genome_size: Option<u64>` and pass it to `build_cdbg_with_contigs_ext`.
   *Gate:* Verify memory allocation and unitig counts on 2.5 Mb and 4.65 Mb datasets.
3. **Phase 3 (Variable 3A - Conservative Trigger Threshold):**
   In `src/main.rs`: Adjust the subsampling trigger condition to `est_cov > 1.5 * args.depth` (i.e. $> 180.0\times$).
   *Gate:* Verify that across all 10 random seeds on 150X, 100% of runs produce a complete, unbroken contig ($N50 = 199,983\text{ bp}$, contigs = 1).
4. **Phase 4 (Full Regression Suite):**
   Execute full integration benchmarks across all 8 datasets to verify zero misassemblies and close milestone `B2-C1b`.

---

## 4. Stop Rule Compliance

- [x] All compute jobs executed exclusively on Slurm node `vnode1`.
- [x] Production baseline `dac6a57` on `main`/`master` left untouched.
- [x] Audit branch `audit/b2-c1b-3a` synchronized and ready for final push.
- [x] No automatic advance to `B2-C1b-3B` implementation or Unicycler `U-series`.
- [x] Execution halted awaiting user review.
