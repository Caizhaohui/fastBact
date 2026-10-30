# F26 Stage D — Cross-Species Baseline Benchmark Report

**Document ID:** `F26-REP-001`  
**Dispatch:** `DISP-F26-018`  
**Date:** 2026-10-30  
**Phase:** F26 — Cross-Species Assembly Benchmarking, Robustness Hardening & Competitive Validation  
**Stage:** D — Preliminary Cross-Species Baseline Evaluation  
**Author:** Developer Worker (Cursor Agent)  
**Governing contract:** [`F26_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md`](../../F26_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md)  
**Metric source:** `benchmarks/f26_results/summary.json` (`schema_version` `f26-benchmark-1`, `n_runs` 16, collated by `scripts/collate_benchmark_results.py`)  
**Companion table:** `benchmarks/f26_results/summary.tsv`  
**Dataset registry:** `F26_DATASET_MANIFEST.tsv` (read only; not modified)

This report records the completed replicate-1 matrix. It does not certify Stage E coverage sensitivity, Stage G three-replicate determinism, Stage H independent acceptance, or re-execution of the frozen F25 *E. coli* oracles.

---

## 1. Document Metadata and Executive Summary

Stage D is complete for the pre-registered cross-species panel. The collated matrix contains **16/16 runs** with `status` `PASS`, assembler `exit_code` **0**, and `quast_exit_code` **0**.

Assemblers, as recorded in `tool_version_text`:

| Tool key | Recorded version | Role |
|---|---|---|
| `fastbact` | `fastBact 0.1.0` | Frozen F25 baseline. Source commit `02872f81811402d2bd603ff7c0a2ecf8b3dbe528` (`02872f8`). `git diff 02872f8 -- src/` is empty. All four fastBact runs used binary SHA-256 `517402653fcefc625cb79912ade7e7e382b0d05343bbeb483d7c82dd807947d1`. |
| `spades` | `SPAdes genome assembler v4.2.0` | Competitor |
| `skesa` | `SKESA 2.5.1` | Competitor |
| `unicycler` | `Unicycler v0.5.1` | Competitor |

Evaluation is Illumina paired-end only, 8 threads, QUAST invoked as `/hpcfs/fhome/caizhh/.local/bin/quast.py` with `--min-contig 200`. QUAST v5.2.0 is the version pinned by `F26_COMPETITOR_BENCHMARK_CONTRACT.md` for that executable. `summary.json` does not store a separate QUAST version string.

The four evaluated taxa are the F26 cross-species rows of `F26_DATASET_MANIFEST.tsv`. GC percentages below are the manifest values.

| Dataset | Organism | Gram | GC (%) | Genome (bp) | SRA | Reference |
|---|---|---|---:|---:|---|---|
| CS-BS01 | *Bacillus subtilis* subsp. subtilis str. 168 | Positive | 43.51 | 4,215,606 | SRR40739147 | `NC_000964.3` |
| CS-PA01 | *Pseudomonas aeruginosa* PAO1 | Negative | 66.56 | 6,264,404 | SRR41244734 | `NC_002516.2` |
| CS-SA01 | *Staphylococcus aureus* MW2 | Positive | 32.83 | 2,820,462 | SRR33682404 | `NC_003923.1` |
| CS-SE01 | *Salmonella enterica* subsp. enterica serovar Typhimurium str. LT2 | Negative | 52.22 | 4,857,450 | SRR41010285 | `NC_003197.2` |

At one decimal place the manifest GC values are 43.5%, 66.6%, 32.8%, and 52.2%. The panel spans four genera and both Gram categories. Every run in this matrix is **replicate 1**. The `determinism` array has 16 entries, each with `checked: false`, `pass: null`, and `reason: "not enough replicates: 1 < 3"`.

Collation provenance in `summary.json` `context.git_head` is `3cc7e27188150e7c1b9304bf435f8af285b598ba`. That hash is the worktree HEAD at collation time. It is not the frozen assembler source. fastBact source identity for this baseline is `02872f8`, as stated above.

---

## 2. Complete 16-Run Cross-Species Benchmark Matrix

Primary measurements are the JSON literals in `benchmarks/f26_results/summary.json`. Two wall-clock values are stored with binary-float expansions: CS-BS01 fastBact `153.07999999999998` and CS-SE01 fastBact `106.27000000000001`.

Derived columns are not stored in the JSON. They are computed from those literals:

- **Speedup vs fastBact** = `wall_clock_s(tool) / wall_clock_s(fastBact)` on the same dataset. The fastBact row is 1.000. A larger value means that tool took longer; fastBact was that many times faster.
- **Memory reduction vs fastBact (%)** = `100 × (peak_rss_mb(tool) − peak_rss_mb(fastBact)) / peak_rss_mb(tool)`. The fastBact row is 0.00. A positive value means fastBact used less RSS than that tool, as a percentage of the tool's own peak RSS. A negative value means fastBact used more RSS.

Speedup is rounded to 3 decimal places and memory reduction to 2 decimal places. Rounding is half-up from the exact float ratio.

Extensive misassemblies are `metrics.misassemblies`. Local misassemblies are `metrics.local_misassemblies`. Base-error rates are per 100 kbp.

| Dataset | Tool | Status | Wall-clock (s) | Speedup vs fastBact | Peak RSS (MB) | Memory reduction vs fastBact (%) | Genome fraction (%) | N50 (bp) | Largest contig (bp) | Extensive misassemblies | Local misassemblies | Mismatches/100 kb | Indels/100 kb |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| CS-BS01 | fastBact | PASS | 153.07999999999998 | 1.000 | 1425.3828125 | 0.00 | 98.576 | 78260 | 354932 | 0 | 0 | 0.94 | 0.82 |
| CS-BS01 | SPAdes | PASS | 552.21 | 3.607 | 5239.390625 | 72.79 | 98.946 | 624246 | 1046381 | 0 | 1 | 6.46 | 0.45 |
| CS-BS01 | SKESA | PASS | 278.44 | 1.819 | 10447.70703125 | 86.36 | 98.927 | 926769 | 1295672 | 1 | 0 | 0.22 | 0.36 |
| CS-BS01 | Unicycler | PASS | 2683.36 | 17.529 | 5245.00390625 | 72.82 | 98.903 | 2220390 | 2220390 | 0 | 1 | 0.24 | 0.41 |
| CS-PA01 | fastBact | PASS | 180.57 | 1.000 | 1527.0703125 | 0.00 | 98.74 | 23514 | 90388 | 2 | 0 | 2.34 | 0.63 |
| CS-PA01 | SPAdes | PASS | 1016.62 | 5.630 | 5245.16015625 | 70.89 | 99.516 | 571752 | 1344826 | 4 | 2 | 4.11 | 0.67 |
| CS-PA01 | SKESA | PASS | 240.12 | 1.330 | 8846.60546875 | 82.74 | 99.334 | 359176 | 850497 | 3 | 1 | 2.2 | 0.58 |
| CS-PA01 | Unicycler | PASS | 3318.3 | 18.377 | 5249.26171875 | 70.91 | 99.313 | 546789 | 1275111 | 3 | 2 | 1.7 | 0.61 |
| CS-SA01 | fastBact | PASS | 148.55 | 1.000 | 1196.7421875 | 0.00 | 98.683 | 25381 | 77704 | 2 | 0 | 2.11 | 0.25 |
| CS-SA01 | SPAdes | PASS | 287.78 | 1.937 | 4395.28515625 | 72.77 | 99.008 | 172807 | 473664 | 0 | 0 | 1.79 | 0.29 |
| CS-SA01 | SKESA | PASS | 1110.08 | 7.473 | 10136.2265625 | 88.19 | 98.655 | 175451 | 372593 | 1 | 0 | 0.32 | 0.0 |
| CS-SA01 | Unicycler | PASS | 1191.97 | 8.024 | 4410.3359375 | 72.87 | 98.785 | 343472 | 721378 | 0 | 0 | 0.32 | 0.07 |
| CS-SE01 | fastBact | PASS | 106.27000000000001 | 1.000 | 7854.46875 | 0.00 | 98.61 | 47464 | 178426 | 1 | 1 | 0.48 | 0.06 |
| CS-SE01 | SPAdes | PASS | 413.0 | 3.886 | 5237.25390625 | −49.97 | 98.895 | 192527 | 627401 | 1 | 3 | 7.5 | 0.17 |
| CS-SE01 | SKESA | PASS | 255.29 | 2.402 | 9003.1328125 | 12.76 | 98.771 | 302057 | 736780 | 1 | 2 | 0.58 | 0.1 |
| CS-SE01 | Unicycler | PASS | 2073.25 | 19.509 | 5245.00390625 | −49.75 | 98.737 | 298595 | 820290 | 1 | 2 | 0.48 | 0.08 |

NGA50 is collected by F26-P2 and is present in the same JSON object. It is omitted from the required columns above and reported in §3.4.

---

## 3. Detailed Comparative Performance Analysis

### 3.1 Runtime and throughput

fastBact wall-clock time on the four species is `106.27000000000001` s (CS-SE01) through `180.57` s (CS-PA01). That is the 106 s–180 s span of the four fastBact runs.

fastBact was faster than every competitor on every dataset. Exact competitor/fastBact wall-clock ratios, from the JSON floats:

| Dataset | Fastest competitor ratio | Slowest competitor ratio |
|---|---|---|
| CS-BS01 (*B. subtilis*) | SKESA 1.819× | Unicycler 17.529× |
| CS-PA01 (*P. aeruginosa*) | SKESA 1.330× | Unicycler 18.377× |
| CS-SA01 (*S. aureus*) | SPAdes 1.937× | Unicycler 8.024× |
| CS-SE01 (*S. enterica*) | SKESA 2.402× | Unicycler 19.509× |

Rounded to one decimal place, those ranges are 1.8×–17.5× on *B. subtilis*, 1.3×–18.4× on *P. aeruginosa*, 1.9×–8.0× on *S. aureus*, and 2.4×–19.5× on *S. enterica*. The one-decimal statement matches the exact ratios above. Unicycler is the slowest tool on all four datasets (1,191.97 s to 3,318.3 s). SKESA is the fastest competitor on three datasets and is slower than SPAdes only on CS-SA01 (1,110.08 s versus 287.78 s).

### 3.2 Memory consumption

Peak RSS is `peak_rss_mb` from GNU time, as stored in the JSON. Dividing by 1024 converts that figure to GiB.

On three of the four species, fastBact peak RSS stayed between 1,196.7421875 MB and 1,527.0703125 MB:

| Dataset | fastBact peak RSS (MB) | GiB (MB / 1024) |
|---|---:|---:|
| CS-SA01 | 1196.7421875 | 1.169 |
| CS-BS01 | 1425.3828125 | 1.392 |
| CS-PA01 | 1527.0703125 | 1.491 |

That is the 1.17–1.49 GiB band. Against SPAdes, SKESA, and Unicycler on these three datasets, fastBact used **70.89% to 88.19%** less RSS (the minimum is CS-PA01 versus SPAdes; the maximum is CS-SA01 versus SKESA). A 70%–86% summary covers the SPAdes, Unicycler, and most SKESA comparisons in this band. The CS-SA01 SKESA comparison is higher, at 88.19%.

CS-SE01 is the exception. fastBact peak RSS is **7854.46875 MB** (7.670 GiB). That is 49.97% higher than SPAdes (5237.25390625 MB) and 49.75% higher than Unicycler (5245.00390625 MB), and 12.76% lower than SKESA (9003.1328125 MB). The manifest lists CS-SE01 as the largest read set in the panel (1,006,879,350 bases, estimated full coverage about 207×). This report records the RSS observation. It does not assign a stage-level cause.

SKESA is the highest-RSS tool on every dataset (8,846.60546875 MB to 10,447.70703125 MB).

### 3.3 Structural accuracy and misassemblies

Extensive misassemblies (`metrics.misassemblies`):

| Dataset | fastBact | SPAdes | SKESA | Unicycler |
|---|---:|---:|---:|---:|
| CS-BS01 | 0 | 0 | 1 | 0 |
| CS-PA01 | 2 | 4 | 3 | 3 |
| CS-SA01 | 2 | 0 | 1 | 0 |
| CS-SE01 | 1 | 1 | 1 | 1 |

- *B. subtilis*: fastBact has **0** extensive misassemblies, tied with SPAdes and Unicycler. SKESA has 1.
- *P. aeruginosa*: fastBact has **2**, fewer than SKESA (3), Unicycler (3), and SPAdes (4).
- *S. aureus*: fastBact has **2**, which is more than SPAdes (0), Unicycler (0), and SKESA (1).
- *S. enterica*: fastBact has **1**, tied with SPAdes, SKESA, and Unicycler.

Local misassemblies for fastBact are 0, 0, 0, and 1 on CS-BS01, CS-PA01, CS-SA01, and CS-SE01. On every dataset that count is less than or equal to each competitor (competitors reach 1, 2, and 3 local events; the maximum is SPAdes on CS-SE01, with 3).

The *S. aureus* extensive-misassembly excess is a single-replicate observation. Stage D does not isolate the responsible assembly stage and does not repeat the run. It is recorded here and is not reclassified as an algorithm defect.

### 3.4 Base accuracy and contiguity

Genome fraction for all 16 runs is above 98.5%. The minimum in the matrix is fastBact on CS-BS01 at **98.576%**. fastBact genome fraction on the four species is 98.576% (CS-BS01), 98.74% (CS-PA01), 98.683% (CS-SA01), and 98.61% (CS-SE01). Competitor genome fractions are higher on most rows. The lowest competitor value is SKESA on CS-SA01 at 98.655%.

Mismatch rates versus SPAdes (`mismatches_per_100kbp`):

| Dataset | fastBact | SPAdes | SPAdes / fastBact |
|---|---:|---:|---:|
| CS-BS01 | 0.94 | 6.46 | 6.872 |
| CS-PA01 | 2.34 | 4.11 | 1.756 |
| CS-SA01 | 2.11 | 1.79 | 0.848 |
| CS-SE01 | 0.48 | 7.5 | 15.625 |

The largest fastBact advantage versus SPAdes is **15.625×** on CS-SE01 (7.5 / 0.48), which rounds to 15.6×. That ratio is the maximum, not a panel-wide result. On CS-SA01, fastBact mismatches are higher than SPAdes (2.11 versus 1.79). SKESA and Unicycler mismatch rates are lower than fastBact on CS-BS01 (0.22 and 0.24 versus 0.94) and on CS-SA01 (both 0.32 versus 2.11). On CS-SE01, fastBact (0.48) matches Unicycler (0.48) and is below SKESA (0.58) and SPAdes (7.5).

Indels per 100 kbp are likewise mixed. fastBact is lowest on CS-SE01 (0.06 versus 0.17, 0.1, and 0.08) and highest on CS-BS01 (0.82 versus 0.45, 0.36, and 0.41). Pre-registered biological differences (evolved *B. subtilis* lineage, *S. aureus* phage co-culture, *P. aeruginosa* GFP reporter, *S. enterica* ALE lineage) can contribute to the absolute mismatch and indel floor. Relative ranking on identical inputs remains the F26-P2 comparison. Those interpretation rules are in `F26_DATASET_PROVENANCE_AND_SELECTION.md` §6.1.

Contiguity is the clear fastBact tradeoff. N50 (bp):

| Dataset | fastBact | SPAdes | SKESA | Unicycler |
|---|---:|---:|---:|---:|
| CS-BS01 | 78260 | 624246 | 926769 | 2220390 |
| CS-PA01 | 23514 | 571752 | 359176 | 546789 |
| CS-SA01 | 25381 | 172807 | 175451 | 343472 |
| CS-SE01 | 47464 | 192527 | 302057 | 298595 |

NGA50 from the same JSON follows the same ordering, except where QUAST's aligned-block NGA50 differs slightly from N50: fastBact NGA50 is 78260, 23265, 25462, and 48812 on CS-BS01, CS-PA01, CS-SA01, and CS-SE01. Unicycler NGA50 on CS-PA01 is 479845, below its N50 of 546789. Competitor N50 is 4.06× to 28.37× the fastBact N50 on the same dataset. This matrix does not treat that gap as a crash or as a genome-recovery failure. Genome fraction remains above 98.5% on every fastBact run.

---

## 4. Stage F Algorithmic Hardening Evaluation

F26-P4 forbids benchmark-driven speculative patching. An algorithmic repair is allowed only when all five conditions hold:

1. A reproducible failure case is identified.
2. The responsible assembly stage is isolated.
3. A minimal failing unit or regression test is formulated.
4. Conservative evidence isolation (`p > b + H`) and loop re-entrancy stay intact.
5. An independent audit precedes merge.

Stage D evidence against opening Stage F:

- **No assembler crash.** All four fastBact runs, and all 12 competitor runs, have `status` `PASS` and `exit_code` 0. QUAST also exited 0 on every run.
- **Genome recovery stays high.** fastBact genome fraction is 98.576%, 98.74%, 98.683%, and 98.61%, each above 98.5%.
- **Misassembly profile, stated exactly.** Extensive counts are 0 on *B. subtilis* (tied with the best competitors), 2 on *P. aeruginosa* (lower than SKESA 3, Unicycler 3, and SPAdes 4), 1 on *S. enterica* (tied with all three competitors), and 2 on *S. aureus* (higher than SPAdes 0, Unicycler 0, and SKESA 1). Local counts are less than or equal to every competitor on all four datasets. The *S. aureus* extensive excess is one replicate, with no isolated stage and no failing test. It does not meet condition 1–3.
- **Frozen logic was not modified.** `src/` is unchanged relative to `02872f8`. The four fastBact runs share one binary digest. The accepted F25 implementation of `p > b + H` is therefore the code that ran. These QUAST tables do not add a new dynamic trace of that inequality. They also show no evidence that the invariant was loosened: no source edit was introduced to obtain these results.
- **Contiguity and the CS-SE01 RSS outlier are not repair tickets.** Lower N50 is the observed cost of the frozen conservative command (`--min-len 500`, `--kmers auto`, `--depth 120`, subsampling off). Raising N50 by retuning those parameters, or by changing resolver thresholds, would be speculative. The CS-SE01 RSS increase is an observation without a localized failing stage.

**Conclusion:** speculative algorithmic changes and parameter retuning are not necessary and are not scientifically justified for Stage F on this Stage D matrix. Stage F stays closed unless a later replicate set produces a reproducible failure that can be isolated and encoded as a minimal regression test under F26-P4. The *S. aureus* extensive-misassembly count (2 versus 0/0/1) and the panel-wide N50 gap remain open observations for Stage G replication. They are not authorization to patch the frozen baseline.

---

## 5. Verification and Acceptance Criteria Mapping

Stage D in contract §3 is: "Preliminary Cross-Species Baseline Evaluation (Evaluating frozen baseline on new species)."

| Criterion | Contract location | Evidence in this matrix | Disposition |
|---|---|---|---|
| Evaluate the frozen baseline on the new species before any code change | Stage D; F26-P4 ordering | `src/` diff against `02872f8` is empty. Binary SHA-256 `517402653fcefc625cb79912ade7e7e382b0d05343bbeb483d7c82dd807947d1` on all four fastBact runs. | **Met** |
| At least three non-*E. coli* species, at least two genera, both Gram categories | F26-P1 | Four species, four genera, Gram-positive CS-BS01 and CS-SA01, Gram-negative CS-PA01 and CS-SE01. | **Met** by the evaluated panel |
| SRA/ENA accession, BioSample, Illumina paired-end layout, complete reference | F26-P1 | Manifest columns match the four runs. This report does not recompute reference or FASTQ hashes. | **Identity matched; hashes not re-audited here** |
| Head-to-head versus SPAdes v4.2.0, SKESA v2.5.1, Unicycler v0.5.1 | F26-P2 | `tool_version_text` matches those versions on every competitor run. fastBact reports `fastBact 0.1.0`. | **Met** |
| Short-read only, 8 threads, memory cap | F26-P2; benchmark contract §2 | Recorded argv uses `-t 8` or `--cores 8`. SPAdes and Unicycler receive `-m 14`. SKESA receives `--memory 14`. | **Met** as recorded in `command` |
| QUAST `--min-contig 200` and the F26-P2 metric set | F26-P2 | Every `quast_command` includes `--min-contig` `200`. Each run has N50, NGA50, genome fraction, extensive misassemblies, local misassemblies, mismatches/100 kb, indels/100 kb, wall time, and peak RSS. `quast_exit_code` is 0 on 16/16 runs. | **Met** |
| Complete preliminary matrix | Stage D | 4 datasets × 4 tools × replicate 1 = 16 runs, all `PASS`. | **Met** |
| No speculative hardening from this baseline | F26-P4; Stage F | §4. Conditions 1–3 for a repair are not met. | **Stage F not indicated** |
| Coverage gradients (30×, 50×, 100×, 150×, full) | F26-P3; Stage E | Not present in `summary.json`. | **Out of scope** |
| Three-replicate determinism, `σ² = 0` | Stage G; benchmark contract | `determinism[].checked` is false because `1 < 3`. | **Not claimed** |
| F25 *E. coli* oracle non-regression (`Δ = 0` on R1-001, R1-003, R1-001-Val, R1-003-Val) | Contract §2 | Those four datasets are absent from this 16-run file. | **Not re-executed in Stage D** |
| Independent scientific acceptance | Stage H | Not this document. | **Not claimed** |

**Stage D disposition:** `STAGE_D_REPLICATE1_MATRIX_COMPLETE`. Sixteen of sixteen collated runs passed. The frozen `02872f8` baseline recovered more than 98.5% genome fraction on all four species, with no algorithmic crash and no source change. Stage F intervention is not justified by these results. Determinism, coverage sensitivity, *E. coli* oracle replay, and final acceptance remain later gates.
