# F26 Stage E & Stage G — Coverage Robustness and Determinism Report

**Document ID:** `F26-REP-002`  
**Date:** 2026-10-30  
**Phase:** F26 — Cross-Species Assembly Benchmarking, Robustness Hardening & Competitive Validation  
**Stages:** E — Coverage Robustness and Sensitivity; G — fastBact three-replicate determinism on the full-depth cross-species panel  
**Author:** Implementation and scientific reporting Worker (Cursor Agent)  
**Governing contract:** [`F26_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md`](../../F26_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md)  
**Stage D comparator:** [`F26_CROSS_SPECIES_BASELINE_REPORT.md`](./F26_CROSS_SPECIES_BASELINE_REPORT.md)  
**Metric source:** `benchmarks/f26_results/summary.json` (`schema_version` `f26-benchmark-1`, `context.n_runs` 40, collated by `scripts/collate_benchmark_results.py`)  
**Companion table:** `benchmarks/f26_results/summary.tsv`  
**Depth manifest:** `F26_STAGE_E_MANIFEST.tsv`  
**Dataset registry:** `F26_DATASET_MANIFEST.tsv`  
**Subsample log:** `logs/f26_stage_e_subsamples_3075001.err` (Slurm job `3075001`)

This report records the fastBact coverage gradient and the fastBact N=3 full-depth determinism check. It does not replace the Stage D competitor matrix, it does not assemble the historical *E. coli* depth FASTQs, and it does not issue Stage H acceptance.

---

## 1. Executive Summary and Verification Scope

### 1.1 Scope

The collated file contains **40/40** runs, all with `status` `PASS`, assembler `exit_code` **0**, and `quast_exit_code` **0**:

| Block | Contents | Runs |
|---|---|---:|
| Stage D | 4 cross-species datasets × 4 tools × replicate 1 | 16 |
| Stage E | 4 cross-species datasets × depths 30×, 50×, 100×, 150× × fastBact replicate 1 | 16 |
| Stage G | fastBact replicates 2 and 3 on the four full-depth datasets | 8 |

Stage G certification uses all three full-depth fastBact replicates (the Stage D replicate 1 plus the eight additional runs). Competitor tools remain replicate 1. Their Stage D metrics are unchanged and are not repeated here.

The evaluated taxa are the four F26 cross-species rows. GC percentages are the `F26_DATASET_MANIFEST.tsv` values. One-decimal rounding, as used in Stage D, is 43.5%, 66.6%, 32.8%, and 52.2%.

| Dataset | Organism | Gram | GC (%) | Genome (bp) | SRA | Reference |
|---|---|---|---:|---:|---|---|
| CS-BS01 | *Bacillus subtilis* subsp. subtilis str. 168 | Positive | 43.51 | 4,215,606 | SRR40739147 | `NC_000964.3` |
| CS-PA01 | *Pseudomonas aeruginosa* PAO1 | Negative | 66.56 | 6,264,404 | SRR41244734 | `NC_002516.2` |
| CS-SA01 | *Staphylococcus aureus* MW2 | Positive | 32.83 | 2,820,462 | SRR33682404 | `NC_003923.1` |
| CS-SE01 | *Salmonella enterica* subsp. enterica serovar Typhimurium str. LT2 | Negative | 52.22 | 4,857,450 | SRR41010285 | `NC_003197.2` |

`F26_STAGE_E_MANIFEST.tsv` also contains 16 *E. coli* depth rows (R1-001, R1-003, R1-001-Val, R1-003-Val at 30×, 50×, 100×, and 150×). Job `3075001` generated those FASTQ pairs. None of them appear in `summary.json`. This report does not assign them assembly metrics and does not replay the frozen F25 oracles.

Collation `context.git_head` is `c23a6f96e00ece7bdc7389ed5a756700afb06238`. That is the worktree HEAD at collation. It is not the frozen assembler source. `git diff 02872f81811402d2bd603ff7c0a2ecf8b3dbe528 -- src/` is empty. Every fastBact run in the per-run raw summaries records binary SHA-256 `517402653fcefc625cb79912ade7e7e382b0d05343bbeb483d7c82dd807947d1` and `tool_version_text` `fastBact 0.1.0`.

### 1.2 Findings that the 40-run file supports

- **Coverage behavior is organism-specific, and it is bounded.** On *B. subtilis*, genome fraction stays above 98.2% at 30×, 50×, 100×, 150×, and full depth (minimum `98.272` at 30×), with **0** extensive and **0** local misassemblies at every one of those five inputs. On *S. enterica*, extensive misassemblies stay at **1** and local misassemblies stay at **1** at every depth, including full.
- **Low-GC threshold on *S. aureus*.** At 30×, genome fraction is `64.582` (64.58% at two decimal places), N50 is 1,266 bp, and QUAST counts 1,713 contigs ≥200 bp. At 50×, genome fraction is `95.199` (95.20% at two decimal places). At 150× it is `98.445` (98.44% at two decimal places). Extensive and local misassemblies are **0** at 30×, 50×, 100×, and 150×. The full-depth baseline is not in that zero class: CS-SA01 has **2** extensive misassemblies and **0** local misassemblies, identical across replicates 1–3.
- **High-GC trajectory on *P. aeruginosa*.** Extensive misassemblies are 4 at 30×, 3 at 50×, and 2 at 100×. Local misassemblies are 1 at 30× and 50×, then 0 at 100×. The nominal 150× row is not an independent 150× experiment: the subsampler clamped it to the entire library (`pairs_out` = `pairs_in` = 1,518,403). Its contig SHA-256 equals the full-depth assembly, and its QUAST metrics equal the full-depth replicate 1 metrics, including **2** extensive and **0** local misassemblies.
- **Stage E runtime and memory.** The 16 depth runs finish in `35.5`–`167.36` s, all under 3 minutes. Peak RSS is `660.1015625`–`4072.14453125` MiB (0.645–3.977 GiB). Thirteen of the 16 runs are under 1.5 GiB (1,536 MiB). The three above that line are CS-BS01_150× (`1554.4453125` MiB), CS-SE01_100× (`3787.87109375` MiB), and CS-SE01_150× (`4072.14453125` MiB).
- **Full-depth memory is not inside the Stage E ceiling.** CS-SE01 full-depth peak RSS is `7854.46875`, `7777.0703125`, and `7858.8125` MiB across replicates 1–3 (7.675 GiB at replicate 3). A bound of 4.1 GiB describes the 16 depth runs. It does not describe the full-depth *S. enterica* baseline. Every *P. aeruginosa* fastBact run, depth and full, stays at or below `1537.9296875` MiB (1.502 GiB).
- **Full-depth wall clock is not uniformly under 3 minutes.** Three full-depth replicates exceed 180 s: CS-PA01 replicate 1 `180.57` s, CS-PA01 replicate 2 `185.91` s, and CS-BS01 replicate 3 `183.21` s. The maximum is `185.91` s.
- **Contig determinism.** For each of the four species, replicates 1, 2, and 3 share one `assembly_sha256`. The harness marks all four groups `checked: true`, `pass: true`, and `mismatches: []`. Stored `variance_max` is `0.0` for CS-PA01, CS-SA01, and CS-SE01. For CS-BS01 it is `2.0194839173657902e-28`. That value is the population variance of three identical `genome_fraction_pct` literals `98.576` after division by 3 in IEEE-754 arithmetic (`(x+x+x)/3` is `98.57599999999998`). The three stored metric values compare equal. Wall clock and peak RSS are outside `METRIC_FIELDS` and do vary between replicates.

No fastBact run crashed. No run in the 40-run file has a non-zero assembler or QUAST exit. `src/` is unchanged relative to `02872f8`. These results do not meet the F26-P4 conditions for an algorithmic patch.

---

## 2. Methodology

### 2.1 Deterministic paired subsampling

Stage E FASTQs were produced by Slurm job `3075001` on the `qfree` partition (`logs/f26_stage_e_subsamples_3075001.err`). The generator is `scripts/generate_stage_e_subsamples.py`, calling `scripts/subsample_fastq_pairs.py`. The log ends with `rows=32 generated=32 skipped=0`.

The selector is `random.Random(seed).sample`, then the chosen pair indices are sorted, so mate pairing and original pair order are preserved. The seed recorded in every manifest `notes` field and in every subsample command is **42**. Pair count for a target depth is

```text
k = round((target_depth * genome_size_bp) / (2 * read_len))
```

`read_len` is the first R1 sequence length when `--read-len` is omitted. When R1 lengths are not uniform, the log warns and still uses that first length. *B. subtilis* and *S. enterica* produced no such warning (151 bp and 150 bp). *S. aureus* warned and used 151 bp. *P. aeruginosa* warned and used 251 bp. Formula depth below is therefore the subsampler's depth, not a recount of every base.

When `k` exceeds the library, the tool keeps every pair and logs `clamped=yes`. That happened once in the cross-species panel: CS-PA01 at the 150× target asked for 1,871,834 pairs and the library contains 1,518,403. The keep-all path returns the pairs in file order.

`F26_STAGE_E_MANIFEST.tsv` records depth, seed, source dataset, genome size, and reference SHA-256. It does not have columns for subsampled pair counts or FASTQ SHA-256. Those values are the `[subsample]` lines of job `3075001`. `total_spots` in the manifest is the registered parent SRA spot count, not `pairs_out`.

Scanned input pairs (`pairs_in`) equal the manifest `total_spots` for CS-SA01 (2,239,653), CS-BS01 (2,253,807), and CS-PA01 (1,518,403). CS-SE01 was scanned as **3,353,233** pairs; the manifest SRA spot count is **3,359,296**. The report uses the scanned pair count as the FASTQ size that was actually subsampled. It does not reconcile the 6,063-spot difference.

Formula depth is `pairs_out * 2 * read_len / genome_size_bp`.

| Dataset | Read length used (bp) | Pairs in | Target | Pairs out | Formula depth | Clamped | R1 SHA-256 | R2 SHA-256 |
|---|---:|---:|---:|---:|---:|---|---|---|
| CS-SA01_30x | 151 | 2239653 | 30× | 280178 | 30.000× | no | `c571c4210171909b133abb01b706395dc66e3593b341e0dc48383a4be0bd9aef` | `b038cf0c0b784628783999aa458751f06d96d7c78058fee90f33e1ce9fe37dc3` |
| CS-SA01_50x | 151 | 2239653 | 50× | 466964 | 50.000× | no | `1aff39b5c4e6222ec58ac3da527e8bd67007ac053207f65a71e082bf58b5e4f5` | `1f0a7e0ddd30a2e2981d1b51778f414cc4658afa508971ca2f00d84d2d9ad604` |
| CS-SA01_100x | 151 | 2239653 | 100× | 933928 | 100.000× | no | `995805f8a2208bf8ff4d333624e4c8c2582c23c863c693314d471f1a91a1beaf` | `5c58c3c01ac5c732ee6a32294edecf84b9a2a44feb8261f9a749cd3d7bb0f31e` |
| CS-SA01_150x | 151 | 2239653 | 150× | 1400892 | 150.000× | no | `5371a04ae22444a2f513de9d0c37266a11c1cd7167c8e8c730d909a42b2ce900` | `0e44917753852311304f1f878f689ae711fd137d724f6da01f9ce302e3cecbd0` |
| CS-BS01_30x | 151 | 2253807 | 30× | 418769 | 30.000× | no | `07992363126bb11ad404713accda0bbdcd382030724d922d20641215b1163230` | `280597eb62220a3b9434ddc672936d4b9d490bcadfcc224778128e71b5a9c064` |
| CS-BS01_50x | 151 | 2253807 | 50× | 697948 | 50.000× | no | `a23be408101e0934a8ccd019bc7712789845b5fc09ce8edf226dd36c9e52757c` | `c15742f8b47460864a4dbfead372cc67b5b3143f188919b99d797a701b855f6b` |
| CS-BS01_100x | 151 | 2253807 | 100× | 1395896 | 100.000× | no | `ccce5ec23d27c147de4654269f324b53a7df4dfddc18fb2d101ad85e4ad35f47` | `550fb7762c891c0631ec656ac8adbcb6ab4e4c00b76ef74095460e15b62add93` |
| CS-BS01_150x | 151 | 2253807 | 150× | 2093844 | 150.000× | no | `479e63dd63f5dbecf23f4b7c686640c626c77e95910f102e4276b0c891375595` | `173fa4c0618bb532a705b3347764a31541ad477bab8519fd66332ed80f362845` |
| CS-PA01_30x | 251 | 1518403 | 30× | 374367 | 30.000× | no | `9c6836ae34d5e69173fda8c666742805dbee29fd64395b9d5ae4dcbd04d198b7` | `9805a6607975c55c72c4a8d9a6a42b2a2073c1a3f9575b190365280b883878f9` |
| CS-PA01_50x | 251 | 1518403 | 50× | 623945 | 50.000× | no | `93d10c4ad6e88d5d8c868ea7ef977aad038c3274926277797d15b497e1a5c76f` | `707034701a6592f9bc3fc2532566f1a7bf7f96a595ff27a03f82f56da8c723c8` |
| CS-PA01_100x | 251 | 1518403 | 100× | 1247889 | 100.000× | no | `e9c6765815301ded466e863529d990ad063c3996d451821ad0d53ab5976a191a` | `0c8ba020f0a2a29706ec51624b63d769841062a283114799ffcd80f31de195be` |
| CS-PA01_150x | 251 | 1518403 | 150× | 1518403 | 121.678× | **yes** | `b3bbcae57004d6df658a8cf433aba10c68be5831ea33c08daa0d0e11bacd59ba` | `c2652aaabbfc31da1e054de8ac2a9e8a39e9ec737f193fbc4257952c284a54ea` |
| CS-SE01_30x | 150 | 3353233 | 30× | 485745 | 30.000× | no | `fbc748dc5017b1ece10edaf062fc63778613c2785b711f9ac2b8e94379d307f6` | `089c8a15bd1c1effc8f0e206a37ce13390df957bc3867f35c448793bf74c1edc` |
| CS-SE01_50x | 150 | 3353233 | 50× | 809575 | 50.000× | no | `03fef28f07c1c7efd3bc0f701833c136491baaea12adf64d8628391968b33328` | `960b17684f01a5d18a69bff6e5c5ff2f6715af768e5beefa261f477d0446c3a8` |
| CS-SE01_100x | 150 | 3353233 | 100× | 1619150 | 100.000× | no | `06211ad17a2f7f11ba216e9aad1e70daa4c367ba19b841a605a1e3870bd158e8` | `08c12a9ca3bafd7371fb3f14a53ec7d9e05c8acb0028af6db8fe96db374212f3` |
| CS-SE01_150x | 150 | 3353233 | 150× | 2428725 | 150.000× | no | `16fa663a4e4affc6cef37b3509107ffef9f2067e11095e721ade4e77bc0298bd` | `05d91693944d467157f820b4b57c57fba54716b883ece38cdf4db73c0857b1ba` |

Full-library formula depth, using the same read length, is 239.810× (CS-SA01), 161.460× (CS-BS01), 121.678× (CS-PA01), and 207.098× (CS-SE01). The *S. aureus* and *P. aeruginosa* figures use a non-uniform first-read length.

"Actual read pairs" in §3 are these `pairs_out` values. Individual reads are twice that number because both mates are retained. Full-depth rows use `pairs_in`.

### 2.2 Benchmark harness

Each assembly is one Slurm job submitted through `scripts/submit_stage_e_benchmarks.sh` and `scripts/slurm_benchmark_harness.sbatch` on partition `qfree`, 8 CPUs, 14G memory, wall limit `03:00:00`. The harness is `scripts/run_competitor_benchmark.py`.

Every fastBact command in the 28-run set uses the same flags after the input paths:

```text
-t 8 --min-len 500 --kmers auto --depth 120
```

`--enable-auto-subsample` is absent. Under the benchmark contract and `src/main.rs`, `--depth 120` without that flag is a prior and does not discard reads. The coverage gradient is the external FASTQ subsample.

QUAST is `/hpcfs/fhome/caizhh/.local/bin/quast.py` with `-R` set to the dataset reference, `-t 8`, and `--min-contig 200`. QUAST v5.2.0 is the version pinned by `F26_COMPETITOR_BENCHMARK_CONTRACT.md`. `summary.json` does not store a separate QUAST version string. Contig counts in this report are `metrics.contigs_ge_200bp`. The fastBact writer cutoff is 500 bp, so contigs shorter than 500 bp are absent from `contigs.fa`.

Extensive misassemblies are `metrics.misassemblies`. Local misassemblies are `metrics.local_misassemblies`. Indels are `metrics.indels_per_100kbp`. Peak RSS is `peak_rss_mb` (GNU time maximum resident set size, in MiB). GiB figures in the prose are that value divided by 1024.

### 2.3 Slurm provenance for the 28 fastBact runs

| Job | Node | Dataset | Replicate |
|---:|---|---|---:|
| 3074506 | vnode1 | CS-SA01 | 1 |
| 3074507 | vnode2 | CS-BS01 | 1 |
| 3074508 | vnode1 | CS-PA01 | 1 |
| 3074509 | vnode2 | CS-SE01 | 1 |
| 3076128 | vnode1 | CS-SA01_30x | 1 |
| 3076129 | vnode2 | CS-SA01_50x | 1 |
| 3076130 | vnode2 | CS-SA01_100x | 1 |
| 3076131 | vnode1 | CS-SA01_150x | 1 |
| 3076132 | vnode2 | CS-BS01_30x | 1 |
| 3076133 | vnode1 | CS-BS01_50x | 1 |
| 3076134 | vnode2 | CS-BS01_100x | 1 |
| 3076135 | vnode1 | CS-BS01_150x | 1 |
| 3076136 | vnode2 | CS-PA01_30x | 1 |
| 3076137 | vnode2 | CS-PA01_50x | 1 |
| 3076138 | vnode1 | CS-PA01_100x | 1 |
| 3076139 | vnode2 | CS-PA01_150x | 1 |
| 3076140 | vnode1 | CS-SE01_30x | 1 |
| 3076141 | vnode2 | CS-SE01_50x | 1 |
| 3076142 | vnode1 | CS-SE01_100x | 1 |
| 3076143 | vnode2 | CS-SE01_150x | 1 |
| 3076145 | vnode1 | CS-SA01 | 2 |
| 3076146 | vnode2 | CS-BS01 | 2 |
| 3076147 | vnode1 | CS-PA01 | 2 |
| 3076148 | vnode2 | CS-SE01 | 2 |
| 3076150 | vnode2 | CS-SA01 | 3 |
| 3076151 | vnode1 | CS-BS01 | 3 |
| 3076152 | vnode2 | CS-PA01 | 3 |
| 3076153 | vnode1 | CS-SE01 | 3 |

Jobs `3074506`–`3074509` are the Stage D fastBact replicate 1 runs. Jobs `3076128`–`3076143` are Stage E. Jobs `3076145`–`3076153` are the Stage G replicates. Job ids `3076144` and `3076149` are not in this fastBact map.

---

## 3. Quantitative Results Matrix

Primary measurements are the JSON literals in `summary.json` for fastBact replicate 1. Replicates 2 and 3 have the same QUAST metrics; their wall clock and RSS are in §5. Two wall-clock values are stored with binary-float expansions: CS-BS01 `153.07999999999998` and CS-SE01 `106.27000000000001`. The same class of expansion appears on CS-PA01_150× (`163.07999999999998`) and CS-SA01_150× (`75.28999999999999`).

Read pairs come from job `3075001`, not from QUAST. The full-depth CS-SE01 pair count is the scanned FASTQ count, 3,353,233.

| Dataset | Organism | GC (%) | Target depth | Actual read pairs | Wall time (s) | Peak RSS (MiB) | Genome fraction (%) | Contigs (≥200 bp) | N50 (bp) | Largest contig (bp) | Misassemblies | Local misassemblies | Indels / 100 kbp |
|---|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| CS-BS01_30x | *B. subtilis* | 43.51 | 30× | 418769 | 36.46 | 1096.85546875 | 98.272 | 119 | 66973 | 204588 | 0 | 0 | 0.51 |
| CS-BS01_50x | *B. subtilis* | 43.51 | 50× | 697948 | 56.81 | 1033.9765625 | 98.346 | 170 | 39818 | 234253 | 0 | 0 | 0.84 |
| CS-BS01_100x | *B. subtilis* | 43.51 | 100× | 1395896 | 88.59 | 1152.92578125 | 98.549 | 172 | 47159 | 128251 | 0 | 0 | 0.6 |
| CS-BS01_150x | *B. subtilis* | 43.51 | 150× | 2093844 | 167.36 | 1554.4453125 | 98.575 | 111 | 91207 | 214696 | 0 | 0 | 0.55 |
| CS-BS01 | *B. subtilis* | 43.51 | full | 2253807 | 153.07999999999998 | 1425.3828125 | 98.576 | 110 | 78260 | 354932 | 0 | 0 | 0.82 |
| CS-PA01_30x | *P. aeruginosa* | 66.56 | 30× | 374367 | 57.44 | 1104.46484375 | 99.036 | 217 | 53277 | 223412 | 4 | 1 | 0.77 |
| CS-PA01_50x | *P. aeruginosa* | 66.56 | 50× | 623945 | 76.65 | 1081.39453125 | 99.078 | 191 | 59072 | 308504 | 3 | 1 | 0.58 |
| CS-PA01_100x | *P. aeruginosa* | 66.56 | 100× | 1247889 | 165.96 | 1454.51171875 | 98.861 | 424 | 30828 | 201253 | 2 | 0 | 0.67 |
| CS-PA01_150x | *P. aeruginosa* | 66.56 | 150×, clamped to library | 1518403 | 163.07999999999998 | 1533.12109375 | 98.74 | 576 | 23514 | 90388 | 2 | 0 | 0.63 |
| CS-PA01 | *P. aeruginosa* | 66.56 | full | 1518403 | 180.57 | 1527.0703125 | 98.74 | 576 | 23514 | 90388 | 2 | 0 | 0.63 |
| CS-SA01_30x | *S. aureus* | 32.83 | 30× | 280178 | 37.91 | 815.08984375 | 64.582 | 1713 | 1266 | 105012 | 0 | 0 | 0.05 |
| CS-SA01_50x | *S. aureus* | 32.83 | 50× | 466964 | 38.39 | 660.1015625 | 95.199 | 859 | 4927 | 63184 | 0 | 0 | 0.0 |
| CS-SA01_100x | *S. aureus* | 32.83 | 100× | 933928 | 46.87 | 710.5390625 | 98.194 | 217 | 32244 | 96025 | 0 | 0 | 0.25 |
| CS-SA01_150x | *S. aureus* | 32.83 | 150× | 1400892 | 75.28999999999999 | 789.5859375 | 98.445 | 253 | 28035 | 87734 | 0 | 0 | 0.0 |
| CS-SA01 | *S. aureus* | 32.83 | full | 2239653 | 148.55 | 1196.7421875 | 98.683 | 239 | 25381 | 77704 | 2 | 0 | 0.25 |
| CS-SE01_30x | *S. enterica* | 52.22 | 30× | 485745 | 35.5 | 1381.6640625 | 98.11 | 225 | 46607 | 132906 | 1 | 1 | 0.04 |
| CS-SE01_50x | *S. enterica* | 52.22 | 50× | 809575 | 38.76 | 1428.65625 | 98.401 | 322 | 27151 | 91340 | 1 | 1 | 0.15 |
| CS-SE01_100x | *S. enterica* | 52.22 | 100× | 1619150 | 69.25 | 3787.87109375 | 98.647 | 220 | 39823 | 208121 | 1 | 1 | 0.23 |
| CS-SE01_150x | *S. enterica* | 52.22 | 150× | 2428725 | 76.68 | 4072.14453125 | 98.669 | 155 | 58899 | 199626 | 1 | 1 | 0.1 |
| CS-SE01 | *S. enterica* | 52.22 | full | 3353233 | 106.27000000000001 | 7854.46875 | 98.61 | 198 | 47464 | 178426 | 1 | 1 | 0.06 |

### 3.1 Comparison with the full-depth fastBact baseline

The baseline is fastBact replicate 1 on the unsampled dataset. Differences of genome fraction are the stored percentages subtracted and shown to three decimal places. N50 and misassembly differences are exact integer subtraction. CS-PA01_150× matches the baseline on every QUAST field in the table because the 150× target was clamped to the full library; the contig SHA-256 is the same string (`b15ebd2289fb7f7ca39a49698f9049aa281d83a00982199f69da0722a068a9ff`). Wall clock and peak RSS of that row still differ, because it is a separate execution.

| Dataset | Genome fraction (%) | Δ vs full (pp) | Misassemblies | Δ misassemblies | Local misassemblies | N50 (bp) | Δ N50 (bp) |
|---|---:|---:|---:|---:|---:|---:|---:|
| CS-BS01_30x | 98.272 | −0.304 | 0 | 0 | 0 | 66973 | −11287 |
| CS-BS01_50x | 98.346 | −0.230 | 0 | 0 | 0 | 39818 | −38442 |
| CS-BS01_100x | 98.549 | −0.027 | 0 | 0 | 0 | 47159 | −31101 |
| CS-BS01_150x | 98.575 | −0.001 | 0 | 0 | 0 | 91207 | +12947 |
| CS-PA01_30x | 99.036 | +0.296 | 4 | +2 | 1 | 53277 | +29763 |
| CS-PA01_50x | 99.078 | +0.338 | 3 | +1 | 1 | 59072 | +35558 |
| CS-PA01_100x | 98.861 | +0.121 | 2 | 0 | 0 | 30828 | +7314 |
| CS-PA01_150x | 98.74 | 0.000 | 2 | 0 | 0 | 23514 | 0 |
| CS-SA01_30x | 64.582 | −34.101 | 0 | −2 | 0 | 1266 | −24115 |
| CS-SA01_50x | 95.199 | −3.484 | 0 | −2 | 0 | 4927 | −20454 |
| CS-SA01_100x | 98.194 | −0.489 | 0 | −2 | 0 | 32244 | +6863 |
| CS-SA01_150x | 98.445 | −0.238 | 0 | −2 | 0 | 28035 | +2654 |
| CS-SE01_30x | 98.11 | −0.500 | 1 | 0 | 1 | 46607 | −857 |
| CS-SE01_50x | 98.401 | −0.209 | 1 | 0 | 1 | 27151 | −20313 |
| CS-SE01_100x | 98.647 | +0.037 | 1 | 0 | 1 | 39823 | −7641 |
| CS-SE01_150x | 98.669 | +0.059 | 1 | 0 | 1 | 58899 | +11435 |

N50 is not monotonic in target depth. On *B. subtilis*, 50× N50 (39,818 bp) is below 30× (66,973 bp), and 150× N50 (91,207 bp) is above the full-depth N50 (78,260 bp). On *P. aeruginosa*, 30× and 50× have higher N50 than 100× and the full library, together with more misassemblies. Genome fraction above the full-depth value at low *P. aeruginosa* depth is a QUAST observation. It is not evidence that the low-depth assembly is structurally cleaner: extensive misassemblies are higher there (4 and 3 versus 2).

### 3.2 Depth-run contig SHA-256

These are replicate-1 `assembly_sha256` values. Only CS-PA01_150× matches its full-depth hash.

| Dataset | Contig FASTA SHA-256 |
|---|---|
| CS-BS01_30x | `8862338aa29dddc6293e5ba511edd475bdda93a38e85b8c78ca39d26d228366c` |
| CS-BS01_50x | `c702745d2e30220d1e20cd53a1601c668fe4e616c08a6f9fbd8817353aefabcc` |
| CS-BS01_100x | `8068a60e694b405edab0f979e6f0d942d4b8a3b0fffb772abe2c8e9567319d67` |
| CS-BS01_150x | `d8323d7b4105c0660b2d55f50aa3fb9e42c8c7c1d9fdba3b5e288f5c3df090b8` |
| CS-PA01_30x | `4ef99ab10098a4ac67113a320e8f6db953b107b4653e8ff1fb59f5c39c0b346b` |
| CS-PA01_50x | `55322f134cb4b4d5ba8e2035bf8374a5c9dab8cab1872de47cb82aa64597f94f` |
| CS-PA01_100x | `5859f0390096d1fdaf394dc23fd018800b276b45a0a12abff948ced58861f93d` |
| CS-PA01_150x | `b15ebd2289fb7f7ca39a49698f9049aa281d83a00982199f69da0722a068a9ff` |
| CS-SA01_30x | `ec61418b0a171c2afb29f3ef1309e81bdb4897e5e432c6f41d1b6a3a013bbbdd` |
| CS-SA01_50x | `36b2d5677b32a6a92f5e080733582279c569eff743c76d7a04a29ad31f6168c0` |
| CS-SA01_100x | `743103b0d36eb39090a0a50016d5629b2989b84cb470c28fc259d7aee407684c` |
| CS-SA01_150x | `b848a59cc35e575d096eaa89d24e0e8a18259af66c22db48f01647c21c0be81d` |
| CS-SE01_30x | `45a23fe57c1a594cd649690ab8501669ebcdbda974faa9d200a6bed9af0d77b4` |
| CS-SE01_50x | `3d2d668ec48fde6a7044a0f452ed656d582b9b4c5be2ec9d87f193deba9b4418` |
| CS-SE01_100x | `adcdd6d8076aa54a009645749d2a267199be3d877f24dfcf75eae2b689c51f30` |
| CS-SE01_150x | `c162c33bb1adb193b8dee7c4b993c0d1de86a342298e9bb75c7f4eff73dfdee3` |

---

## 4. Organism Trajectory and GC-Stress Analysis

GC classes below use the manifest percentages. The trajectories are the five fastBact replicate-1 rows in §3.

### 4.1 Low GC — *S. aureus* MW2, 32.83%

The 30× input is the panel's coverage threshold. Genome fraction falls to `64.582`, N50 to 1,266 bp, NGA50 to 777 bp, and the ≥200 bp contig count rises to 1,713. Indels are `0.05` per 100 kbp and mismatches are `0.33` per 100 kbp. Extensive and local misassemblies are both 0. The run still exits 0. This is an incomplete recovery of the reference, not an assembler crash.

Recovery at the next planned depth is sharp. At 50×, genome fraction is `95.199`, contigs fall to 859, and N50 rises to 4,927 bp, still with 0 extensive and 0 local misassemblies and `0.0` indels per 100 kbp. At 100×, genome fraction is `98.194`, N50 is 32,244 bp, and contigs are 217. At 150×, genome fraction is `98.445`, N50 is 28,035 bp, contigs are 253, and indels return to `0.0`. Misassemblies stay 0 across all four gradient depths.

The full library (2,239,653 pairs; formula depth about 239.8× at 151 bp) reaches genome fraction `98.683` and then records **2** extensive misassemblies, with local misassemblies still 0, indels `0.25`, and mismatches `2.11`. Those two events are absent at 150× (1,400,892 pairs). Stage D already recorded the full-depth count of 2 against SPAdes 0, Unicycler 0, and SKESA 1. The gradient shows that this count is depth-associated on this dataset: it is not present at 30–150× and it is present, and replicate-stable, at full depth. This report does not assign the events to a resolver stage.

Peak RSS on the four *S. aureus* depth runs stays between `660.1015625` and `815.08984375` MiB. Full-depth replicates stay between `1196.7421875` and `1250.46484375` MiB.

### 4.2 High GC — *P. aeruginosa* PAO1, 66.56%

The library cannot supply a true 150× subsample under the 251 bp first-read formula. Achieved formula depth of the clamped row, and of the full library, is **121.678×**. Comparisons that treat CS-PA01_150× as deeper than CS-PA01 are comparisons of two executions of the same pair count. Their contig hashes match. Their QUAST metrics match. Their wall times do not (`163.07999999999998` s versus `180.57` s).

On the depths that are real subsamples, extensive misassemblies fall as depth rises: **4** at 30×, **3** at 50×, **2** at 100×. Local misassemblies are **1** at 30× and 50× and **0** at 100×. The clamped 150× row and the full library remain at 2 and 0. Genome fraction moves the other way: `99.036` at 30× and `99.078` at 50×, then `98.861` at 100× and `98.74` at the full library. N50 is 53,277 bp at 30× and 59,072 bp at 50×, then 30,828 bp at 100× and 23,514 bp at full depth. Higher N50 at low depth coincides with more misassemblies. The matrix does not support a claim that added depth only improves every QUAST column.

Memory on the four depth rows is `1104.46484375`, `1081.39453125`, `1454.51171875`, and `1533.12109375` MiB. The full-depth replicates are `1527.0703125`, `1537.9296875`, and `1496.4765625` MiB. The maximum *P. aeruginosa* peak RSS in the fastBact set is therefore `1537.9296875` MiB, which is 1.502 GiB.

### 4.3 Intermediate GC — *B. subtilis* 43.51% and *S. enterica* 52.22%

*B. subtilis* keeps 0 extensive and 0 local misassemblies at 30×, 50×, 100×, 150×, and full depth. Genome fraction is `98.272`, `98.346`, `98.549`, `98.575`, and `98.576`. The minimum is above 98.2%. Contig counts stay between 110 and 172. N50 does not rise smoothly: 66,973, 39,818, 47,159, 91,207, and 78,260 bp. The 150× assembly is more contiguous by N50 than the full library (91,207 versus 78,260) and slightly less complete by genome fraction (`98.575` versus `98.576`). Largest contig is highest at full depth (354,932 bp). Peak RSS stays near 1.0–1.5 GiB except CS-BS01_150× at `1554.4453125` MiB (1.518 GiB). Full-depth replicate 3 reaches `1534.1171875` MiB and `183.21` s.

*S. enterica* is stable on the structural counts that this matrix measures. Every depth, including full, has exactly **1** extensive misassembly and **1** local misassembly. Genome fraction stays high: `98.11`, `98.401`, `98.647`, `98.669`, and `98.61`. The 100× and 150× rows are 0.037 and 0.059 percentage points above the full-depth genome fraction. N50 is again non-monotonic: 46,607, 27,151, 39,823, 58,899, and 47,464 bp. Indels per 100 kbp are `0.04`, `0.15`, `0.23`, `0.1`, and `0.06`.

Memory on *S. enterica* is the panel's scaling outlier. Depth-run peak RSS is `1381.6640625` MiB at 30×, `1428.65625` MiB at 50×, `3787.87109375` MiB at 100× (3.699 GiB), and `4072.14453125` MiB at 150× (3.977 GiB). Full depth then reaches `7858.8125` MiB (7.675 GiB) on replicate 3. Wall clock does not follow that jump: the five replicate-1 times are `35.5`, `38.76`, `69.25`, `76.68`, and `106.27000000000001` s. Stage D already recorded the full-depth RSS outlier and declined to assign it to a stage. The gradient localizes the step: RSS is about 1.4 GiB at 30–50× and rises above 3.6 GiB from 100× upward, while misassembly counts do not change.

---

## 5. Determinism Certification (Stage G)

The harness groups runs by dataset and tool. A group is checked only when it has at least three replicates. It passes only when every replicate is `PASS`, every `assembly_sha256` is the same non-empty value, and every field in `METRIC_FIELDS` compares equal. `variance_max` is the maximum population variance of those numeric fields, `sum((value - mean)^2) / n`, with `n = 3`. Wall clock, CPU time, and peak RSS are not in that set.

Depth rows and competitor rows have one replicate. Their determinism entries stay `checked: false` with reason `not enough replicates: 1 < 3`. This section certifies the four full-depth fastBact groups only.

### 5.1 Contig FASTA SHA-256

| Dataset | Replicate 1 | Replicate 2 | Replicate 3 | Identical |
|---|---|---|---|---|
| CS-BS01 | `2bb0ecdd10f09b2320f4ba557b15a79bc5069f02f89352eceb9ec053b43a8861` | `2bb0ecdd10f09b2320f4ba557b15a79bc5069f02f89352eceb9ec053b43a8861` | `2bb0ecdd10f09b2320f4ba557b15a79bc5069f02f89352eceb9ec053b43a8861` | yes |
| CS-PA01 | `b15ebd2289fb7f7ca39a49698f9049aa281d83a00982199f69da0722a068a9ff` | `b15ebd2289fb7f7ca39a49698f9049aa281d83a00982199f69da0722a068a9ff` | `b15ebd2289fb7f7ca39a49698f9049aa281d83a00982199f69da0722a068a9ff` | yes |
| CS-SA01 | `86081b422341266a5a32feef424080bc0595f8ce599d996817b3e7d6c26a1962` | `86081b422341266a5a32feef424080bc0595f8ce599d996817b3e7d6c26a1962` | `86081b422341266a5a32feef424080bc0595f8ce599d996817b3e7d6c26a1962` | yes |
| CS-SE01 | `4b72cfaa1b4fa8406be3f9d29c790047a8cae2ce2cb996cea9c9f31859314bfb` | `4b72cfaa1b4fa8406be3f9d29c790047a8cae2ce2cb996cea9c9f31859314bfb` | `4b72cfaa1b4fa8406be3f9d29c790047a8cae2ce2cb996cea9c9f31859314bfb` | yes |

The same four strings are the `determinism[].assembly_sha256` values. Each group has `mismatches: []`.

### 5.2 Metric variance

| Dataset | Replicates | checked | pass | variance_max | Mismatched metric fields |
|---|---:|---|---|---:|---|
| CS-BS01 | 3 | true | true | `2.0194839173657902e-28` | none |
| CS-PA01 | 3 | true | true | `0.0` | none |
| CS-SA01 | 3 | true | true | `0.0` | none |
| CS-SE01 | 3 | true | true | `0.0` | none |

For CS-BS01, all three `genome_fraction_pct` values are the JSON number `98.576`, and every other QUAST field is likewise equal across the three replicates, including misassemblies 0, local misassemblies 0, N50 78,260, and indels `0.82`. Recomputing the harness formula on three copies of `98.576` reproduces `variance_max` exactly, because the mean `(98.576 + 98.576 + 98.576) / 3` is `98.57599999999998`. This is not a contig difference and not a QUAST-table difference. The other three species store `variance_max` as `0.0`.

Contig identity is bit-for-bit on all four species. The statement `variance_max = 0.0` is true for CS-PA01, CS-SA01, and CS-SE01. It is not the literal stored value for CS-BS01.

### 5.3 Resource measurements are not deterministic

QUAST metrics match. Runtime and RSS do not. Reporting σ² = 0 for the assembly does not extend to these counters.

| Dataset | Rep | Wall time (s) | Peak RSS (MiB) | CPU time (s) |
|---|---:|---:|---:|---:|
| CS-BS01 | 1 | 153.07999999999998 | 1425.3828125 | 411.8 |
| CS-BS01 | 2 | 153.64 | 1418.22265625 | 414.63 |
| CS-BS01 | 3 | 183.21 | 1534.1171875 | 543.48 |
| CS-PA01 | 1 | 180.57 | 1527.0703125 | 625.06 |
| CS-PA01 | 2 | 185.91 | 1537.9296875 | 636.63 |
| CS-PA01 | 3 | 162.67000000000002 | 1496.4765625 | 535.07 |
| CS-SA01 | 1 | 148.55 | 1196.7421875 | 413.07 |
| CS-SA01 | 2 | 138.31 | 1250.46484375 | 372.68 |
| CS-SA01 | 3 | 123.31 | 1248.37890625 | 313.78 |
| CS-SE01 | 1 | 106.27000000000001 | 7854.46875 | 267.87 |
| CS-SE01 | 2 | 105.50999999999999 | 7777.0703125 | 265.26 |
| CS-SE01 | 3 | 114.41 | 7858.8125 | 294.8 |

---

## 6. Robustness Hardening and Invariant Assessment (Stage F)

F26-P4 allows an algorithmic repair only when all five conditions hold:

1. A reproducible failure case is identified.
2. The responsible assembly stage is isolated.
3. A minimal failing unit or regression test is formulated.
4. Conservative evidence isolation (`p > b + H`) and loop re-entrancy stay intact.
5. An independent audit precedes merge.

### 6.1 What the depth matrix shows about the running code

The 16 depth runs and the 12 full-depth fastBact replicates all completed with assembler exit 0 and QUAST exit 0. The other 12 Stage D competitor runs in the same file also have exit 0. There is no crash, no non-zero assembler status, and no QUAST failure to use as condition 1.

`git diff 02872f81811402d2bd603ff7c0a2ecf8b3dbe528 -- src/` is empty at the collation HEAD. The 28 fastBact executions share one binary digest, `517402653fcefc625cb79912ade7e7e382b0d05343bbeb483d7c82dd807947d1`. The accepted F25 implementation of `p > b + H`, and the loop re-entrancy safeguards in that same source, are the code that ran. The QUAST tables are not a dynamic trace of the inequality. They also show no source edit introduced to obtain these depth results. Loop re-entrancy is assessed here as an empirical completion record: every submitted fastBact job returned, and the longest fastBact wall time in the file is `185.91` s, inside the 3-hour Slurm limit. That is not a proof that the invariant holds for every depth.

### 6.2 Observations that do not authorize a patch

- ***S. aureus* 30× genome fraction `64.582`.** The assembler finished. The same binary recovers `95.199` at 50× and `98.445` at 150× with 0 misassemblies. The low genome fraction tracks the planned input depth. No stage was isolated, and no minimal regression test was written. Changing k-mer schedules, `--depth`, or resolver thresholds to lift 30× recovery would be speculative.
- ***S. aureus* full-depth extensive misassemblies = 2.** The count is reproducible across three replicates, and it is absent at 30–150×. Reproducibility satisfies part of condition 1. Conditions 2 and 3 are not met: the matrix does not name the stage that creates the two events, and it does not supply a minimal test. Stage D already left this count as an open observation. The gradient narrows the input range. It still does not justify a code change.
- ***P. aeruginosa* 4 → 3 → 2 extensive misassemblies.** The count declines as the real subsample depth rises, and it matches the full-library count of 2 once the input is the whole library. Local misassemblies fall from 1 to 0 by 100×. This is not a new defect relative to the Stage D baseline.
- ***S. enterica* RSS step above 100×.** Misassemblies stay 1 + 1 while RSS moves from about 1.4 GiB to 3.977 GiB at 150× and 7.675 GiB at full depth. The observation is reproducible on the full-depth replicates (all three are above 7.5 GiB). The responsible stage is not isolated. Raising or lowering a memory threshold without that isolation would be speculative.
- **Non-monotonic N50.** Several 150× or low-depth N50 values beat the full-depth N50. Contiguity is not a crash. Retuning `--min-len`, `--kmers`, or `--depth 120` to chase N50 remains outside F26-P4.

**Disposition:** speculative code modifications are rejected. Stage F stays closed. The existing binary is the one that produced 28/28 fastBact passes under this depth and replicate stress.

---

## 7. Conclusion and Handoff to Stage H

**Disposition:** `STAGE_E_GRADIENT_REPORTED` and `STAGE_G_FASTBACT_DETERMINISM_REPORTED`.

The frozen `02872f8` binary completed the 16-point cross-species depth grid and three full-depth replicates on each of the four taxa. Structural behavior splits by organism. *B. subtilis* holds 0 misassemblies and genome fraction above 98.2% from 30× through full depth. *S. enterica* holds exactly one extensive and one local misassembly across the same grid. *S. aureus* drops to `64.582` genome fraction at 30×, recovers to `95.199` at 50× and `98.445` at 150× with 0 misassemblies on the gradient, and carries 2 extensive misassemblies at full depth. *P. aeruginosa* extensive misassemblies fall from 4 at 30× to 2 at 100×, with the nominal 150× row clamped to the full library at formula depth 121.678×. Stage E wall times are `35.5`–`167.36` s. Stage E peak RSS peaks at `4072.14453125` MiB. Full-depth *S. enterica* peak RSS reaches `7858.8125` MiB. Contig SHA-256 matches across N=3 for every species. Stored metric `variance_max` is `0.0` on three species and `2.0194839173657902e-28` on *B. subtilis* from identical `98.576` values.

Stage H is the independent scientific acceptance audit. This document does not perform it and does not request a merge. Items that the audit should check against `summary.json`, rather than against a rounded summary, are:

1. CS-PA01_150× clamp: `pairs_out` 1,518,403, formula depth 121.678×, contig SHA-256 equal to full-depth CS-PA01.
2. CS-BS01 `variance_max` `2.0194839173657902e-28` with equal contig hashes and equal QUAST fields.
3. CS-SA01 extensive misassemblies: 0 on the four gradient rows, 2 on all three full-depth replicates.
4. CS-SE01 peak RSS: about 1.4 GiB at 30–50×, 3.699 GiB at 100×, 3.977 GiB at 150×, and 7.595–7.675 GiB on the three full-depth replicates, with misassemblies unchanged at 1 + 1.
5. Empty `src/` diff against `02872f8`, and one fastBact binary digest on all 28 runs.
6. Absence of the 16 *E. coli* depth rows from this 40-run file. F25 oracle non-regression remains unreplayed.

Stage H's auditor, under the F26 contract, is Codex Sol. Acceptance, merge, push, and release stay with the user as `FINAL_AUTHORIZATION_AUTHORITY`.
