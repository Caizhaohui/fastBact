The analysis is complete, but I could not save `assembly_comparison.md` because this session enforces read-only filesystem access. The complete Markdown report follows.

# fastBact Phase F20: Assembly Comparison

## 1. Executive Summary and Baseline Provenance

**Overall verdict: `REGRESSION_DETECTED`. C0-02R remains strictly NOT ACCEPTED.**

Candidate B shows dataset-dependent results. R1-001 regresses in six of seven acceptance metrics: N50, NGA50, genome fraction, misassemblies, mismatches, and indels. R1-003 improves N50, NGA50, misassemblies, mismatches, and indels, preserves local misassemblies, but regresses in genome fraction.

Across both datasets, `acceptance.tsv` records **7 PASS and 7 REGRESSION results**. Improvements on R1-003 do not offset failures on R1-001 or the completeness regressions on both datasets.

### Frozen identities

| Version | Role | Commit identity |
|---|---|---|
| A | Frozen baseline | `b9727c3ab372f7b995e439d92b72c7a4bfd051cf` |
| B | Candidate benchmark HEAD | `0606be423d187771f5925a83da2cc377c40c2a3d` |
| B | Candidate code commit | `a9930471a73612159a83694d3d661f060042efc8` |

- **Baseline A binary SHA-256:** `2b184db0f108cb508ca7e591240eddb551bea8ebad90ac5eb8ee360bac78ea38`
- **Candidate B binary SHA-256:** `4fbafeeeffe94abeef42efaa5710f75dbd0af8b234edd2bffba43cc4f50198f0`
- **Candidate B production `src/` tree:** `30967bf4940844aeb6b744732bcad6d2a72daab4`

These identities are recorded in [binary_checksums.tsv](binary_checksums.tsv). F20 baseline A corresponds to version B in the historical three-version benchmark; historical labels must not be interpreted as F20 labels.

### Evidence and execution

The copies of `assembly_metrics.tsv` and `acceptance.tsv` under the F20 evidence directory and `benchmark-results/C0_02R_f20_selective_20261022/` were compared and are byte-identical.

Slurm job **3001561** ran on **vnode1**, partition **qfree**, with total elapsed time **00:13:13 (793 seconds)**. All eight numbered sub-steps completed with exit code `0:0`; none failed, were cancelled, or encountered OOM. Successful execution establishes benchmark completion, not assembly acceptance.

The accounting record is authoritative for completion; the submission-era ledger still contains `RUNNING`.

Each dataset/version combination has **one replicate**. All assemblies were evaluated with a recorded minimum contig threshold of **200 bp**. Comparisons are descriptive; performance significance and run-to-run variability cannot be estimated.

## 2. Quantitative Comparison

For each metric:

\[
\Delta=B-A,\qquad \%\Delta=100(B-A)/A
\]

Genome-fraction absolute differences are expressed in **percentage points (pp)**; their relative changes use the formula above. Relative changes are undefined when A is zero and are shown as `N/A`.

Higher N50, NGA50, and genome fraction are favorable. Lower error counts, error rates, wall time, and peak RSS are favorable.

For graph telemetry, classifications use an explicit **simplification proxy**: more valid joins and fewer residual active unitigs are favorable. These classifications describe graph simplification, not independently established assembly-quality improvements or acceptance gates.

### R1-001

| Metric | A | B | Δ (B − A) | Relative change | Classification |
|---|---:|---:|---:|---:|---|
| N50 (bp) | 58,318 | 57,462 | −856 | −1.468% | REGRESSED |
| NGA50 (bp) | 58,318 | 57,127 | −1,191 | −2.042% | REGRESSED |
| Genome fraction (%) | 98.364 | 98.281 | −0.083 pp | −0.084% | REGRESSED |
| Misassemblies | 4 | 5 | +1 | +25.000% | REGRESSED |
| Local misassemblies | 0 | 0 | 0 | N/A | UNCHANGED |
| Mismatches/100 kbp | 2.14 | 2.38 | +0.24 | +11.215% | REGRESSED |
| Indels/100 kbp | 0.48 | 1.03 | +0.55 | +114.583% | REGRESSED |
| Assembly wall time (s) | 150.30 | 148.87 | −1.43 | −0.951% | IMPROVED |
| Assembly peak RSS (KiB) | 1,543,760 | 1,479,280 | −64,480 | −4.177% | IMPROVED |
| Valid k joins | 987 | 947 | −40 | −4.053% | REGRESSED |
| Valid k−1 joins | 18 | 10 | −8 | −44.444% | REGRESSED |
| Rejected ambiguous junctions | 0 | 0 | 0 | N/A | UNCHANGED |
| Rejected invalid junctions | 0 | 0 | 0 | N/A | UNCHANGED |
| Active unitigs after exSPAnder | 315 | 1,728 | +1,413 | +448.571% | REGRESSED |

**Acceptance:** 1 PASS, 6 REGRESSION. The sole passing gate is unchanged local misassemblies.

### R1-003

| Metric | A | B | Δ (B − A) | Relative change | Classification |
|---|---:|---:|---:|---:|---|
| N50 (bp) | 67,452 | 70,868 | +3,416 | +5.064% | IMPROVED |
| NGA50 (bp) | 63,648 | 70,010 | +6,362 | +9.996% | IMPROVED |
| Genome fraction (%) | 98.330 | 98.233 | −0.097 pp | −0.099% | REGRESSED |
| Misassemblies | 2 | 1 | −1 | −50.000% | IMPROVED |
| Local misassemblies | 1 | 1 | 0 | 0.000% | UNCHANGED |
| Mismatches/100 kbp | 1.51 | 1.08 | −0.43 | −28.477% | IMPROVED |
| Indels/100 kbp | 0.04 | 0.00 | −0.04 | −100.000% | IMPROVED |
| Assembly wall time (s) | 219.60 | 221.48 | +1.88 | +0.856% | REGRESSED |
| Assembly peak RSS (KiB) | 1,861,880 | 1,757,880 | −104,000 | −5.586% | IMPROVED |
| Valid k joins | 564 | 497 | −67 | −11.879% | REGRESSED |
| Valid k−1 joins | 7 | 2 | −5 | −71.429% | REGRESSED |
| Rejected ambiguous junctions | 0 | 0 | 0 | N/A | UNCHANGED |
| Rejected invalid junctions | 0 | 0 | 0 | N/A | UNCHANGED |
| Active unitigs after exSPAnder | 300 | 704 | +404 | +134.667% | REGRESSED |

**Acceptance:** 6 PASS, 1 REGRESSION. Genome fraction is the failing gate.

## 3. Multi-Domain Metric Analysis

### Structural correctness

R1-001 misassemblies increase from **4 to 5 (+25%)**, while local misassemblies remain zero. R1-003 misassemblies decrease from **2 to 1 (−50%)**, while local misassemblies remain one.

The summed misassembly count is six for both versions, but that aggregate conceals the R1-001 regression. Acceptance must remain dataset-specific.

### Base-level accuracy

R1-001 mismatches increase **11.215%**, and indels increase **114.583%**. The candidate therefore more than doubles the reported indel rate on this dataset.

R1-003 mismatches decrease **28.477%**, and the reported indel rate decreases from **0.04 to 0.00 per 100 kbp**. The latter is a 100% reduction in the reported metric; it should not be generalized into a claim that the entire assembly is error-free.

These normalized rates should not be averaged across datasets without their underlying aligned-base denominators.

### Assembly completeness and contiguity

R1-001 loses both contiguity and aligned contiguity: N50 falls **856 bp**, and NGA50 falls **1,191 bp**. Genome fraction also decreases **0.083 pp**.

R1-003 gains **3,416 bp in N50** and **6,362 bp in NGA50**, accompanied by fewer structural and base-level errors. Nevertheless, genome fraction decreases **0.097 pp**. Better contiguity therefore does not establish improved completeness.

Both genome-fraction decreases are small in absolute terms but are explicitly classified as regressions by the supplied acceptance results.

### Graph simplification and junction telemetry

Candidate B records fewer valid k and k−1 joins on both datasets:

- R1-001: valid k joins decrease **40**, and valid k−1 joins decrease **8**.
- R1-003: valid k joins decrease **67**, and valid k−1 joins decrease **5**.

Residual active unitigs increase from **315 to 1,728** on R1-001, approximately **5.49×**, and from **300 to 704** on R1-003, approximately **2.35×**.

Together, these measurements indicate reduced simplification under the stated proxy. They do not prove a causal relationship with the final quality metrics: R1-003 improves several assembly outcomes despite the same telemetry direction.

Both rejection counters remain zero throughout. This establishes that no rejections were recorded in those categories; it does not establish that all joins are biologically correct.

### Computing cost

R1-001 is **1.43 seconds faster (−0.951%)**, while R1-003 is **1.88 seconds slower (+0.856%)**. Summed assembly wall time changes from **369.90 to 370.35 seconds**, an increase of **0.45 seconds (+0.122%)**.

Peak RSS decreases on both datasets:

- R1-001: **64,480 KiB lower (−4.177%)**.
- R1-003: **104,000 KiB lower (−5.586%)**.

These are single-run observations. They support a measured memory reduction, but no statistical performance claim.

The report uses `assembly_peak_rss_kib` from the requested metrics TSV. Slurm `MaxRSS` values differ and are not substituted. Total job elapsed time includes assembly, evaluation, and orchestration, so it is not equivalent to summed assembly wall time.

## 4. Historical Comparison with Phase D

The historical candidate is commit **`f241d9d3c9d21470f63bd85b5e1bd6eb0bfc8d30`**, labeled **C** in the earlier three-version results.

Historical values below come from:

`benchmark-results/C0_02R_pe_fix_20261019/results/C0_02R_ASSEMBLY_AB.tsv`

### R1-001: historical candidate versus F20 B

| Metric | Phase D old candidate | F20 B | Δ | Relative change | Classification |
|---|---:|---:|---:|---:|---|
| N50 (bp) | 58,318 | 57,462 | −856 | −1.468% | REGRESSED |
| NGA50 (bp) | 58,318 | 57,127 | −1,191 | −2.042% | REGRESSED |
| Genome fraction (%) | 98.375 | 98.281 | −0.094 pp | −0.096% | REGRESSED |
| Misassemblies | 4 | 5 | +1 | +25.000% | REGRESSED |
| Local misassemblies | 0 | 0 | 0 | N/A | UNCHANGED |
| Mismatches/100 kbp | 2.86 | 2.38 | −0.48 | −16.783% | IMPROVED |
| Indels/100 kbp | 0.96 | 1.03 | +0.07 | +7.292% | REGRESSED |

R1-001 mismatches improve from **2.86 to 2.38**, but remain above baseline A’s **2.14**, by **0.24 (+11.215%)**. This is a partial recovery, not restoration of baseline accuracy. Historical indels, structural correctness, completeness, and contiguity all worsen.

### R1-003: historical candidate versus F20 B

| Metric | Phase D old candidate | F20 B | Δ | Relative change | Classification |
|---|---:|---:|---:|---:|---|
| N50 (bp) | 67,452 | 70,868 | +3,416 | +5.064% | IMPROVED |
| NGA50 (bp) | 63,648 | 70,010 | +6,362 | +9.996% | IMPROVED |
| Genome fraction (%) | 98.330 | 98.233 | −0.097 pp | −0.099% | REGRESSED |
| Misassemblies | 2 | 1 | −1 | −50.000% | IMPROVED |
| Local misassemblies | 1 | 1 | 0 | 0.000% | UNCHANGED |
| Mismatches/100 kbp | 1.67 | 1.08 | −0.59 | −35.329% | IMPROVED |
| Indels/100 kbp | 0.04 | 0.00 | −0.04 | −100.000% | IMPROVED |

**Historical attribution correction:** the supplied task describes an R1-003 mismatch reduction from **1.51 to 1.08**. That comparison is correct for **baseline A versus F20 B**. The archived row for historical candidate `f241d9d` reports **1.67**, making the direct historical-candidate reduction **1.67 to 1.08**. Both comparisons show improvement, but their reference versions must remain distinct.

R1-003 misassemblies decrease from **2 to 1** against both reference versions.

Historical comparisons provide context across benchmark phases; they do not replace the current F20 A/B acceptance evidence.

## 5. Disposition

**C0-02R remains strictly NOT ACCEPTED.**

The completed benchmark establishes:

- **R1-001:** six acceptance regressions, including structural correctness, base-level accuracy, contiguity, and completeness.
- **R1-003:** five improved quality metrics and unchanged local misassemblies, but a genome-fraction regression.
- **Overall:** seven failed metric comparisons out of fourteen, consistent with the recorded `REGRESSION_DETECTED` status.

Successful Slurm execution, reduced peak memory, and partial historical recovery do not satisfy the failed assembly-quality gates.