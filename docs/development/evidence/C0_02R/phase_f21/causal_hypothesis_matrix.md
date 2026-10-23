# fastBact Phase F21-A: Causal Hypothesis Evaluation Matrix

This matrix formalizes the evaluation of competing causal hypotheses explaining the Phase F20 benchmark results across `R1-001` and `R1-003`.

In strict adherence to Phase F21-A guidelines, all findings are categorized into three evidentiary tiers:
- **`OBSERVED`**: Direct empirical facts established by exact measurement in Phase F20 / F21.
- **`SUPPORTED`**: Mechanistic attributions verified by intermediate graph states, sequence alignments, and boundary accounting.
- **`UNVERIFIED`**: Hypotheses lacking sufficient locus-level evidence.

---

## 1. Summary of Evaluated Hypotheses

| Hypothesis | Description | Evidentiary Status | Primary Finding & Mechanism |
| :--- | :--- | :---: | :--- |
| **H1** | **exSPAnder accepted a new chimeric join on R1-001** | **`OBSERVED (NEW JOIN)` / `SUPPORTED (HYPOTHESIS)`** | **NEW_JOIN_ERROR_CONFIRMED**. `contig_16` in Candidate B is a true chimeric relocation (-2,095,535 bp) connecting reference 3.66 Mb to 1.57 Mb across 7 short intermediate nodes. It exists bit-identically in Stage 5 `unitig_1631`. The specific mechanism involves multi-hop repeat search and cross-locus repeat homology. |
| **H2** | **Terminal contig boundaries drive Stage 6 polishing errors** | **`OBSERVED (BOUNDARY ERROR)` / `SUPPORTED (HYPOTHESIS)`** | **POLISHING_BOUNDARY_REGRESSION_SUPPORTED**. At Stage 5 (pre-polish), Candidate B had *fewer* mismatches (24 vs 29) and *identical* indels (4 vs 4) compared to Baseline A. 84.6% of Candidate B's Stage 6 errors (132/156) are localized within $\le 350\text{ bp}$ of contig termini. In deep internal contig sequence ($> 350\text{ bp}$), Candidate B has fewer errors than Baseline A (20 vs 22). |
| **H3** | **Conservative repeat withholding left short unitigs that were excluded by output length thresholds** | **`OBSERVED (EXCLUSION)` / `SUPPORTED (COVERAGE FACTOR)`** | **OUTPUT_FILTER_SEQUENCE_LOSS_SUPPORTED**. In Candidate B, Stage 5 active unitigs rose to 1,728 (R1-001) and 704 (R1-003). 1,561 unitigs totalling 235 kb in R1-001 fall below the 500 bp threshold and were not emitted into `contigs.fa`. At Stage 5 (all contigs evaluated), Candidate B had higher genome fraction (+0.38 to +0.41 pp) than Baseline A. |
| **H4** | **QUAST alignment artificially split or recounted misassemblies at contig breaks** | **`REFUTED (FOR R1-001 NEW EVENT)`** | Both versions share exactly 3 common contig groups (accounting for 4 misassemblies). The 5th misassembly in Candidate B is a single, distinct, continuous contig (`contig_16`) with an authentic 2.1 Mb relocation breakpoint. It is not an artifact of contig fragmentation or QUAST recounting. |
| **H5** | **Observed divergence reflects cumulative A/B version drift rather than search admissibility alone** | **`OBSERVED (MULTIPLE MODULE DIFFS)`** | Baseline A is commit `b9727c3` (Phase FB-P0-02). Candidate B includes P1 fragment span verification, single-use lifecycle guards, repeat completeness probes, and short prefix reuse across 1,500+ lines of diff. Cumulative changes are an acknowledged factor. |

---

## 2. Detailed Evidentiary Breakdown

### Hypothesis H1: New Chimeric Join Accepted by exSPAnder
- **Level**: **`SUPPORTED`** (Label: `NEW_JOIN_ERROR_CONFIRMED`)
- **Direct Observations (`OBSERVED`)**:
  - QUAST reports 4 extensive misassemblies in Baseline A and 5 in Candidate B on `R1-001`.
  - The 5th misassembly is exclusively on Candidate B's `contig_16` (71,975 bp).
  - Breakpoint: Contig coordinate 8,704 / 8,705 bp.
  - Left flank: `NC_000913.3: 3,657,472..3,666,175` (8,704 bp, 99.53% identity).
  - Right flank: `NC_000913.3: 1,570,641..1,633,911` (63,271 bp, 99.99% identity).
  - Inconsistency: $-2,095,535\text{ bp}$.
- **Intermediate Evidence (`SUPPORTED`)**:
  - In `stage_04_multik_127_simplified.fa`, the two flanks exist as `unitig_1927` (5,073 bp) and `unitig_169` (12,633 bp).
  - In `stage_05_exspander.fa`, Candidate B formed `unitig_1631` (71,975 bp), containing the chimeric join connecting `unitig_1927` to `unitig_169` through 6 short nodes: `unitig_1926` (267 bp) $\to$ `unitig_3275` (184 bp) $\to$ `unitig_1475` (197 bp) $\to$ `unitig_1476` (255 bp) $\to$ `unitig_4349` (479 bp) $\to$ `unitig_4465` (312 bp) $\to$ `unitig_3919` (398 bp).
  - Baseline A did not traverse this path because its search was bounded by `MAX_INTERMEDIATE_NODES = 4`.
  - The reference sequence at both flanks shares an identical 31-mer repeat motif (`AACGTTATCAGGT[A/G]TGTTTAAAGCTGTTCTG`). Paired-end read pairs from the repeat family provided false-positive support that satisfied the fragment span verification filter.

---

### Hypothesis H2: Polishing Boundary Read Recruitment Regression
- **Level**: **`SUPPORTED`** (Label: `POLISHING_BOUNDARY_REGRESSION_SUPPORTED`)
- **Direct Observations (`OBSERVED`)**:
  - In `R1-001`, Candidate B's mismatches rose from $2.14$ to $2.38$ per 100 kbp (+11.2%).
  - Candidate B's indels rose from $0.48$ to $1.03$ per 100 kbp (+114.6%, raw count: $22 \to 47$).
- **Intermediate Evidence (`SUPPORTED`)**:
  - In Stage 5 (pre-polish, `stage_05_exspander.fa`), Candidate B had **24 raw mismatches** vs Baseline A's **29 raw mismatches** (fewer errors!).
  - In Stage 5, Candidate B had **4 raw indels** vs Baseline A's **4 raw indels** (identical!).
  - **100% of the mismatch excess (85 introduced in B vs 69 in A) and indel excess (43 introduced in B vs 18 in A) occurred during Stage 6 consensus polishing.**
  - **Boundary Localization**:
    - In `Bin 101-350 bp` from contig ends, Candidate B contains **88 mismatches** and **44 indels** (132 total errors, accounting for **84.6% of all errors in the assembly**).
    - 130 of these 132 boundary errors (98.5%) were introduced during polishing.
    - In `Bin >350 bp` (deep internal contig sequence), Candidate B has only **18 mismatches** and **2 indels** (total 20 errors), which is **FEWER** than Baseline A's 22 internal errors (18 mismatches, 4 indels).
  - Forensic conclusion: Candidate B's base sequence generated by exSPAnder is cleaner than Baseline A's. The elevated mismatch and indel rates in the final FASTA are an artifact of read-mapping recruitment at newly created terminal contig boundaries during Stage 6 consensus polishing.

---

### Hypothesis H3: Sequence Loss from Length Threshold Filtering
- **Level**: **`SUPPORTED`** (Label: `OUTPUT_FILTER_SEQUENCE_LOSS_SUPPORTED`)
- **Direct Observations (`OBSERVED`)**:
  - QUAST genome fraction on $\ge 500\text{ bp}$ contigs dropped slightly in Candidate B:
    - R1-001: $98.364\% \to 98.281\%$ ($-0.083$ pp).
    - R1-003: $98.330\% \to 98.233\%$ ($-0.097$ pp).
- **Intermediate Evidence (`SUPPORTED`)**:
  - Active unitigs post-exSPAnder rose from 315 to 1,728 in R1-001 (5.49×) and from 300 to 704 in R1-003 (2.35×).
  - In R1-001 Candidate B, **1,370 unitigs (185,579 bp)** are $< 200\text{ bp}$ and were filtered by `--min-contig-len 200`.
  - Another **191 unitigs (49,383 bp)** are between 200 and 499 bp.
  - In total, **234,962 bp** of assembled sequence in 1,561 unitigs fell below the 500 bp threshold.
  - In Stage 5 (where all contigs $\ge 200\text{ bp}$ were evaluated by QUAST):
    - Candidate B achieved **98.835% genome fraction** on R1-001 (vs 98.456% for Baseline A, **+0.379 pp higher**).
    - Candidate B achieved **98.857% genome fraction** on R1-003 (vs 98.444% for Baseline A, **+0.413 pp higher**).
  - The drop in Stage 6 genome fraction is strictly an artifact of length threshold filtering applied to unmerged short unitigs.

---

## 3. Final Multi-Mechanism Synthesis

The forensic investigation establishes that three complementary mechanisms explain the Phase F20 benchmark results:
1. **`NEW_JOIN_ERROR_CONFIRMED`**: A single false-positive join occurred on `R1-001` (`contig_16`), caused by multi-hop search expansion traversing a 6-short-node repeat path whose homologous repeat units satisfied the paired-end fragment span filter.
2. **`POLISHING_BOUNDARY_REGRESSION_SUPPORTED`**: The entire increase in mismatches and indels was introduced during Stage 6 consensus polishing, concentrated within $\le 350\text{ bp}$ of contig termini.
3. **`OUTPUT_FILTER_SEQUENCE_LOSS_SUPPORTED`**: The apparent decrease in genome fraction was caused by length threshold filtering of unmerged short unitigs, whereas the raw graph at Stage 5 actually captured greater genome coverage than Baseline A.
