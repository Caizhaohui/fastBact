# Phase F21-A Independent Forensic Audit Report

**Audit ID**: AUDIT-20261023-016  
**Auditor**: Codex CLI (gpt-6.1-sol, medium reasoning effort)  
**Overall Forensic Verdict**: **FAIL for an unqualified PASS on forensic accuracy and causal attribution**  
**Stage C0-02R Acceptance**: **Strictly NOT ACCEPTED**

The principal measurements are reproducible, including the additional R1-001 misassembly and terminal concentration of final base errors. However, incorrect execution parameters, erroneous path accounting, a nonexistent baseline search-cap attribution, and unsupported locus-specific causal claims prevent the requested PASS verdict.

This audit was strictly read-only. No files were created or modified.

## 1. Manifest, provenance, checksums, and Slurm

**Verified:**

- Baseline commit: `b9727c3ab372f7b995e439d92b72c7a4bfd051cf`.
- Candidate evidence commit: `0606be423d187771f5925a83da2cc377c40c2a3d`.
- Candidate code commit: `a9930471a73612159a83694d3d661f060042efc8`.
- Forensic base: `4b956a46935248779e17b8e3453f4bcb188f8e1f`.
- Candidate `src/` tree: `30967bf4940844aeb6b744732bcad6d2a72daab4`.
- Both binary SHA-256 digests, both Cargo.lock digests, four FASTQ digests, and two reference digests passed fresh verification: **10/10 entries**.
- All **242 F20 artifact checksums** and **140 Stage 5 QUAST checksums** passed.
- The installed evaluator and archived QUAST logs identify **QUAST v5.2.0**.
- An archived compiler record confirms `rustc 1.89.0 (29483883e 2025-08-04)`; this corroborates the recorded toolchain without independently rebuilding either binary.

**Material manifest error:** [forensic_manifest.tsv:20](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21-forensics-worktree/docs/development/evidence/C0_02R/phase_f21/forensic_manifest.tsv:20) claims:

```text
-k 21,33,55,77 -t 8 --min-contig-len 200
```

The actual benchmark command supplied neither that k-mer schedule nor that length option. It used automatic k-mer selection and the default assembler minimum length, **500 bp**. Execution logs report **[31, 63, 127]**. QUAST separately used **`--min-contig 200`**.

Scheduler verification has these limits:

| Job | Evidence verified | Limitation |
|---|---|---|
| 3001561 | Archived accounting: `COMPLETED`, exit `0:0`, elapsed `00:13:13`; eight numbered steps completed; `qfree`, `vnode1`, 8 CPUs, 12G | Fresh scheduler query unavailable |
| 3007620 | Allocation confirms `qfree`, `vnode1`, 8 CPUs, 12G; logs show four QUAST evaluations and successful script completion | Saved allocation was captured while `RUNNING`; no completed accounting record was found among the inspected F21 logs |

A fresh `sacct` query failed with **Operation not permitted**. Job 3007620’s script completion is corroborated, but completed scheduler accounting cannot be independently certified.

## 2. Misassembly events and graph-path reconstruction

**The final R1-001 event reconciliation passes.** Both raw `mis_contigs.info` and `all_alignments_contigs.tsv` confirm:

| Event | Baseline A | Candidate B |
|---|---|---|
| Two inversions near 1.20 Mb | `contig_6` | `contig_7` |
| Relocation near 1.87 Mb ↔ 291 kb | `contig_149` | `contig_154` |
| Relocation near 1.87 Mb ↔ 1.97 Mb | `contig_163` | `contig_166` |
| Additional relocation near 3.66 Mb ↔ 1.57 Mb | Absent | `contig_16` |

The three common contig pairs are **sequence-identical**, establishing that the four shared events are identical beyond their reported coordinates.

Candidate B’s additional event is verified:

- Contig length: **71,975 bp**.
- Breakpoint: **8,704/8,705**.
- Left reference alignment: **3,657,472–3,666,175**, 99.53% identity.
- Right reference alignment: **1,570,641–1,633,911**, 99.99% identity.
- QUAST inconsistency: **−2,095,535 bp**.

The Stage 5 sequence `unitig_1631` is **exactly identical** to final `contig_16`. Thus, polishing did not create this join. **`NEW_JOIN_ERROR_CONFIRMED` is justified as a reference-relative assembly defect**, rather than a QUAST recounting artifact.

### Path-accounting corrections

The claimed “six intermediate nodes totaling 1,825 bp” is incorrect. Exact sequence tiling gives these positions within Stage 5 `unitig_1631`:

| Unitig | Length | Orientation | Zero-based start |
|---|---:|---|---:|
| 1927 | 5,073 | Reverse complement | 3,590 |
| 1926 | 267 | Forward | 8,536 |
| 3275 | 184 | Reverse complement | 8,676 |
| 1475 | 197 | Reverse complement | 8,733 |
| 1476 | 255 | Reverse complement | 8,803 |
| 4349 | 479 | Reverse complement | 8,931 |
| 4465 | 312 | Reverse complement | 9,283 |
| 3919 | 398 | Reverse complement | 9,468 |
| 169 | 12,633 | Forward | 9,739 |

There are **seven** listed short nodes, totaling **2,092 bp**. The reported **1,825 bp** excludes `unitig_1926`. Their sequence overlaps also need to be accounted for explicitly when calculating traversal separation.

### Causal attribution does not pass

The [locus analysis](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21-forensics-worktree/docs/development/evidence/C0_02R/phase_f21/misassembly_locus_analysis.md:88) overstates the mechanism:

- The actual baseline `b9727c3:src/resolve/exspander.rs` contains **no `MAX_INTERMEDIATE_NODES = 4` search cap**. It implements direct/local continuation and repeat-exit decisions.
- Candidate R1-001 logs report median insert size **374 bp**, MAD **53 bp**. The mapper’s formula yields a maximum near **559 bp**, which does not substantiate the document’s claim that its calculated **1,069 bp** separation was within the insert window.
- No inspected locus-specific decision record establishes the supporting pair identities, competing scores, or alleged empty/truncated alternative set.
- Candidate code explicitly withholds decisions when search status is `Truncated`; a particular bypass or alternative explanation requires direct evidence.

The new Stage 5 join is confirmed. Its attribution specifically to cap relaxation, prefix reuse, and false paired-end support remains **unverified**.

## 3. Eight-condition Stage 5/Stage 6 comparison

Nine principal report fields were checked across all eight conditions, giving **72 matching comparisons**: misassemblies, local misassemblies, mismatch/indel rates, genome fraction, N50, NGA50, largest contig, and aligned length.

| Isolate | Version | Stage | Extensive misassemblies | Genome fraction |
|---|---|---|---:|---:|
| R1-001 | A | 5 | 4 | 98.456% |
| R1-001 | A | 6 | 4 | 98.364% |
| R1-001 | B | 5 | 6 | 98.835% |
| R1-001 | B | 6 | 5 | 98.281% |
| R1-003 | A | 5 | 2 | 98.444% |
| R1-003 | A | 6 | 2 | 98.330% |
| R1-003 | B | 5 | 2 | 98.857% |
| R1-003 | B | 6 | 1 | 98.233% |

For R1-001, the requested base-error counts are verified:

| Measurement | A | B |
|---|---:|---:|
| Stage 5 mismatches | 29 | 24 |
| Stage 5 indel events | 4 | 4 |
| Stage 5 indel length | 176 bp | 176 bp |
| Stage 6 mismatches | 98 | 109 |
| Stage 6 indel events | 22 | 47 |
| Net mismatch increase | 69 | 85 |
| Net indel-event increase | 18 | 43 |

The final excess therefore emerges during the **Stage 5-to-final-output transition**. However, **85/69 and 43/18 are net changes, not counts of newly introduced loci**.

The locus table itself classifies:

- A: **78 new mismatch loci and 19 new indel loci**.
- B: **93 new mismatch loci and 44 new indel loci**.

Furthermore, Stage 6 includes polishing **and curation**, with different evaluated contig populations. Assigning the entire change specifically to consensus read recruitment requires matched-contig analysis. The polisher changes bases without inserting or deleting sequence, so increased QUAST indel calls also require alignment-level interpretation.

**Table schema error:** [stage5_stage6_comparison.tsv:1](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21-forensics-worktree/docs/development/evidence/C0_02R/phase_f21/stage5_stage6_comparison.tsv:1) labels raw all-length totals as `TotalBpGte500` and `ContigsGte500`. For R1-001 B:

- Raw Stage 5: **1,728 sequences, 4,808,057 bp**.
- QUAST-eligible ≥200 bp: **358 sequences, 4,622,478 bp**.
- Actual ≥500 bp: **167 sequences, 4,573,095 bp**.

These categories must be separated.

## 4. Error loci and terminal-boundary distribution

**The coordinate and distribution measurements pass.**

All **396 final error records** matched raw QUAST records by contig, contig position, reference position, and alleles. Their recorded contig lengths matched the final FASTA.

Every recorded distance reproduced the table’s one-based convention:

```text
min(position, contig_length − position + 1)
```

For R1-001:

| Distance bin | A mismatches / indels | B mismatches / indels |
|---|---:|---:|
| 0–100 bp | 10 / 1 | 3 / 1 |
| 101–350 bp | 70 / 17 | 88 / 44 |
| >350 bp | 18 / 4 | 18 / 2 |

Thus:

- **132/156 = 84.615%** of B’s final error events are in the **101–350 bp** zone.
- B has **20 internal errors**, compared with A’s **22**.
- The table marks **130/132** B errors in that terminal zone as absent from Stage 5.

The pre-existing/new tags reproduce reference-coordinate-and-allele matching against Stage 5 calls. This supports temporal association but does not independently establish the read-recruitment mechanism. These are **error-event counts**, not affected-base totals or exposure-normalized regional error rates.

## 5. Unitig fate accounting

The complete R1-001 B Stage 5 length accounting passes:

| Length class | Count | Total sequence length |
|---|---:|---:|
| <200 bp | 1,370 | 185,579 bp |
| 200–499 bp | 191 | 49,383 bp |
| ≥500 bp | 167 | 4,573,095 bp |
| **Total** | **1,728** | **4,808,057 bp** |

Below 500 bp, there are **1,561 unitigs totaling 234,962 bp**:

- **90.34% by unitig count**.
- **4.89% by total sequence length**.

The report’s output-filter explanation needs correction:

1. The assembler used its **500 bp default output threshold**.
2. QUAST used **200 bp**, not 500 bp.
3. The final B assembly contains only **167 contigs**. The 191 sequences of 200–499 bp were **not emitted into final `contigs.fa`**, contrary to the contrast report.
4. Summed unitig lengths do not measure lost unique reference coverage.

Short-sequence exclusion is strongly supported. Claiming it **fully explains** the genome-fraction reduction requires reference-coordinate coverage accounting that separates filtering, polishing, trimming, and deduplication.

## 6. R1-003 contrast

The final improvements are verified:

- Extensive misassemblies: **2 → 1**.
- Indel events: **2 → 0**.
- NGA50: **63,648 → 70,010 bp**, approximately **+10.0%**.
- Genome fraction nevertheless declines: **98.330% → 98.233%**.

The **+590,017 bp** relocation on A’s `contig_4` is absent from B. B’s Stage 5 misassembly report already lacks that event, placing the structural difference before polishing.

However, the assertion that **β(curr) and `SearchStatus::Truncated` specifically caused its elimination** is unsupported by the inspected saved logs. The source contains the conservative rule; that does not demonstrate its execution at this locus. [r1_003_contrast.md:49](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21-forensics-worktree/docs/development/evidence/C0_02R/phase_f21/r1_003_contrast.md:49)

The shared repeat-associated locus is strongly supported:

- R1-003’s A `contig_16` and B `contig_19` are sequence-identical.
- Corresponding breakpoint-centered reference windows are identical between isolates.
- Within each isolate, the two locus windows share **847 distinct 31-mers**.

This supports homologous repeat context across isolates, but it does not establish identical supporting read-pair histories or decision mechanisms.

## 7. H1–H5 evidentiary assessment

| Hypothesis | Independent assessment |
|---|---|
| **H1: New chimeric join** | **OBSERVED** for the new reference-relative Stage 5 join. Specific paired-end/search explanation remains **UNVERIFIED**. Correct the path count and remove the nonexistent baseline-cap explanation. |
| **H2: Boundary-associated polishing regression** | **SUPPORTED** as a boundary-associated Stage 5-to-final regression. Specific read-recruitment causality remains **UNVERIFIED** without matched sequence and mapping evidence. |
| **H3: Short-unitig filtering** | **OBSERVED** for counts and actual output exclusion; **SUPPORTED** as a coverage-loss contributor. “Strictly explained” is too strong, and the QUAST-threshold explanation is incorrect. |
| **H4: QUAST recounting artifact** | **REFUTED for the additional final R1-001 event**: a distinct continuous chimeric sequence exists before polishing. This does not eliminate alignment effects from other metrics. |
| **H5: Cumulative version differences** | **OBSERVED** that A/B differ across multiple source modules; cumulative differences are a valid confounder. The claimed four-hop-cap relaxation is unsupported and contradicted by the stated baseline source. |

The [hypothesis matrix](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21-forensics-worktree/docs/development/evidence/C0_02R/phase_f21/causal_hypothesis_matrix.md:14) repeatedly combines verified observations with unverified mechanisms under one `SUPPORTED` label. Those claims should be separated.

## Production-source integrity and recommendations

**Production `src/` is untouched by F21-A.** Staged and working-tree source diffs are empty; no untracked source files were found. HEAD, F19 candidate, F19 code commit, and F18 audited commit share `src/` tree:

```text
30967bf4940844aeb6b744732bcad6d2a72daab4
```

Only the pre-existing untracked F21 evidence directory appeared in repository status. Stage 4 simplified FASTAs are also byte-identical between A/B for both isolates, strengthening localization of the structural divergence to subsequent processing.

Recommended corrections before forensic approval:

1. Correct execution parameters, length thresholds, table headers, node counts, and path-length arithmetic.
2. Preserve `NEW_JOIN_ERROR_CONFIRMED`, while downgrading unsupported decision-level explanations.
3. Obtain completed accounting for Job 3007620.
4. Compare matched retained contigs before/after polishing and curation; distinguish new, resolved, and alignment-dependent error calls.
5. Measure unique reference coverage removed by output filtering.
6. Capture locus-specific candidate paths, pair evidence, scores, β values, and search status before selecting a repair.

**C0-02R remains strictly NOT ACCEPTED.** The measurements establish a real regression; the current forensic narrative does not yet establish all of its claimed causes.