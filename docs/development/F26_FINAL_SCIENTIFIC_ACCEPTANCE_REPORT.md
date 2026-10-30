# Phase F26 Final Independent Scientific Acceptance Report

**Document ID:** F26-SOL-FINAL-AUDIT-001  
**Date:** 2026-10-30  
**Role:** Sol, independent scientific reviewer and acceptance auditor under Multi-AI Collaboration Policy V3.2  
**Audited worktree:** `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f26-robustness-worktree`  
**Audited HEAD:** `998948fb25f535249ac7603ecfe7d7d24ee2275b`  
**Frozen comparator:** `02872f81811402d2bd603ff7c0a2ecf8b3dbe528`  
**Disposition:** `PARTIALLY_VALIDATED / F26_ACCEPTANCE_WITHHELD`

## 1. Executive Summary & Audit Mandate

F26 evaluates the frozen fastBact implementation beyond the historical *E. coli* panel, compares it with established short-read assemblers, characterizes coverage sensitivity, and permits hardening only after localized failure evidence. This Stage H audit independently assesses Stages D, E, F and G against both governing contracts. It does not modify the algorithm, waive the scientific contract, or grant integration authorization.

The measurements are credible within their recorded scope: **40/40 assemblies and QUAST evaluations passed**, including **28/28 fastBact executions**. Fresh checks reproduce all 40 assembly hashes, all 520 QUAST metric values, the JSON/TSV correspondence, and all 44 distinct benchmark input/reference file hashes. All 40 recorded benchmark jobs are COMPLETED on `qfree`, exit `0:0`, eight CPUs and 14G requested memory. A fresh `cargo test --release -- --nocapture` audit job passes **213/213 tests**, with zero failures or ignored tests. `src/` is unchanged against the frozen baseline.

**Unconditional scientific acceptance is withheld.** The requested labels `F26_SCIENTIFICALLY_VALIDATED` and `READY_FOR_USER_ACCEPTANCE` are not supported as whole-phase conclusions under the unchanged contract:

| Finding | Severity | Independent determination |
| --- | --- | --- |
| F26-H01: invalid formal replicate pooling | P1, acceptance blocker | Collated identity fields are empty. Actual replicate-1 Git HEAD differs from replicates 2–3 on all four taxa. This violates benchmark contract §10.7 despite identical binary, inputs and outputs. |
| F26-H02: F25 oracle non-regression not replayed | P1, acceptance evidence gap | None of the four frozen full-depth *E. coli* oracles occurs in the 40-run collection. Zero source diff supports preservation but does not measure contract §2 metric deltas. |
| F26-H03: nominal depth is not bacterial chromosome depth | P2, interpretation limitation | Variable read lengths, a phage-dominated MW2 library, and PAO1's keep-all 150× clamp prevent an unqualified bacterial-depth or GC-threshold claim. |
| F26-H04: biological attribution remains incomplete | P2, interpretation limitation | Raw QUAST counts are valid observations; worker reports do not provide Stage D/G per-contig Staph1N attribution or resolve the MW2 extensive events. No biological exemption is granted. |

The empirical bitwise-output result is positive for all four full-depth fastBact triplets. It must be distinguished from formal contract-compliant Stage G certification. Runtime advantages are **1.330×–19.509×**, not a general 10×–50× claim. Successful process completion is not proof of mathematical invariants.

## 2. Dataset Provenance & Experimental Rigor Audit

### 2.1 Registered panel and primary metadata

The eight-row parent manifest contains four historical oracle rows and the following four new species. The panel spans four genera, two Gram-positive and two Gram-negative taxa, reference GC **32.83%–66.56%** (32.8%–66.6% rounded), and chromosome lengths 2.82–6.26 Mb. Gram classification is species metadata, not a phenotype measured from these reads.

| Dataset | Species / reference lineage | Gram | GC (%) | Chromosome bp | SRA | BioProject / BioSample | RefSeq |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CS-SA01 | Staphylococcus aureus MW2 | Positive | 32.83 | 2820462 | SRR33682404 | PRJNA1263016 / SAMN48627418 | NC_003923.1 |
| CS-BS01 | Bacillus subtilis subsp. subtilis str. 168 | Positive | 43.51 | 4215606 | SRR40739147 | PRJNA1530750 / SAMN63192558 | NC_000964.3 |
| CS-PA01 | Pseudomonas aeruginosa PAO1 | Negative | 66.56 | 6264404 | SRR41244734 | PRJNA1545662 / SAMN64277765 | NC_002516.2 |
| CS-SE01 | Salmonella enterica subsp. enterica serovar Typhimurium str. LT2 | Negative | 52.22 | 4857450 | SRR41010285 | PRJNA1539015 / SAMN63983458 | NC_003197.2 |

The retained NCBI primary responses `/hpcfs/fhome/caizhh/.local/state/f26-audit-001/sra.xml` and `ncbi_reference_summary.json` were inspected directly. They contain the four requested public SRA runs and the six distinct versioned chromosome records, with matching accessions, BioSample identifiers, spot counts, reference titles and lengths. The retained Stage B/C audit documents the BioProject/layout checks and their original retrieval endpoints. This audit relies on those archived primary responses and local files; it does not claim a fresh database fetch or exact identity between an evolved isolate and its ancestral reference.

The latest section F26-AUDIT-003 of [the Stage B/C audit](F26_INDEPENDENT_AUDIT_STAGE_BC.md) supersedes its historical blocked sections. Its MW2-lineage conclusion is accepted with its stated limits: the SRA title conflicts with the MW2 attributes, but retained alignment/discriminatory-marker evidence supports MW2-lineage bacterial content plus Staph1N co-culture. The reference is chromosome-only `NC_003923.1`, not an exact isolate reconstruction. Other accepted contexts include chloramphenicol-evolved *B. subtilis*, a PAO1-GFP mutant, and an ALE-derived LT2 sample. These confound causal GC comparisons and absolute reference-discrepancy interpretation. Identical inputs do not eliminate every bias in whole-assembly rankings when viral recovery differs.

### 2.2 File integrity and read accounting

Fresh raw and normalized-sequence SHA-256 calculations match **8/8 parent-manifest rows** (six distinct references). Normalization concatenates uppercase sequence lines, excluding FASTA headers and whitespace. Length and rounded GC also match every row. The four cross-species digest pairs are:

| Dataset | Raw reference SHA-256 | Normalized sequence SHA-256 |
| --- | --- | --- |
| CS-SA01 | `19ba62a339065c9c498efb72bb2491e372344e4f23a4ac1e3e34e7345a25c63c` | `7c4947b9106af48ffee6bdf8dc29753333d97d2824092fb8228b57fc753d315d` |
| CS-BS01 | `a334e891ffc0e307f23f48842775d3383177a9d9cb5d5075b552a2cccddfe139` | `71572ae0b98f44a6382b29e194e56c1f7d87477be8f53915fb06816f886b17b7` |
| CS-PA01 | `902fa552c8c25aaaf82eccfed000e359d6a04c8744d9ec0a8134699cb6ba53e6` | `f969830bbe6aa58beccd852ae6e129a9489a92557ad082283f639b6e8539a340` |
| CS-SE01 | `59fceaea5ca57836e02b6c4ad3582847623c87268f7e5511874909eedb4acafe` | `37854e9ce0ff2b3aa4a3394bfab5a63f4922d5e9479acb930558979ffca81326` |

Using the richer retained per-job `raw/*/summary.json` identities, fresh hashes match **44/44 distinct benchmark input/reference paths**: eight original FASTQ mates, 32 cross-species subsample mates, and four chromosome FASTAs. Original-mate identities also correspond to [F26_FASTQ_CHECKSUMS.tsv](F26_FASTQ_CHECKSUMS.tsv). Stage B/C had already checked complete gzip framing and mate-ID synchronization; those full-stream checks were not rerun here. This audit verifies the present bytes, not independent provenance merely from a hash.

LT2 has **3,359,296 SRA spots**, **3,353,233 paired records per mate**, and **6,063 singleton records**, reconciled in Stage B/C. All assemblers use only the paired files. The parent spot total is not the paired benchmark count. No reads were deduplicated for this audit.

### 2.3 Stage E reproducibility and achieved-depth limits

[F26_STAGE_E_MANIFEST.tsv](../../F26_STAGE_E_MANIFEST.tsv) has 32 planned rows: 16 cross-species and 16 historical *E. coli* depth inputs. Every notes field records **seed=42**, source, target and genome size. The retained `logs/f26_stage_e_subsamples_3075001.err` ends `rows=32 generated=32 skipped=0`. Only the 16 cross-species rows have assembly results in the supplied collection.

Source inspection of the subsampler confirms sorted `random.Random(seed).sample(range(N), k)` pair indices, paired streaming, preservation of original order, and deterministic gzip headers (`mtime=0`, no filename). The selector uses `k=round(target_depth × genome_size / (2 × first_R1_length))`, clamped to N. Fresh subsample hashes match retained job identities; generation itself was not independently replayed. One seed and one assembly per gradient cell do not establish sampling uncertainty or depth-cell determinism.

Manifest `total_spots` remains the parent SRA field; achieved pair counts and subsample hashes must be read from the generation log and per-job summaries. R1 length used is 151 bp for BS01, 151 bp with variable-length warning for SA01, 251 bp with warning for PA01, and 150 bp for SE01. PA01's nominal 150× requested 1,871,834 pairs but retained **all 1,518,403 pairs**; its formula depth is **121.678×**, not 150×. The corresponding contig hash and all QUAST fields equal the full-library run. Formula depths for variable-length libraries are not actual counted-base depths.

MW2 is especially important: the retained strain-discrimination evidence measures **69.9473× bacterial chromosome depth** and approximately **2,984.21× phage depth** in the full library. The generator's full-library formula is 239.810× because it treats all sequence as bacterial and uses the first read length. Scaling the retained chromosome depth by sampling fractions predicts approximately 8.75×, 14.58×, 29.17× and 43.75× bacterial depth at nominal 30×, 50×, 100× and 150×. These are expectations, not direct mapping measurements of the subsamples. The low-depth recovery transition therefore cannot be assigned to AT-rich sequence or a true 30× bacterial threshold alone.

## 3. Stage D Cross-Species Benchmark Audit (vs Competitors)

### 3.1 Execution and exact comparative measurements

The 16-run Stage D matrix is four taxa × four tools, replicate 1. Recorded versions are fastBact 0.1.0, **SPAdes v4.2.0**, **SKESA 2.5.1**, and **Unicycler v0.5.1**. Raw QUAST stdout confirms **v5.2.0** on every run, beyond the worker reports' version-pin statements. QUAST argv uses the correct chromosome reference, eight threads and `--min-contig 200`. GNU Time elapsed and maximum-RSS fields independently match all 40 collected run records. Resources below measure assembly only, excluding QUAST.

All tools use paired short reads and eight threads, with the same 14G Slurm allocation. Commands match the documented modes: SPAdes `--isolate`, SKESA comma-joined paired `--reads`, Unicycler short-read normal mode, and fastBact `--min-len 500 --kmers auto --depth 120` without auto-subsampling. The fastBact 500-bp output cutoff versus competitor 200-bp products is a disclosed asymmetry. N50/largest-contig/total-length are whole-assembly statistics, whereas genome fraction and error rates use chromosome alignments. In particular, viral contigs can affect the former.

The following values come directly from the current JSON, cross-checked against each retained QUAST TSV. Error rates are per 100 kbp. Wall times are displayed to 0.01 s; MiB is the actual unit of the misleadingly named `peak_rss_mb` field.

| Dataset | Tool | GF % | N50 bp | NGA50 bp | Largest bp | Extensive | Local | Mismatch/100kb | Indel/100kb | Wall s | RSS MiB |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CS-BS01 | fastbact | 98.576 | 78260 | 78260 | 354932 | 0 | 0 | 0.94 | 0.82 | 153.08 | 1425.3828125 |
| CS-BS01 | spades | 98.946 | 624246 | 624246 | 1046381 | 0 | 1 | 6.46 | 0.45 | 552.21 | 5239.390625 |
| CS-BS01 | skesa | 98.927 | 926769 | 926769 | 1295672 | 1 | 0 | 0.22 | 0.36 | 278.44 | 10447.70703125 |
| CS-BS01 | unicycler | 98.903 | 2220390 | 2220390 | 2220390 | 0 | 1 | 0.24 | 0.41 | 2683.36 | 5245.00390625 |
| CS-PA01 | fastbact | 98.74 | 23514 | 23265 | 90388 | 2 | 0 | 2.34 | 0.63 | 180.57 | 1527.0703125 |
| CS-PA01 | spades | 99.516 | 571752 | 571752 | 1344826 | 4 | 2 | 4.11 | 0.67 | 1016.62 | 5245.16015625 |
| CS-PA01 | skesa | 99.334 | 359176 | 359176 | 850497 | 3 | 1 | 2.2 | 0.58 | 240.12 | 8846.60546875 |
| CS-PA01 | unicycler | 99.313 | 546789 | 479845 | 1275111 | 3 | 2 | 1.7 | 0.61 | 3318.30 | 5249.26171875 |
| CS-SA01 | fastbact | 98.683 | 25381 | 25462 | 77704 | 2 | 0 | 2.11 | 0.25 | 148.55 | 1196.7421875 |
| CS-SA01 | spades | 99.008 | 172807 | 175250 | 473664 | 0 | 0 | 1.79 | 0.29 | 287.78 | 4395.28515625 |
| CS-SA01 | skesa | 98.655 | 175451 | 175451 | 372593 | 1 | 0 | 0.32 | 0.0 | 1110.08 | 10136.2265625 |
| CS-SA01 | unicycler | 98.785 | 343472 | 343472 | 721378 | 0 | 0 | 0.32 | 0.07 | 1191.97 | 4410.3359375 |
| CS-SE01 | fastbact | 98.61 | 47464 | 48812 | 178426 | 1 | 1 | 0.48 | 0.06 | 106.27 | 7854.46875 |
| CS-SE01 | spades | 98.895 | 192527 | 192527 | 627401 | 1 | 3 | 7.5 | 0.17 | 413.00 | 5237.25390625 |
| CS-SE01 | skesa | 98.771 | 302057 | 302057 | 736780 | 1 | 2 | 0.58 | 0.1 | 255.29 | 9003.1328125 |
| CS-SE01 | unicycler | 98.737 | 298595 | 298595 | 820290 | 1 | 2 | 0.48 | 0.08 | 2073.25 | 5245.00390625 |

### 3.2 Scientific interpretation

fastBact is faster on all 12 matched competitor comparisons, but the requested **10×–50×** range is unsupported. Correct competitor/fastBact wall-time ratios, rounded to three decimals, are:

| Dataset | SPAdes / fastBact | SKESA / fastBact | Unicycler / fastBact |
| --- | --- | --- | --- |
| CS-BS01 | 3.607× | 1.819× | 17.529× |
| CS-PA01 | 5.630× | 1.330× | 18.377× |
| CS-SA01 | 1.937× | 7.473× | 8.024× |
| CS-SE01 | 3.886× | 2.402× | 19.509× |

These are single-run ratios, not confidence intervals. fastBact uses 70.89%–88.19% less peak RSS than the competitors on BS01/PA01/SA01. LT2 reverses part of that advantage: fastBact is 7,854.46875 MiB, approximately 49.97% above SPAdes and 49.75% above Unicycler, though 12.76% below SKESA. A uniformly low-memory claim is rejected.

Structural accuracy is mixed. BS01 has zero extensive and local events for fastBact. PA01 fastBact has fewer extensive events (2 versus SPAdes 4, SKESA 3, Unicycler 3). SA01 has more (2 versus 0, 1, 0); Stage G shows the two events are repeatable, so they cannot be dismissed as single-run noise. SE01 has one extensive event for all tools, and fastBact's one local event is below each competitor's count. No cross-species raw event is reclassified as a biological artifact without sequence-level evidence.

Genome recovery exceeds 98.5% in all full-depth Stage D runs, but fastBact generally recovers less chromosome than competitors. Its N50 is lower than every competitor on every taxon; competitor N50 is approximately 4.06×–28.37× larger. NGA50 also exposes a substantial contiguity gap. These are material tradeoffs, not grounds for reporting superiority from speed alone.

Base-error ranking is organism-dependent. fastBact's indel rate is highest on BS01 (0.82) and lowest on SE01 (0.06). On SA01 its mismatch rate is 2.11 versus SPAdes 1.79 and SKESA/Unicycler 0.32. Its SE01 mismatch rate 0.48 ties Unicycler and is below SPAdes 7.50. Reference lineage differences limit attribution of absolute discrepancy counts, but do not erase the measured comparisons.

## 4. Stage E Coverage Robustness & Sensitivity Audit

This table contains the **16 actual cross-species depth assemblies**. Targets are nominal library formula depths, with the PA01 clamp qualified above. Each row is one sampled input and one assembly, not a biological replicate series.

| Dataset / target | GF % | Contigs ≥200 bp | N50 bp | NGA50 bp | Largest bp | Extensive | Local | Mismatch/100kb | Indel/100kb | Wall s | RSS MiB |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CS-BS01_30x | 98.272 | 119 | 66973 | 66458 | 204588 | 0 | 0 | 0.68 | 0.51 | 36.46 | 1096.85546875 |
| CS-BS01_50x | 98.346 | 170 | 39818 | 39195 | 234253 | 0 | 0 | 0.55 | 0.84 | 56.81 | 1033.9765625 |
| CS-BS01_100x | 98.549 | 172 | 47159 | 45013 | 128251 | 0 | 0 | 3.33 | 0.6 | 88.59 | 1152.92578125 |
| CS-BS01_150x | 98.575 | 111 | 91207 | 91207 | 214696 | 0 | 0 | 3.6 | 0.55 | 167.36 | 1554.4453125 |
| CS-PA01_30x | 99.036 | 217 | 53277 | 52412 | 223412 | 4 | 1 | 3.54 | 0.77 | 57.44 | 1104.46484375 |
| CS-PA01_50x | 99.078 | 191 | 59072 | 59067 | 308504 | 3 | 1 | 2.33 | 0.58 | 76.65 | 1081.39453125 |
| CS-PA01_100x | 98.861 | 424 | 30828 | 30428 | 201253 | 2 | 0 | 2.26 | 0.67 | 165.96 | 1454.51171875 |
| CS-PA01_150x | 98.74 | 576 | 23514 | 23265 | 90388 | 2 | 0 | 2.34 | 0.63 | 163.08 | 1533.12109375 |
| CS-SA01_30x | 64.582 | 1713 | 1266 | 777 | 105012 | 0 | 0 | 0.33 | 0.05 | 37.91 | 815.08984375 |
| CS-SA01_50x | 95.199 | 859 | 4927 | 4577 | 63184 | 0 | 0 | 0.48 | 0.0 | 38.39 | 660.1015625 |
| CS-SA01_100x | 98.194 | 217 | 32244 | 32900 | 96025 | 0 | 0 | 0.76 | 0.25 | 46.87 | 710.5390625 |
| CS-SA01_150x | 98.445 | 253 | 28035 | 28326 | 87734 | 0 | 0 | 0.68 | 0.0 | 75.29 | 789.5859375 |
| CS-SE01_30x | 98.11 | 225 | 46607 | 46607 | 132906 | 1 | 1 | 0.52 | 0.04 | 35.50 | 1381.6640625 |
| CS-SE01_50x | 98.401 | 322 | 27151 | 27162 | 91340 | 1 | 1 | 0.77 | 0.15 | 38.76 | 1428.65625 |
| CS-SE01_100x | 98.647 | 220 | 39823 | 39823 | 208121 | 1 | 1 | 1.5 | 0.23 | 69.25 | 3787.87109375 |
| CS-SE01_150x | 98.669 | 155 | 58899 | 58899 | 199626 | 1 | 1 | 0.58 | 0.1 | 76.68 | 4072.14453125 |

- **BS01:** zero extensive/local events at every target and full depth. Genome fraction progresses 98.272 → 98.346 → 98.549 → 98.575%, then 98.576% at full depth. N50 is non-monotonic; 150× (91,207 bp) exceeds full depth (78,260 bp).
- **SA01:** nominal 30× yields 64.582% genome fraction, 1,713 contigs, N50 1,266 bp and NGA50 777 bp. Recovery rises to 95.199%, 98.194% and 98.445% at the next targets. All four gradient rows have zero extensive/local events; full depth has two extensive events and 98.683% genome fraction. This is a real input-fraction sensitivity, strongly confounded by phage dilution; it is not a demonstrated GC-specific mechanism or a measured bacterial 30× limit.
- **PA01:** extensive counts decline 4 → 3 → 2 as nominal targets increase 30× → 50× → 100×; local counts decline 1 → 1 → 0. The clamped top row remains 2/0. N50 and genome fraction do not monotonically improve: lower-depth higher N50 accompanies more extensive events. There is no genuine 150× measurement.
- **SE01:** every gradient and full-depth input has one extensive and one local event. Genome fraction remains 98.110%–98.669%. RSS rises from 1,381.6640625/1,428.65625 MiB at 30×/50× to 3,787.87109375/4,072.14453125 MiB at 100×/150×, then up to 7,858.8125 MiB over full-depth replicates. The responsible stage is not localized.

All 16 gradient runs complete in **35.50–167.36 s**, peak RSS **660.1015625–4,072.14453125 MiB**. Full-depth runs must not be included in that ceiling: maximum fastBact time is **185.91 s** and maximum RSS is **7,858.8125 MiB**. More depth does not consistently improve contiguity or accuracy. Biological GC, read length, phage abundance and lineage differ simultaneously, so this design characterizes observed taxa rather than estimating an isolated GC effect.

## 5. Stage G Determinism Certification Audit

### 5.1 Empirical bitwise and metric equality

Fresh SHA-256 calculations confirm the following complete FASTA bytes are identical across replicates 1, 2 and 3 for each taxon. Every one of the 13 QUAST metric fields compares exactly equal within each triplet.

| Dataset | Replicates | Common contig FASTA SHA-256 | Stored variance_max |
| --- | --- | --- | --- |
| CS-BS01 | 1, 2, 3 | `2bb0ecdd10f09b2320f4ba557b15a79bc5069f02f89352eceb9ec053b43a8861` | 2.0194839173657902e-28 |
| CS-PA01 | 1, 2, 3 | `b15ebd2289fb7f7ca39a49698f9049aa281d83a00982199f69da0722a068a9ff` | 0.0 |
| CS-SA01 | 1, 2, 3 | `86081b422341266a5a32feef424080bc0595f8ce599d996817b3e7d6c26a1962` | 0.0 |
| CS-SE01 | 1, 2, 3 | `4b72cfaa1b4fa8406be3f9d29c790047a8cae2ce2cb996cea9c9f31859314bfb` | 0.0 |

The exact population variance of equal decimal metric values is **σ²=0** on all four taxa. Literal serialized `variance_max` is **2.0194839173657902e-28** for BS01, not zero: binary floating-point averaging of three copies of 98.576 gives 98.57599999999998. This numerical residue does not indicate output variation. The other three store `0.0`. Runtime, CPU and RSS vary and are outside the determinism metric set.

### 5.2 F26-H01: formal certification fails the frozen identity rule

Benchmark contract §10.7 requires the same dataset/tool, read paths and hashes, reference path/hash, binary hash, **Git HEAD**, threads and memory to pool replicates. Direct inspection finds:

1. All 40 `rep_*/metrics.json` files lack identity fields such as `git_head`, read hashes and binary hash. The collator reads these records. Every full-depth collated `determinism.identity` therefore contains empty strings for all ten non-dataset/tool identity fields. Empty unknown identities are being treated as equal, yielding `pass=true` without enforcing the required evidence.
2. The original per-job `raw/{dataset}_fastbact_r{n}/summary.json` retains the identities. In every species, replicate 1 records **`968b06f177c90ccf47dd3ffa42b55a83b1e9aaec`**, whereas replicates 2 and 3 record **`c23a6f96e00ece7bdc7389ed5a756700afb06238`**. All other listed identity fields agree within each taxon, including the binary hash. Under the frozen grouping rule these are groups of **one and two**, neither N=3.
3. Source inspection explains the missing identity propagation: `run_tool` writes `metrics.json` before the main execution path calls `stamp_run_identity`; the enriched records are subsequently written to the per-job summary, not back to those metrics files. `collate_benchmark_results.py` consumes the earlier files. The identity tokenizer maps absent fields to empty strings. No protected artifact was regenerated or repaired in this audit.

Thus **empirical bitwise reproducibility is verified; formal Stage G contract certification is withheld**. Same binary and unchanged `src/` are strong scientific supporting evidence but do not silently waive an explicit Git-HEAD requirement. Repair the evidence producer, reject incomplete identities, preserve historical records, and acquire a third full-depth replicate at the same complete identity as replicates 2–3 (or a new common three-run identity). Any alternative pooling criterion requires an authorized contract decision. Competitor and gradient runs have N=1 and receive no determinism certification.

## 6. Stage F Robustness Hardening & Invariant Assessment

The engineering decision to refrain from speculative algorithm edits is appropriate under F26-P4. All 28 fastBact runs and their QUAST evaluations exit zero, and the source diff is empty. No run crashed, timed out or exhausted its allocated memory. This establishes observed completion on this panel, not error-free assemblies: QUAST reports residual extensive/local discrepancies as tabulated.

Static inspection of `src/resolve/exspander.rs` confirms the inherited safeguards remain present. `score_evidence_rows` and its scoring counterpart reject starts `p <= b + H` when the entry-shadow condition applies; arithmetic uses widened values. Traversal retains insert-distance budgets, short-repeat re-entry probing, long-node revisit pruning, circular-node withholding and a bounded search-state cap. This is preservation of the frozen conditional safeguards, not a new universal mathematical proof or exhaustive re-audit of every algorithm path.

The fresh release suite includes the 42-test state-semantics gate and existing provenance/attribution regressions. No F26 dynamic per-decision invariant trace is supplied. Therefore the requested inference that 28 successful runs **confirm the mathematical invariants across all regimes** is too strong and is not certified.

The reproducible SA01 two-event excess warrants forensic work even though the process succeeds. No responsible stage or minimal counterexample has been supplied. Likewise LT2 memory scaling and low-depth MW2 fragmentation are observations that require localization before a repair. Stage F may properly close as **no algorithm modification justified by currently localized evidence**, with these scientific limitations retained; success exit codes must not close the observations as resolved defects. The evidence-producer defect F26-H01 is a separate engineering finding for the implementation owner.

## 7. Code and Test Regression Audit

The exact required command was run:

```sh
git diff 02872f81811402d2bd603ff7c0a2ecf8b3dbe528 -- src/
```

It returns exit 0 and **zero output**. The fresh test job additionally uses `git diff --exit-code` for the same comparison. Audited `HEAD:src` tree object is **`fe98d216f11b2b579cd358de2418773d9e4ff683`**. All 28 retained fastBact job summaries identify binary SHA-256 **`517402653fcefc625cb79912ade7e7e382b0d05343bbeb483d7c82dd807947d1`**; the current release executable has this digest before and after testing.

**Fresh independent regression evidence:** Slurm job **3076232**, partition **qfree**, node **vnode1**, eight CPUs, 12G, requested 20-minute limit, command **`cargo test --release -- --nocapture`**. Scheduler accounting records COMPLETED and exit **`0:0`**. Summing the actual target results gives **213 passed, 0 failed, 0 ignored**. The log records audited HEAD, source tree, command output, binary hashes and `F26_SOL_RELEASE_TEST_EXIT=0`:

- Log: `logs/f26_sol_release_audit_3076232.log`
- Log SHA-256: `9eb653edca0e843f1f37f5870a7dab1bfdecea54d40151b63a02c3df618386a7`

Tests ran only on the Slurm compute node. No real-data assembly or QUAST benchmark was rerun in this audit; retained benchmark products were read and hashed. The release suite exercised its existing synthetic/integration scenarios. No new Rust source, tests or algorithm parameters were written.

**F26-H02 remains open:** contract §2 freezes full-depth R1-001, R1-003, R1-001-Val and R1-003-Val metric deltas. The supplied 40 runs contain none of those datasets. Generating 16 historical depth FASTQ pairs is not assembly verification and is not a full-depth oracle replay. The unchanged source/binary evidence supports non-regression but cannot substitute for observed paired oracle deltas under the stated input/environment conditions. Historical F25 results remain frozen; they are not newly measured F26 results. No waiver is issued here.

## 8. Residual Limitations, Risks & Recommendations for Phase F27

Acceptance closure work belongs to **F26**, rather than being relabeled as optional F27 work:

1. Close F26-H01 with identity-complete collation and valid N=3 groups; preserve all original run evidence and distinguish missing identity from equal known identity.
2. Supply F26-H02 full-depth historical-oracle replay and metric-by-metric comparison, or obtain explicit authorized adjudication of the missing gate. Do not adjust the frozen values or infer zero deltas from source identity alone.
3. Retain raw chromosome-reference QUAST results and qualify the MW2 depth claim. Under the Stage B/C binding interpretation, any viral exemption requires per-contig sequence attribution, with partial/chimeric contigs investigated and unresolved sequence left unresolved. The Stage D/G reports do not supply that attribution. This audit grants no exemption and makes no claim that all nonchromosomal contigs are Staph1N.

Recommended F27 investigations after acceptance closure are: stage-localized forensics of the SA01 extensive events and persistent PA01/SE01 discrepancies; chromosome-mapped base-depth calibration and targeted AT-rich tests using cleaner independent libraries; profiling LT2 memory allocation by assembly stage; and conservative contiguity recovery only with minimal failing examples and invariant-preserving tests. The current contiguity gap, especially against Unicycler on BS01 and SPAdes on PA01, should remain explicit in product claims.

One SRA library per species, one seed per depth and N=1 competitor runtimes do not estimate biological variability, sampling variance or robust performance distributions. The GC range is broad but not a factorial GC experiment. Chromosome-only references omit plasmids and phage. Evolved lineages limit interpretation of error rates as assembler defects. A zero-error-output or universal determinism claim across versions, platforms, thread counts or coverage levels is not supported.

The older F23 blocked text in AGENTS.md conflicts with later accepted historical records; HANDOFF already records that discrepancy. This F26 audit does not adjudicate C0-02R or change that scientific status. HANDOFF/ACTIVE_TASK still lead with F25 and must not be read as an F26 acceptance declaration. The Controller should carry these audit findings into the existing state records; this assignment writes only this report plus the fresh verification log.

## 9. Final Scientific Acceptance Determination

**Official verdict: `PARTIALLY_VALIDATED / F26_ACCEPTANCE_WITHHELD`.**

| Audit area | Determination |
| --- | --- |
| Dataset diversity, archived accession provenance, present reference/input integrity | Verified, with lineage/co-culture limitations |
| Stage D 16-run comparative measurements | Verified; speed, memory and contiguity tradeoffs qualified |
| Stage E 16-run library-fraction response | Verified within nominal-depth scope; not true 150× PAO1 or a causal GC threshold |
| Stage F preservation / no speculative algorithm changes | Supported; no universal invariant proof claimed |
| Stage G bitwise assembly and exact metric equality | Observed on all four triplets |
| Stage G frozen same-identity N=3 certification | **Not satisfied: F26-H01** |
| Zero `src/` diff and 213/213 release tests | Independently verified |
| Frozen F25 full-depth oracle replay / zero metric deltas | **Not demonstrated: F26-H02** |
| Whole-phase `F26_SCIENTIFICALLY_VALIDATED / READY_FOR_USER_ACCEPTANCE` | **Withheld pending closure or authorized contract adjudication** |

The requested positive verdict is retained here as an **unachieved target**, not granted merely because it appears in the report brief. Empirical results are not invalidated by the acceptance gaps, and the gaps must not be concealed to obtain the target label.

**Local main is NOT merged with F26 and the F26 branch is NOT pushed.** Read-only `git show-ref`/`git worktree list` show local `main` at **`02872f81811402d2bd603ff7c0a2ecf8b3dbe528`** while the F26 branch remains at the audited HEAD. Fresh `git ls-remote --heads origin` also returns that baseline for remote main and lists no F26 branch. This audit creates no Git commit, merge, push, tag or release. Human authorization remains required for scientific-contract changes and any later merge, push or release; no such operation is authorized by this report.

### Audit evidence fingerprint

Primary raw evidence is under `benchmarks/f26_results/raw/`, including each run's `metrics.json`, enriched per-job `summary.json`, `contigs.fa`, `quast/report.tsv`, command/version files and GNU Time logs. All 40 source records match the collated values and all TSV cells. The following digests freeze the inspected contract/manifest/summary snapshots:

| File | SHA-256 |
| --- | --- |
| F26_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md | `f2913504f8d0c533a7012b9ae22d5b69842f5522ca247a4ce32be3e40e4d1f0d` |
| F26_COMPETITOR_BENCHMARK_CONTRACT.md | `bfe46ff5e4caacd06399fc0fcabd65627b54a45271dc8be25740b121407be80a` |
| F26_DATASET_MANIFEST.tsv | `e23dd244b42d7d565b46b68483c2b68ab1cedfd024b7f28228b3f8e31ac70dc9` |
| F26_STAGE_E_MANIFEST.tsv | `fa9aace0e360e7eb1581611294631aaf86b2be1c7b9b3317c365662cae305d69` |
| benchmarks/f26_results/summary.json | `c00ba3c41697bc4d13cc435c3390fa7769ea6be64ed08bd471369b97388a0f9c` |
| benchmarks/f26_results/summary.tsv | `c92935ac83bdcdc5900cf5bb0b717d4db46c4944e77900c5abb4c1b7965b5a1e` |

Reproduction of the audit does not require overwriting benchmark evidence: hash the retained FASTAs and resolved inputs; compare all 13 metric fields against each QUAST TSV and flat summary; inspect enriched per-job identities separately from the un-enriched metrics files; query recorded job IDs with `sacct`; and run release tests only on `qfree`. This report's verification is independent of the worker reports' PASS labels.
