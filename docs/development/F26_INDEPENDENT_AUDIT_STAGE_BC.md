# F26 independent audit of Stages B and C

- Task: `F26-AUDIT-001`.
- Role: independent auditor, Codex Sol assignment; no implementation changes.
- Audit date: 2026-10-30 (UTC).
- Audited source: `075ecb91f26fa102c08c0bd5d0dac894dda384e5`.
- Branch: `feature/f26-cross-species-robustness`.
- Worktree: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f26-robustness-worktree`.
- Verdict: **STAGE_BC_AUDIT_BLOCKED**.

## 1. Findings, in severity order

### BC-01 / P1 / blocking: CS-SA01 reference-strain identity is not conclusively verified

The fresh NCBI SRA response for `SRR33682404` identifies sample `SRS25134113` / BioSample `SAMN48627418`. Its sample alias is `MW2-cotreat2`, its run alias begins `MW2_A_Staph1N`, and its strain attribute is `MW2`; however, its sample TITLE is **MRSA252-cotreat2**. A separate BioSample fetch independently reproduces that title/strain conflict. The selected reference `NC_003923.1` is MW2. This is contradictory submitter metadata, not proof that the reads are MRSA252, and not something an auditor may silently resolve in favor of either field.

`F26_DATASET_PROVENANCE_AND_SELECTION.md:115` reports the supporting alias/strain attributes but omits the contradictory title; its unconditional strain-correspondence PASS is therefore too strong. Its documented mixed-culture, phage and oxacillin context further limits a claim of an exact isolate-specific reference. GC diversity and matching reference-file hashes cannot settle read identity.

**Required closure:** resolve the title discrepancy with authoritative sample documentation/submitter clarification or a discriminatory sequence-based identity analysis with retained evidence. If unresolved, have the responsible scientific authority decide replacement or explicitly qualified use under the frozen contract. Do not change the reference or dataset based on this audit alone. Other derived strains can be described as reference-strain-lineage comparators, not exact finished genomes of the sequenced evolved/engineered isolates.

Evidence: newly retrieved `sra.xml` and `biosample_sa.xml` in the evidence directory in §7; [NCBI BioSample SAMN48627418](https://www.ncbi.nlm.nih.gov/biosample/SAMN48627418), [SRA SRR33682404](https://www.ncbi.nlm.nih.gov/sra/SRR33682404). F26-P1 explicitly requires a complete reference for the sequenced isolate or verified reference strain.

### BC-02 / P1 / blocking: advertised benchmark wall time conflicts with qfree

`scripts/slurm_benchmark_harness.sbatch:18` requests `36:00:00`; benchmark contract §2 repeats 36 hours and sets an 8-hour per-assembler timeout. Fresh `scontrol show partition qfree` reports `MaxTime=03:00:00`, so the checked-in default exceeds the advertised partition limit by 12-fold. Four assemblers plus QUAST cannot in general consume their advertised individual timeout budgets within three hours.

Important qualification: `scontrol show config` reports `EnforcePartLimits=NO`, and an actual `sbatch --test-only scripts/slurm_benchmark_harness.sbatch --dataset CS-SA01 --tool all --outdir /tmp/f26-audit-unused` returned exit 0 with a start estimate. **Submission rejection was not observed**, and no 36-hour benchmark was submitted. This finding is a configuration/cluster-policy inconsistency, not a fabricated failed submission or measured benchmark timeout.

**Required closure:** align the wrapper and contract with the authorized three-hour qfree limit, choose appropriate per-run timeouts and job granularity, or document an authorized cluster exception backed by administrator/scheduler evidence. Verify the resulting wrapper on qfree. CPU and memory defaults themselves pass (§5).

### BC-03 / P2 / blocking Stage C contract compliance: two replicates can receive a determinism PASS

Benchmark contract §10 requires **three** replicates for its F26 variance claim. `determinism_report()` at `scripts/run_competitor_benchmark.py:810` declines checking only when there are fewer than two records; two identical PASS records receive `checked=true`, `pass=true`, and zero numeric variance. The shipped self-test explicitly supplies two records and asserts this PASS; the independent run reproduced `SELF_TEST_PASS`.

**Required closure:** distinguish a two-replicate comparison from the contract's three-replicate acceptance gate, and require the prescribed replicate count for the latter. Prevent pooling different experiment identities: current grouping uses only dataset/tool, not read hashes, reference hash, source revision or binary hash. This is a pre-benchmark harness finding; no real F26 determinism result has been rejected or certified here.

### BC-04 / P2 / advisory: release binary exists, but frozen-baseline identity remains unproven

On vnode2, `target/release/fastBact --version` returned `fastBact 0.1.0`, exit 0. Its SHA-256 is `517402653fcefc625cb79912ade7e7e382b0d05343bbeb483d7c82dd807947d1`. HANDOFF's F25 post-integration release hash is `7c70e3572996ec777b62b8be0080cd8a0a1528ab01f03d0d4b9f3f56640a8112`. These bytes differ; version text and the current checkout SHA alone cannot establish that this executable is the frozen F25 comparator. A rebuild can legitimately differ, so this is **not** evidence of an algorithm change. `git diff 02872f81811402d2bd603ff7c0a2ecf8b3dbe528 HEAD -- src` is empty.

Before Stage D, retain a build/source/toolchain binding or use the verified frozen artifact. The harness records the fastBact binary hash in a per-run file, but accepts any existing explicit/default release binary without an expected-hash gate. This audit did not rebuild or replace it.

### BC-05 / P2 / advisory: download/validator safeguards have bounded coverage

The download script performs full gzip checks, equal nonzero mate counts, retries and persistent checksums; the standalone validator additionally checks every normalized header ID. Current files pass these checks. Remaining static concerns are:

- `scripts/download_f26_datasets.sh:150` sends the uncompressed singleton `.fastq` to `count_lines()`, which invokes `pigz -dc`; gzip tools may reject or pass through such input with an error status. The successful historical job does not prove correct error propagation.
- `process_accession` is invoked under `if ! process_accession ...`; Bash disables ordinary `errexit` propagation within such a function. Several compression/count/hash operations have no explicit failure return. Do not rely on `set -e` alone for those failures.
- Resume checks return before checksum recording, so an existing intact pair need not gain a missing/stale ledger entry on resume.
- The validator checks four-line framing and header markers, not plus-line syntax or sequence/quality length equality. Default discovery can omit an entirely missing accession; explicit expected accessions should be provided for a completeness gate. Manifest spot enforcement is deliberately optional and must not confuse SRA spots with paired reads when singletons exist.

These are implementation-hardening observations from source inspection, not reproduced corruption of this dataset. Repairs belong to the implementation owner; protected and non-owned files were not changed.

## 2. Dataset provenance and selection

Fresh NCBI EFetch retrieved exactly four requested SRA experiment packages, all public and load-complete. Each has WGS/GENOMIC library metadata, a PAIRED layout, and an Illumina instrument. All four BioProject/BioSample mappings and total spot counts match the manifest. Nominal read lengths are labels; total bases show variable/truncated reads and do not prove every read has the nominal length.

| Dataset / run | BioProject / BioSample | Declared strain; sequencing context | Instrument; nominal PE | Reference; measured GC; Gram |
|---|---|---|---|---|
| CS-SA01 / SRR33682404 | PRJNA1263016 / SAMN48627418 | MW2 attribute, conflicting MRSA252 title; phage/oxacillin mixed culture | NovaSeq X Plus; 2×150 | NC_003923.1; 32.83%; positive |
| CS-BS01 / SRR40739147 | PRJNA1530750 / SAMN63192558 | 168; chloramphenicol laboratory-evolved culture | NextSeq 550; 2×151 | NC_000964.3; 43.51%; positive |
| CS-PA01 / SRR41244734 | PRJNA1545662 / SAMN64277765 | PAO1; PAO1-GFP isolate, mexR/aroB mutant, trimmed submission | NovaSeq 6000; 2×240 | NC_002516.2; 66.56%; negative |
| CS-SE01 / SRR41010285 | PRJNA1539015 / SAMN63983458 | LT2; antimicrobial-resistance ALE, BAM submission, RANDOM PCR | NovaSeq 6000; 2×150 | NC_003197.2; 52.22%; negative |

Four non-E. coli species, four genera, two Gram-positive and two Gram-negative taxa, and reference GC from 32.83% to 66.56% provide a scientifically sensible diversity panel. Gram assignments are conventional species-level phenotypes; they were not measured from these FASTQs. The panel does not isolate GC as a causal factor because genome size, library preparation, read length and treatment also vary.

NCBI nuccore summaries independently confirm all six unique versioned references' strain names, complete-genome/sequence titles and lengths. These are valid named-strain chromosomes. They do not establish exact sequence identity of derived isolates, and BC-01 prevents unconditional verification of all four run/reference links. Chromosome-only evaluation omits plasmids and phage content; the provenance report appropriately discloses these confounds. Its numerical estimates of mutation/cassette/phage sizes are expectations, not measurements certified by this audit. Follow the report's later forensic obligations before attributing QUAST differences to assembler defects.

## 3. Reference physical verification

Shared data root is `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data`. Manifest `data/` prefixes are stripped when resolving under that root, as documented and implemented by the harness.

**All eight manifest rows resolve to existing FASTAs; there are six distinct physical reference files, not eight unique files.** R1-001 and R1-001-Val share one file; R1-003 and R1-003-Val share another. Independently recomputed raw and normalized-sequence SHA-256 values match every row (8/8 raw, 8/8 normalized). Every file has one FASTA record, zero N bases, and manifest-matching length and rounded GC.

| Dataset rows | Path below shared data root | Length | Raw SHA-256 (matched) |
|---|---|---:|---|
| R1-001, R1-001-Val | references/m6_identity/NC_000913.3.fasta | 4,641,652 | `6b195feda4c66140f6762742eb8b30c2652f02b45878b174f5b00ef85ecc95d7` |
| R1-003, R1-003-Val | references/m6_identity/CP009273.1.fasta | 4,631,469 | `37eade01311f1f81659fad53bf907a93139075a75d115fd3aec2ef09beac2e90` |
| CS-SA01 | f26_cross_species/references/NC_003923.1.fasta | 2,820,462 | `19ba62a339065c9c498efb72bb2491e372344e4f23a4ac1e3e34e7345a25c63c` |
| CS-BS01 | f26_cross_species/references/NC_000964.3.fasta | 4,215,606 | `a334e891ffc0e307f23f48842775d3383177a9d9cb5d5075b552a2cccddfe139` |
| CS-PA01 | f26_cross_species/references/NC_002516.2.fasta | 6,264,404 | `902fa552c8c25aaaf82eccfed000e359d6a04c8744d9ec0a8134699cb6ba53e6` |
| CS-SE01 | f26_cross_species/references/NC_003197.2.fasta | 4,857,450 | `59fceaea5ca57836e02b6c4ad3582847623c87268f7e5511874909eedb4acafe` |

Normalized hashing concatenated uppercase sequence lines with headers/newlines excluded. NCBI summary lengths matched independently. A fresh bulk NCBI FASTA fetch timed out after 60 seconds with a partial response; it was not used as a complete sequence-comparison oracle. Full online sequence equivalence beyond the frozen local manifest hashes is not claimed.

## 4. Physical FASTQ verification

The supplied `validation_report.json` was inspected and then independently rerun against current files on qfree, without overwriting the supplied report. Fresh results agree with it. Every mate exists and is readable; full-stream decompression succeeds; all line counts divide by four; mate counts are equal and nonzero; all compared normalized read IDs match; header-malformation counters are zero.

| Accession | R1 lines | R2 lines | Reads per mate / compared pairs | ID mismatches | gzip / framing / parity |
|---|---:|---:|---:|---:|---|
| SRR33682404 | 8,958,612 | 8,958,612 | 2,239,653 | 0 | PASS |
| SRR40739147 | 9,015,228 | 9,015,228 | 2,253,807 | 0 | PASS |
| SRR41244734 | 6,073,612 | 6,073,612 | 1,518,403 | 0 | PASS |
| SRR41010285 | 13,412,932 | 13,412,932 | 3,353,233 | 0 | PASS |

Every one of the **nine** entries in `docs/development/F26_FASTQ_CHECKSUMS.tsv` was independently hashed from disk; byte sizes and SHA-256 digests all match. That includes eight mates and the Salmonella singleton gzip. A separate `pigz -t` over all nine gzip files returned 0. The immutable ledger supplies the exact digests; the fresh verification log records all observed values and `match True` for each.

For SRR41010285, the manifest's 3,359,296 SRA spots differ from 3,353,233 paired reads by **6,063**. The singleton file contains 24,252 lines / 6,063 records, exactly reconciling that difference under split-3 extraction. This is not observed data loss or mate desynchronization. The harness selects only the paired files for all four assemblers, consistently excluding those singletons. Its paired input therefore must not be described as all SRA reads or assigned depth from the unadjusted SRA base total.

The provenance report §6.1 requests a post-download duplicate assessment before Stage C. A fresh exact paired-sequence duplicate count on vnode3 found 3,007,915 distinct ordered sequence pairs among 3,353,233 pairs: 345,318 excess duplicate pairs, **10.298062%**. This is below the report's 20% flag threshold under the stated `(total - distinct)/total` definition. Reads in duplicate groups total 647,768 pairs (19.3177% approximately). Sequence identity is a descriptive statistic, not proof of PCR origin or a mapping-based duplicate estimate. No deduplication or input transformation was performed.

## 5. Competitor harness and resource audit

| Check | Observed result |
|---|---|
| qfree partition | Wrapper explicitly requests qfree; independent verification jobs used qfree |
| CPU | One task, one node, eight requested CPUs; assembler and QUAST flags use eight; BLAS/MKL use one |
| Memory | Wrapper requests 14G; SPAdes/SKESA internal cap 14; Unicycler forwards `-m 14` to SPAdes |
| Physical ceiling | `scontrol show node vnode2`: RealMemory=15500; AllocMem=14336 for 14G. Request is below reported capacity |
| Wall time | FAIL: 36-hour wrapper/contract versus qfree MaxTime three hours; see BC-02 and its enforcement qualification |
| SPAdes | Installed version v4.2.0, exit 0; `--isolate -1 R1 -2 R2 -t 8 -m 14`; contigs.fasta |
| SKESA | Installed version 2.5.1, exit 0; single comma-joined `--reads R1,R2`, `--cores 8 --memory 14 --min_contig 200`; contigs.fa |
| Unicycler | Installed version v0.5.1, exit 0; paired short reads, normal mode, min_fasta_length 200, pinned SPAdes path and `--spades_options '-m 14'`; assembly.fasta |
| fastBact | Release-path executable present; version 0.1.0, exit 0; `-t 8 --min-len 500 --kmers auto --depth 120 --dump-stages ...`; no auto-subsampling flag; see BC-04 |
| QUAST | Installed version v5.2.0, exit 0; correct strain-reference path, `-t 8 --min-contig 200`, explicit contig thresholds, no --fast/--no-snps |
| Dry-run | Eight dataset invocations × four tools = 32 DRY_RUN records; all exit 0, no missing reads/references or reference-hash mismatches |
| Harness self-test | SELF_TEST_PASS, exit 0; parsing, command construction, timed subprocess success/failure/timeout and hash-based replicate checks exercised |
| Shell syntax | `bash -n scripts/download_f26_datasets.sh scripts/slurm_benchmark_harness.sbatch`: exit 0 |

The Python harness executes assemblers sequentially and uses GNU time separately from QUAST. It stores stdout/stderr, argv, exit codes, assembly hash, metrics, wall/CPU/RSS data, and reference/input paths. It checks reference hashes and executable versions for real runs; dry-run does not execute or verify binaries, so the separate compute-node version checks above are essential. The version checks are substring-based, not cryptographic software provenance checks.

The Slurm wrapper refuses invocation without SLURM_JOB_ID; the Python CLI itself has no compute-node guard. Resource flags can be overridden and merely produce a deviation note, so the approved wrapper/defaults and external Slurm allocation remain necessary. Actual cgroup enforcement was not independently stress-tested; node allocation alone is not proof of its behavior under OOM.

The fastBact 500-bp output cutoff differs from competitors' 200-bp outputs. This inherited F25 asymmetry is explicitly disclosed in the benchmark contract and reproduced correctly in the harness; interpret small-contig/count comparisons accordingly. No assembler pipelines or real QUAST benchmark executions were required or performed in this preparatory audit. A dry-run does not establish end-to-end assembler compatibility, benchmark quality, runtime, OOM handling, or scientific robustness.

## 6. Gate disposition and required next action

| Gate | Disposition |
|---|---|
| Species/Gram/reference-GC diversity | PASS, with experimental confounds disclosed |
| Four public Illumina paired-end WGS accessions and manifest identifiers | PASS |
| All four isolate/reference assignments conclusively verified | BLOCKED by BC-01 |
| Eight reference-row checks, six distinct files | PASS |
| Four physical FASTQ pairs, gzip, framing, count parity, full ID sync | PASS |
| Nine recorded FASTQ byte-size and SHA-256 entries | PASS |
| Salmonella singleton reconciliation and duplicate screening | PASS within the stated measurement scope |
| Installed tool versions and command construction | PASS; frozen fastBact binary binding remains advisory |
| All-manifest dry-run without missing inputs | PASS, 32/32 |
| Benchmark cluster wall-time consistency | BLOCKED by BC-02 |
| Three-replicate acceptance semantics | BLOCKED by BC-03 |

Return BC-01 to the provenance/scientific owner and BC-02/BC-03 to the harness owner; independently recheck their evidence before advancing Stage B/C acceptance. Preserve all current data, manifests and historical artifacts. This audit changes no scientific contract, grants no Stage D/G acceptance, and leaves the repository's existing F23 status/history discrepancy outside scope. HANDOFF and ACTIVE_TASK were read; only this explicitly owned audit report is modified.

**Definitive verdict: STAGE_BC_AUDIT_BLOCKED.** Passing storage-integrity checks and command previews do not close the three findings above.

## 7. Reproducibility and evidence inventory

Evidence directory (outside the Git worktree):

`/hpcfs/fhome/caizhh/.local/state/f26-audit-001/`

- `verify.sh`, submitted with `sbatch`: job **3073629**, qfree/vnode2, eight CPUs, 14G, 20-minute requested wall limit; COMPLETED, Slurm exit **0:0**, elapsed 1m29s. Log `verify-3073629.log` contains independent FASTA/hash derivations, nine FASTQ hashes, explicit command exit codes, dry-run exits and installed versions. Each inner check's result was inspected; batch completion alone was not used as proof.
- `independent_validation.json`: fresh complete four-accession FASTQ validation, validator exit **0**; original supplied JSON preserved.
- `harness-self-test.log`: all checks passed, `SELF_TEST_PASS`; `dry-{dataset}.log` and `dry/{dataset}/summary.json` / `.tsv`: all eight dataset previews, four tools each.
- `duplicate.sh`: exact ordered paired-sequence count using streamed mate sequence lines, paste, `LC_ALL=C sort -S 2G --parallel=4`, `uniq -c` and count aggregation. Job **3073637**, qfree/vnode3, eight CPUs, 14G, ten-minute requested limit; COMPLETED, exit **0:0**, elapsed 12s; `duplicate-3073637.log` also records singleton line count.
- `sra.xml`: fresh NCBI SRA EFetch; targeted four accessions, one response, four packages, no record filters or identifier conversion. `biosample_sa.xml`: separate targeted BioSample fetch; one record. `ncbi_reference_summary.json`: nuccore ESummary; six requested versioned accessions, six returned records. Access date 2026-10-30; unauthenticated GET, bounded calls, no pagination needed.

Reproducible database endpoints and parameters:

```text
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
  db=sra
  id=SRR33682404,SRR40739147,SRR41244734,SRR41010285
  rettype=full&retmode=xml&tool=fastbact_f26_audit
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
  db=biosample&id=SAMN48627418&retmode=xml&tool=fastbact_f26_audit
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi
  db=nuccore
  id=NC_000913.3,CP009273.1,NC_003923.1,NC_000964.3,NC_002516.2,NC_003197.2
  retmode=json&tool=fastbact_f26_audit
```

Local checks can be repeated using the retained `verify.sh` and `duplicate.sh` on qfree with new output paths to preserve audit evidence. Scheduler facts were retrieved with `scontrol show partition qfree`, `scontrol show config`, `scontrol show node vnode2`, and `sacct -j 3073629,3073637 --format=JobID,Partition,State,ExitCode,Elapsed,AllocCPUS,ReqMem,NodeList -P`. The original 36-hour wrapper was tested only with `sbatch --test-only`, not submitted for execution.

Evidence SHA-256:

| Artifact | SHA-256 |
|---|---|
| sra.xml | `f0be679d16bf5afc7a8b3d0c598b492f3c5cc4f80b8358c4e81a33dea9871077` |
| biosample_sa.xml | `9d3d157a78a5fa106468e7e2677484fc9e1ac0658ed9a9fed5ee2ef087f2d334` |
| ncbi_reference_summary.json | `807f530c241f5e51440a6a460e9f5c413075fe938fcc392490a14fb120b0944c` |
| verify-3073629.log | `9a4baf3f9a9cda9347d18e29d5755ce6f317afab3cb8f67ccb4e71b4fc78d6fa` |
| independent_validation.json | `3dc5d4e355ecf1a5658b18e062964ea1e518e323c20a9bbee181ce2f5bc5a8c7` |
| duplicate-3073637.log | `697f0cc6cda176c2b15000e21bea17569c2bb0da6d96b76f8c486f356a9fdef3` |

Audited input SHA-256:

| Repository/input artifact | SHA-256 |
|---|---|
| F26_DATASET_MANIFEST.tsv | `e23dd244b42d7d565b46b68483c2b68ab1cedfd024b7f28228b3f8e31ac70dc9` |
| F26_DATASET_PROVENANCE_AND_SELECTION.md | `9e8a167520615da678e0c89316a4b2c5f2ae79b77baf9a33a46dcb7320cd2a0d` |
| F26_COMPETITOR_BENCHMARK_CONTRACT.md | `c7147e9075354e8b98673e719030bab13d7857d4d77fe7a5445edec4858aa7a5` |
| scripts/download_f26_datasets.sh | `c4476637695a423e819ae89e09a806d5069c85207755c7697d27ddfd2e26b581` |
| scripts/validate_fastq_pairs.py | `bb85cd25d3ce2819e4dd13b3c2a68b69ca3884dceb1964827a8213d5c7cceaeb` |
| scripts/run_competitor_benchmark.py | `bfa4d807f1d5d5b759b8b684c28512dbfd54f01758c06759b773fb70953dfddd` |
| scripts/slurm_benchmark_harness.sbatch | `51cc1af42f65b62d30d14353ebe03878a3a7552eb70ce79de7db8364dcb272e5` |
| docs/development/F26_FASTQ_CHECKSUMS.tsv | `51b0aa1417e7a192c78b5c104229ca5bb8fecd45a6033233eff8d2b6e4d87afc` |
| supplied validation_report.json | `94a51c931bd189df304d398e926ed8e1b6317d862aeed30f980e1ff0dcbf05ac` |

No push, merge, release, source edit, manifest edit, reference edit, FASTQ edit or scientific-status acceptance was performed. This report is the sole repository deliverable.
