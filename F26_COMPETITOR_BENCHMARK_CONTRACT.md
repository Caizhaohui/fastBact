# fastBact Phase F26 — Unified Competitor Benchmark Contract

**Document ID:** `FASTBACT-F26-BM-CONTRACT-01`  
**Task:** `F26-BM-001`  
**Revision:** `F26-BM-002` (2026-10-30) — qfree wall clock is 03:00:00; one job is one dataset and one tool; determinism PASS requires three replicates of the same run identity.  
**Date:** 2026-10-30  
**Phase:** F26 — Cross-Species Assembly Benchmarking, Robustness Hardening & Competitive Validation  
**Harness:** `scripts/run_competitor_benchmark.py`  
**Slurm entry:** `scripts/slurm_benchmark_harness.sbatch`  
**Governing contract:** `F26_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md` (F26-P2)  
**Base commit:** `02872f81811402d2bd603ff7c0a2ecf8b3dbe528`

This document freezes the execution protocol for head-to-head short-read assemblies. It does not change assembler source, dataset bytes, or `F26_DATASET_MANIFEST.tsv`.

---

## 1. Scope

One invocation evaluates one manifest row with one or all of:

| Tool key | Software | Pinned binary |
|---|---|---|
| `fastbact` | fastBact (Cargo package `fastBact`, binary name `fastBact`) | `{worktree}/target/release/fastBact` |
| `spades` | SPAdes v4.2.0 | `/hpcfs/fhome/caizhh/.conda/envs/BactGenome/bin/spades.py` |
| `skesa` | SKESA v2.5.1 | `/hpcfs/fhome/caizhh/.conda/envs/LineaBact/bin/skesa` |
| `unicycler` | Unicycler v0.5.1 | `/hpcfs/fhome/caizhh/.conda/envs/BactGenome/bin/unicycler` |

Evaluation is QUAST v5.2.0 at `/hpcfs/fhome/caizhh/.local/bin/quast.py`.

Mode is Illumina paired-end only. No long reads, no `--meta`, no `--careful`, no plasmid-only mode, and no fastBact subsampling flag is passed.

---

## 2. Resource limits

| Resource | Limit | Enforcement |
|---|---|---|
| Threads | 8 | Tool flags plus `OMP_NUM_THREADS=8`, `RAYON_NUM_THREADS=8`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1` |
| Memory | 14 GB | Slurm `--mem=14G` on partition `qfree` (cgroup is the hard cap). SPAdes and SKESA also receive an internal 14 GB cap. |
| CPUs on node | 16 available, 8 requested | One task, one node (`vnode1`–`vnode4`) |
| Assembler wall timeout | 8100 s (2 h 15 min) default | One assembler inside one qfree job. Harness status `TIMEOUT`, exit code 124. Override with `--timeout-seconds`. |
| QUAST wall timeout | 1500 s (25 min) default | Separate from assembler time. Override with `--quast-timeout-seconds`. |
| Slurm wall clock | 03:00:00 (3 h) | `#SBATCH -t 03:00:00` on partition `qfree`. This matches the partition's advertised MaxTime. Do not request a longer wall clock. |
| Job granularity | one dataset and one tool | `--dataset` selects one manifest row. `--tool` selects one assembler. Default timeouts fit that pair once and leave 1200 s before the 10800 s qfree limit (8100 + 1500). A real run whose worst-case product `tools × replicates × (assembler timeout + QUAST timeout)` exceeds 10800 s exits 2 before any assembler starts. `--dry-run` prints the same warning and continues. N=3 determinism is three jobs (`--replicate 1`, `--replicate 2`, `--replicate 3`) that share `--outdir`. |

`ulimit -v` is not set. SPAdes reserves a large virtual address space, and a virtual-memory ulimit would abort it before the RSS cap is reached. Peak RSS is the measured memory quantity.

Tool-internal memory flags (`-m 14`, `--memory 14`, Unicycler `--spades_options "-m 14"`) are decimal-GB limits as defined by each tool. The Slurm cgroup remains the allocation that the job cannot exceed. An out-of-memory kill (exit 137) is recorded as `OOM_KILLED` and does not stop the remaining tools.

---

## 3. Inputs

### 3.1 Manifest

The harness reads `F26_DATASET_MANIFEST.tsv` (or `--manifest`). Required columns:

`dataset_id`, `sra_accession`, `bioproject`, `reference_accession`, `reference_fasta_path`, `reference_raw_sha256`

Optional columns `fastq_r1` and `fastq_r2`, when non-empty, override path discovery. This harness does not write the manifest.

### 3.2 Read and reference roots

Default data root:

`/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data`

Historical *E. coli* runs (`bioproject=PRJNA1088182`):

```text
{data_root}/prjna1088182/fastq/{SRA}/{SRA}_1.fastq.gz
{data_root}/prjna1088182/fastq/{SRA}/{SRA}_2.fastq.gz
```

Cross-species runs:

```text
{data_root}/f26_cross_species/fastq/{SRA}_1.fastq.gz
{data_root}/f26_cross_species/fastq/{SRA}_2.fastq.gz
```

Reference paths in the manifest are relative (`data/...`). The harness resolves them under the data root with the leading `data/` segment removed, then under the worktree. Before assembly it checks `reference_raw_sha256`. A mismatch yields `REFERENCE_HASH_MISMATCH` and that dataset is not assembled.

### 3.3 fastBact binary

Cargo `[[bin]] name` is `fastBact`, so the release path is `target/release/fastBact`. The harness also accepts `target/release/fastbact` and `--fastbact-bin`. The binary is produced on a `qfree` compute node by `cargo build --release`. This task does not build it. Each run records the binary SHA-256 and `git rev-parse HEAD`.

---

## 4. Exact assembler commands

Placeholders: `{R1}`, `{R2}`, `{TOOL_OUT}`, `{THREADS}=8`, `{MEM_GB}=14`.

The process environment always includes:

```text
OMP_NUM_THREADS=8
RAYON_NUM_THREADS=8
OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1
LC_ALL=C
TZ=UTC
```

`PATH` is prefixed with the BactGenome and LineaBact `bin` directories so Unicycler can find `spades.py`, `bowtie2`, `samtools`, BLAST, and Pilon without a conda activate.

### 4.1 fastBact

Matches the F25 frozen invocation in `run_f25_benchmark.sbatch`, with defaults written out so the command is auditable. Subsampling stays off: `--enable-auto-subsample` is not passed, so `--depth 120` does not discard reads.

```text
{FASTBACT} \
  -1 {R1} \
  -2 {R2} \
  -o {TOOL_OUT} \
  -t 8 \
  --min-len 500 \
  --kmers auto \
  --depth 120 \
  --dump-stages {TOOL_OUT}/stages
```

`--min-len 500` is the F25 default. QUAST still scores that FASTA at 200 bp. Contigs shorter than 500 bp are absent from the fastBact FASTA, so they cannot contribute to the 200 bp QUAST tables. SKESA, SPAdes, and Unicycler are not given a 500 bp writer cutoff. This asymmetry is inherited from the frozen F25 command and is part of the comparison, not a silent filter added by this harness.

Evaluated product: `{TOOL_OUT}/contigs.fa`.

### 4.2 SPAdes v4.2.0

Isolate Illumina mode. `--isolate` is the SPAdes-recommended mode for high-coverage bacterial isolates. It is not combined with `--meta` or `--careful`.

```text
{SPADES} --isolate \
  -1 {R1} \
  -2 {R2} \
  -o {TOOL_OUT} \
  -t 8 \
  -m 14
```

Evaluated product: `{TOOL_OUT}/contigs.fasta`. `scaffolds.fasta` is not scored.

### 4.3 SKESA v2.5.1

Paired reads are one comma-joined `--reads` argument. Two separate `--reads` flags would be treated as independent single-end inputs and are forbidden.

```text
{SKESA} \
  --reads {R1},{R2} \
  --cores 8 \
  --memory 14 \
  --min_contig 200 \
  --contigs_out {TOOL_OUT}/contigs.fa
```

`--min_contig 200` equals the SKESA 2.5.1 default and the QUAST floor. `--memory` is gigabytes and applies to the sorted k-mer counter.

Evaluated product: `{TOOL_OUT}/contigs.fa`.

### 4.4 Unicycler v0.5.1

Short-read only: `-l` is omitted, so the long-read bridge is not used. `--mode normal` is the Unicycler default and is passed explicitly. `--spades_path` pins the same SPAdes v4.2.0 binary. `--min_fasta_length 200` makes the emitted FASTA match the QUAST length floor.

```text
{UNICYCLER} \
  -1 {R1} \
  -2 {R2} \
  -o {TOOL_OUT} \
  -t 8 \
  --mode normal \
  --min_fasta_length 200 \
  --spades_path {SPADES} \
  --spades_options "-m 14"
```

Evaluated product: `{TOOL_OUT}/assembly.fasta`.

### 4.5 Canonical FASTA

The harness copies the evaluated product to `{run_dir}/contigs.fa` and hashes that copy. Preference order if several files exist: `contigs.fa`, `contigs.fasta`, `assembly.fasta`. Scaffold files are ignored.

---

## 5. Resource measurement

Assembler and QUAST are timed separately. Assembler wall time, CPU time, and peak RSS are the published resource numbers. QUAST time is stored beside them and is not added to assembler wall time.

Wrapper, without a shell:

```text
/usr/bin/time -v -o {time.txt} -- {argv...}
```

GNU Time 1.9 fields:

| GNU Time line | Recorded field |
|---|---|
| `User time (seconds)` + `System time (seconds)` | `cpu_time_s` |
| `Elapsed (wall clock) time` | `wall_clock_s` |
| `Maximum resident set size (kbytes)` | `peak_rss_kib`, and `peak_rss_mb = kib / 1024` |

Elapsed values `m:ss.ss`, `h:mm:ss`, and `h:mm:ss.ss` are accepted. The harness also records a Python monotonic clock. If the GNU Time file cannot be parsed, wall time falls back to that clock and peak RSS is null, with `resource_parse_ok=false`.

The timed command is started in its own session. On timeout the harness signals that process group, then SIGKILL after 30 s. Sibling tools still run.

---

## 6. QUAST protocol

QUAST runs only after a non-empty `contigs.fa` exists. Command:

```text
{QUAST} \
  -R {REFERENCE_FASTA} \
  -o {run_dir}/quast \
  -t 8 \
  --min-contig 200 \
  --contig-thresholds 0,200,1000,5000,10000,25000,50000 \
  --labels {dataset}-{tool}-rep{n} \
  {run_dir}/contigs.fa
```

`-R` and `-r` are the same QUAST 5.2.0 reference option. `-R` is used because the F25 oracle used `-R`.

Rules:

- Reference comparison is on. `--no-snps` and `--fast` are forbidden: mismatches and indels per 100 kbp must be computed.
- `--min-contig 200` is the scoring filter for every tool.
- `--contig-thresholds` adds an explicit `# contigs (>= 200 bp)` row. It does not change the set of contigs QUAST analyzes under `--min-contig 200`.
- Bacterial chromosomes in this manifest are below 10 Mb, so `--large` is not used.
- One assembly is passed per QUAST process. Tools are not co-loaded into one QUAST report.

---

## 7. Metric definitions

Values come from `{run_dir}/quast/report.tsv` (first value column). A cell of `-` is stored as JSON `null` and means QUAST left the statistic undefined. It is not coerced to zero.

| JSON / TSV field | QUAST row | Definition under this protocol |
|---|---|---|
| `contigs_ge_200bp` | `# contigs (>= 200 bp)`, else `# contigs` | Contigs at least 200 bp. With `--min-contig 200`, `# contigs` is that count. When both rows exist they must be equal; a disagreement is `METRICS_INCOMPLETE`. |
| `n50` | `N50` | Contig length such that half of the assembly length is in contigs of that length or longer. Computed on contigs ≥ 200 bp. |
| `nga50` | `NGA50` | N50 after breaking contigs at extensive misassemblies, relative to the reference length. Null when no aligned block covers half of the reference. Required by F26-P2. Null NGA50 alone is not a parse failure. |
| `l50` | `L50` | Number of contigs that comprise half of the assembly length. |
| `largest_contig` | `Largest contig` | Longest contig length in bp among contigs ≥ 200 bp. |
| `total_length` | `Total length` | Sum of contig lengths ≥ 200 bp. |
| `gc_pct` | `GC (%)` | GC content of those contigs. |
| `misassemblies` | `# misassemblies` | Extensive misassemblies (relocations, translocations, inversions). This is the F26-P2 raw extensive misassembly count. Local events are not included. |
| `local_misassemblies` | `# local misassemblies` | Local misassemblies, reported separately. |
| `mismatches_per_100kbp` | `# mismatches per 100 kbp` | Mismatch bases per 100 kbp of aligned sequence. |
| `indels_per_100kbp` | `# indels per 100 kbp` | Indel bases per 100 kbp of aligned sequence. |
| `genome_fraction_pct` | `Genome fraction (%)` | Percent of the reference covered by aligned contigs. |
| `duplication_ratio` | `Duplication ratio` | Aligned assembly length divided by the aligned reference length. |

Resource fields stored with each run: `wall_clock_s`, `cpu_time_s`, `peak_rss_mb`, `peak_rss_kib`, assembler `exit_code`, QUAST `quast_exit_code`, `assembly_sha256`.

A missing required row other than the NGA50 exception sets status `METRICS_INCOMPLETE`.

---

## 8. Output layout

```text
{outdir}/
  summary.json
  summary.tsv
  {dataset_id}/{tool}/rep_{n}/
    command.txt
    stdout.log
    stderr.log
    time.txt
    tool_out/...
    contigs.fa
    quast/report.tsv
    quast.stdout.log
    quast.stderr.log
    quast.time.txt
    metrics.json
```

`summary.json` schema version is `f26-benchmark-1`. Re-running into the same `--outdir` replaces rows with the same `dataset_id`, `tool`, and `replicate`, and keeps other rows. `--force` deletes only that replicate directory before rerunning. Without `--force`, an existing replicate directory is `EXISTS` and is not overwritten.

`summary.tsv` is the flat table of every retained run. `determinism` is always embedded in `summary.json`. A group can receive `pass=true` only when it has at least three replicates and those replicates share one run identity.

---

## 9. Failure behavior

Each tool is isolated. An exception, non-zero exit, timeout, missing input, hash mismatch, version mismatch, or QUAST failure is written into that run's `metrics.json` and into the summary. The harness then continues with the next tool. The process exits 1 if any requested run is not `PASS` or `DRY_RUN`. Usage errors (unknown dataset, bad manifest, missing `--outdir`) and a real schedule that exceeds the qfree 10800 s worst-case budget exit 2 before any assembler starts.

| Status | Meaning |
|---|---|
| `PASS` | Assembler exit 0, QUAST exit 0, required metrics present |
| `DRY_RUN` | Commands and paths resolved; nothing executed |
| `FAIL` | Non-zero assembler exit |
| `TIMEOUT` | Assembler or QUAST exceeded its wall limit (exit 124) |
| `OOM_KILLED` | Exit 137 |
| `MISSING_BINARY` | Pinned executable absent |
| `VERSION_MISMATCH` | Executable version is not the pinned version |
| `MISSING_INPUT` | Read pair absent |
| `MISSING_REFERENCE` | Reference FASTA absent |
| `REFERENCE_HASH_MISMATCH` | File SHA-256 ≠ manifest `reference_raw_sha256` |
| `MISSING_ASSEMBLY` | Assembler exited 0 but no contig FASTA was produced |
| `QUAST_FAIL` | QUAST non-zero exit or missing `report.tsv` |
| `METRICS_INCOMPLETE` | Report parsed but a required row was missing or the 200 bp contig rows disagreed |
| `EXISTS` | Replicate directory already present and `--force` was not set |

Version pins checked before assembly: SPAdes `v4.2.0`, SKESA `2.5.1`, Unicycler `v0.5.1`, QUAST `v5.2.0`. fastBact has no separate version gate beyond the recorded binary SHA-256.

---

## 10. Determinism protocol

Use this protocol for the F26 variance claim. Do not compare wall time or peak RSS.

1. Build fastBact once (`cargo build --release` on `qfree`) and record `sha256sum target/release/fastBact`.
2. On one node class (partition `qfree`), run three replicates as three jobs. Default timeouts fit one tool and one replicate inside the 3-hour limit, so `--replicates 3` in a single job is rejected unless both timeouts are lowered enough that `3 × (assembler timeout + QUAST timeout) ≤ 10800`.

```text
sbatch scripts/slurm_benchmark_harness.sbatch \
  --manifest F26_DATASET_MANIFEST.tsv \
  --dataset {DATASET} \
  --tool {fastbact|spades|skesa|unicycler} \
  --outdir {OUT} \
  --replicate {1|2|3}
```

3. Each job uses `scripts/slurm_benchmark_harness.sbatch`, so the cgroup is 8 CPUs, 14 GB, and `#SBATCH -t 03:00:00`. The three jobs share `{OUT}`. Re-running the same `dataset_id`, `tool`, and `replicate` replaces that row and keeps the other replicates.
4. Pass condition, all three replicates in one identity group:
   - status `PASS`
   - identical `assembly_sha256`
   - exact equality of every metric in section 7, including null NGA50
5. `wall_clock_s`, `cpu_time_s`, and `peak_rss_mb` are excluded. They are scheduler measurements, not assembly outputs.
6. The harness writes `determinism.pass=true` only when the group has at least three replicates and those equalities hold, and records `variance_max=0` for numeric metrics. Fewer than three replicates yields `checked=false`, `pass=null`, and reason `not enough replicates`. Any hash or metric difference is `pass=false`. Replicates must not be averaged to hide a mismatch.
7. Replicates pool only when dataset, tool, read paths, read SHA-256, reference path, reference SHA-256, tool-binary SHA-256, git HEAD, threads, and memory all match. A change to the fastBact binary, the read files, or the reference is a different experiment and cannot be pooled with the previous replicates. Read SHA-256 is recorded on real runs and left empty on `--dry-run`.

The F25 oracle statement `σ² = 0` on identical inputs is the same predicate: bitwise assembly identity and exact metric identity.

---

## 11. How to run

Self-test (no assembler, safe on a login node):

```text
python3 scripts/run_competitor_benchmark.py --self-test
```

Command preview of one dataset and one tool:

```text
python3 scripts/run_competitor_benchmark.py \
  --manifest F26_DATASET_MANIFEST.tsv \
  --dataset R1-001 \
  --tool spades \
  --outdir /tmp/f26-dry-run \
  --dry-run
```

Compute-node job, one dataset and one tool:

```text
mkdir -p logs
sbatch scripts/slurm_benchmark_harness.sbatch \
  --manifest F26_DATASET_MANIFEST.tsv \
  --dataset R1-001 \
  --tool spades \
  --replicate 1 \
  --outdir benchmarks/f26_results
```

Repeat with `--tool fastbact`, `--tool skesa`, and `--tool unicycler`. Repeat `--replicate 2` and `--replicate 3` into the same `--outdir` for the N=3 determinism gate. The batch script requests partition `qfree`, 8 CPUs, 14 GB, wall clock `03:00:00`, and refuses to start when `SLURM_JOB_ID` is unset.

---

## 12. Non-goals

- No edits under `src/`, `data/`, or `F26_DATASET_MANIFEST.tsv`.
- No coverage subsampling (that is F26-P3, not this harness).
- No claim that a dry-run or self-test is a scientific benchmark result.
