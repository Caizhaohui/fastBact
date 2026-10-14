# `fastBact` — Automated Benchmark & Validation Execution Plan (FB-BENCH-001-EXEC)

## 0. Agent Execution Directives & Core Rules (STRICTLY ENFORCED)

You are an expert bioinformatics QA and systems benchmarking engineer executing the official benchmark and validation suite (`FB-BENCH-001`) for **`fastBact`**.

### 0.1 Core Philosophy & Prohibitions
1. **Priority Order**: Correctness → Biological Assembly Accuracy → Reproducibility & Resource Limits → Speed. A faster assembly with unsupported structural joins (misassemblies) is a regression and an immediate `FAIL`.
2. **HPC Login Node Ban (CRITICAL)**:
   - **NEVER** run `cargo`, assemblers (`fastBact`, `spades.py`, `shovill`, `skesa`), read simulators, `fasterq-dump`, `pigz`, or `quast.py` on the login node.
   - The login node may ONLY be used for writing scripts/configs, lightweight environment inspection, and submitting jobs via `./scripts/hpc_run.sh` or `./benchmarks/scripts/bench_submit.sh` (`sbatch --wait` to partition `qfree` on `vnode1`).
   - **DO NOT create any new Conda environments**. Discover and reuse existing Conda environments on the server for `spades.py`, `shovill`, `skesa`, `quast.py`, `seqtk`, and `python3`.
3. **Strict Anti-Cheating Rules**:
   - Do NOT modify reference truth sequences, datasets, or QUAST thresholds to make `fastBact` pass.
   - Do NOT drop failed, timed-out, or OOM runs from summaries; never backfill failed runs with `0` elapsed time or `0` RSS (use `null` in JSON and `NA` in TSV).
   - Commit only scripts, schemas, tiny `S0/N` fixtures (< 100 KB), and summary markdown/TSV reports to Git. Keep all simulated/real FASTQ files and `benchmark-results/` outside Git (`.gitignore`).

---

## 1. Directory Architecture & Data Contracts

Create and maintain the following exact directory structure:

```text
fastBact/
├── BENCHMARK_PLAN.md
├── benchmarks/
│   ├── README.md
│   ├── datasets.tsv                    # 28-column immutable dataset manifest
│   ├── configs/
│   │   ├── tools.json                  # Discovered paths & versions of fastBact, SPAdes, Shovill, SKESA, QUAST
│   │   └── suite_matrix.json           # Parameters for S0, S1, S2, R1, R2, N suites
│   ├── scripts/
│   │   ├── validate_manifest.py        # Validates datasets.tsv schema, SHA256 & admission_status
│   │   ├── build_fixtures.py           # Generates deterministic S0, N, S1, S2 datasets + truth_events.json
│   │   ├── curate_real_data.sh         # Fetches/links & validates R1/R2 real bacterial FASTQ + reference
│   │   ├── bench_submit.sh             # Login-node sbatch --wait wrapper
│   │   ├── bench_compute.sh            # Compute-node single-run execution harness (/usr/bin/time -v)
│   │   ├── evaluate_truth.py           # Independent junction TP/FP/FN & plasmid circularity validator
│   │   ├── collect_metrics.py          # Parses time.txt, slurm sacct, QUAST (500bp + 200bp passes) -> metrics.json
│   │   └── generate_report.py          # Aggregates all runs into scorecard tables & PASS/WARN/FAIL verdicts
│   ├── fixtures/                       # Tiny deterministic S0 (1-10kb) & N (negative) fixtures
│   └── schemas/
│       ├── run.schema.json
│       └── metrics.schema.json
└── benchmark-results/                  # Ignored by git (.gitignore)
    ├── data_store/                     # Immutable generated S1/S2 & downloaded R1/R2 FASTQ.gz + references
    └── <campaign>/<dataset_id>/<tool>/<commit>/<rep>/
        ├── command.txt
        ├── environment.json
        ├── run.json
        ├── time.txt
        ├── slurm.txt
        ├── stdout.log
        ├── stderr.log
        ├── quast_500/                  # Standard headline QUAST (--min-contig 500)
        ├── quast_200/                  # Small-contig/plasmid QUAST pass (--min-contig 200)
        └── metrics.json
```

### 1.1 `benchmarks/datasets.tsv` Schema (28 Mandatory Columns)
Every registered test input MUST have one row with these exact tab-separated columns:
`dataset_id	kind	organism	source_url	study_accession	run_accession	biosample	reference_accession_version	reference_confidence	replicons	platform	layout	read_length	fragment_mean	fragment_sd	coverage_nominal	coverage_measured	simulator_version	random_seed	fastq_r1	fastq_r2	reference_fasta	truth_events_file	r1_sha256	r2_sha256	reference_sha256	admission_status	notes`

* **Rules**: Missing numeric values MUST be `NA` (never `0`). `kind` $\in \{\text{S0}, \text{S1}, \text{S2}, \text{R1}, \text{R2}, \text{N}\}$. `reference_confidence` $\in \{\text{verified}, \text{qualified}, \text{unavailable}\}$. Only rows with `admission_status=APPROVED` are executed in official campaigns.

### 1.2 `run.json` & `metrics.json` Contracts
Each execution directory `<campaign>/<dataset_id>/<tool>/<commit>/<rep>/` must emit:

**`run.json`**:
```json
{
  "campaign_id": "FB-BENCH-v0.1.0",
  "dataset_id": "S2-001-ecoli-100x",
  "track": "E2E",
  "tool": "fastBact",
  "tool_version": "0.1.0",
  "git_commit": "<full SHA>",
  "git_dirty": false,
  "threads": 8,
  "memory_limit_gib": 32,
  "node_class": "<hostname + CPU model>",
  "slurm_job_id": "<job_id>",
  "repetition": 1,
  "seed": 101,
  "command": "<exact command executed>",
  "exit_code": 0,
  "status": "PASS",
  "wall_s": 17.34,
  "cpu_user_s": 88.20,
  "cpu_system_s": 2.15,
  "peak_rss_kib": 579072,
  "scratch_peak_bytes": 0,
  "reference_confidence": "verified",
  "assembly_sha256": "<SHA256 contigs.fa of>",
  "quast_metrics_file": "quast_500/report.tsv",
  "evaluation_version": "QUAST 5.2.0"
}
```

**`metrics.json`**:
```json
{
  "schema_version": "1.0",
  "extensive_misassemblies": 0,
  "local_misassemblies": 0,
  "genome_fraction_pct": 99.201,
  "nga50": 198925,
  "ng50": 198925,
  "n50": 198925,
  "largest_contig": 298924,
  "total_length": 2503193,
  "contig_count_500bp": 15,
  "contig_count_200bp": 15,
  "duplication_ratio": 1.001,
  "unaligned_length": 0,
  "mismatches_per_100kbp": 0.12,
  "indels_per_100kbp": 0.04,
  "junction_precision": 1.0,
  "junction_recall": 0.92,
  "plasmid_recovery_pct": 100.0,
  "plasmid_circularity_correct": true,
  "plasmid_chimera_count": 0
}
```

---

## 2. Two Comparison Tracks (Do Not Mix)

1. **Track 1 — E2E Product Comparison (`track=E2E`)**:
   - Supply identical **raw PE FASTQ** to each tool with `-t 8` (8 CPUs) and 32 GiB RAM limit:
     - `fastBact`: `./target/release/fastBact --r1 ${R1} --r2 ${R2} --outdir${OUT} --threads 8 --depth 120 --min-len 500`
     - `SPAdes`: `spades.py --isolate -1 ${R1} -2 ${R2} -o${OUT} -t 8 -m 32`
     - `Shovill`: `shovill --R1 ${R1} --R2 ${R2} --outdir${OUT} --cpus 8 --ram 32`
     - `SKESA` (if installed on server): `skesa --reads ${R1},${R2} --cores 8 --memory 32 >${OUT}/contigs.fa`
   - Includes each tool's internal decompression, subsampling, trimming, error correction, assembly, polishing, and temporary disk I/O in its measured wall time and scratch usage.
2. **Track 2 — Controlled Core Comparison (`track=CORE`)**:
   - Pre-generate an immutable, adapter-trimmed, pair-synchronized FASTQ dataset fixed at exactly `60x` and `120x` depth outside the timed block.
   - Run `fastBact` (`--depth 1000` to bypass further subsampling), `spades.py --only-assembler --isolate -t 8`, and `shovill --depth 0 --noreadcorr --nostitch --nocorr --cpus 8` (verify supported flags on the installed Shovill binary) to measure pure graph construction and repeat resolution speed on identical effective input bases.

---

## 3. Phased Implementation Workflow (Phase B0 → Phase B5)

Execute the benchmark suite construction and validation strictly in order:

### Phase B0: HPC Tool Discovery & Harness Bootstrap
**Tasks**:
1. Inspect all existing Conda environments on the server (`conda info --envs`) to locate exact binary paths and versions for `fastBact` (`./target/release/fastBact`), `spades.py`, `shovill`, `skesa` (mark `UNAVAILABLE` in `tools.json` if not installed), `quast.py`, `seqtk`, and `python3`. Save to `benchmarks/configs/tools.json`.
2. Record the current `git rev-parse HEAD` of `fastBact` and ensure `benchmark-results/` is added to `.gitignore`.
3. Implement `benchmarks/scripts/validate_manifest.py`, `benchmarks/scripts/bench_submit.sh` (submits to `qfree` on `vnode1` via `sbatch --wait` requesting `--cpus-per-task=8 --mem=32G`), and `benchmarks/scripts/bench_compute.sh` (captures `/usr/bin/time -v`, monitors scratch directory peak bytes, runs `sacct`, and writes `run.json`).

**Phase B0 Verification Gate**:
- Run `python3 benchmarks/scripts/validate_manifest.py` and test `./benchmarks/scripts/bench_submit.sh` on a dummy command to verify `environment.json`, `time.txt`, and `slurm.txt` are accurately captured on `vnode1`.

---

### Phase B1: Negative (`N`) & Micro (`S0`) Invariant Suite (Level L0)
**Tasks**:
1. Implement in `benchmarks/scripts/build_fixtures.py` deterministic generators for **Tier `S0`** and **Tier `N`**:
   - `S0-001`: 5 kb linear genome, no repeats, 50× PE.
   - `S0-002`: 4 kb isolated circular plasmid (`circular=true` truth), 80× PE.
   - `S0-003`: 8 kb sequence with synthetic sequencing error tips and a 1-SNP bubble, 60× PE.
   - `S0-004`: 6 kb sequence containing palindromic k-mers and reverse-complement junction branches.
   - `N-001`: Truncated/corrupted gzip FASTQ file.
   - `N-002`: Unequal mate counts (R1 has 1,000 reads, R2 has 850 reads).
   - `N-003`: Reads consisting of >50% `N` bases or all-low Phred quality (`#`).
2. Register all `S0` and `N` datasets in `benchmarks/datasets.tsv` with SHA256 checksums and `admission_status=APPROVED`.
3. Verify hard gate **H1 & H2**: `fastBact` must assemble `S0-001..004` with 100% accuracy and must cleanly reject or safely handle `N-001..003` with clear error messages (zero panics, zero infinite loops, zero silent pair corruption). If `fastBact` panics on `N-001` or `N-002`, fix input validation in `src/ingest/` before proceeding!

**Phase B1 Verification Gate**:
- Submit the `S0` + `N` test run to `vnode1` via `./benchmarks/scripts/bench_submit.sh`. All 7 cases must pass their expected outcomes (`PASS` for `S0`, graceful non-crash rejection/handling for `N`).

---

### Phase B2: Synthetic Repeat, Plasmid, Coverage & GC Matrix (`S1` & `S2`, Level L1)
**Tasks**:
1. Extend `benchmarks/scripts/build_fixtures.py` (with fixed random seeds and recorded truth junction coordinates in `truth_events.json`) to generate the mandatory synthetic matrix in `benchmark-results/data_store/`:
   - **`S1-REP-001..003` (Repeat Length & Evidence Limit Matrix, 150 kb genomes, fragment 450±30 bp)**:
     - `S1-REP-001`: Exact repeats of `35 bp`, `65 bp`, `120 bp` (resolvable by Multi-$k$).
     - `S1-REP-002`: Exact repeat of `350 bp` (3 copies, resolvable by `exSPAnder` PE links).
     - `S1-REP-003`: Exact repeats of `800 bp` and `1500 bp` (exceeds 450 bp insert size; tests **conservatism** — assembler MUST abstain and split cleanly at repeat boundaries with `0` misassemblies rather than guessing false joins!).
   - **`S1-PLAS-001..002` (Plasmid & Shared IS Matrix, 200 kb chromosome)**:
     - `S1-PLAS-001`: Three circular plasmids (`2.2 kb` at 25× copy ratio, `5.5 kb` at 5× copy ratio, `45 kb` at 0.5× copy ratio).
     - `S1-PLAS-002`: Chromosome + `8 kb` circular plasmid sharing an identical `1.0 kb` IS element (tests hard gate **H4**: zero false chromosome–plasmid chimeras).
   - **`S1-COV-001..003` (Coverage Matrix, 200 kb genome)**: `15x` (low depth), `60x` (standard), `300x` (extreme depth testing subsampling & memory cap).
   - **`S1-GC-001..002` (Extreme GC Matrix, 200 kb genome)**: `30% GC` (AT-rich) and `70% GC` (GC-rich with GC-correlated coverage dropout).
   - **`S1-EDIT-001` (Structural Variant Truth)**: 200 kb parental genome with a known `1 kb` deletion, `1.2 kb` insertion, and `5 kb` inversion.
   - **`S2-REAL-001` (Full-Scale 4.6 Mb Bacterial Isolate)**: Complete 4.6 Mb chromosome (modeled on *E. coli* K-12 scale with 7 copies of 5 kb rRNA operons, 10 copies of 1.2 kb IS elements) + 15 kb circular plasmid (8× copy number), simulated at `100x` depth ($2 \times 150\text{ bp}$, fragment $420 \pm 35\text{ bp}$, $0.5\%$ substitution + indel error profile with 3' quality decay).
2. Implement `benchmarks/scripts/evaluate_truth.py` and `benchmarks/scripts/collect_metrics.py` to run both QUAST passes (`--min-contig 500` and `--min-contig 200`) and compute junction precision/recall and plasmid circularity accuracy.

**Phase B2 Verification Gate**:
- Execute all `S1` and `S2` datasets on `vnode1` for `fastBact`, `SPAdes`, and `Shovill`. Confirm `fastBact` achieves **0 extensive misassemblies** across all `S1/S2` cases (including `S1-REP-003` and `S1-PLAS-002`).

---

### Phase B3: Real-World Isolate Curation & Evaluation (`R1` & `R2`, Level L2)
**Tasks**:
1. Implement `benchmarks/scripts/curate_real_data.sh` to discover local real bacterial paired-end FASTQ datasets on the server OR retrieve/prepare at least **3 real-world bacterial isolate datasets** (e.g., from the candidate pool: *E. coli* K-12 MG1655 `PRJEB2323`, *S. aureus* `SRR569301`, *R. sphaeroides* `SRR522246` high-GC, or existing real bacterial sequencing datasets available on the server filesystem).
2. Verify read pairing, compute SHA256 hashes, assign `reference_confidence` (`verified` if matched reference FASTA is available, or `qualified`/`unavailable` if evaluated reference-free or against a close strain), and register them as `APPROVED` in `benchmarks/datasets.tsv`.
3. Run the E2E comparison (`fastBact`, `SPAdes --isolate`, `Shovill`, and `SKESA` if available) on all curated `R1/R2` datasets on `vnode1`.

**Phase B3 Verification Gate**:
- Verify all `R1/R2` runs produce complete `run.json` and `metrics.json` bundles and pass Hard Gate **H1** (0 crashes/OOMs).

---

### Phase B4: Reproducibility, ABBA Repetitions & Thread Scaling (`L3` Performance Matrix)
**Tasks**:
1. **Track 1 (E2E) & Track 2 (Controlled Core) Repeated Benchmarks**:
   - On the full-scale `S2-REAL-001` (4.6 Mb) and primary `R1` real isolate dataset, run **3 repetitions** in alternating tool order (`fastBact -> SPAdes -> Shovill`, then reversed) on `vnode1` (8 threads, 32 GiB limit).
   - Compute median, min, max, and MAD for `wall_s` and `peak_rss_kib`.
2. **Thread Scaling Matrix**:
   - Run `fastBact` on `S2-REAL-001` with `threads = 1, 2, 4, 8` on `vnode1`. Compute speedup $S_n = T_1 / T_n$ and parallel efficiency $E_n = T_1 / (n \times T_n)$.
3. **Deterministic Reproducibility Check**:
   - Verify that repeated runs of `fastBact` with the same seed and thread count produce semantically identical assemblies.

**Phase B4 Verification Gate**:
- Confirm all repetition and scaling runs are logged in `benchmark-results/` with valid `run.json` and `metrics.json` files.

---

### Phase B5: Automated Scorecard & Final Benchmark Report Generation
**Tasks**:
1. Implement `benchmarks/scripts/generate_report.py` to scan `benchmark-results/`, validate all JSON outputs against `schemas/`, and generate `benchmarks/BENCHMARK_REPORT.md` containing:
   - **Executive Gate Summary**: Explicit `PASS / WARN / FAIL` status for Hard Gates `H1–H5`.
   - **Table 1 — L0 (`S0` & `N`) Invariant & Negative Suite Scorecard**.
   - **Table 2 — L1 (`S1` & `S2`) Synthetic Truth Matrix Comparison (`fastBact` vs `SPAdes` vs `Shovill`)**: Reporting in strict priority order: `(1) Fatal/Extensive Misassemblies -> (2) Local Misassemblies -> (3) Genome Fraction % & Plasmid Recovery -> (4) Mismatches & Indels / 100 kbp -> (5) NGA50 & Contig Count -> (6) Median Wall Time (s) & Peak RSS (MB)`.
   - **Table 3 — L2 (`R1` & `R2`) Real-World Isolate Comparison**.
   - **Table 4 — Track 1 (E2E) vs Track 2 (Controlled Core) & Thread Scaling (1, 2, 4, 8 Threads) Analysis**.
   - **Triage & Regression Analysis**: Highlight any edge cases where `fastBact` showed warnings or differences vs `SPAdes`, with concrete recommendations.

**Phase B5 Final Verification Gate**:
- Run `python3 benchmarks/scripts/generate_report.py` and verify `benchmarks/BENCHMARK_REPORT.md` is complete, accurate, and contains zero fabricated or placeholder numbers.