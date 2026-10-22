# Phase F20 Independent Benchmark Audit

**Audit ID:** AUDIT-20261022-015  
**Auditor:** Codex CLI  
**Date:** 2026-10-22  
**Scope:** Selective real-data A/B benchmark for R1-001 and R1-003  
**Overall verdict:** **FAIL (REGRESSION DETECTED)**  
**C0-02R disposition:** **Strictly NOT ACCEPTED**

## 1. Findings

1. **Acceptance blocker:** R1-001 fails six quality criteria: extensive misassemblies, mismatches, indels, genome fraction, N50, and NGA50.
2. **Acceptance blocker:** R1-003 fails genome fraction despite improvements in five other criteria.
3. **Documentation inconsistency:** `slurm_job_ledger.tsv` retains `RUNNING`; archived accounting and the execution trace establish completion.
4. **Qualification:** Baseline assembly quality and telemetry reproduce historical results exactly. Timing and memory measurements are new observations and do not match historical measurements bit-for-bit.
5. **Provenance limitation:** Live scheduler verification was blocked by sandbox networking. Execution findings therefore rely on archived scheduler evidence cross-checked against allocation, host, command, and timing records.

No discrepancy was found in the audited metric extraction or acceptance calculations.

## 2. Provenance and execution integrity

### Inputs and binaries

Independently reran SHA-256 verification against `logs/expected.sha256`. All ten entries passed:

- R1-001: both SRR28349173 FASTQ files and reference `NC_000913.3.fasta`.
- R1-003: both SRR28349164 FASTQ files and reference `CP009273.1.fasta`.
- Baseline A and Candidate B binaries.
- Both recorded Cargo.lock files.

The six dataset digests agree with `dataset_checksums.tsv` and the benchmark command script.

| Binary | Independently verified SHA-256 |
|---|---|
| Baseline A | `2b184db0f108cb508ca7e591240eddb551bea8ebad90ac5eb8ee360bac78ea38` |
| Candidate B | `4fbafeeeffe94abeef42efaa5710f75dbd0af8b234edd2bffba43cc4f50198f0` |

Saved before/after checksum records report successful verification of these same files. The command trace records both verification passes.

The benchmark script’s independently computed SHA-256 matches `benchmark_commands.sha256`:

```text
771bc0a25aeae97a2e7d4c8b8924fba4b3bd8550fa7845185825215a23f75512
```

The QUAST launcher also passes its recorded checksum verification; the recorded version is **QUAST v5.2.0**.

### Artifact integrity

Rechecked all **242 entries** in `artifact_checksums.sha256`: **242 passed, zero failed**. No duplicate manifest paths were found.

The manifest covers 240 files under assemblies, QUAST outputs, and eligible logs, plus the metrics and acceptance tables. The evidence-directory copies of these three files are byte-identical to their benchmark-results counterparts:

- `assembly_metrics.tsv`
- `acceptance.tsv`
- `artifact_checksums.sha256`

The manifest intentionally excludes Slurm output/error logs and command traces. It also does not seal every surrounding evidence document or the final status file. Checksum success establishes consistency with the supplied manifest; it does not provide an independently signed historical record.

### Slurm execution

Archived accounting records:

| Property | Recorded value |
|---|---|
| Job | 3001561 |
| Partition | qfree |
| Node | vnode1 |
| Allocated CPUs | 8 |
| Requested/allocated memory | 12G |
| Final state | COMPLETED |
| Exit code | 0:0 |
| Elapsed | 00:13:13 |

All eight numbered steps, `3001561.0` through `3001561.7`, report `COMPLETED`, exit `0:0`, on `vnode1`. Batch and extern records also completed successfully.

Allocation and host records corroborate the node and resources. The trace records four assembly and four QUAST invocations using eight CPUs. All eight `/usr/bin/time` records report exit zero.

This supports execution entirely on `vnode1`. It does **not** establish dedicated possession of the node: the allocation reports `OverSubscribe=OK`.

A fresh `sacct` query failed with `Operation not permitted` when opening the scheduler connection. Archived completion could therefore be corroborated locally but not independently refreshed.

### Production-source immutability

The benchmarked candidate worktree is:

```text
/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f19-codex-worktree
```

Its HEAD is `0606be423d187771f5925a83da2cc377c40c2a3d`, and its production source tree is:

```text
30967bf4940844aeb6b744732bcad6d2a72daab4
```

That same tree object is present at:

- F18 audited commit `c37ef6e03475230dc7f8596e0fc6788ba4fc6355`.
- F19 baseline `39f11ca90a5de13a5b3d5572d29c191c40d9c03d`.
- F19 audited commit `a9930471a73612159a83694d3d661f060042efc8`.
- F20 candidate HEAD `0606be423d187771f5925a83da2cc377c40c2a3d`.

All **28 candidate source files** were independently hashed against their committed Git blob identities without discrepancies. Candidate `src/`, Cargo.toml, and Cargo.lock are clean; no untracked source files were reported. Execution traces show matching identity and cleanliness checks before and after the benchmark.

**Conclusion:** Candidate production contents match the F18 audited tree, with no recorded endpoint drift. These checks cannot exclude an unobserved transient modification during execution.

The current main workspace has a different `src/` tree, `2cab60dd172daaaf247e49f3dc744eb91479dfc3`. The F18 identity finding applies specifically to the worktree used for Candidate B.

## 3. Metric authenticity and historical reproduction

Independently reconciled **56 recorded fields** across four runs against raw QUAST reports, assembly timing records, and junction/node log entries. Result: **zero mismatches**.

Source inspection confirms that extraction selects exact QUAST metric labels and preserves their reported values. Wall time and peak RSS come from `/usr/bin/time` using `%e`, `%M`, and `%x`. Slurm MaxRSS is a separate measurement and was not substituted.

Direct FASTA length calculations independently reproduced QUAST contig counts, total lengths, and N50:

| Dataset/version | Contigs | Total bp | N50 bp |
|---|---:|---:|---:|
| R1-001 A | 164 | 4,576,828 | 58,318 |
| R1-001 B | 167 | 4,573,069 | 57,462 |
| R1-003 A | 135 | 4,563,880 | 67,452 |
| R1-003 B | 144 | 4,558,923 | 70,868 |

For both datasets, F20 A’s FASTA, stage telemetry, and complete QUAST `report.tsv` are **byte-identical** to historical version B under:

```text
benchmark-results/C0_02R_pe_fix_20261019/results/assemblies/
```

Historical version B is the frozen baseline used as F20 A. Historical version A is a different binary.

The reproduced QUAST report hashes are:

| Dataset | SHA-256 |
|---|---|
| R1-001 | `b0f537aad84adab06cbbbacc7bfcecc0ff46bd32e595d08b66e5db846a73de33` |
| R1-003 | `ac5c40eaae1de65b5ce69976b4a7491180b426b33c4e9e96d40083b080fda515` |

The **quality and telemetry reproduction claim passes**. A claim that the entire metrics row reproduces bit-for-bit would be incorrect because runtime and memory measurements differ.

## 4. Independent acceptance evaluation

Applied the specified zero-regression rules independently: higher-or-equal N50, NGA50, and genome fraction; lower-or-equal structural and base-error metrics.

| Criterion | R1-001 A → B | Result | R1-003 A → B | Result |
|---|---:|---|---:|---|
| Extensive misassemblies | 4 → 5 | FAIL | 2 → 1 | PASS |
| Local misassemblies | 0 → 0 | PASS | 1 → 1 | PASS |
| Mismatches/100 kbp | 2.14 → 2.38 | FAIL | 1.51 → 1.08 | PASS |
| Indels/100 kbp | 0.48 → 1.03 | FAIL | 0.04 → 0.00 | PASS |
| Genome fraction | 98.364% → 98.281% | FAIL | 98.330% → 98.233% | FAIL |
| N50 bp | 58,318 → 57,462 | FAIL | 67,452 → 70,868 | PASS |
| NGA50 bp | 58,318 → 57,127 | FAIL | 63,648 → 70,010 | PASS |

**Total: 7 PASS and 7 REGRESSION comparisons.** All independently calculated outcomes agree with `acceptance.tsv`.

R1-001 loses 0.083 percentage points of genome fraction; its indel rate rises approximately 114.58%. R1-003 reduces extensive misassemblies by 50% and mismatches by approximately 28.48%, but loses 0.097 percentage points of genome fraction.

The final status:

```text
REGRESSION_DETECTED	0	3001561
```

is correct. Exit zero denotes successful benchmark execution; it does not denote quality acceptance.

## 5. Verdict and recommendation

**Overall verdict: FAIL (REGRESSION DETECTED).**  
**Stage C0-02R: Strictly NOT ACCEPTED.**

Stop benchmark expansion and preserve all evidence. Conduct focused investigation of R1-001’s additional misassembly and increased indels by linking QUAST breakpoints and error loci to matched A/B graph paths and pre-/post-polishing sequences. Investigate completeness losses on both datasets.

The existing scientific interpretation appropriately treats causal explanations as hypotheses. This audit confirms the measured regressions but does not establish their locus-level mechanisms. One replicate per dataset/version also limits performance and variability conclusions.

