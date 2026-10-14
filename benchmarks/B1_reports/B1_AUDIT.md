# `fastBact` Benchmark Reliability Forensic Audit (B1-1)

**Document ID**: FB-BENCH-B1-AUDIT  
**Date**: 2026-10-14  
**Audit Target**: `fastBact v0.1.0` (Baseline Commit `4259fe3ffe077f67946700e8092ccf46fbbb9bc6`)  
**Audit Scope**: Manifest reconciliation, comparator versions, failure attribution, isolate provenance, and plasmid accounting.

---

## 1. Dataset Count Reconciliation (22 vs 19) — `[VERIFIED]`

### Forensic Finding:
- In `benchmarks/datasets.tsv`, there are exactly **22 registered datasets**:
  - **Level L0 (`S0`) Micro Invariants**: 4 datasets (`S0-001`, `S0-002`, `S0-003`, `S0-004`).
  - **Level L0 (`N`) Negative Controls**: 3 datasets (`N-001`, `N-002`, `N-003`).
  - **Level L1 (`S1`) Synthetic Matrix**: 11 datasets (`S1-REP-001..003`, `S1-PLAS-001..002`, `S1-COV-001..003`, `S1-GC-001..002`, `S1-EDIT-001`).
  - **Level L1 (`S2`) Full Bacterial Genome**: 1 dataset (`S2-REAL-001`).
  - **Level L2 (`R1`) Real Isolate Cohort**: 3 datasets (`R1-001`, `R1-002`, `R1-003`).
  - **Total Unique Datasets**: $4 + 3 + 11 + 1 + 3 = 22$.
- **Discrepancy Explained**:
  - The number of unique positive dataset IDs is $4 (\text{S0}) + 11 (\text{S1}) + 1 (\text{S2}) + 3 (\text{R1}) = 19$.
  - In the initial Gate H1 report, "22 positive tests verified" referred to the **22 completed execution runs** of `fastBact`:
    - 19 unique positive dataset runs + 3 independent reproducibility runs on `S2-REAL-001` (`rep_1`, `rep_2`, `rep_3`).
  - Each dataset and run is 100% accounted for and mapped to `benchmarks/datasets.tsv`.

---

## 2. Shovill Failure Root Cause Verification — `[VERIFIED]`

### Forensic Finding:
Shovill v1.1.0 failed on 6 synthetic datasets: `S1-REP-001`, `S1-REP-003`, `S1-PLAS-001`, `S1-COV-001`, `S1-COV-002`, and `S1-EDIT-001`. Forensic inspection of `tool_out/shovill.log` and compute node verification experiments identified two distinct technical failure mechanisms:

1. **Root Cause A: FLASH Empty Merged Reads Bug in Shovill (`[VERIFIED]`)**:
   - In datasets where insert size ($450\text{ bp}$) $> 2 \times$ read length ($300\text{ bp}$), FLASH cannot merge any reads (`Combined pairs: 0`, `Percent combined: 0.00%`).
   - FLASH writes an empty gzip file `flash.extendedFrags.fastq.gz` (20 bytes).
   - In `/hpcfs/fhome/caizhh/.conda/envs/LineaBact/bin/shovill` (lines 287–288):
     ```perl
     $opts .= " $stitch_arg $SE" if $SE;
     ```
     Shovill checks whether the string `$SE` is non-empty, but fails to check whether `-s $SE` (file size $> 0$).
   - Consequently, Shovill invokes `spades.py` with `--merged flash.extendedFrags.fastq.gz`.
   - `spades.py v4.2.0` crashes with:
     ```text
     == Error == file is empty: .../flash.extendedFrags.fastq.gz
     ```
   - **Verification Experiment**: When Shovill was executed on `vnode1` with `--nostitch` on `S1-REP-001`, it completed successfully in 17 seconds, assembling a single 149,986 bp contig.

2. **Root Cause B: KMC `-ci10` Low-Coverage Dropout (`[VERIFIED]`)**:
   - On `S1-COV-001` (15x depth), Shovill executes `kmc -k21 -ci10` to estimate genome size.
   - At 15x nominal coverage, true genomic $k$-mers have mean depth $\sim 15$ with a Poisson tail dropping below 10. `kmc -ci10` discarded $> 95\%$ of genomic $k$-mers, estimating genome size as 16,119 bp instead of 200,000 bp.
   - This caused incorrect coverage calculations and erroneous read subsampling, compounded by the empty `--merged` crash.

---

## 3. Comparator Version Audit — `[VERIFIED]`

All comparator tools were probed directly on compute node `vnode1` within the `LineaBact` Conda environment:

| Tool Component | Executable Path | Version Output | Role in Benchmark |
| :--- | :--- | :--- | :--- |
| **fastBact** | `./target/release/fastBact` | `v0.1.0` (commit `4259fe3`) | Under evaluation |
| **SPAdes** | `/hpcfs/fhome/caizhh/.conda/envs/LineaBact/bin/spades.py` | `v4.2.0` | Standalone de novo assembler |
| **Shovill** | `/hpcfs/fhome/caizhh/.conda/envs/LineaBact/bin/shovill` | `1.1.0` (Perl v5.32.1) | Assembly pipeline wrapper |
| **Internal SPAdes** | `/hpcfs/fhome/caizhh/.conda/envs/LineaBact/bin/spades.py` | `v4.2.0` | Invoked by Shovill engine |
| **KMC** | `/hpcfs/fhome/caizhh/.conda/envs/LineaBact/bin/kmc` | `ver. 3.2.4 (2024-02-09)` | K-mer counting in Shovill |
| **FLASH** | `/hpcfs/fhome/caizhh/.conda/envs/LineaBact/bin/flash` | `FLASH v1.2.11` | Overlap merging in Shovill |
| **Lighter** | `/hpcfs/fhome/caizhh/.conda/envs/LineaBact/bin/lighter` | `Lighter v1.1.3` | Read error correction in Shovill |
| **BWA** | `/hpcfs/fhome/caizhh/.conda/envs/LineaBact/bin/bwa` | `Version: 0.7.19-r1273` | Read mapping in Shovill |
| **Pilon** | `/hpcfs/fhome/caizhh/.conda/envs/LineaBact/bin/pilon` | `Version 1.24` | Assembly polishing in Shovill |
| **SKESA** | `/hpcfs/fhome/caizhh/.conda/envs/LineaBact/bin/skesa` | `SKESA 2.5.1` | Strategic alternative engine |
| **QUAST** | `/hpcfs/fhome/caizhh/.local/bin/quast.py` | `QUAST v5.2.0` | Independent metrics evaluation |

---

## 4. Real Dataset Provenance & Reference Audit (`R1-001..003`) — `[VERIFIED]`

Audit of NCBI BioProject `PRJNA1088182` SRA metadata and publication truth tables:

1. **Sequencing Platform & Library Architecture**:
   - Instrument: **Illumina NovaSeq X Plus** (paired-end 2x150 bp).
   - Library strategy: `WGS`, Selection: `RANDOM`, Source: `GENOMIC`.
   - Read characteristics: Mean read length 150 bp, insert size $\sim 300\text{--}450\text{ bp}$, forward-reverse (`FR`) orientation.
2. **Biological Cohort Details**:
   - `R1-001` (`SRR28349173`): BioSample `SAMN40463664`. *E. coli* K-12 MG1655 transformed with Cas9 targeting vector (`Cas9_110_111`). Total reads: 2,409,728 pairs ($\sim 78\times$ coverage). Validated against reference `NC_000913.3`.
   - `R1-002` (`SRR28349221`): BioSample `SAMN40463696`. *E. coli* BW25113 wild-type control (`control_17_18_19_20`). Total reads: 1,993,895 pairs ($\sim 65\times$ coverage). Validated against reference `CP009273.1`.
   - `R1-003` (`SRR28349164`): BioSample `SAMN40463672`. *E. coli* BW25113 recombineering strain with $\lambda\text{Red}$ helper plasmid (`lambdaRed_9_10`). Total reads: 4,093,992 pairs ($\sim 132\times$ coverage). Validated against reference `CP009273.1`.
3. **Plasmids / Vector Contigs in Real Reads (`[LIKELY]`)**:
   - Unaligned contigs in `R1-001` (2,251 bp unaligned length) and `R1-003` represent exogenous helper plasmid sequences (e.g., pKD46, Cas9/sgRNA vectors) that do not align to the circular chromosomal reference genome.

---

## 5. Accounting & Plasmid Audit — `[VERIFIED]`

1. **Short Contig Accounting ($< 500\text{ bp}$ vs $< 200\text{ bp}$)**:
   - In all tested datasets (`S1-COV-001`, `S2-REAL-001`, `R1-001`), QUAST reports identical contig counts at `--min-contig 200` and `--min-contig 500` (e.g., exactly 144 contigs for `S2-REAL-001`).
   - Root cause: `fastBact`'s default `--min-len 500` filter discards unitigs $< 500\text{ bp}$ during `curate_contigs` before writing `contigs.fa`.
2. **Plasmid Recovery Audit Across Test Cases**:
   - `S0-002` (4.0 kb plasmid): Assembled as a single 4,000 bp contig, 100.0% sequence identity, `circular=true`.
   - `S1-PLAS-001` (3 plasmids: 2.2 kb, 5.5 kb, 45.0 kb):
     - Contig 2: 45,000 bp, `cov=20.1`, `circular=true` (100% recovered).
     - Contig 6: 5,500 bp, `cov=112.6`, `circular=true` (100% recovered).
     - Contig 8: 2,201 bp, `cov=20.0`, `circular=true` (100% recovered).
   - `S1-PLAS-002` (8.0 kb plasmid sharing 1.0 kb IS with chromosome):
     - Contig 3: 8,000 bp, `cov=118.6`, `circular=true` (0 chimeras).
   - `S2-REAL-001` (15.0 kb plasmid):
     - Contig 97: 15,000 bp, `cov=100.8`, `circular=true` (100% recovered).
