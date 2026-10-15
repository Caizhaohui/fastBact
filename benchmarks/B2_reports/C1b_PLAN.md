# fastBact C1b Implementation Plan: Accurate Coverage Estimation

**Document ID:** FB-PLAN-B2-C1b  
**Target Milestone:** `fastBact v0.2.0-alpha.c1b`  
**Author:** AI Pair Programmer & Bioinformatic Lead  
**Date:** 2026-10-15  
**Execution Policy:** Strict 3-Phase Progression. **Only C1b-1 is active in this turn.** C1b-2 and C1b-3 are strictly deferred.

---

## 1. Background & Problem Statement

In the B2-G0 forensics audit (`B2_ESTIMATOR_FORENSICS.md`), we unmasked two compounded failures in `v0.1.0` estimator:
1. **Uncalibrated HLL 32-bit Hash Truncation:** A 32-bit MurmurHash3 truncated to 14 bits caused an estimated genome size of 522 kb on a 4.6 Mb *E. coli* genome (8.8-fold underestimation).
2. **Artificial Read Count Extrapolation via File Bytes:** The estimator scanned only the first 100k reads, recorded file offset `bytes_read`, and extrapolated total pairs via:
   $$\widehat{N}_{pairs} = \frac{\text{file\_size\_bytes}}{\text{bytes\_read\_per\_pair}} = \frac{\text{file\_size\_bytes}}{75.0}$$
   Because `file_size_bytes` is the *compressed gzip stream length* while `bytes_read_per_pair` was measured on *uncompressed ASCII text*, this dimensional inconsistency exaggerated read depth to 626.4x and triggered an aggressive 8.5% subsampling, discarding 91.5% of raw sequencing data and collapsing assembly N50 from 58.3 kb to 3.6 kb.

To eliminate these systematic failures permanently, Phase C1b establishes a mathematically rigorous, 3-stage estimation pipeline.

---

## 2. Three-Stage Implementation Blueprint

```
+-----------------------------------------------------------------------------------+
| Phase C1b-1 (CURRENT TURN ONLY): Exact FASTQ Streaming Counting                   |
| - Abolish approx_bytes_per_pair = 75 and gzip stream extrapolation                |
| - Fast, exact streaming counter: N_pairs, B_total, L_min, L_max, L_mean, L_median  |
| - Dual-end paired header synchronization & truncation integrity validation        |
| - Zero algorithmic change to assembly (retaining 100% reads, bit-for-bit verified)|
+-----------------------------------------------------------------------------------+
                                         |
                                         v
+-----------------------------------------------------------------------------------+
| Phase C1b-2 (NEXT TURN): k-mer Distinct Volume & Accurate Genome Size             |
| - Replace 32-bit MurmurHash3 with 64-bit WyHash/AHash                             |
| - Implement k-mer abundance histogram / multiplicity profiling                    |
| - Correct genome size estimation: G = (B_total - (k-1)*2*N_pairs) / Coverage_peak |
| - High-confidence calibration gate                                                |
+-----------------------------------------------------------------------------------+
                                         |
                                         v
+-----------------------------------------------------------------------------------+
| Phase C1b-3 (FUTURE TURN): Dynamic CMS Sizing & Safe Auto-Subsampling Closure      |
| - Dynamic Count-Min Sketch capacity based on calibrated G                         |
| - Close the loop on auto-subsampling: only enable when G is calibrated            |
| - Full end-to-end benchmark on 10x - 300x depth gradients                         |
+-----------------------------------------------------------------------------------+
```

---

## 3. Detailed Specifications for Phase C1b-1 (Current Turn)

### 3.1 Streaming Exact Counter (`count_fastx_exact`)
- **Input:** Paths to paired FASTQ files (gzipped `.fastq.gz` or uncompressed `.fastq`).
- **Data Structures:**
  - `ExactCountStats`:
    - `exact_read_pairs: usize` (must match R1 and R2 identically)
    - `total_bases_r1: usize`
    - `total_bases_r2: usize`
    - `total_bases: usize` (`total_bases_r1 + total_bases_r2`)
    - `min_read_len: usize`
    - `max_read_len: usize`
    - `mean_read_len: f64`
    - `median_read_len: usize` (via small fixed-size frequency histogram array, e.g. `[usize; 2048]`, avoiding $O(N)$ heap allocations)
    - `elapsed_s: f64`
- **Integrity Validation:**
  - **Pairing Check:** $N_{R1} == N_{R2}$. If one file has more records than the other, raise an explicit fatal error (`FastqPairingMismatch`).
  - **Header Synchronization Check:** For sampled or all records, verify that the read identifier prefix (before `/1`, `/2` or whitespace) matches between R1 and R2.
  - **Truncation Check:** Any truncated 4-line record immediately aborts with an informative IO error.
- **Performance Requirement:**
  - High-throughput streaming via buffered gzip readers (`flate2::read::MultiGzDecoder` or `gzp`).
  - Streaming throughput $\ge 500,000\text{ reads/s}$ per core.
  - Constant $O(1)$ memory overhead ($< 10\text{ MB}$ RSS).

### 3.2 Complete Removal of Approximate Extrapolation
- Delete `approx_bytes_per_pair = 75`.
- Delete all calculations extrapolating read counts from `std::fs::metadata(path)?.len()`.
- The estimator must report `exact_read_pairs` and `total_bases` derived from literal stream consumption.

### 3.3 Safety Architecture & Assembly Invariance
- In Phase C1b-1, `subsampling_enabled` remains `false` by default (safe fallback).
- If `--enable-auto-subsample` is requested without `--genome-size`, the safety guard logs a warning and disables subsampling until Phase C1b-2/C1b-3 calibration is completed.
- Downstream assembly output on all benchmark datasets (`S1-COV-001`, `S2-REAL-001`, `R1-001`) must maintain bit-for-bit SHA-256 identity with `C1a_Default`.

---

## 4. Verification & Acceptance Criteria (Phase C1b-1 Gate)

1. **Count Accuracy Gate:** Across a gradient of 8 datasets spanning 15x to 300x, `fastbact` exact counts must match `zcat | wc -l / 4` and `seqtk` with **$0.00\%$ relative error**.
2. **Assembly Invariance Gate:** A/B assembly comparison against `C1a_Default` must yield **0 misassembly divergence and bit-for-bit identical contig checksums**.
3. **Cluster Execution Gate:** All builds, tests, and benchmarks executed strictly on Slurm compute node `vnode1` via `./scripts/hpc_run.sh`.
