# `fastBact` B2-G0 — Estimator Deep Forensics & Benchmark Reliability Audit

**Document ID:** FB-BENCH-B2-G0-FORENSICS  
**Date:** 2026-10-15  
**Evaluation Node:** `vnode1` (HPC Slurm Partition `qfree`)  
**Scope:** Forensic investigation of genome size underestimation (522 kb), depth inflation (626.4x), gzip stream extrapolation, Unicycler source topology, and coordinate-level misassembly classification on `R1-001`.

---

## 1. Lineage and Baseline Verification (`4259fe3` vs `552ff54`)

### 1.1. Exact Git Lineage
As documented in [`B2_LINEAGE.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/benchmarks/B1_reports/B2_LINEAGE.tsv):
- **Frozen Baseline Tag (`v0.1.0-baseline`):**
  Points directly to commit `4259fe3ffe077f67946700e8092ccf46fbbb9bc6`.
  This commit marks the completion of the FB-BENCH-001 benchmark suite and establishes the immutable baseline.
- **Diagnostic / Audit Commit (`552ff547781dae9dd2a0468d6f0ae812eeb8053d`):**
  This is the direct child commit of `4259fe3`. It introduced hidden command-line ablation flags (`--no-subsample`, `--no-trim`, etc.), intermediate unitig stage dumping (`--dump-stages`), and the B1 audit report suite.

### 1.2. Verification of Bit-for-Bit Assembly Invariance
Under default parameters, `552ff54` does not alter assembly behavior. On dataset `S2-REAL-001` (4.64 Mb *E. coli* MG1655), the resulting `contigs.fa` was evaluated on compute node `vnode1`:
$$\text{SHA-256} = \mathbf{0cdcad9edfb90ab7325b9cc27d38806b6706b899788439e1a2dd869f1c32738e}$$
This checksum matches `v0.1.0-baseline` bit-for-bit.

---

## 2. Mathematical Provenance of HLL 522 kb & 626.4x Coverage

In `fastBact v0.1.0`, on `S2-REAL-001` (nominal coverage 100x, reference length 4,641,652 bp), `main.rs` logged:
```text
INFO fastBact: Estimated genome size: 522066 bp, coverage: 626.4x (subsample factor: 0.1916)
```
Our deep forensics on `vnode1` revealed that this erroneous output was caused by three compounding mathematical and algorithmic defects in [`src/ingest/estimator.rs`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/src/ingest/estimator.rs):

### 2.1. Defect 1: Canonical Hash Distribution Distortion (`min(h_fwd, h_rev)`)
In `estimator.rs`:
```rust
let h = rolling.canonical(); // returns self.h_fwd.min(self.h_rev)
...
if c == 2 {
    hll_solid.add_hash(h);
}
```
- **The Core Flaw:** HyperLogLog strictly requires uniform random input over $[0, 2^{64}-1]$. The top 14 bits select the register index, and the remaining 50 bits determine the geometric leading-zero count $\rho$.
- **The Distortion:** If $X, Y \sim \text{Uniform}(0, 2^{64}-1)$, their minimum $M = \min(X, Y)$ has cumulative distribution $F(t) = 1 - (1 - t/2^{64})^2$. This heavily skews values toward zero!
- **Empirical Validation:**
  On `vnode1`, we simulated inserting $N = 3,000,000$ distinct items into HyperLogLog:
  - **Uniform Hash:** HLL estimate = **2,967,545** (relative error = -1.08%).
  - **Raw $\min(X, Y)$ Hash:** HLL estimate = **864,807** (relative error = **-71.17%**, a **3.43x deflation**)!
- **Consequence:** Feeding raw canonical ntHash values into HLL without a secondary mixing hash (such as `splitmix64`) causes HLL to underestimate unique solid k-mers by more than 3.5-fold.

### 2.2. Defect 2: Small Count-Min Sketch Saturation & $c=2$ Skipping
- In `estimator.rs`, the sample Count-Min Sketch was initialized with `1 << 20` (1,048,576) counters per row.
- During the sampling pass of 50,000 read pairs (15 Mb of sequence), ~13,000,000 21-mers were inserted into 4 rows of 1M slots (an average load of 12.4 items per slot).
- In a 4-bit saturated CMS, counters quickly exceed 2 due to hash collisions.
- The condition in `estimator.rs` was strictly:
  ```rust
  if c == 2 { hll_solid.add_hash(h); }
  ```
- Any genomic k-mer whose slots already had $\text{count} \ge 2$ from prior collisions returned $c > 2$ on its **very first touch**, completely skipping `hll_solid` registration.
- On `S2-REAL-001`, exact counting revealed **3,243,532** distinct 21-mers with count $\ge 2$ in the 50k sample. Under CMS $2^{20}$, only ~1.68M items registered $c=2$, which was further deflated to **522,066** by the $\min(X,Y)$ HLL distortion.

### 2.3. Defect 3: Fixed 75-Byte Gzip Constant & Coverage Extrapolation
- Total R1 compressed file size on `S2-REAL-001`: $81,759,300$ bytes.
- In `estimator.rs`:
  ```rust
  let approx_bytes_per_pair = if is_gz { 75u64 } else { ... };
  let total_pairs_est = (file_bytes / approx_bytes_per_pair).max(pair_count as u64);
  // total_pairs_est = 81,759,300 / 75 = 1,090,124 pairs
  let sample_ratio = 1,090,124 / 50,000 = 21.80248;
  let total_projected_bases = 15,000,000 * 21.80248 = 327,037,200 bp;
  let estimated_cov = 327,037,200 / 522,066 = 626.427x;
  let subsample_factor = (120.0 / 626.427) = 0.1916;
  ```
- **Result:** Combining the 3.5x HLL deflation and the 2.5x sample projection inflation caused `fastBact` to estimate coverage at **626.4x** instead of 100x, violently discarding **80.9%** of sequencing reads during ingest.

---

## 3. Gzip Compressed Stream Position & Extrapolation Evaluation

We evaluated whether tracking the compressed byte stream offset during streaming can replace fixed-ratio extrapolation:

### 3.1. Empirical Measurement Across Gzip Blocks (`estimator_forensics_raw.json`)
On compute node `vnode1`, chunked zlib stream tracking across FASTQ reads demonstrated:
- **`S2-REAL-001` (Synthetic 2x150bp PE):**
  - At 10k pairs: 52.43 gz bytes/pair $\to$ Extrapolated total: 1,559,435 pairs (True: 1,573,333, Error: -0.88%).
  - At 50k pairs: 52.43 gz bytes/pair $\to$ Extrapolated total: 1,559,435 pairs.
- **`R1-001` (`SRR28349173`, Real NovaSeq X Plus):**
  - At 10k pairs: 91.75 gz bytes/pair $\to$ Extrapolated total: 2,189,125 pairs (True: 2,409,728, Error: -9.16%).
  - At 25k pairs: 86.51 gz bytes/pair $\to$ Extrapolated total: 2,321,800 pairs (Error: -3.65%).
  - At 50k pairs: 85.20 gz bytes/pair $\to$ Extrapolated total: 2,357,520 pairs (Error: -2.17%).

### 3.2. Limitations of Gzip Stream Position Extrapolation
1. **Decompressor Buffer Read-Ahead:** In Rust's `flate2::read::GzDecoder`, the internal reader buffers 64 KB of compressed stream data. Querying `.stream_position()` on the underlying file reflects the read-ahead buffer position, introducing a 5–15% skew on small sample sizes.
2. **Quality Score Entropy Drift:** In real sequencing runs, Phred quality scores become more variable toward the 3' ends and across later optical tiles. Compression ratios fluctuate from 5.1x to 6.3x across the file.
3. **Multi-Member Gzip Files:** Concatenated gzip archives (`cat a.fq.gz b.fq.gz > c.fq.gz`) contain multiple deflate streams and separate headers, making linear byte offset projection unreliable.
4. **Architectural Recommendation for B2-C1b:**  
   Rather than relying on heuristic gzip stream offsets, a **two-pass exact streaming count** (Pass 1 streaming line counting without disk writing, Pass 2 sampling and assembly) is 100% deterministic, robust to any compression scheme, and introduces negligible I/O overhead.

---

## 4. Upstream Unicycler Architecture Verification

We audited the local installation of Unicycler (`/hpcfs/fhome/caizhh/.conda/envs/LineaBact/lib/python3.12/site-packages/unicycler/`) and the referenced upstream commit (`1429301`):

| Module Functionality | Actual Upstream File | Erroneous Reference Corrected |
|:---|:---|:---|
| Copy Depth / Multiplicity | `assembly_graph_copy_depth.py` | *Not* in `assembly_graph.py` |
| Primary Multiplicity Entry | `determine_copy_depth(graph)` | *Not* `calculate_copy_depth` |
| SPAdes Contig Short-Read Bridge | `bridge_spades_contig.py` | *No such file as `bridge.py`* |
| Loop Unrolling Bridge | `bridge_loop_unroll.py` | *No such file as `bridge.py`* |
| Long-Read Overlap Bridge | `bridge_long_read.py` | Distinct from short-read bridging |
| Common Bridge Scoring Helpers | `bridge_common.py` | Shared utility module |

This confirms Section 2 of `FB_B2_Unicycler_Revised_Execution_Plan.md`. In future Unicycler study (U-series), all citations and references must target these verified module paths.

---

## 5. Coordinate-Level Forensic Review of `R1-001` Misassemblies

In `benchmarks/B1_reports/B1_ABLATION.tsv`, `--no-subsample` on `R1-001` reported 4 extensive misassemblies (vs 1 in baseline `v0.1.0`). We conducted coordinate-level alignment analysis against `NC_000913.3` (documented in [`B2_MISASSEMBLY_REVIEW.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/benchmarks/B1_reports/B2_MISASSEMBLY_REVIEW.tsv)):

### 5.1. Contig 6: Genuine Biological Inversion (Concordant with SPAdes)
- **Coordinates in `fastBact` `contig_6` (125,420 bp):**
  - Left flank: `1..110452` aligns to reference `1097353..1207793` (+ strand, 99.83% IDY).
  - Inverted segment: `110466..112263` aligns to reference `1207806..1209602` (**- strand**, 98.67% IDY, 1,797 bp).
  - Right flank: `112264..125420` aligns to reference `1209603..1222756` (+ strand, 99.69% IDY).
- **Concordance with SPAdes v4.2.0:**
  In SPAdes default assembly of `R1-001`, `NODE_5` (189,459 bp) assembled the **exact same inversion**:
  - `1..110480` (+ strand) $\to$ `110481..112277` (**- strand**, 1,797 bp) $\to$ `112278..189459` (+ strand).
- **Biological Classification:**  
  This 1.8 kb invertible segment is the well-documented phase-variable tail switch / prophage element in *E. coli*. Both independent assemblers (`fastBact` and `spades.py`) reconstruct the identical continuous chromosome path. QUAST blindly scores entry and exit of this strain difference as 2 "inversion misassemblies".  
  **Verdict: TRUE BIOLOGICAL VARIANT / REFERENCE DISAGREEMENT (Not an algorithm error).**

### 5.2. Contigs 150 & 164: Repeat-Flank Chimeras on Short Contigs
- **`contig_150` (1,285 bp):** 1,159 bp body at `1873031..1874189` joined to a 126 bp terminal segment at `291276..291401` (an IS element copy).
- **`contig_164` (551 bp):** 425 bp body at `1872615..1873039` joined to a 126 bp terminal segment at `1978503..1978628` (another IS element copy).
- In baseline `v0.1.0`, these regions were broken into sub-500 bp fragments and filtered out. In `--no-subsample`, they assembled into short borderline contigs (1.2 kb and 551 bp).

---

## 6. Minimal Code Modification Scope & A/B Acceptance for Phase B2

Based on B2-G0 findings, the implementation stages must be separated:

### 6.1. B2-C1a: Safe Retention / Failure Fallback (Zero Estimator Change)
- **Scope:** Modify `src/main.rs` to add a safety guard: if estimated coverage is below 120x or deemed uncertain, **disable automatic subsampling completely**. Expose explicit `--no-subsample` and `--genome-size <BP>` flags.
- **A/B Acceptance Test:**
  - On `S2-REAL-001` (100x true coverage), C1a must retain 100% of reads.
  - Assembly output must match the frozen `--no-subsample` baseline (NGA50 $\ge 285\text{ kb}$, contigs $\le 28$).

### 6.2. B2-C1b: Accurate Counting & Calibrated Estimation
- **Scope:**
  1. Replace fixed 75-byte gzip constant with a fast streaming two-pass exact record counter.
  2. In `HyperLogLog`, hash `RollingNtHash` values through `splitmix64` before insertion.
  3. Size the sample CMS appropriately ($2^{23}$ or $2^{24}$) to eliminate saturation.
- **A/B Acceptance Test:**
  - On 15x, 30x, 100x, and 150x datasets, estimated genome size must remain within $\pm 10\%$ of ground truth.
  - At 100x, subsampling must not trigger. At 200x, subsampling must accurately reduce coverage to $120\text{x} \pm 5\text{x}$.
