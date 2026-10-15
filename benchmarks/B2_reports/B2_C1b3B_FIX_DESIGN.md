# FB-BENCH-B2-C1b-3B: Specialized Minimal Fix Design Specification

- **Date:** 2026-10-15
- **Branch:** `audit/b2-c1b-3a` (Commit: `f8be9ba`)
- **Baseline Git SHA:** `dac6a57417b390cf7385ff07adf569e550cc8ca9`
- **Execution Target:** HPC Slurm Node `vnode1`

---

## 1. Design Overview & Single-Variable Isolation Principle

This document specifies the experimental design for milestone **`FB-BENCH-B2-C1b-3B`**, directly resolving the gaps identified in the `B2-C1b-3A` audit:
1. Conflation of `--genome-size` with subsampling activation in `src/main.rs`.
2. Omission of explicit genome size propagation to high-k cDBG construction in `src/graph/local_multik.rs`.
3. 150X assembly regression due to premature subsampling across Poisson depth dips.
4. Count-Min Sketch capacity trade-offs under high sequencing depth.

In accordance with strict verification protocol, each of the four proposed interventions is formulated as an **isolated, orthogonal variable**. No combined multi-variable changes will be committed to production without individual verification against baseline `dac6a57`.

---

## 2. Independent Fix Designs

### Variable 1: CLI Parameter Decoupling (`--genome-size` vs. Subsampling)

#### Problem Diagnosis (`VERIFIED`)
In `src/main.rs` (lines 168–176):
```rust
let subsampling_enabled = cli.subsample || cli.genome_size.is_some();
```
Passing `--genome-size 4.65m` implicitly sets `subsampling_enabled = true`. If a user specifies expected genome size purely to optimize memory allocation or telemetry, reads may be discarded unexpectedly when depth exceeds $100\times$.

#### Minimal Fix Design
1. Decouple `subsampling_enabled` from `cli.genome_size`:
   ```rust
   // Only explicit --subsample flag or automated pipeline policy triggers read discarding
   let subsampling_enabled = cli.subsample && !cli.no_subsample;
   ```
2. Retain `cli.genome_size` strictly as an override for `SketchParams` and table capacity.

#### Verification & Metric Gates
- **Control Test:** Run 150X dataset with `--genome-size 200k` (without `--subsample`).
- **Acceptance Criteria:** `total_reads_ingested` must equal 100% of raw reads; `contigs.fa` SHA-256 must be bit-identical to the un-subsampled run (`N50 = 199,983 bp`).
- **Resource Impact:** Peak RSS unchanged ($\pm 0$ MB), wall time unchanged ($\pm 0.05$ s).

---

### Variable 2: Propagating Calibrated Genome Size to High-k Graph Construction

#### Problem Diagnosis (`VERIFIED`)
In `src/graph/local_multik.rs` (lines 82 and 108):
```rust
let mut g63 = build_cdbg_with_contigs::<[u64; 2]>(&active_bank_k63, &prior_contigs_k31, 63, Some(4));
```
`build_cdbg_with_contigs` passes `genome_size = None` to `build_cdbg_with_contigs_ext`. When `None` is passed:
```rust
let est_g = if let Some(g) = genome_size {
    g.max(50_000)
} else if !contigs.is_empty() {
    contigs.iter().map(|c| c.len() as u64).sum::<u64>().max(50_000)
} else { ... };
```
If prior contigs are fragmented (e.g., during complex graph resolution), `contigs.iter().map(...).sum()` can significantly underestimate the true genome size, causing the $k=63$ and $k=127$ CMS tables to be undersized, increasing collision probability.

#### Minimal Fix Design
1. Update `run_local_multik` and `assemble_pipeline` signatures to accept `trusted_genome_size: Option<u64>`.
2. Pass `trusted_genome_size` directly to `build_cdbg_with_contigs_ext` at all multi-k stages ($k=63, 127$).

#### Verification & Metric Gates
- **Control Test:** Run multi-k assembly on 2.5 Mb and 4.65 Mb datasets with and without explicit `trusted_genome_size`.
- **Acceptance Criteria:** Graph building memory conforms strictly to `CountMinSketch::for_genome_size(trusted_genome_size)`; 0 node/edge loss from collision-induced misallocations.
- **Resource Impact:** Peak RSS $+0$ to $+16$ MB depending on genome scale; runtime unchanged.

---

### Variable 3: Conservative Subsampling Activation Thresholds & Target Depth Policy

#### Problem Diagnosis (`VERIFIED`)
Depth tracing in `B2_150X_DEPTH_TRACE.tsv` demonstrated:
- At 150X, raw depth at the Poisson dip (coords 53,371–53,497 bp) is only **5X**.
- Subsampling with trigger $\ge 100\times$ and factor $0.8257$ randomly reduces depth in the dip to $2\times$–$4\times$.
- In $40\%$ of random seeds, depth drops below the $k=127$ solid threshold ($t_{\text{solid}} = 4$), fracturing the assembly into 146,614 bp and 53,490 bp contigs.
- In contrast, fastBact processes 150X un-subsampled reads in **0.40 seconds** using only **130 MB RAM**. Subsampling at 150X provides negligible performance benefit while incurring severe assembly risk.

#### Evaluation of Alternative Candidate Policies
- **Policy 3A (Conservative Trigger Threshold $\ge 180\times$):**
  - Subsampling triggers only when depth exceeds $180\times$ (or $200\times$).
  - Datasets $\le 180\times$ (including 150X) retain 100% of reads, preserving all Poisson dips.
  - Ultra-high depth datasets (300X) are subsampled to target $120\times$ (factor $0.41$).
- **Policy 3B (Target Depth Buffer to $150\times$):**
  - Retain trigger threshold at $\ge 120\times$, but set target coverage to $150\times$.
  - 150X data has factor $1.000$ (no reads dropped). 300X data has factor $0.500$ (target $150\times$, maintaining $\ge 8\times$ in dips).
- **Policy 3C (Two-Pass K-mer Protected Subsampling):**
  - Read pairs containing low-coverage k-mers ($\le 5\times$) are marked as immune from subsampling.
  - High complexity, requires an extra pre-pass scan.

#### Recommendation
**Policy 3A** is strictly recommended. It provides 100% assembly continuity on 150X across all seeds, adds 0 lines of algorithmic overhead, and saves CPU cycles.

---

### Variable 4: Graph Builder CMS Capacity Tuning vs. Exact Solid Verification

#### Problem Diagnosis (`VERIFIED`)
- Ground truth audit established that at 300X un-subsampled coverage, the Dynamic Standard CMS table ($w = 16 \times \widehat{G}$) incurs 1.76M collision false positives, yielding 38.05% precision.
- However, with graph simplification (`clip_tips`, `pop_bubbles`), 100% of collision spurs are pruned, producing a flawless 199,983 bp contig with 0 misassemblies.

#### Evaluation of Alternative Solutions
- **Option 4A (Double Capacity $32 \times \widehat{G}$):**
  - Increase CMS table width multiplier from $16\times$ to $32\times$.
  - Doubles CMS memory (e.g., from 8 MB to 16 MB for 200 kb; from 256 MB to 512 MB for 4.65 Mb).
  - Elevates 300X precision from $38.05\%$ to $74.75\%$, and 4.65 Mb precision from $92.27\%$ to $99.23\%$.
- **Option 4B (Exact Solid Hash Table for Small Genomes $<500\text{ kb}$):**
  - Use exact `FxHashMap<u64, u16>` for small genomes when memory is unconstrained.
  - Yields 100.00% precision and 100.00% recall.
- **Option 4C (Retain $16 \times \widehat{G}$ with Subsampling):**
  - Because automatic subsampling on $>180\times$ reduces read volume by $\sim 60\%$, distinct k-mers drop below table saturation limits, naturally restoring precision to $>90\%$ without increasing RAM.

#### Recommendation
Evaluate Option 4A vs Option 4C in the A/B matrix. Option 4A is safe for bacterial genomes as peak RSS on 4.65 Mb only increases from 1.62 GB to 1.88 GB, well within modern HPC 4GB/thread limits.

---

## 3. Systematic A/B Test Execution Plan

The experimental combinations are cataloged in `B2_C1b3B_AB_MATRIX.tsv`.
All runs will be executed via `./scripts/hpc_run.sh` on Slurm node `vnode1`.
Each test will record:
1. `N50` (bp) & `Contigs Count`
2. `Misassemblies` (against reference)
3. `127-mer Transition Continuity at 53497 bp`
4. `Peak RSS` (MB) & `Wall Time` (s)
5. `Bitwise SHA-256 Consistency` where applicable
