# FB-BENCH-B2-C1b-3B: Specialized Minimal Fix Design Specification (Corrected)

- **Date:** 2026-10-15
- **Branch:** `audit/b2-c1b-3a`
- **Baseline Git SHA:** `dac6a57417b390cf7385ff07adf569e550cc8ca9`
- **Execution Target:** HPC Slurm Node `vnode1`

---

## 1. Design Overview & Single-Variable Isolation Principle

This document specifies the technical design for milestone **`FB-BENCH-B2-C1b-3B`**, directly resolving the gaps identified in the `B2-C1b-3A-R` audit:
1. Conflation of `args.genome_size` with subsampling activation in `src/main.rs`.
2. Omission of explicit genome size propagation to high-k cDBG construction in `upgrade_multik()`.
3. 150X assembly regression due to premature subsampling across Poisson depth dips.
4. Count-Min Sketch capacity trade-offs under high sequencing depth.

In accordance with strict verification protocol, each of the four proposed interventions is formulated as an **isolated, orthogonal variable** using **actual, verified production code interfaces**:
- CLI fields: `args.enable_auto_subsample: bool`, `args.genome_size: Option<u64>`, `args.depth: f64`, `args.no_subsample: bool`.
- Graph functions: `upgrade_multik<K2: KmerInt>(&mut CdbgGraph, &ReadBank, usize)` and `build_cdbg_with_contigs_ext<K: KmerInt>(&ReadBank, &[Vec<u8>], usize, Option<u8>, Option<u64>)`.

---

## 2. Independent Fix Designs

### Variable 1: CLI Parameter Decoupling (`args.genome_size` vs. Subsampling)

#### Problem Diagnosis (`VERIFIED`)
In `src/main.rs` (lines 75–95):
```rust
} else if let Some(explicit_g) = args.genome_size {
    if explicit_g == 0 {
        (1.0, false, "EXPLICIT_GENOME_SIZE_ZERO_INVALID".to_string())
    } else {
        let explicit_cov = estimate.total_bases as f64 / (explicit_g as f64);
        if args.depth > 0.0 && explicit_cov > 1.1 * args.depth {
            let factor = (args.depth / explicit_cov).clamp(0.01, 1.0);
            (
                factor,
                true,
                format!("EXPLICIT_GENOME_SIZE_CAP (cov={:.1}x > target={:.1}x, factor={:.4})", explicit_cov, args.depth, factor),
            )
        } else { ... }
    }
}
```
If a user specifies `--genome-size 4.65m` solely to configure memory sketches or progress tracking, `explicit_cov > 1.1 * args.depth` (i.e. $> 132.0\times$) automatically activates read subsampling and drops reads, even when `--enable-auto-subsample` was NOT set.

#### Minimal Fix Design
1. Decouple `args.genome_size` from read discarding:
   ```rust
   // Only explicit --enable-auto-subsample flag triggers read subsampling
   let (subsample_factor, subsampling_enabled, subsampling_reason) = if args.no_subsample {
       (1.0, false, "EXPLICIT_NO_SUBSAMPLE_FLAG".to_string())
   } else if args.enable_auto_subsample {
       // Evaluate confidence and est_cov / explicit_cov against target
       ...
   } else {
       // Safe Default: 100% data retention
       (1.0, false, "DISABLED_BY_DEFAULT_SAFE_FALLBACK".to_string())
   };
   ```
2. Retain `args.genome_size` strictly as an override for `effective_g` (used in `CountMinSketch::for_genome_size(effective_g)`).

#### Verification & Metric Gates
- **Control Test:** Run 150X dataset with `--genome-size 200000` (without `--enable-auto-subsample`).
- **Baseline Measured Value:** 150X with `--genome-size 200000` dropped reads (factor $0.8257$), causing 2 contigs (N50 = 146,614 bp).
- **Experimental Expectation:** Retains 100% of reads ($N50 = 199,983\text{ bp}$, 1 contig).
- **Resource Impact:** Peak RSS $\pm 0$ MB, wall time $\pm 0.02$ s.

---

### Variable 2: Propagating Calibrated Genome Size to High-k Graph Construction

#### Problem Diagnosis (`VERIFIED`)
In `src/graph/local_multik.rs` (line 50):
```rust
let mut new_graph = build_cdbg_with_contigs::<K2>(&active_bank, &prior_contigs, k2, Some(4));
```
`build_cdbg_with_contigs` calls `build_cdbg_with_contigs_ext` with `genome_size = None`. In `build.rs`:
```rust
let est_g = if let Some(g) = genome_size {
    g.max(50_000)
} else if !contigs.is_empty() {
    contigs.iter().map(|c| c.len() as u64).sum::<u64>().max(50_000)
} else { ... };
```
If prior contigs are fragmented, `contigs.iter().map(...).sum()` underestimates true genome size, undersizing the CMS table for $k=63$ and $k=127$.

#### Minimal Fix Design
1. Provide `upgrade_multik_ext<K2: KmerInt>(prev_graph: &mut CdbgGraph, bank: &ReadBank, k2: usize, trusted_genome_size: Option<u64>)`.
2. Pass `trusted_genome_size` directly to `build_cdbg_with_contigs_ext` at all multi-k stages ($k=63, 127$).

#### Verification & Metric Gates
- **Control Test:** Run multi-k assembly on 2.5 Mb and 4.65 Mb datasets with and without explicit `trusted_genome_size`.
- **Baseline Measured Value:** S2-REAL-001 $k=63, 127$ sketches sized from prior contigs ($\sim 4.6\text{ Mb}$).
- **Experimental Expectation:** Guaranteed $256\text{ MB}$ CMS allocation; zero collision-induced edge dropout in high-k graphs.
- **Resource Impact:** Peak RSS $+0$ to $+16$ MB; runtime unchanged.

---

### Variable 3: Conservative Subsampling Activation Thresholds & Target Depth Policy

#### Problem Diagnosis (`VERIFIED`)
In `src/main.rs`:
```rust
if args.depth > 0.0 && est_cov > 1.1 * args.depth { ... }
```
With default `--depth 120.0`, the trigger threshold is $1.1 \times 120 = 132.0\times$.
For `S1-COV-150X` (estimated coverage $\sim 145.3\times$), subsampling triggers with factor $120 / 145.3 = 0.8257$.
However:
- In `S1-COV-150X`, the Poisson dip at coordinates 53,371–53,375 bp has raw depth of only **5X**.
- Subsampling with factor $0.8257$ drops depth across the 5 consecutive 127-mers to $2\times$–$3\times$.
- Because $2 < t_{\text{solid}} = 4$, the junction breaks into 146,614 bp and 53,490 bp contigs in 40% of random seeds.
- Subsampling 150X saves only $0.05$ s of runtime and $<20\text{ MB}$ RAM, but causes severe assembly fragmentation.

#### Evaluation of Alternative Candidate Policies
- **Policy 3A (Conservative Trigger Threshold $\ge 1.5 \times \text{target}$ / $180\times$):**
  - Trigger condition updated to `est_cov > 1.5 * args.depth` (i.e. $> 180.0\times$).
  - Datasets $\le 180\times$ (such as 150X) are fully assembled with 100% of reads, preserving Poisson dips.
  - Ultra-high depth datasets (300X) trigger subsampling to target $120\times$.
- **Policy 3B (Target Depth Buffer to $150\times$):**
  - Retain trigger condition at $1.1 \times \text{target}$, but set default target depth to $150\times$.
  - 150X data has factor $1.000$ (no reads dropped). 300X subsamples to $150\times$ (retaining $\ge 7\times$ in dips).

#### Recommendation
**Policy 3A** (`est_cov > 1.5 * args.depth`) is strictly recommended.

---

### Variable 4: Graph Builder CMS Capacity Tuning vs. Exact Solid Verification

#### Problem Diagnosis (`VERIFIED`)
- Ground truth audit established that at 300X un-subsampled coverage, the Dynamic Standard CMS table ($w = 16 \times \widehat{G}$) incurs 1.76M collision false positives, yielding 38.05% precision.
- However, graph simplification prunes 100% of collision spurs, producing 0 misassemblies.

#### Evaluation of Alternative Solutions
- **Option 4A (Double Capacity $32 \times \widehat{G}$):**
  - Doubles CMS memory (16 MB for 200 kb; 512 MB for 4.65 Mb).
  - Elevates 300X precision from $38.05\%$ to $74.75\%$, and 4.65 Mb precision from $92.27\%$ to $99.23\%$.
- **Option 4C (Retain $16 \times \widehat{G}$ with Subsampling):**
  - Automatic subsampling on $>180\times$ reduces read volume by $\sim 60\%$, naturally restoring precision to $>90\%$ without increasing RAM.

---

## 3. Reference Genome Reality Check (S2-REAL-001)

| Metric | S1-COV-150X (200 kb Synthetic) | S2-REAL-001 (4.65 Mb Real E. coli) |
| :--- | :--- | :--- |
| **True Genome Size** | 200,000 bp | 4,655,000 bp |
| **Known Repeat Complexity** | 0 repeats $> 60\text{ bp}$ | **Seven 5kb `rrn` ribosomal operon repeats** |
| **Baseline Contigs Count** | 1 (un-subsampled) / 2 (seed 0 subsampled) | **2 to 3 contigs** (bounded by 5kb `rrn` repeats) |
| **Baseline N50** | 199,983 bp | **4,082,514 bp** (~4.08 Mb) |
| **Baseline NGA50** | 199,983 bp | **2,897,224 bp** (~2.90 Mb) |
| **Baseline Misassemblies** | 0 | **1** (5kb `rrn` repeat inversion documented in `B2_MISASSEMBLY_REVIEW.tsv`) |
| **Target for B2-C1b-3B** | N50 = 199,983 bp across 100% of seeds | N50 $\ge 4.08\text{ Mb}$, NGA50 $\ge 2.89\text{ Mb}$, 0 regression |
