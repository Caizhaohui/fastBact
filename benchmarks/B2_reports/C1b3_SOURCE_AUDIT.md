# Source Code & Architecture Audit: FB-BENCH-B2-C1b-3

**Document ID:** FB-AUDIT-B2-C1b-3-SOURCE  
**Author:** AI Bioinformatic Systems Engineer  
**Date:** 2026-10-15  
**Audit Target:** CLI Argument Semantics, Dynamic CMS Fallback, Estimator Confidence Edge Cases, and Subsampling Policy Thresholds  
**Scope:** Phase C1b-3 Final Acceptance Gap Audit  
**Status:** Completed (`VERIFIED`)

---

## 1. Finding 1: `--genome-size` Triggers Subsampling Without `--enable-auto-subsample` (`VERIFIED`)

### Code Location:
[`src/main.rs:71-91`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/src/main.rs#L71-L91):
```rust
let (subsample_factor, subsampling_enabled, subsampling_reason) = if args.no_subsample {
    (1.0, false, "EXPLICIT_NO_SUBSAMPLE_FLAG".to_string())
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

### Analysis:
- `Cli::depth` defaults to `120.0` ([`src/cli.rs:31`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/src/cli.rs#L31)).
- If a user passes `--genome-size 4650000` on a 250x library *without* passing `--enable-auto-subsample`, `explicit_cov = 250.0 > 132.0`.
- FastBact **unconditionally downsamples** reads to ~120x, discarding ~50% of sequencing reads.
- **Root Cause:** In Phase C1a, explicit `--genome-size` was introduced as a manual escape hatch to allow safe downsampling when heuristic estimation was broken. But because `--depth` has a non-zero default, providing `--genome-size` implicitly acts as an opt-in to subsampling.
- **User Impact:** A user who passes `--genome-size` simply to assist Count-Min Sketch capacity allocation or inform assembly heuristics unexpectedly suffers read downsampling unless they remember to also supply `--no-subsample` or `--depth 0`.

### Minimal Recommendation:
In the next minor update, require explicit user intent for subsampling:
- Subsample only if `--enable-auto-subsample` is true, OR if `--depth` is explicitly provided alongside `--genome-size`.

---

## 2. Finding 2: `build_cdbg_with_contigs_ext` Fallback Capacity Formula (`VERIFIED`)

### Code Location:
[`src/graph/build.rs:62-74`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/src/graph/build.rs#L62-L74):
```rust
pub fn build_cdbg_with_contigs_ext<K: KmerInt>(
    bank: &ReadBank,
    contigs: &[Vec<u8>],
    k: usize,
    target_solid_cov: Option<u8>,
    genome_size: Option<u64>,
) -> CdbgGraph {
    let est_g = if let Some(g) = genome_size {
        g.max(50_000)
    } else if !contigs.is_empty() {
        contigs.iter().map(|c| c.len() as u64).sum::<u64>().max(50_000)
    } else {
        (bank.total_bases() as u64 / 60).max(50_000)
    };
    let cms = CountMinSketch::for_genome_size(est_g);
```

### Analysis:
- In production (`src/main.rs`), `build_cdbg_ext` is invoked with `Some(effective_g)`, so `est_g` correctly receives the calibrated genome size.
- However, the legacy helper function:
  ```rust
  pub fn build_cdbg<K: KmerInt>(bank: &ReadBank, k: usize, target_solid_cov: Option<u8>) -> CdbgGraph {
      build_cdbg_with_contigs_ext::<K>(bank, &[], k, target_solid_cov, None)
  }
  ```
  passes `genome_size = None`.
- If an external integration or unit test calls `build_cdbg`, it falls back to `bank.total_bases() / 60`.
- **User Impact:** Zero impact on production CLI (which uses `build_cdbg_ext`), but creates internal API duality.

### Minimal Recommendation:
- Retain `(bank.total_bases() / 60)` as a fallback when no estimate is available.
- Mark `build_cdbg` with `#[deprecated(note = "Use build_cdbg_ext with explicit genome size")]`.

---

## 3. Finding 3: EstimatorConfidence Edge Cases & False-Positive Risk (`VERIFIED`)

### Code Location:
[`src/ingest/estimator.rs:406-413`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/src/ingest/estimator.rs#L406-L413):
```rust
let confidence = if d_k < 5.0 || solid_ratio < 0.05 || u_ge_2 < 1_000 {
    EstimatorConfidence::Unreliable
} else if d_k < 20.0 || solid_ratio < 0.15 {
    EstimatorConfidence::LowCoverageWarning
} else {
    EstimatorConfidence::High
};
```

### Edge Case Analysis:
1. **Low-depth Interception Gate ($D_k < 20\text{x}$):**
   Confirmed to correctly protect `S1-COV-001` (15x, $D_k = 12.7\text{x}$), safely preventing any downsampling.
2. **High-depth Multi-strain Contamination:**
   If a sequencing library contains a mixture of two bacterial species (e.g., 50% *E. coli* 4.6 Mb and 50% *K. pneumoniae* 5.3 Mb), HLL observes distinct solid k-mers from both genomes ($U_{\ge 2} \approx 9.9\text{ Mb}$).
   Estimated genome size becomes $\approx 9.9\text{ Mb}$, which *underestimates* library depth as $\sim 50\text{x}$ instead of $100\text{x}$.
   **Outcome:** Auto-subsampling is safely *bypassed* (retaining 100% data). Contamination does not cause aggressive data loss.
3. **High-error Sequencing Chemistry:**
   If sequencing error rate rises above 2%, singleton error k-mers explode. However, because $U_{\ge 2}$ uses an idempotent insertion filter ($c \ge 2$), error k-mers do not inflate $U_{\ge 2}$, while $U_{\ge 1}$ expands, driving `solid_ratio` below 0.15.
   **Outcome:** Successfully triggers `LowCoverageWarning` or `Unreliable`, intercepting downsampling.
4. **Verdict on Estimator Confidence:**
   The state machine is conservative and directionally safe: errors in estimation consistently fail toward **keeping more data**, never toward accidental catastrophic deletion.

---

## 4. Finding 4: Conservative Policy Recommendation for Moderate-High Depth (130x–170x) (`VERIFIED`)

### Empirical Evidence from Forensics Report (`C1b3_150X_REGRESSION.md`):
- At 150x depth (`S1-COV-150X`), downsampling to 120x:
  - Runtime reduction: **0.48s** (4.22s -> 3.74s, negligible).
  - Memory reduction: **23.9 MB** (109.7 MB -> 85.8 MB, negligible).
  - Assembly regression: N50 drops from **199,983 bp** to **146,614 bp** in 40% of random seeds due to Poisson dip depletion at $k=127$ ($t_{\text{solid}} = 4$).
- At 300x depth (`S1-COV-003`), downsampling to 120x:
  - Runtime reduction: **4.45s** (8.06s -> 3.61s, 2.23x speedup).
  - Memory reduction: **153.8 MB** (241.2 MB -> 87.4 MB, 63.8% reduction).
  - Assembly quality: **100% preserved** across all 10 tested seeds (N50 = 199,993 bp, 0 misassemblies).
- On Real Isolate `R1-003` (~132x nominal, 263.7x total library depth):
  - Downsampling from 263.7x to 120x saved **105.7s** (48% speedup), cut RAM by **467 MB**, improved N50 by **+7.7 kb**, and eliminated **2 misassemblies** caused by depth-induced chimeric repeat collisions.

### Recommendation:
Set the auto-subsampling activation threshold to a more conservative margin:
- **Option A (Higher Activation Threshold):** Only trigger downsampling when estimated coverage $> 1.5 \times \text{target\_cov}$ ($> 180\text{x}$ for default 120x target). Libraries between 120x and 180x are kept at 100% data retention.
- **Option B (Raised Target Depth):** Increase default target depth to $150\text{x}$ (subsampling only libraries $> 165\text{x}$).
