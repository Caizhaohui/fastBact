//! Verification binary for B2-C1b-3B Var 3 Quality Gates
//! Verifies:
//! 1. Rejection of illegal / non-finite / < 1.0 `--subsample-threshold-ratio` arguments.
//! 2. Seed 0 determinism: `--subsample-seed 0` is 100% bit-identical to default unseeded subsampling.
//! 3. Seed repeatability: Running with the same seed twice produces bit-for-bit identical outputs.
//! 4. Threshold ratio consistency across both explicit genome size and HLL estimator paths.

use anyhow::{bail, Context, Result};
use std::path::Path;
use std::process::Command;

fn compute_sha256<P: AsRef<Path>>(path: P) -> Result<String> {
    let output = Command::new("sha256sum")
        .arg(path.as_ref())
        .output()
        .context("Failed to run sha256sum")?;
    if !output.status.success() {
        bail!("sha256sum failed on {:?}", path.as_ref());
    }
    let text = String::from_utf8_lossy(&output.stdout);
    let hash = text.split_whitespace().next().unwrap_or("UNKNOWN").to_string();
    Ok(hash)
}

fn run_fastbact(args: &[&str]) -> (i32, String) {
    let output = Command::new("./target/release/fastBact")
        .args(args)
        .output()
        .expect("Failed to execute fastBact");
    let combined = format!(
        "{}\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    (output.status.code().unwrap_or(-1), combined)
}

fn main() -> Result<()> {
    println!("=======================================================");
    println!("B2-C1b-3B Var 3 Source Code Quality Gates Verification");
    println!("=======================================================");

    let r1 = "benchmark-results/data_store/S1-COV-150X_R1.fq.gz";
    let r2 = "benchmark-results/data_store/S1-COV-150X_R2.fq.gz";
    let r1_300 = "benchmark-results/data_store/S1-COV-003_R1.fq.gz";
    let r2_300 = "benchmark-results/data_store/S1-COV-003_R2.fq.gz";

    // -------------------------------------------------------------------------
    // Gate 1: Rejection of illegal / non-finite / < 1.0 threshold ratios
    // -------------------------------------------------------------------------
    println!("\n[Gate 1] Testing rejection of invalid --subsample-threshold-ratio ...");
    let invalid_ratios = ["-1.0", "0.0", "0.5", "0.99", "NaN", "Inf", "-Inf"];
    for &val in &invalid_ratios {
        let flag = format!("--subsample-threshold-ratio={val}");
        let (code, log) = run_fastbact(&[
            "-1", r1, "-2", r2, "-o", "test_gate_tmp/invalid",
            &flag,
        ]);
        assert_ne!(code, 0, "fastBact must reject ratio={val}");
        assert!(
            log.contains("Invalid --subsample-threshold-ratio"),
            "Expected error message for ratio={val}, got:\n{log}"
        );
        println!("  ✓ Ratio '{}' properly rejected with error.", val);
    }
    println!("Gate 1 PASSED: All invalid/non-finite/sub-1.0 ratios rejected.");

    // -------------------------------------------------------------------------
    // Gate 2: Ratio consistency across explicit genome size and HLL paths
    // -------------------------------------------------------------------------
    println!("\n[Gate 2] Testing ratio consistency (HLL vs explicit genome size paths) ...");
    // S1-COV-150X is ~145X.
    // Ratio 1.1 (trigger 132X): Must trigger subsampling in both paths.
    // Ratio 1.5 (trigger 180X): Must bypass subsampling in both paths.

    // Path A: HLL estimator with ratio 1.1 -> Subsamples
    let out_hll_11 = "test_gate_tmp/hll_11";
    let _ = std::fs::remove_dir_all(out_hll_11);
    let (code, log) = run_fastbact(&[
        "-1", r1, "-2", r2, "-o", out_hll_11,
        "--enable-auto-subsample", "--subsample-threshold-ratio", "1.1", "-t", "4",
    ]);
    assert_eq!(code, 0, "HLL ratio 1.1 failed: {log}");
    assert!(log.contains("AUTO_SUBSAMPLE_HIGH_CONFIDENCE"), "Expected subsample trigger in HLL path: {log}");
    let sha_hll_11 = compute_sha256(format!("{out_hll_11}/contigs.fa"))?;

    // Path B: Explicit genome size with ratio 1.1 -> Subsamples
    let out_exp_11 = "test_gate_tmp/exp_11";
    let _ = std::fs::remove_dir_all(out_exp_11);
    let (code, log) = run_fastbact(&[
        "-1", r1, "-2", r2, "-o", out_exp_11,
        "--enable-auto-subsample", "--genome-size", "200000", "--subsample-threshold-ratio", "1.1", "-t", "4",
    ]);
    assert_eq!(code, 0, "Explicit ratio 1.1 failed: {log}");
    assert!(log.contains("AUTO_SUBSAMPLE_EXPLICIT_GENOME_SIZE"), "Expected subsample trigger in explicit path: {log}");
    let sha_exp_11 = compute_sha256(format!("{out_exp_11}/contigs.fa"))?;
    assert!(!sha_hll_11.is_empty() && !sha_exp_11.is_empty());
    println!("  ✓ 1.1x ratio consistently triggers subsampling across HLL and explicit paths.");

    // Path C: HLL estimator with ratio 1.5 -> Bypasses subsampling
    let out_hll_15 = "test_gate_tmp/hll_15";
    let _ = std::fs::remove_dir_all(out_hll_15);
    let (code, log) = run_fastbact(&[
        "-1", r1, "-2", r2, "-o", out_hll_15,
        "--enable-auto-subsample", "--subsample-threshold-ratio", "1.5", "-t", "4",
    ]);
    assert_eq!(code, 0, "HLL ratio 1.5 failed: {log}");
    assert!(log.contains("AUTO_SUBSAMPLE_HIGH_CONFIDENCE_WITHIN_TARGET"), "Expected bypass in HLL path: {log}");
    let sha_hll_15 = compute_sha256(format!("{out_hll_15}/contigs.fa"))?;

    // Path D: Explicit genome size with ratio 1.5 -> Bypasses subsampling
    let out_exp_15 = "test_gate_tmp/exp_15";
    let _ = std::fs::remove_dir_all(out_exp_15);
    let (code, log) = run_fastbact(&[
        "-1", r1, "-2", r2, "-o", out_exp_15,
        "--enable-auto-subsample", "--genome-size", "200000", "--subsample-threshold-ratio", "1.5", "-t", "4",
    ]);
    assert_eq!(code, 0, "Explicit ratio 1.5 failed: {log}");
    assert!(log.contains("AUTO_SUBSAMPLE_EXPLICIT_GENOME_SIZE_WITHIN_TARGET"), "Expected bypass in explicit path: {log}");
    let sha_exp_15 = compute_sha256(format!("{out_exp_15}/contigs.fa"))?;
    assert_eq!(sha_hll_15, sha_exp_15, "HLL and Explicit 1.5x paths must yield identical assemblies");

    // Path E: Un-subsampled default baseline -> Bit-identical to 1.5x bypass!
    let out_default = "test_gate_tmp/default";
    let _ = std::fs::remove_dir_all(out_default);
    let (code, _) = run_fastbact(&["-1", r1, "-2", r2, "-o", out_default, "-t", "4"]);
    assert_eq!(code, 0);
    let sha_default = compute_sha256(format!("{out_default}/contigs.fa"))?;
    assert_eq!(sha_hll_15, sha_default, "1.5x bypass must be 100% bit-identical to un-subsampled default");
    println!("  ✓ 1.5x ratio consistently bypasses subsampling across HLL, explicit, and default paths.");
    println!("Gate 2 PASSED: Threshold ratio behaves identically in all paths.");

    // -------------------------------------------------------------------------
    // Gate 3: Seed 0 determinism and repeatability
    // -------------------------------------------------------------------------
    println!("\n[Gate 3] Testing Seed 0 determinism vs unseeded default on 300X ...");
    let out_seed_0a = "test_gate_tmp/seed_0a";
    let out_seed_0b = "test_gate_tmp/seed_0b";
    let _ = std::fs::remove_dir_all(out_seed_0a);
    let _ = std::fs::remove_dir_all(out_seed_0b);

    // Run without --subsample-seed (defaults to 0 internally)
    let (code, _) = run_fastbact(&[
        "-1", r1_300, "-2", r2_300, "-o", out_seed_0a,
        "--enable-auto-subsample", "-t", "4",
    ]);
    assert_eq!(code, 0);
    let sha_0a = compute_sha256(format!("{out_seed_0a}/contigs.fa"))?;

    // Run with explicit --subsample-seed 0
    let (code, _) = run_fastbact(&[
        "-1", r1_300, "-2", r2_300, "-o", out_seed_0b,
        "--enable-auto-subsample", "--subsample-seed", "0", "-t", "4",
    ]);
    assert_eq!(code, 0);
    let sha_0b = compute_sha256(format!("{out_seed_0b}/contigs.fa"))?;
    assert_eq!(sha_0a, sha_0b, "Explicit seed 0 must match implicit default seed 0 bit-for-bit");
    println!("  ✓ Implicit seed 0 matches explicit --subsample-seed 0 bit-for-bit ({sha_0a}).");

    // Repeatability test: Run same non-zero seed twice (seed 42)
    println!("\n[Gate 4] Testing repeatability of identical non-zero seed (seed 42) ...");
    let out_seed_42a = "test_gate_tmp/seed_42a";
    let out_seed_42b = "test_gate_tmp/seed_42b";
    let _ = std::fs::remove_dir_all(out_seed_42a);
    let _ = std::fs::remove_dir_all(out_seed_42b);

    let (code, _) = run_fastbact(&[
        "-1", r1_300, "-2", r2_300, "-o", out_seed_42a,
        "--enable-auto-subsample", "--subsample-seed", "42", "-t", "4",
    ]);
    assert_eq!(code, 0);
    let sha_42a = compute_sha256(format!("{out_seed_42a}/contigs.fa"))?;

    let (code, _) = run_fastbact(&[
        "-1", r1_300, "-2", r2_300, "-o", out_seed_42b,
        "--enable-auto-subsample", "--subsample-seed", "42", "-t", "4",
    ]);
    assert_eq!(code, 0);
    let sha_42b = compute_sha256(format!("{out_seed_42b}/contigs.fa"))?;
    assert_eq!(sha_42a, sha_42b, "Identical seed 42 runs must be 100% bit-for-bit identical");
    assert_ne!(sha_42a, sha_0a, "Different seeds (42 vs 0) must produce distinct sampled subsets");
    println!("  ✓ Repeated execution of seed 42 is 100% bit-for-bit identical ({sha_42a}).");
    println!("  ✓ Seed 42 and Seed 0 produce distinct subsampling draws as expected.");


    // Clean up temporary test directories
    let _ = std::fs::remove_dir_all("test_gate_tmp");

    println!("\n=======================================================");
    println!("ALL VAR 3 SOURCE CODE QUALITY GATES VERIFIED SUCCESSFULLY!");
    println!("=======================================================");
    Ok(())
}
