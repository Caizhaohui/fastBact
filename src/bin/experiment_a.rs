//! Experiment A: Multi-Scale, Multi-GC HLL Pure Cardinality & Genome Estimation Generalization.
//!
//! Evaluates HLL cardinality accuracy vs exact FxHashSet truth across:
//! - Genome sizes: 50 kb, 200 kb, 1.0 Mb, 2.5 Mb, 4.65 Mb, 8.0 Mb
//! - GC content: 30%, 50%, 70%
//! - Verifies |HLL_fixed - Exact| / Exact <= 2.5% vs raw min(X,Y) collapse (-71%).

use fastbact::common::kmer::RollingNtHash;
use fastbact::common::sketch::{splitmix64, HyperLogLog};
use rustc_hash::FxHashSet;
use std::fs::File;
use std::io::Write;

fn make_dna_seq(len: usize, gc: f64, seed: u64) -> Vec<u8> {
    let mut seq = Vec::with_capacity(len);
    let mut rng = seed;
    let threshold = (gc * (u64::MAX as f64)) as u64;

    for _ in 0..len {
        rng = splitmix64(rng);
        let is_gc = rng < threshold;
        rng = splitmix64(rng);
        let b = if is_gc {
            if (rng & 1) == 0 { b'C' } else { b'G' }
        } else {
            if (rng & 1) == 0 { b'A' } else { b'T' }
        };
        seq.push(b);
    }
    seq
}

fn main() {
    println!("==========================================================================================");
    println!(" EXPERIMENT A: Multi-Scale (50kb - 8Mb) & Multi-GC (30% - 70%) HLL Cardinality Benchmark ");
    println!("==========================================================================================");

    let sizes = [
        (50_000, "50 kb (Plasmid/Phage)"),
        (200_000, "200 kb (S1 Scale)"),
        (1_000_000, "1.0 Mb (Mycoplasma)"),
        (2_500_000, "2.5 Mb (S. aureus)"),
        (4_650_000, "4.65 Mb (E. coli)"),
        (8_000_000, "8.0 Mb (Streptomyces)"),
    ];

    let gcs = [
        (0.30, "30% GC (AT-rich)"),
        (0.50, "50% GC (Balanced)"),
        (0.70, "70% GC (GC-rich)"),
    ];

    let k = 21;
    let mut tsv_lines = Vec::new();
    tsv_lines.push(
        "genome_scale\tgc_target\ttrue_len_bp\texact_distinct_kmers\traw_min_hll_est\traw_error_pct\tsplitmix64_hll_est\tfixed_error_pct\tgate_verdict".to_string()
    );

    println!(
        "{:<22} {:<18} {:>10} {:>12} {:>12} {:>12} {:>12} {:>12} {:<6}",
        "Scale", "GC", "True BP", "Exact Kmers", "Raw HLL", "Raw Err%", "Fixed HLL", "Fixed Err%", "Gate"
    );
    println!("{}", "-".repeat(115));

    let mut all_pass = true;

    for &(len, scale_name) in &sizes {
        for &(gc, gc_name) in &gcs {
            let seed = (len as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ ((gc * 1000.0) as u64);
            let seq = make_dna_seq(len, gc, seed);

            let mut exact_set = FxHashSet::default();
            let mut hll_raw = HyperLogLog::new();
            let mut hll_fixed = HyperLogLog::new();

            if seq.len() >= k {
                let mut rolling = RollingNtHash::new(&seq[..k]);
                let raw_h = rolling.canonical();
                exact_set.insert(raw_h);
                hll_raw.add_hash(raw_h);
                hll_fixed.add_canonical_hash(raw_h);

                for i in k..seq.len() {
                    let raw_h = rolling.roll(seq[i - k], seq[i]);
                    exact_set.insert(raw_h);
                    hll_raw.add_hash(raw_h);
                    hll_fixed.add_canonical_hash(raw_h);
                }
            }

            let exact = exact_set.len();
            let raw_est = hll_raw.estimate().round() as usize;
            let raw_err = (raw_est as f64 - exact as f64) / exact as f64 * 100.0;

            let fixed_est = hll_fixed.estimate().round() as usize;
            let fixed_err = (fixed_est as f64 - exact as f64) / exact as f64 * 100.0;
            let pass = fixed_err.abs() <= 2.50;

            if !pass {
                all_pass = false;
            }

            let verdict = if pass { "PASS" } else { "FAIL" };

            println!(
                "{:<22} {:<18} {:>10} {:>12} {:>12} {:>11.2}% {:>12} {:>11.2}% {:<6}",
                scale_name, gc_name, len, exact, raw_est, raw_err, fixed_est, fixed_err, verdict
            );

            tsv_lines.push(format!(
                "{}\t{}\t{}\t{}\t{}\t{:.2}%\t{}\t{:.2}%\t{}",
                scale_name, gc_name, len, exact, raw_est, raw_err, fixed_est, fixed_err, verdict
            ));
        }
    }

    println!("{}", "-".repeat(115));
    if all_pass {
        println!(">>> ALL 18 CONDITIONS PASSED GATE: |SplitMix64 HLL Error| <= 2.50% (vs -71% raw collapse) <<<");
    } else {
        eprintln!(">>> GATE FAILED: one or more conditions exceeded 2.50% error! <<<");
        std::process::exit(1);
    }

    let out_path = "benchmarks/B2_reports/C1b_EXPERIMENT_A_HLL_MULTISCALE.tsv";
    let mut file = File::create(out_path).expect("Failed to create Experiment A TSV");
    for line in tsv_lines {
        writeln!(file, "{}", line).unwrap();
    }
    println!("Saved Experiment A report to {}", out_path);
}
