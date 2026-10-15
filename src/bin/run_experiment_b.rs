//! Experiment B: Comprehensive Evaluation across 13 Benchmarking Datasets.
//!
//! Evaluates C1b-2 genome size estimation, sequencing depth, confidence classification,
//! and compares with v0.1.0 uncalibrated baseline.
//! Outputs `benchmarks/B2_reports/C1b_ESTIMATOR_ACCURACY.tsv`.

use fastbact::common::kmer::RollingNtHash;
use fastbact::common::sketch::{CountMinSketch, HyperLogLog};
use fastbact::ingest::estimator::estimate_genome;
use std::fs::File;
use std::io::Write;
use std::path::Path;

struct DatasetInfo {
    id: &'static str,
    r1: &'static str,
    r2: &'static str,
    true_genome_bp: u64,
    true_cov: f64,
}

fn simulate_v010_estimator<P: AsRef<Path>>(r1_path: P, r2_path: P) -> (u64, f64) {
    let mut r1_reader = fastbact::ingest::trim_merge::open_fastx(&r1_path).unwrap();
    let mut r2_reader = fastbact::ingest::trim_merge::open_fastx(&r2_path).unwrap();

    let k = 21;
    let mut hll_solid = HyperLogLog::new();
    let cms = CountMinSketch::new(1 << 20); // 1M counters per row

    let mut pair_count = 0usize;
    let sample_target = 50_000usize;
    let mut sampled_bases = 0u64;

    while pair_count < sample_target {
        let rec1 = r1_reader.next();
        let rec2 = r2_reader.next();
        match (rec1, rec2) {
            (Some(Ok(r1)), Some(Ok(r2))) => {
                pair_count += 1;
                let s1 = r1.seq();
                let s2 = r2.seq();
                sampled_bases += s1.len() as u64 + s2.len() as u64;

                if s1.len() >= k {
                    let mut rolling = RollingNtHash::new(&s1[..k]);
                    let h = rolling.canonical(); // min(h_fwd, h_rev) WITHOUT splitmix64
                    let c = cms.add_and_estimate(h);
                    if c == 2 {
                        hll_solid.add_hash(h);
                    }
                    for i in k..s1.len() {
                        let h = rolling.roll(s1[i - k], s1[i]);
                        let c = cms.add_and_estimate(h);
                        if c == 2 {
                            hll_solid.add_hash(h);
                        }
                    }
                }

                if s2.len() >= k {
                    let mut rolling = RollingNtHash::new(&s2[..k]);
                    let h = rolling.canonical();
                    let c = cms.add_and_estimate(h);
                    if c == 2 {
                        hll_solid.add_hash(h);
                    }
                    for i in k..s2.len() {
                        let h = rolling.roll(s2[i - k], s2[i]);
                        let c = cms.add_and_estimate(h);
                        if c == 2 {
                            hll_solid.add_hash(h);
                        }
                    }
                }
            }
            _ => break,
        }
    }

    let mut est_g = hll_solid.estimate() as u64;
    if est_g < 10_000 {
        est_g = est_g.max((sampled_bases / 30).max(50_000));
    }

    // v0.1.0 extrapolation via R1 file size / 75
    let file_bytes = std::fs::metadata(r1_path.as_ref()).map(|m| m.len()).unwrap_or(sampled_bases);
    let approx_bytes_per_pair = 75u64;
    let total_pairs_est = (file_bytes / approx_bytes_per_pair).max(pair_count as u64);
    let sample_ratio = (total_pairs_est as f64) / (pair_count as f64);
    let total_projected_bases = (sampled_bases as f64 * sample_ratio).max(1.0);
    let est_cov = total_projected_bases / (est_g as f64);

    (est_g, est_cov)
}

fn main() {
    println!("==========================================================================================");
    println!(" EXPERIMENT B: Benchmark Evaluation across 13 Standard Datasets (C1b-2 Verification)      ");
    println!("==========================================================================================");

    let datasets = [
        DatasetInfo {
            id: "S1-COV-001",
            r1: "benchmark-results/data_store/S1-COV-001_R1.fq.gz",
            r2: "benchmark-results/data_store/S1-COV-001_R2.fq.gz",
            true_genome_bp: 200_000,
            true_cov: 15.0,
        },
        DatasetInfo {
            id: "S1-COV-030X",
            r1: "benchmark-results/data_store/S1-COV-030X_R1.fq.gz",
            r2: "benchmark-results/data_store/S1-COV-030X_R2.fq.gz",
            true_genome_bp: 200_000,
            true_cov: 30.0,
        },
        DatasetInfo {
            id: "S1-COV-002",
            r1: "benchmark-results/data_store/S1-COV-002_R1.fq.gz",
            r2: "benchmark-results/data_store/S1-COV-002_R2.fq.gz",
            true_genome_bp: 200_000,
            true_cov: 60.0,
        },
        DatasetInfo {
            id: "S1-COV-100X",
            r1: "benchmark-results/data_store/S1-COV-100X_R1.fq.gz",
            r2: "benchmark-results/data_store/S1-COV-100X_R2.fq.gz",
            true_genome_bp: 200_000,
            true_cov: 100.0,
        },
        DatasetInfo {
            id: "S1-COV-150X",
            r1: "benchmark-results/data_store/S1-COV-150X_R1.fq.gz",
            r2: "benchmark-results/data_store/S1-COV-150X_R2.fq.gz",
            true_genome_bp: 200_000,
            true_cov: 150.0,
        },
        DatasetInfo {
            id: "S1-COV-003",
            r1: "benchmark-results/data_store/S1-COV-003_R1.fq.gz",
            r2: "benchmark-results/data_store/S1-COV-003_R2.fq.gz",
            true_genome_bp: 200_000,
            true_cov: 300.0,
        },
        DatasetInfo {
            id: "S1-GC-001",
            r1: "benchmark-results/data_store/S1-GC-001_R1.fq.gz",
            r2: "benchmark-results/data_store/S1-GC-001_R2.fq.gz",
            true_genome_bp: 200_000,
            true_cov: 60.0,
        },
        DatasetInfo {
            id: "S1-GC-002",
            r1: "benchmark-results/data_store/S1-GC-002_R1.fq.gz",
            r2: "benchmark-results/data_store/S1-GC-002_R2.fq.gz",
            true_genome_bp: 200_000,
            true_cov: 60.0,
        },
        DatasetInfo {
            id: "SYN-2.5M",
            r1: "benchmark-results/data_store/SYN-2.5M_R1.fq.gz",
            r2: "benchmark-results/data_store/SYN-2.5M_R2.fq.gz",
            true_genome_bp: 2_500_000,
            true_cov: 100.0,
        },
        DatasetInfo {
            id: "S2-REAL-001",
            r1: "benchmark-results/data_store/S2-REAL-001_R1.fq.gz",
            r2: "benchmark-results/data_store/S2-REAL-001_R2.fq.gz",
            true_genome_bp: 4_655_000,
            true_cov: 100.0,
        },
        DatasetInfo {
            id: "R1-001",
            r1: "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349173/SRR28349173_1.fastq.gz",
            r2: "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349173/SRR28349173_2.fastq.gz",
            true_genome_bp: 4_641_652,
            true_cov: 77.8,
        },
        DatasetInfo {
            id: "R1-002",
            r1: "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349221/SRR28349221_1.fastq.gz",
            r2: "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349221/SRR28349221_2.fastq.gz",
            true_genome_bp: 4_631_469,
            true_cov: 64.5,
        },
        DatasetInfo {
            id: "R1-003",
            r1: "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349164/SRR28349164_1.fastq.gz",
            r2: "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349164/SRR28349164_2.fastq.gz",
            true_genome_bp: 4_631_469,
            true_cov: 132.3,
        },
    ];

    let header = "dataset_id\ttrue_genome_bp\tv0.1.0_est_genome_bp\tv0.1.0_error_pct\tc1b2_est_genome_bp\tc1b2_error_pct\ttrue_cov\tv0.1.0_est_cov\tc1b2_est_cov\tc1b2_cov_error_pct\tu_ge_1_hll\tu_ge_2_hll\tsample_kmer_depth\testimator_confidence";
    let mut tsv_lines = vec![header.to_string()];

    println!(
        "{:<12} {:>10} {:>10} {:>9} {:>10} {:>9} {:>6} {:>8} {:>8} {:>9} {:>8} {:<18}",
        "Dataset", "True G", "v0.1.0 G", "v0.1 Err", "C1b2 G", "C1b2 Err", "TrueCov", "v0.1 Cov", "C1b2 Cov", "Cov Err", "SampleD", "Confidence"
    );
    println!("{}", "-".repeat(128));

    let mut all_pass = true;

    for d in &datasets {
        // 1. Simulate v0.1.0 estimator
        let (v010_g, v010_cov) = simulate_v010_estimator(d.r1, d.r2);
        let v010_err_pct = (v010_g as f64 - d.true_genome_bp as f64) / (d.true_genome_bp as f64) * 100.0;

        // 2. Run C1b-2 estimator
        let c1b2 = estimate_genome(d.r1, d.r2, 250_000, 0.0).unwrap();
        let c1b2_err_pct = (c1b2.genome_size as f64 - d.true_genome_bp as f64) / (d.true_genome_bp as f64) * 100.0;
        let c1b2_cov_err_pct = (c1b2.estimated_cov - d.true_cov) / d.true_cov * 100.0;

        // Verify criteria
        let error_pass = if d.true_cov >= 30.0 {
            c1b2_err_pct.abs() <= 10.0
        } else {
            // For 15x, allow up to 10% as well, and must be LOW_COVERAGE_WARNING
            c1b2_err_pct.abs() <= 10.0
        };

        let conf_pass = if d.id == "S1-COV-001" {
            c1b2.confidence == fastbact::ingest::EstimatorConfidence::LowCoverageWarning
        } else {
            c1b2.confidence == fastbact::ingest::EstimatorConfidence::High
        };

        if !error_pass || !conf_pass {
            all_pass = false;
        }

        println!(
            "{:<12} {:>10} {:>10} {:>8.1}% {:>10} {:>8.2}% {:>6.1}x {:>7.1}x {:>7.1}x {:>8.1}% {:>7.2}x {:<18}",
            d.id,
            d.true_genome_bp,
            v010_g,
            v010_err_pct,
            c1b2.genome_size,
            c1b2_err_pct,
            d.true_cov,
            v010_cov,
            c1b2.estimated_cov,
            c1b2_cov_err_pct,
            c1b2.sample_kmer_depth,
            c1b2.confidence.to_string()
        );

        tsv_lines.push(format!(
            "{}\t{}\t{}\t{:.2}%\t{}\t{:.2}%\t{:.1}x\t{:.1}x\t{:.1}x\t{:.2}%\t{}\t{}\t{:.2}x\t{}",
            d.id,
            d.true_genome_bp,
            v010_g,
            v010_err_pct,
            c1b2.genome_size,
            c1b2_err_pct,
            d.true_cov,
            v010_cov,
            c1b2.estimated_cov,
            c1b2_cov_err_pct,
            c1b2.u_ge_1,
            c1b2.u_ge_2,
            c1b2.sample_kmer_depth,
            c1b2.confidence
        ));
    }

    println!("{}", "-".repeat(128));
    if all_pass {
        println!(">>> ALL 13 DATASETS PASSED ACCEPTANCE GATES (G-error <= 10%, Confidence states strictly verified) <<<");
    } else {
        eprintln!(">>> ACCEPTANCE GATES FAILED! <<<");
        std::process::exit(1);
    }

    let out_path = "benchmarks/B2_reports/C1b_ESTIMATOR_ACCURACY.tsv";
    let mut file = File::create(out_path).expect("Failed to create C1b_ESTIMATOR_ACCURACY.tsv");
    for line in tsv_lines {
        writeln!(file, "{}", line).unwrap();
    }
    println!("Saved Experiment B accuracy report to {}", out_path);
}
