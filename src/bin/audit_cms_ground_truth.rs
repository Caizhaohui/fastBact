//! Audit tool for FB-BENCH-B2-C1b-3A:
//! Exact Ground Truth validation of Count-Min Sketch across multi-scale genomes & coverage gradients.

use anyhow::Result;
use fastbact::common::bitseq::ASCII_TO_2BIT;
use fastbact::common::kmer::KmerInt;
use fastbact::common::sketch::{splitmix64, CountMinSketch};
use fastbact::ingest::trim_merge::open_fastx;
use hashbrown::HashMap;
use rayon::prelude::*;
use rustc_hash::FxHasher;
use std::hash::BuildHasherDefault;
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::time::Instant;

type FxHashMap<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;

const NUM_SHARDS: usize = 64;

#[inline(always)]
fn shard_idx(kmer: u64) -> usize {
    splitmix64(kmer) as usize & (NUM_SHARDS - 1)
}

struct ExactKmerCounter {
    shards: Vec<Mutex<FxHashMap<u64, u16>>>,
}

impl ExactKmerCounter {
    fn new() -> Self {
        let mut shards = Vec::with_capacity(NUM_SHARDS);
        for _ in 0..NUM_SHARDS {
            shards.push(Mutex::new(FxHashMap::default()));
        }
        Self { shards }
    }

    fn add(&self, kmer: u64) {
        let s = shard_idx(kmer);
        let mut map = self.shards[s].lock().unwrap();
        let entry = map.entry(kmer).or_insert(0);
        *entry = (*entry).saturating_add(1).min(15);
    }
}

fn read_peak_rss_mb() -> f64 {
    #[cfg(target_os = "linux")]
    {
        if let Ok(status) = std::fs::read_to_string("/proc/self/status") {
            for line in status.lines() {
                if line.starts_with("VmHWM:") {
                    let parts: Vec<&str> = line.split_whitespace().collect();
                    if parts.len() >= 2 {
                        if let Ok(kb) = parts[1].parse::<f64>() {
                            return kb / 1024.0;
                        }
                    }
                }
            }
        }
    }
    0.0
}

struct DatasetSpec {
    id: &'static str,
    r1: &'static str,
    r2: &'static str,
    true_genome_bp: u64,
    true_cov: f64,
}

fn stream_kmers<P: AsRef<Path>, F: Fn(u64) + Send + Sync>(
    r1_path: P,
    r2_path: P,
    k: usize,
    f: F,
) -> Result<usize> {
    let mut r1_reader = open_fastx(r1_path)?;
    let mut r2_reader = open_fastx(r2_path)?;
    let mut pairs = 0;

    let chunk_size = 8192;
    let mut chunk_r1 = Vec::with_capacity(chunk_size);
    let mut chunk_r2 = Vec::with_capacity(chunk_size);

    loop {
        chunk_r1.clear();
        chunk_r2.clear();

        while chunk_r1.len() < chunk_size {
            match (r1_reader.next(), r2_reader.next()) {
                (Some(Ok(rec1)), Some(Ok(rec2))) => {
                    pairs += 1;
                    chunk_r1.push(rec1.seq().to_vec());
                    chunk_r2.push(rec2.seq().to_vec());
                }
                _ => break,
            }
        }

        if chunk_r1.is_empty() {
            break;
        }

        // Parallel processing of chunks
        chunk_r1.par_iter().zip(chunk_r2.par_iter()).for_each(|(s1, s2)| {
            if s1.len() >= k {
                let mut fwd = <u64 as KmerInt>::from_ascii(&s1[0..k]);
                let mut rev = fwd.rev_comp(k);
                let can = if rev < fwd { rev } else { fwd };
                f(can);
                for i in k..s1.len() {
                    let base = ASCII_TO_2BIT[s1[i] as usize];
                    let (next_fwd, next_rev) = u64::roll_next(fwd, rev, base, k);
                    fwd = next_fwd;
                    rev = next_rev;
                    let can = if rev < fwd { rev } else { fwd };
                    f(can);
                }
            }
            if s2.len() >= k {
                let mut fwd = <u64 as KmerInt>::from_ascii(&s2[0..k]);
                let mut rev = fwd.rev_comp(k);
                let can = if rev < fwd { rev } else { fwd };
                f(can);
                for i in k..s2.len() {
                    let base = ASCII_TO_2BIT[s2[i] as usize];
                    let (next_fwd, next_rev) = u64::roll_next(fwd, rev, base, k);
                    fwd = next_fwd;
                    rev = next_rev;
                    let can = if rev < fwd { rev } else { fwd };
                    f(can);
                }
            }
        });
    }

    Ok(pairs)
}

struct CmsMetrics {
    config_name: String,
    cms_width: usize,
    cms_memory_mb: f64,
    total_distinct: usize,
    exact_solid: usize,
    exact_singleton: usize,
    cms_solid: usize,
    false_positives: usize,
    false_negatives: usize,
    solid_precision_pct: f64,
    solid_recall_pct: f64,
    total_overcount: u64,
    mean_overcount: f64,
    max_overcount: u8,
    saturation_kmer_count: usize,
    saturation_kmer_pct: f64,
    table_sat_pct: f64,
    wall_time_s: f64,
    peak_rss_mb: f64,
}

fn evaluate_cms_on_exact(
    exact: &ExactKmerCounter,
    dataset: &DatasetSpec,
    k: usize,
    width_mult: f64,
    config_label: &str,
) -> Result<CmsMetrics> {
    let t0 = Instant::now();
    let base_width = CountMinSketch::for_genome_size(dataset.true_genome_bp).width;
    let target_width = ((base_width as f64 * width_mult).round() as usize)
        .max(1 << 20)
        .min(1 << 28);
    let cms = CountMinSketch::new(target_width);
    let cms_mem_mb = (4 * cms.width * 1) as f64 / (2.0 * 1024.0 * 1024.0);

    // Populate CMS with same reads
    stream_kmers(dataset.r1, dataset.r2, k, |kmer| {
        cms.add(kmer);
    })?;

    let fill_time = t0.elapsed().as_secs_f64();
    let peak_rss = read_peak_rss_mb();

    // Now compare CMS against exact counts across all shards in parallel
    let mut total_distinct = 0;
    let mut exact_solid = 0;
    let mut exact_singleton = 0;
    let mut cms_solid = 0;
    let mut false_positives = 0;
    let mut false_negatives = 0;
    let mut total_overcount = 0u64;
    let mut max_overcount = 0u8;
    let mut saturation_kmer_count = 0;

    for shard in &exact.shards {
        let map = shard.lock().unwrap();
        for (&kmer, &true_cnt) in map.iter() {
            total_distinct += 1;
            if true_cnt >= 2 {
                exact_solid += 1;
            } else {
                exact_singleton += 1;
            }

            let cms_cnt = cms.estimate(kmer);
            if cms_cnt >= 2 {
                cms_solid += 1;
            }

            if true_cnt == 1 && cms_cnt >= 2 {
                false_positives += 1;
            }
            if true_cnt >= 2 && cms_cnt < 2 {
                false_negatives += 1;
            }

            let over = cms_cnt.saturating_sub(true_cnt as u8);
            total_overcount += over as u64;
            if over > max_overcount {
                max_overcount = over;
            }

            if cms_cnt >= 15 {
                saturation_kmer_count += 1;
            }
        }
    }

    let precision = if cms_solid > 0 {
        ((exact_solid - false_negatives) as f64 / cms_solid as f64) * 100.0
    } else {
        100.0
    };
    let recall = if exact_solid > 0 {
        ((exact_solid - false_negatives) as f64 / exact_solid as f64) * 100.0
    } else {
        100.0
    };
    let mean_overcount = if total_distinct > 0 {
        total_overcount as f64 / total_distinct as f64
    } else {
        0.0
    };
    let sat_kmer_pct = if total_distinct > 0 {
        (saturation_kmer_count as f64 / total_distinct as f64) * 100.0
    } else {
        0.0
    };

    // Calculate table saturation rate: fraction of 4-bit counters == 15
    let mut sat_counters = 0usize;
    let total_table_counters = cms.width * 4;
    for byte in &cms.table {
        let val = byte.load(std::sync::atomic::Ordering::Relaxed);
        let low = val & 0x0F;
        let high = (val >> 4) & 0x0F;
        if low >= 15 {
            sat_counters += 1;
        }
        if high >= 15 {
            sat_counters += 1;
        }
    }
    let table_sat_pct = (sat_counters as f64 / total_table_counters as f64) * 100.0;

    Ok(CmsMetrics {
        config_name: config_label.to_string(),
        cms_width: cms.width,
        cms_memory_mb: cms_mem_mb,
        total_distinct,
        exact_solid,
        exact_singleton,
        cms_solid,
        false_positives,
        false_negatives,
        solid_precision_pct: precision,
        solid_recall_pct: recall,
        total_overcount,
        mean_overcount,
        max_overcount,
        saturation_kmer_count,
        saturation_kmer_pct: sat_kmer_pct,
        table_sat_pct,
        wall_time_s: fill_time,
        peak_rss_mb: peak_rss,
    })
}

fn main() -> Result<()> {
    let datasets = [
        DatasetSpec {
            id: "S1-COV-001",
            r1: "benchmark-results/data_store/S1-COV-001_R1.fq.gz",
            r2: "benchmark-results/data_store/S1-COV-001_R2.fq.gz",
            true_genome_bp: 200_000,
            true_cov: 15.0,
        },
        DatasetSpec {
            id: "S1-COV-030X",
            r1: "benchmark-results/data_store/S1-COV-030X_R1.fq.gz",
            r2: "benchmark-results/data_store/S1-COV-030X_R2.fq.gz",
            true_genome_bp: 200_000,
            true_cov: 30.0,
        },
        DatasetSpec {
            id: "S1-COV-002",
            r1: "benchmark-results/data_store/S1-COV-002_R1.fq.gz",
            r2: "benchmark-results/data_store/S1-COV-002_R2.fq.gz",
            true_genome_bp: 200_000,
            true_cov: 60.0,
        },
        DatasetSpec {
            id: "S1-COV-100X",
            r1: "benchmark-results/data_store/S1-COV-100X_R1.fq.gz",
            r2: "benchmark-results/data_store/S1-COV-100X_R2.fq.gz",
            true_genome_bp: 200_000,
            true_cov: 100.0,
        },
        DatasetSpec {
            id: "S1-COV-150X",
            r1: "benchmark-results/data_store/S1-COV-150X_R1.fq.gz",
            r2: "benchmark-results/data_store/S1-COV-150X_R2.fq.gz",
            true_genome_bp: 200_000,
            true_cov: 150.0,
        },
        DatasetSpec {
            id: "S1-COV-003",
            r1: "benchmark-results/data_store/S1-COV-003_R1.fq.gz",
            r2: "benchmark-results/data_store/S1-COV-003_R2.fq.gz",
            true_genome_bp: 200_000,
            true_cov: 300.0,
        },
        DatasetSpec {
            id: "SYN-2.5M",
            r1: "benchmark-results/data_store/SYN-2.5M_R1.fq.gz",
            r2: "benchmark-results/data_store/SYN-2.5M_R2.fq.gz",
            true_genome_bp: 2_500_000,
            true_cov: 100.0,
        },
        DatasetSpec {
            id: "S2-REAL-001",
            r1: "benchmark-results/data_store/S2-REAL-001_R1.fq.gz",
            r2: "benchmark-results/data_store/S2-REAL-001_R2.fq.gz",
            true_genome_bp: 4_655_000,
            true_cov: 100.0,
        },
    ];

    let k = 31;
    let out_tsv = PathBuf::from("benchmarks/B2_reports/C1b3_CMS_EXACT_VALIDATION.tsv");
    std::fs::create_dir_all(out_tsv.parent().unwrap())?;

    let header = [
        "dataset_id",
        "true_genome_bp",
        "true_cov",
        "cms_config",
        "cms_width",
        "cms_mem_mb",
        "total_distinct_kmers",
        "exact_solid_kmers",
        "exact_singleton_kmers",
        "cms_solid_kmers",
        "false_positives",
        "false_negatives",
        "solid_precision_pct",
        "solid_recall_pct",
        "total_overcount",
        "mean_overcount",
        "max_overcount",
        "sat_kmer_count",
        "sat_kmer_pct",
        "table_sat_pct",
        "wall_time_s",
        "peak_rss_mb",
    ];

    let mut rows: Vec<String> = vec![header.join("\t")];

    for d in &datasets {
        eprintln!("\n=======================================================");
        eprintln!("Computing Exact Ground Truth for {} ({}bp, {}x)...", d.id, d.true_genome_bp, d.true_cov);
        eprintln!("=======================================================");

        let t_exact = Instant::now();
        let exact = ExactKmerCounter::new();
        let pairs = stream_kmers(d.r1, d.r2, k, |kmer| {
            exact.add(kmer);
        })?;
        eprintln!("  Ingested {} pairs in {:.2?} (Peak RSS: {:.1} MB)", pairs, t_exact.elapsed(), read_peak_rss_mb());

        // Test 3 CMS configurations: Half (0.5x), Dynamic_Standard (1.0x), Double (2.0x)
        let configs = [
            ("Half_Capacity", 0.5),
            ("Dynamic_Standard", 1.0),
            ("Double_Capacity", 2.0),
        ];

        for (cfg_name, mult) in configs {
            eprintln!("  Testing CMS config {} (mult={})...", cfg_name, mult);
            let m = evaluate_cms_on_exact(&exact, d, k, mult, cfg_name)?;
            eprintln!(
                "    Result: Precision={:.4}%, Recall={:.4}%, FP={}, FN={}, MeanOvercount={:.4}, TableSat={:.2}%, Mem={:.1}MB",
                m.solid_precision_pct, m.solid_recall_pct, m.false_positives, m.false_negatives, m.mean_overcount, m.table_sat_pct, m.cms_memory_mb
            );

            let row = [
                d.id.to_string(),
                d.true_genome_bp.to_string(),
                format!("{:.1}", d.true_cov),
                m.config_name,
                m.cms_width.to_string(),
                format!("{:.1}", m.cms_memory_mb),
                m.total_distinct.to_string(),
                m.exact_solid.to_string(),
                m.exact_singleton.to_string(),
                m.cms_solid.to_string(),
                m.false_positives.to_string(),
                m.false_negatives.to_string(),
                format!("{:.4}", m.solid_precision_pct),
                format!("{:.4}", m.solid_recall_pct),
                m.total_overcount.to_string(),
                format!("{:.4}", m.mean_overcount),
                m.max_overcount.to_string(),
                m.saturation_kmer_count.to_string(),
                format!("{:.2}", m.saturation_kmer_pct),
                format!("{:.2}", m.table_sat_pct),
                format!("{:.2}", m.wall_time_s),
                format!("{:.1}", m.peak_rss_mb),
            ];
            rows.push(row.join("\t"));
        }
    }

    std::fs::write(&out_tsv, rows.join("\n") + "\n")?;
    eprintln!("\nSuccessfully written exact CMS validation report to {:?}", out_tsv);

    Ok(())
}
