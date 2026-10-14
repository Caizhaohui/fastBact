//! Genome size estimation, coverage calculation, and dynamic multi-k schedule.
//!
//! Replaces KMC/Mash by streaming the initial read pairs, populating HyperLogLog
//! and Count-Min Sketch to estimate genome size and sequencing depth.

use crate::common::kmer::RollingNtHash;
use crate::common::sketch::{CountMinSketch, HyperLogLog};
use anyhow::Result;
use std::path::Path;

/// Genome estimation summary statistics.
#[derive(Clone, Debug, PartialEq)]
pub struct GenomeEstimate {
    /// Estimated unique genome size in bp (excluding singleton error k-mers)
    pub genome_size: u64,
    /// Median read length observed in the sample
    pub median_read_len: usize,
    /// Maximum read length observed
    pub max_read_len: usize,
    /// Estimated sequencing depth
    pub estimated_cov: f64,
    /// Subsampling keep probability in (0.0, 1.0]
    pub subsample_factor: f64,
    /// Dynamic odd multi-k schedule
    pub kmers: Vec<usize>,
}

/// Dynamically computes the odd multi-k schedule based on median read length and pairing.
///
/// Ensures all returned k-values are odd integers within [21, 127] and bounded by read length.
pub fn compute_dynamic_kmers(median_len: usize, target_max_k: Option<usize>) -> Vec<usize> {
    let min_k = if median_len < 75 { 21 } else { 31 };

    // For paired-end reads, merged reads reach ~1.6x read length
    let effective_len = (median_len * 16) / 10;
    let auto_max = ((effective_len as f64) * 0.85) as usize;
    let mut max_k = target_max_k.unwrap_or(auto_max).clamp(min_k, 127);

    // Ensure max_k is odd
    if max_k % 2 == 0 {
        max_k -= 1;
    }

    if max_k <= min_k {
        return vec![min_k];
    }

    if median_len >= 100 {
        // High-performance 3-step ladder aligned with hardware registers:
        // k=31 (u64 native), k=63 (u128 native), k=max_k (up to 127 via [u64; 4])
        let mid_k = 63;
        if max_k > mid_k + 10 {
            vec![min_k, mid_k, max_k]
        } else if max_k > min_k {
            vec![min_k, max_k]
        } else {
            vec![min_k]
        }
    } else {
        // Generate up to 3 odd steps for short reads
        let steps = 2;
        let mut step_size = (max_k - min_k) / steps;
        if step_size % 2 == 1 {
            step_size += 1;
        }
        if step_size < 6 {
            step_size = 6;
        }

        let mut kmers = Vec::new();
        let mut k = min_k;
        while k <= max_k {
            kmers.push(k);
            k += step_size;
        }

        if *kmers.last().unwrap() != max_k && max_k > *kmers.last().unwrap() {
            if max_k - kmers.last().unwrap() >= 10 {
                kmers.push(max_k);
            } else {
                *kmers.last_mut().unwrap() = max_k;
            }
        }

        kmers.dedup();
        kmers
    }
}

/// Estimates genome size and coverage by streaming up to `max_pairs` read pairs from R1 and R2.
pub fn estimate_genome<P: AsRef<Path>>(
    r1_path: P,
    r2_path: P,
    max_pairs: usize,
    target_cov: f64,
) -> Result<GenomeEstimate> {
    let mut r1_reader = crate::ingest::trim_merge::open_fastx(&r1_path)?;
    let mut r2_reader = crate::ingest::trim_merge::open_fastx(&r2_path)?;

    let k = 21;
    let mut hll_solid = HyperLogLog::new();
    let cms = CountMinSketch::new(1 << 20); // 1M counters per row (2MB RAM for sample)

    let mut read_lens = Vec::with_capacity(max_pairs.min(100_000));
    let mut sampled_bases: u64 = 0;
    let mut pair_count = 0usize;
    let sample_target = max_pairs.min(50_000);

    while pair_count < sample_target {
        let rec1_opt = r1_reader.next();
        let rec2_opt = r2_reader.next();

        match (rec1_opt, rec2_opt) {
            (Some(Ok(rec1)), Some(Ok(rec2))) => {
                pair_count += 1;
                let s1 = rec1.seq();
                let s2 = rec2.seq();

                sampled_bases += s1.len() as u64 + s2.len() as u64;
                if read_lens.len() < 100_000 {
                    read_lens.push(s1.len());
                    read_lens.push(s2.len());
                }

                // Process k-mers in R1 using rolling ntHash
                if s1.len() >= k {
                    let mut rolling = RollingNtHash::new(&s1[..k]);
                    let h = rolling.canonical();
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

                // Process k-mers in R2 using rolling ntHash
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

    if pair_count == 0 {
        anyhow::bail!("Input FASTQ files contain 0 valid read pairs");
    }

    read_lens.sort_unstable();
    let median_read_len = read_lens[read_lens.len() / 2];
    let max_read_len = *read_lens.last().unwrap_or(&150);

    // Estimated unique solid 21-mers
    let mut est_g = hll_solid.estimate() as u64;
    if est_g < 10_000 {
        est_g = est_g.max((sampled_bases / 30).max(50_000));
    }

    // Project total bases from file size if we only sampled a subset
    let sample_ratio = if pair_count >= sample_target {
        if let Ok(meta) = std::fs::metadata(r1_path.as_ref()) {
            let file_bytes = meta.len();
            let is_gz = r1_path.as_ref().to_string_lossy().ends_with(".gz");
            let approx_bytes_per_pair = if is_gz { 75u64 } else { (median_read_len * 4 + 40) as u64 };
            let total_pairs_est = (file_bytes / approx_bytes_per_pair).max(pair_count as u64);
            (total_pairs_est as f64) / (pair_count as f64)
        } else {
            1.0
        }
    } else {
        1.0
    };

    let total_projected_bases = (sampled_bases as f64 * sample_ratio).max(1.0);
    let estimated_cov = total_projected_bases / (est_g as f64);

    let subsample_factor = if target_cov > 0.0 && estimated_cov > 1.1 * target_cov {
        (target_cov / estimated_cov).clamp(0.01, 1.0)
    } else {
        1.0
    };


    let kmers = compute_dynamic_kmers(median_read_len, None);

    Ok(GenomeEstimate {
        genome_size: est_g,
        median_read_len,
        max_read_len,
        estimated_cov,
        subsample_factor,
        kmers,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_dynamic_kmers_standard_150bp() {
        let kmers = compute_dynamic_kmers(150, None);
        assert_eq!(*kmers.first().unwrap(), 31);
        assert_eq!(*kmers.last().unwrap(), 127);
        for &k in &kmers {
            assert_eq!(k % 2, 1, "k-mer {} must be odd", k);
            assert!(k >= 31 && k <= 127);
        }
    }

    #[test]
    fn test_compute_dynamic_kmers_short_75bp() {
        let kmers = compute_dynamic_kmers(50, None);
        assert_eq!(*kmers.first().unwrap(), 21);
        for &k in &kmers {
            assert_eq!(k % 2, 1);
        }
    }
}
