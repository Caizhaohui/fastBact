//! Genome size estimation, coverage calculation, and dynamic multi-k schedule.
//!
//! Provides exact streaming FASTQ pair/base counting, dual-end pairing integrity checks,
//! and initial k-mer cardinality profiling via HyperLogLog and Count-Min Sketch.

use crate::common::kmer::RollingNtHash;
use crate::common::sketch::{splitmix64, CountMinSketch, HyperLogLog};
use anyhow::Result;
use std::path::Path;

/// Confidence state of genome size and depth estimation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EstimatorConfidence {
    /// High confidence: ample solid depth (D_k >= 20x) and clean solid k-mer spectrum.
    High,
    /// Low coverage warning: marginal solid depth (5x <= D_k < 20x); higher Poisson truncation variance.
    LowCoverageWarning,
    /// Unreliable: insufficient coverage (D_k < 5x), excessive error rate, or tiny dataset.
    Unreliable,
}

impl std::fmt::Display for EstimatorConfidence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EstimatorConfidence::High => write!(f, "HIGH"),
            EstimatorConfidence::LowCoverageWarning => write!(f, "LOW_COVERAGE_WARNING"),
            EstimatorConfidence::Unreliable => write!(f, "UNRELIABLE"),
        }
    }
}

/// Statistics from exact streaming FASTQ parsing.
#[derive(Clone, Debug, PartialEq)]
pub struct ExactCountStats {
    /// Exact count of valid read pairs
    pub exact_read_pairs: usize,
    /// Total bases observed in R1
    pub total_bases_r1: u64,
    /// Total bases observed in R2
    pub total_bases_r2: u64,
    /// Exact sum of all sequenced bases (R1 + R2)
    pub total_bases: u64,
    /// Minimum read length observed
    pub min_read_len: usize,
    /// Maximum read length observed
    pub max_read_len: usize,
    /// Mean read length across all reads
    pub mean_read_len: f64,
    /// Median read length
    pub median_read_len: usize,
}

/// Genome estimation summary statistics.
#[derive(Clone, Debug, PartialEq)]
pub struct GenomeEstimate {
    /// Estimated unique genome size in bp (excluding singleton error k-mers)
    pub genome_size: u64,
    /// Exact count of read pairs
    pub exact_read_pairs: usize,
    /// Exact total bases in R1 and R2
    pub total_bases: u64,
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
    /// Confidence classification of the estimation
    pub confidence: EstimatorConfidence,
    /// HLL cardinality of all distinct k-mers (c >= 1)
    pub u_ge_1: u64,
    /// HLL cardinality of solid distinct k-mers (c >= 2)
    pub u_ge_2: u64,
    /// Average k-mer depth across solid k-mers in the sample window
    pub sample_kmer_depth: f64,
}

/// Extracts read identifier prefix for paired-end synchronization checks.
///
/// Strips Casava/Illumina trailing whitespace comments and `/1` or `/2` suffixes.
#[inline]
pub fn extract_read_id_prefix(id: &[u8]) -> &[u8] {
    let non_ws = match id.iter().position(|&b| b.is_ascii_whitespace()) {
        Some(pos) => &id[..pos],
        None => id,
    };
    if non_ws.ends_with(b"/1") || non_ws.ends_with(b"/2") {
        &non_ws[..non_ws.len() - 2]
    } else {
        non_ws
    }
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

/// Performs exact, zero-bias streaming counting of paired-end FASTQ reads.
///
/// Validates pairing synchronization, checks header prefixes, computes base totals,
/// and calculates exact length statistics in O(1) memory.
pub fn count_fastx_exact<P: AsRef<Path>>(
    r1_path: P,
    r2_path: P,
) -> Result<ExactCountStats> {
    let mut r1_reader = crate::ingest::trim_merge::open_fastx(&r1_path)?;
    let mut r2_reader = crate::ingest::trim_merge::open_fastx(&r2_path)?;

    let mut pair_count = 0usize;
    let mut total_bases_r1 = 0u64;
    let mut total_bases_r2 = 0u64;
    let mut min_read_len = usize::MAX;
    let mut max_read_len = 0usize;

    // Fixed-size frequency histogram for read lengths [0..4095] (O(1) 32 KB RAM)
    let mut len_hist = [0usize; 4096];
    let mut overflow_lens = Vec::new();

    loop {
        let rec1_opt = r1_reader.next();
        let rec2_opt = r2_reader.next();

        match (rec1_opt, rec2_opt) {
            (Some(Ok(rec1)), Some(Ok(rec2))) => {
                pair_count += 1;
                let id1_p = extract_read_id_prefix(rec1.id());
                let id2_p = extract_read_id_prefix(rec2.id());

                if id1_p != id2_p {
                    anyhow::bail!(
                        "FASTQ paired-end synchronization error at pair #{}: R1 ID '{:?}' does not match R2 ID '{:?}'",
                        pair_count,
                        String::from_utf8_lossy(rec1.id()),
                        String::from_utf8_lossy(rec2.id())
                    );
                }

                let l1 = rec1.seq().len();
                let l2 = rec2.seq().len();
                total_bases_r1 += l1 as u64;
                total_bases_r2 += l2 as u64;
                min_read_len = min_read_len.min(l1).min(l2);
                max_read_len = max_read_len.max(l1).max(l2);

                if l1 < 4096 {
                    len_hist[l1] += 1;
                } else {
                    overflow_lens.push(l1);
                }
                if l2 < 4096 {
                    len_hist[l2] += 1;
                } else {
                    overflow_lens.push(l2);
                }
            }
            (None, None) => break,
            (Some(Ok(_)), None) => {
                anyhow::bail!(
                    "FASTQ pairing mismatch: R1 has more reads than R2 (R2 ended prematurely after {} pairs)",
                    pair_count
                );
            }
            (None, Some(Ok(_))) => {
                anyhow::bail!(
                    "FASTQ pairing mismatch: R2 has more reads than R1 (R1 ended prematurely after {} pairs)",
                    pair_count
                );
            }
            (Some(Err(e)), _) => anyhow::bail!(
                "Corrupted/truncated FASTQ record in R1 at pair #{}: {}",
                pair_count + 1,
                e
            ),
            (_, Some(Err(e))) => anyhow::bail!(
                "Corrupted/truncated FASTQ record in R2 at pair #{}: {}",
                pair_count + 1,
                e
            ),
        }
    }

    if pair_count == 0 {
        anyhow::bail!("Input FASTQ files contain 0 valid read pairs");
    }

    let total_bases = total_bases_r1 + total_bases_r2;
    let total_reads = pair_count * 2;
    let mean_read_len = total_bases as f64 / total_reads as f64;

    // Compute median length from frequency histogram
    let median_target = (total_reads + 1) / 2;
    let mut cum = 0usize;
    let mut median_read_len = 0usize;

    for (len, &count) in len_hist.iter().enumerate() {
        cum += count;
        if cum >= median_target {
            median_read_len = len;
            break;
        }
    }
    if median_read_len == 0 && !overflow_lens.is_empty() {
        overflow_lens.sort_unstable();
        let idx = (overflow_lens.len() / 2).min(overflow_lens.len() - 1);
        median_read_len = overflow_lens[idx];
    }

    Ok(ExactCountStats {
        exact_read_pairs: pair_count,
        total_bases_r1,
        total_bases_r2,
        total_bases,
        min_read_len,
        max_read_len,
        mean_read_len,
        median_read_len,
    })
}

/// Estimates genome size and coverage by streaming initial read pairs and utilizing exact counting.
pub fn estimate_genome<P: AsRef<Path>>(
    r1_path: P,
    r2_path: P,
    max_pairs: usize,
    target_cov: f64,
) -> Result<GenomeEstimate> {
    // 1. Perform exact, zero-bias streaming counting of all read pairs and bases
    let count_stats = count_fastx_exact(&r1_path, &r2_path)?;

    // 2. Sample initial reads for dual HyperLogLog (total U>=1 and solid U>=2) and Count-Min Sketch
    let mut r1_reader = crate::ingest::trim_merge::open_fastx(&r1_path)?;
    let mut r2_reader = crate::ingest::trim_merge::open_fastx(&r2_path)?;

    let k = 21;
    let mut hll_total = HyperLogLog::new();
    let mut hll_solid = HyperLogLog::new();
    // 1 << 23 width = 8,388,608 slots per row (16 MB RAM total), ensuring minimal collision inflation
    let cms = CountMinSketch::new(1 << 23);

    let sample_target = max_pairs.min(150_000).min(count_stats.exact_read_pairs);
    let mut m_sample: u64 = 0;
    let mut pair_count = 0usize;

    while pair_count < sample_target {
        let rec1_opt = r1_reader.next();
        let rec2_opt = r2_reader.next();

        match (rec1_opt, rec2_opt) {
            (Some(Ok(rec1)), Some(Ok(rec2))) => {
                pair_count += 1;
                let s1 = rec1.seq();
                let s2 = rec2.seq();

                // Process k-mers in R1 using rolling ntHash with splitmix64 mixing
                if s1.len() >= k {
                    m_sample += (s1.len() - k + 1) as u64;
                    let mut rolling = RollingNtHash::new(&s1[..k]);
                    let raw = rolling.canonical();
                    let c = cms.add_and_estimate(raw);
                    let mixed = splitmix64(raw);
                    hll_total.add_hash(mixed);
                    if c >= 2 {
                        hll_solid.add_hash(mixed);
                    }
                    for i in k..s1.len() {
                        let raw = rolling.roll(s1[i - k], s1[i]);
                        let c = cms.add_and_estimate(raw);
                        let mixed = splitmix64(raw);
                        hll_total.add_hash(mixed);
                        if c >= 2 {
                            hll_solid.add_hash(mixed);
                        }
                    }
                }

                // Process k-mers in R2 using rolling ntHash with splitmix64 mixing
                if s2.len() >= k {
                    m_sample += (s2.len() - k + 1) as u64;
                    let mut rolling = RollingNtHash::new(&s2[..k]);
                    let raw = rolling.canonical();
                    let c = cms.add_and_estimate(raw);
                    let mixed = splitmix64(raw);
                    hll_total.add_hash(mixed);
                    if c >= 2 {
                        hll_solid.add_hash(mixed);
                    }
                    for i in k..s2.len() {
                        let raw = rolling.roll(s2[i - k], s2[i]);
                        let c = cms.add_and_estimate(raw);
                        let mixed = splitmix64(raw);
                        hll_total.add_hash(mixed);
                        if c >= 2 {
                            hll_solid.add_hash(mixed);
                        }
                    }
                }

                // Adaptive sampling: if solid k-mer sample depth lambda_s >= 15.0,
                // the sample has reached sufficient Poisson saturation (>99.99% of genome is solid c>=2)
                // Stop sampling early to avoid accumulating recurrent sequencing errors into c>=2.
                if pair_count >= 5_000 && pair_count % 2_000 == 0 {
                    let u1 = hll_total.estimate().round() as u64;
                    let u2 = hll_solid.estimate().round() as u64;
                    let e1 = u1.saturating_sub(u2);
                    let m_tr = m_sample.saturating_sub(e1);
                    let lam = if u2 > 0 { m_tr as f64 / u2 as f64 } else { 0.0 };
                    if lam >= 15.0 {
                        break;
                    }
                }
            }
            _ => break,
        }
    }

    if pair_count == 0 {
        anyhow::bail!("Input FASTQ files contain 0 valid read pairs");
    }

    // 3. Compute distinct cardinality and zero/one-truncated Poisson genome size
    let u_ge_1 = hll_total.estimate().round() as u64;
    let u_ge_2 = hll_solid.estimate().round() as u64;
    let e_1 = u_ge_1.saturating_sub(u_ge_2);
    let m_true = m_sample.saturating_sub(e_1);
    let lambda_s = if u_ge_2 > 0 { m_true as f64 / u_ge_2 as f64 } else { 0.0 };

    // Zero/one-truncated Poisson correction:
    // P(X >= 2) = 1 - (1 + lambda_s) * exp(-lambda_s)
    let p_ge_2 = 1.0 - (1.0 + lambda_s) * (-lambda_s).exp();
    let mut est_g = if p_ge_2 > 0.05 {
        (u_ge_2 as f64 / p_ge_2).round() as u64
    } else {
        u_ge_2.max(10_000)
    };
    if est_g < 10_000 {
        est_g = 10_000;
    }

    // Full library estimated sequencing coverage
    let estimated_cov = count_stats.total_bases as f64 / (est_g as f64);
    let d_k = estimated_cov * ((count_stats.mean_read_len - k as f64 + 1.0) / count_stats.mean_read_len);
    let solid_ratio = if u_ge_1 > 0 { u_ge_2 as f64 / u_ge_1 as f64 } else { 0.0 };

    // Confidence state machine
    let confidence = if d_k < 5.0 || solid_ratio < 0.05 || u_ge_2 < 1_000 {
        EstimatorConfidence::Unreliable
    } else if d_k < 20.0 || solid_ratio < 0.15 {
        EstimatorConfidence::LowCoverageWarning
    } else {
        EstimatorConfidence::High
    };

    let subsample_factor = if target_cov > 0.0 && estimated_cov > 1.5 * target_cov {
        (target_cov / estimated_cov).clamp(0.01, 1.0)
    } else {
        1.0
    };

    let kmers = compute_dynamic_kmers(count_stats.median_read_len, None);

    Ok(GenomeEstimate {
        genome_size: est_g,
        exact_read_pairs: count_stats.exact_read_pairs,
        total_bases: count_stats.total_bases,
        median_read_len: count_stats.median_read_len,
        max_read_len: count_stats.max_read_len,
        estimated_cov,
        subsample_factor,
        kmers,
        confidence,
        u_ge_1,
        u_ge_2,
        sample_kmer_depth: lambda_s,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_extract_read_id_prefix() {
        assert_eq!(extract_read_id_prefix(b"@read1/1"), b"@read1");
        assert_eq!(extract_read_id_prefix(b"@read1/2"), b"@read1");
        assert_eq!(
            extract_read_id_prefix(b"@SRR28349173.1 1 length=150"),
            b"@SRR28349173.1"
        );
        assert_eq!(
            extract_read_id_prefix(b"@SRR28349173.1 2 length=150"),
            b"@SRR28349173.1"
        );
        assert_eq!(extract_read_id_prefix(b"@read_only"), b"@read_only");
    }

    #[test]
    fn test_count_fastx_exact_success() -> Result<()> {
        let mut f1 = NamedTempFile::new()?;
        let mut f2 = NamedTempFile::new()?;

        for i in 1..=10 {
            writeln!(f1, "@pair_{}/1\nACGTACGTACGT\n+\nIIIIIIIIIIII", i)?;
            writeln!(f2, "@pair_{}/2\nTGCATGCATGCA\n+\nIIIIIIIIIIII", i)?;
        }
        f1.flush()?;
        f2.flush()?;

        let stats = count_fastx_exact(f1.path(), f2.path())?;
        assert_eq!(stats.exact_read_pairs, 10);
        assert_eq!(stats.total_bases_r1, 120);
        assert_eq!(stats.total_bases_r2, 120);
        assert_eq!(stats.total_bases, 240);
        assert_eq!(stats.min_read_len, 12);
        assert_eq!(stats.max_read_len, 12);
        assert_eq!(stats.median_read_len, 12);
        assert!((stats.mean_read_len - 12.0).abs() < 1e-6);

        Ok(())
    }

    #[test]
    fn test_count_fastx_exact_mismatched_pair_count() -> Result<()> {
        let mut f1 = NamedTempFile::new()?;
        let mut f2 = NamedTempFile::new()?;

        for i in 1..=5 {
            writeln!(f1, "@pair_{}/1\nACGT\n+\nIIII", i)?;
        }
        for i in 1..=4 {
            writeln!(f2, "@pair_{}/2\nACGT\n+\nIIII", i)?;
        }
        f1.flush()?;
        f2.flush()?;

        let res = count_fastx_exact(f1.path(), f2.path());
        assert!(res.is_err());
        let err_msg = res.unwrap_err().to_string();
        assert!(err_msg.contains("R1 has more reads than R2"));

        Ok(())
    }

    #[test]
    fn test_count_fastx_exact_header_desync() -> Result<()> {
        let mut f1 = NamedTempFile::new()?;
        let mut f2 = NamedTempFile::new()?;

        writeln!(f1, "@pair_A/1\nACGT\n+\nIIII")?;
        writeln!(f2, "@pair_B/2\nACGT\n+\nIIII")?;
        f1.flush()?;
        f2.flush()?;

        let res = count_fastx_exact(f1.path(), f2.path());
        assert!(res.is_err());
        let err_msg = res.unwrap_err().to_string();
        assert!(err_msg.contains("synchronization error"));

        Ok(())
    }

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
