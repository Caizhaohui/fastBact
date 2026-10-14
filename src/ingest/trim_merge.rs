//! Quality trimming, adapter clipping, and FLASH-style paired-end overlap merging.
//!
//! Merges overlapping paired-end reads (insert size < 2 * read_len) into long single-end reads
//! in `ReadBank::merged_reads`, while retaining unmerged pairs in `pe_forward` / `pe_reverse`.

use crate::common::bitseq::{reverse_complement_ascii, ReadBank};
use crate::common::sketch::splitmix64;
use anyhow::{Context, Result};
use std::path::Path;

/// Standard adapter prefixes for 3' clipping.
pub const ADAPTER_TRUSEQ: &[u8] = b"AGATCGGAAGAGC";
pub const ADAPTER_NEXTERA: &[u8] = b"AGATGTGTATAAGAGACAG";
pub const ADAPTER_NEXTERA_RC: &[u8] = b"CTGTCTCTTATACACATCT";

/// Trims 3' adapter sequences if detected with high confidence near the read end.
pub fn trim_adapter<'a>(seq: &'a [u8], qual: &'a [u8]) -> (&'a [u8], &'a [u8]) {
    if seq.len() < 10 {
        return (seq, qual);
    }
    let mut cut_pos = seq.len();
    // Search adapter towards the 3' end (insert size < read length)
    let start = if seq.len() > 100 { seq.len() - 100 } else { 0 };

    for i in start..=(seq.len() - 10) {
        let b0 = seq[i].to_ascii_uppercase();
        let b1 = seq[i + 1].to_ascii_uppercase();

        // Check if prefix could match TruSeq/Nextera (AGAT) or Nextera_RC (CTGT)
        let is_agat = (b0 == b'A' && b1 == b'G') || (b0 == b'A') || (b1 == b'G');
        let is_ctgt = (b0 == b'C' && b1 == b'T') || (b0 == b'C') || (b1 == b'T');

        if !is_agat && !is_ctgt {
            continue;
        }

        if is_agat {
            for adapter in [ADAPTER_TRUSEQ, ADAPTER_NEXTERA] {
                let min_match = 10.min(adapter.len());
                let window = &seq[i..i + min_match];
                let mut mismatches = 0;
                for j in 0..min_match {
                    if !window[j].eq_ignore_ascii_case(&adapter[j]) {
                        mismatches += 1;
                        if mismatches > 1 {
                            break;
                        }
                    }
                }
                if mismatches <= 1 {
                    cut_pos = cut_pos.min(i);
                    break;
                }
            }
        }

        if is_ctgt && cut_pos > i {
            let min_match = 10.min(ADAPTER_NEXTERA_RC.len());
            let window = &seq[i..i + min_match];
            let mut mismatches = 0;
            for j in 0..min_match {
                if !window[j].eq_ignore_ascii_case(&ADAPTER_NEXTERA_RC[j]) {
                    mismatches += 1;
                    if mismatches > 1 {
                        break;
                    }
                }
            }
            if mismatches <= 1 {
                cut_pos = cut_pos.min(i);
            }
        }
    }

    (&seq[..cut_pos], &qual[..cut_pos])
}

/// Trims trailing low-quality bases using a sliding window of size 4 with average Phred < min_qual.
pub fn trim_quality<'a>(
    seq: &'a [u8],
    qual: &'a [u8],
    min_qual: u8,
    min_len: usize,
) -> Option<(&'a [u8], &'a [u8])> {
    if seq.len() < min_len || qual.len() != seq.len() {
        return None;
    }

    let mut cut = seq.len();
    let window_size = 4;

    if qual.len() >= window_size {
        // Fast path: if the tail window already meets min_qual, keep entire read
        let win = &qual[cut - window_size..cut];
        let sum_q: usize = win.iter().map(|&q| q.saturating_sub(33) as usize).sum();
        if sum_q >= (min_qual as usize) * window_size {
            return Some((seq, qual));
        }

        while cut >= window_size {
            let win = &qual[cut - window_size..cut];
            let sum_q: usize = win.iter().map(|&q| q.saturating_sub(33) as usize).sum();
            if sum_q / window_size >= min_qual as usize {
                break;
            }
            cut -= 1;
        }
    }

    if cut < min_len {
        None
    } else {
        Some((&seq[..cut], &qual[..cut]))
    }
}

/// Result of paired-end overlap attempt.
#[derive(Debug, PartialEq, Eq)]
pub enum MergeResult {
    /// Read pair successfully merged into a single consensus read
    Merged {
        seq: Vec<u8>,
        overlap_len: usize,
        mismatches: usize,
    },
    /// No high-confidence overlap found; remains paired
    Unmerged,
}

/// Attempts to merge a pair of reads (R1, R2) with Phred qualities (Q1, Q2)
/// by evaluating overlaps between 3' suffix of R1 and 5' prefix of RC(R2).
///
/// Overlap length `d` ranges from `min_overlap` (default 15) up to `min(L1, L2)`.
/// Maximum allowed mismatch density is 10% (`m / d <= 0.10`).
pub fn merge_read_pair(
    r1: &[u8],
    q1: &[u8],
    r2: &[u8],
    q2: &[u8],
    min_overlap: usize,
    max_mismatch_ratio: f64,
) -> MergeResult {
    let l1 = r1.len();
    let l2 = r2.len();

    if l1 < min_overlap || l2 < min_overlap {
        return MergeResult::Unmerged;
    }

    // RC of R2: stack buffer for typical reads <= 512bp to avoid heap allocation
    let mut stack_rc = [0u8; 512];
    let heap_rc;
    let rc2: &[u8] = if l2 <= 512 {
        for i in 0..l2 {
            stack_rc[i] = match r2[l2 - 1 - i] {
                b'A' | b'a' => b'T',
                b'C' | b'c' => b'G',
                b'G' | b'g' => b'C',
                b'T' | b't' => b'A',
                _ => b'A',
            };
        }
        &stack_rc[..l2]
    } else {
        heap_rc = reverse_complement_ascii(r2);
        &heap_rc
    };

    let max_d = l1.min(l2);
    let mut best_d = 0usize;
    let mut best_mismatches = 0usize;
    let mut best_score = f64::NEG_INFINITY;

    for d in (min_overlap..=max_d).rev() {
        let s1 = &r1[l1 - d..l1];
        let s2 = &rc2[..d];

        let max_allowed = ((d as f64) * max_mismatch_ratio) as usize;
        let mut mismatches = 0usize;
        for i in 0..d {
            if !s1[i].eq_ignore_ascii_case(&s2[i]) {
                mismatches += 1;
                if mismatches > max_allowed {
                    break;
                }
            }
        }

        if mismatches <= max_allowed {
            let score = (d as f64) - 2.5 * (mismatches as f64);
            if score > best_score {
                best_score = score;
                best_d = d;
                best_mismatches = mismatches;
            }
        }

        if (d as f64) <= best_score {
            break;
        }
    }

    if best_d >= min_overlap {
        // Construct merged consensus
        let merged_len = l1 + l2 - best_d;
        let mut merged_seq = Vec::with_capacity(merged_len);

        // 1. Prefix from R1 (before overlap)
        merged_seq.extend_from_slice(&r1[..l1 - best_d]);

        // 2. Overlap region: consensus by Phred score
        for i in 0..best_d {
            let b1 = r1[l1 - best_d + i];
            let b2 = rc2[i];
            let qual1 = q1.get(l1 - best_d + i).copied().unwrap_or(b'I');
            let qual2 = if i < l2 { q2.get(l2 - 1 - i).copied().unwrap_or(b'I') } else { b'I' };

            if b1.eq_ignore_ascii_case(&b2) {
                merged_seq.push(b1.to_ascii_uppercase());
            } else if qual1 >= qual2 {
                merged_seq.push(b1.to_ascii_uppercase());
            } else {
                merged_seq.push(b2.to_ascii_uppercase());
            }
        }

        // 3. Suffix from RC(R2) (after overlap)
        merged_seq.extend_from_slice(&rc2[best_d..]);

        MergeResult::Merged {
            seq: merged_seq,
            overlap_len: best_d,
            mismatches: best_mismatches,
        }
    } else {
        MergeResult::Unmerged
    }
}

/// Opens a FASTX file, transparently handling gzip decompression if extension is .gz.
pub fn open_fastx<P: AsRef<Path>>(path: P) -> Result<Box<dyn needletail::FastxReader>> {
    let path_ref = path.as_ref();
    let file = std::fs::File::open(path_ref)
        .with_context(|| format!("Failed to open file {:?}", path_ref))?;
    let buf_reader = std::io::BufReader::with_capacity(256 * 1024, file);

    let is_gz = path_ref.extension().map_or(false, |ext| ext == "gz");
    if is_gz {
        let gz = flate2::read::GzDecoder::new(buf_reader);
        let reader = needletail::parse_fastx_reader(gz)
            .map_err(|e| anyhow::anyhow!("Failed to parse gzip FASTX: {}", e))?;
        Ok(reader)
    } else {
        let reader = needletail::parse_fastx_reader(buf_reader)
            .map_err(|e| anyhow::anyhow!("Failed to parse FASTX: {}", e))?;
        Ok(reader)
    }
}

/// Streams paired-end FASTQ files, applies adapter & quality trimming,
/// subsampling, and overlap merging, storing all retained reads directly in `ReadBank`.
pub fn ingest_reads<P: AsRef<Path>>(
    r1_path: P,
    r2_path: P,
    subsample_factor: f64,
    min_len: usize,
    min_qual: u8,
) -> Result<ReadBank> {
    let mut r1_reader = open_fastx(&r1_path)?;
    let mut r2_reader = open_fastx(&r2_path)?;

    let mut bank = ReadBank::with_capacity(500_000, 150);
    let mut pair_idx: u64 = 0;

    let should_subsample = subsample_factor < 0.999;
    let keep_threshold = (subsample_factor * (u64::MAX as f64)) as u64;

    loop {
        let rec1_opt = r1_reader.next();
        let rec2_opt = r2_reader.next();

        match (rec1_opt, rec2_opt) {
            (Some(Ok(rec1)), Some(Ok(rec2))) => {
                pair_idx += 1;

                if should_subsample {
                    let rand_val = splitmix64(pair_idx);
                    if rand_val > keep_threshold {
                        continue;
                    }
                }

                let s1 = rec1.seq();
                let q1 = rec1.qual().unwrap_or(&[]);
                let s2 = rec2.seq();
                let q2 = rec2.qual().unwrap_or(&[]);

                // 1. 3' Adapter Trimming
                let (s1_ad, q1_ad) = trim_adapter(&s1, q1);
                let (s2_ad, q2_ad) = trim_adapter(&s2, q2);

                // 2. Sliding window quality trimming
                let trimmed1 = trim_quality(s1_ad, q1_ad, min_qual, min_len);
                let trimmed2 = trim_quality(s2_ad, q2_ad, min_qual, min_len);

                match (trimmed1, trimmed2) {
                    (Some((t1_seq, t1_q)), Some((t2_seq, t2_q))) => {
                        // 3. FLASH overlap merge attempt
                        match merge_read_pair(t1_seq, t1_q, t2_seq, t2_q, 15, 0.10) {
                            MergeResult::Merged { seq, .. } => {
                                bank.merged_reads.push_seq(&seq);
                            }
                            MergeResult::Unmerged => {
                                bank.pe_forward.push_seq(t1_seq);
                                bank.pe_reverse.push_seq(t2_seq);
                            }
                        }
                    }
                    _ => {
                        // If either read drops below min_len after trimming, discard pair
                    }
                }
            }
            (Some(Err(e)), _) => anyhow::bail!("Corrupted record in R1 FASTQ: {}", e),
            (_, Some(Err(e))) => anyhow::bail!("Corrupted record in R2 FASTQ: {}", e),
            (Some(Ok(_)), None) | (None, Some(Ok(_))) => anyhow::bail!("Mismatched paired-end FASTQ mate counts between R1 and R2"),
            (None, None) => break,
        }
    }

    if r1_reader.next().is_some() || r2_reader.next().is_some() {
        anyhow::bail!("Mismatched paired-end FASTQ mate counts between R1 and R2");
    }

    if bank.total_reads() == 0 {
        anyhow::bail!("No valid reads survived trimming and quality filtering (check input FASTQ files)");
    }

    Ok(bank)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_verification_gate_short_fragment_merge_with_errors() {
        // Verification Gate item 1:
        // Synthetic 2 x 150 bp read pairs simulated from a 230 bp fragment
        // with 2 sequencing errors are accurately merged into a single 230 bp consensus read.
        let mut frag = Vec::with_capacity(230);
        let alphabet = [b'A', b'C', b'G', b'T'];
        for i in 0..230 {
            let r = splitmix64(0x5EED_CAFE_1234_5678 + i as u64);
            frag.push(alphabet[(r % 4) as usize]);
        }

        // R1: first 150 bp of fragment [0..150]
        let mut r1 = frag[..150].to_vec();
        let mut q1 = vec![b'I'; 150]; // High Phred 40 (73 - 33 = 40)

        // Introduce 1 error in R1 overlap region (e.g. at pos 130)
        let orig_b1 = r1[130];
        let err_b1 = if orig_b1 == b'A' { b'C' } else { b'A' };
        r1[130] = err_b1;
        q1[130] = b'!'; // Low Phred score 0 (33 - 33 = 0) so consensus favors R2

        // R2: reverse complement of last 150 bp of fragment [80..230]
        let frag_suffix = &frag[80..230];
        let mut r2 = reverse_complement_ascii(frag_suffix);
        let mut q2 = vec![b'I'; 150];

        // Introduce 1 error in R2 overlap region (corresponding to frag pos 100)
        // In frag_suffix (80..230), frag pos 100 is offset 20.
        // In RC(frag_suffix), offset 20 from 5' end corresponds to pos 150 - 1 - 20 = 129 in R2.
        let orig_b2 = r2[129];
        let err_b2 = if orig_b2 == b'A' { b'T' } else { b'A' };
        r2[129] = err_b2;
        q2[129] = b'!'; // Low Phred score 0

        // Overlap should be 150 + 150 - 230 = 70 bp
        let result = merge_read_pair(&r1, &q1, &r2, &q2, 15, 0.10);
        match result {
            MergeResult::Merged { seq, overlap_len, mismatches } => {
                assert_eq!(seq.len(), 230, "Merged length should be exactly 230 bp");
                assert_eq!(overlap_len, 70, "Overlap should be 70 bp");
                assert_eq!(mismatches, 2, "Should identify 2 mismatches");
                assert_eq!(seq, frag, "Merged consensus must match the ground-truth 230 bp fragment!");
            }
            MergeResult::Unmerged => {
                panic!("Expected 230 bp read pair to be merged!");
            }
        }
    }

    #[test]
    fn test_verification_gate_long_fragment_unmerged() {
        // Verification Gate item 2:
        // Read pairs from a 450 bp fragment remain unmerged in pe_forward / pe_reverse.
        let mut frag = Vec::with_capacity(450);
        let alphabet = [b'A', b'C', b'G', b'T'];
        for i in 0..450 {
            let r = splitmix64(0xABCD_EF01_2345_6789 + i as u64);
            frag.push(alphabet[(r % 4) as usize]);
        }

        let r1 = &frag[..150];
        let q1 = vec![b'I'; 150];

        let frag_suffix = &frag[300..450];
        let r2 = reverse_complement_ascii(frag_suffix);
        let q2 = vec![b'I'; 150];

        let result = merge_read_pair(r1, &q1, &r2, &q2, 15, 0.10);
        assert_eq!(
            result,
            MergeResult::Unmerged,
            "450 bp fragment read pair should NOT merge!"
        );
    }

    #[test]
    fn test_adapter_trimming() {
        let mut seq = b"ACGTACGTACGTACGT".to_vec();
        seq.extend_from_slice(b"AGATCGGAAGAGCextra");
        let qual = vec![b'I'; seq.len()];

        let (trimmed_s, trimmed_q) = trim_adapter(&seq, &qual);
        assert_eq!(trimmed_s, b"ACGTACGTACGTACGT");
        assert_eq!(trimmed_q.len(), trimmed_s.len());
    }
}
