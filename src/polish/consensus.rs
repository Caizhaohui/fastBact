//! In-memory pileup consensus polishing.
//!
//! Replaces BWA + Pilon by anchoring reads with exact 31-mer seeds and accumulating
//! per-position nucleotide base counts, substituting reference alleles with >75% support.

use crate::common::bitseq::{ASCII_TO_2BIT, TWO_BIT_TO_ASCII};
use crate::common::kmer::nthash_canonical;
use crate::graph::cdbg::{CdbgGraph, UnitigNode};
use hashbrown::HashMap;
use rustc_hash::FxHasher;
use std::hash::BuildHasherDefault;

type FxHashMap<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;

fn output_eligible(node: &UnitigNode, min_output_len: usize) -> bool {
    !node.is_deleted && node.seq_len as usize >= min_output_len
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum AnchorOrientation {
    Forward,
    Reverse,
}

#[inline]
fn is_reverse_complement_exact(anchor: &[u8], read_kmer: &[u8]) -> bool {
    if anchor.len() != read_kmer.len() {
        return false;
    }
    let k = anchor.len();
    for j in 0..k {
        let b_anc = anchor[j];
        let b_read = read_kmer[k - 1 - j];
        let comp = match b_read {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => return false,
        };
        if !b_anc.eq_ignore_ascii_case(&comp) {
            return false;
        }
    }
    true
}

/// Statistics collected during consensus polishing.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PolishStats {
    pub bases_modified: usize,
    pub fwd_mapped_reads: usize,
    pub rev_mapped_reads: usize,
    pub fwd_pileup_bases: usize,
    pub rev_pileup_bases: usize,
    pub rejected_anchors: usize,
}

/// Anchors reads from the read pool onto assembled contigs and performs fast in-memory
/// pileup error correction.
///
/// `min_output_len` is the FASTA retention threshold. Shorter residual nodes are kept in
/// the graph for connectivity accounting, but they are not consensus targets. Indexing
/// them lets a repeat copy capture the first read seed and drop coverage on the contig
/// that will actually be emitted.
pub fn polish_contigs(
    graph: &mut CdbgGraph,
    reads: impl Iterator<Item = Vec<u8>>,
    min_output_len: usize,
) -> PolishStats {
    let k = 31;
    let n = graph.nodes.len();
    let mut stats = PolishStats::default();

    let total_active_len: usize = graph
        .nodes
        .iter()
        .filter(|node| output_eligible(node, min_output_len))
        .map(|node| node.seq_len as usize)
        .sum();
    let est_seeds = (total_active_len / 25).max(1024);
    let mut seed_index: FxHashMap<u64, (u32, usize)> =
        FxHashMap::with_capacity_and_hasher(est_seeds, Default::default());
    let mut seed_counts: FxHashMap<u64, u32> =
        FxHashMap::with_capacity_and_hasher(est_seeds, Default::default());

    for id in 0..n as u32 {
        if !output_eligible(&graph.nodes[id as usize], min_output_len) {
            continue;
        }
        let seq = graph.get_unitig_seq(id);
        if seq.len() < k {
            continue;
        }

        // Sample seeds every 25 bp
        let mut i = 0;
        while i + k <= seq.len() {
            let h = nthash_canonical(&seq[i..i + k]);
            *seed_counts.entry(h).or_insert(0) += 1;
            seed_index.insert(h, (id, i));
            i += 25;
        }
    }

    // Keep unique seeds
    seed_index.retain(|h, _| *seed_counts.get(h).unwrap_or(&0) == 1);

    // 2. Pre-extract active contig sequences and initialize pileups
    let mut pileups: Vec<Vec<[u16; 4]>> = Vec::with_capacity(n);
    let mut contig_seqs: Vec<Vec<u8>> = Vec::with_capacity(n);

    for id in 0..n {
        let node = &graph.nodes[id];
        if !output_eligible(node, min_output_len) {
            pileups.push(Vec::new());
            contig_seqs.push(Vec::new());
        } else {
            let seq = graph.get_unitig_seq(id as u32);
            let len = seq.len();
            pileups.push(vec![[0u16; 4]; len]);
            contig_seqs.push(seq);
        }
    }

    // 3. Align reads against seeded contig windows (zero-allocation per read)
    for read in reads {
        if read.len() < k {
            continue;
        }

        // Look for anchor seed in read
        let mut best_anchor: Option<(u32, usize, usize, AnchorOrientation)> = None; // (contig_id, contig_pos, read_pos, orientation)
        let mut had_hash_hit = false;

        let mut r_pos = 0;
        while r_pos + k <= read.len() {
            let h = nthash_canonical(&read[r_pos..r_pos + k]);
            if let Some(&(contig_id, contig_pos)) = seed_index.get(&h) {
                had_hash_hit = true;
                let c_seq = &contig_seqs[contig_id as usize];
                let anchor_kmer = &c_seq[contig_pos..contig_pos + k];
                let read_kmer = &read[r_pos..r_pos + k];

                if anchor_kmer.eq_ignore_ascii_case(read_kmer) {
                    best_anchor = Some((contig_id, contig_pos, r_pos, AnchorOrientation::Forward));
                    break;
                } else if is_reverse_complement_exact(anchor_kmer, read_kmer) {
                    best_anchor = Some((contig_id, contig_pos, r_pos, AnchorOrientation::Reverse));
                    break;
                }
            }
            r_pos += 30;
        }

        if best_anchor.is_none() && had_hash_hit {
            stats.rejected_anchors += 1;
        }

        if let Some((contig_id, c_pos, r_pos, orientation)) = best_anchor {
            let pileup = &mut pileups[contig_id as usize];

            match orientation {
                AnchorOrientation::Forward => {
                    let mut voted = false;
                    let read_start_in_contig = (c_pos as i64) - (r_pos as i64);
                    for (idx, &b) in read.iter().enumerate() {
                        let target_pos = read_start_in_contig + idx as i64;
                        if target_pos >= 0 && (target_pos as usize) < pileup.len() {
                            let code = ASCII_TO_2BIT[b as usize] as usize;
                            pileup[target_pos as usize][code] =
                                pileup[target_pos as usize][code].saturating_add(1);
                            stats.fwd_pileup_bases += 1;
                            voted = true;
                        }
                    }
                    if voted {
                        stats.fwd_mapped_reads += 1;
                    }
                }
                AnchorOrientation::Reverse => {
                    let mut voted = false;
                    // Formula: p(i) = c + k - 1 + r - i
                    let base_coord = (c_pos as i64) + (k as i64) - 1 + (r_pos as i64);
                    for (idx, &b) in read.iter().enumerate() {
                        let target_pos = base_coord - idx as i64;
                        if target_pos >= 0 && (target_pos as usize) < pileup.len() {
                            let raw_code = ASCII_TO_2BIT[b as usize] as usize;
                            let code = crate::common::bitseq::COMPLEMENT_2BIT[raw_code] as usize;
                            pileup[target_pos as usize][code] =
                                pileup[target_pos as usize][code].saturating_add(1);
                            stats.rev_pileup_bases += 1;
                            voted = true;
                        }
                    }
                    if voted {
                        stats.rev_mapped_reads += 1;
                    }
                }
            }
        }
    }

    // 4. Polish contig sequences based on pileup counts
    for id in 0..n as u32 {
        if !output_eligible(&graph.nodes[id as usize], min_output_len) {
            continue;
        }

        let mut seq = graph.get_unitig_seq(id);
        let pileup = &pileups[id as usize];
        let mut modified = false;

        for (pos, base) in seq.iter_mut().enumerate() {
            let counts = pileup[pos];
            let total_cov: u16 = counts.iter().sum();

            if total_cov >= 5 {
                // Find dominant allele
                let mut max_idx = 0;
                let mut max_cnt = counts[0];
                for idx in 1..4 {
                    if counts[idx] > max_cnt {
                        max_cnt = counts[idx];
                        max_idx = idx;
                    }
                }

                // If dominant allele has > 75% support and differs from current base
                if (max_cnt as f64) / (total_cov as f64) > 0.75 {
                    let polished_base = TWO_BIT_TO_ASCII[max_idx];
                    if !base.eq_ignore_ascii_case(&polished_base) {
                        *base = polished_base;
                        modified = true;
                        stats.bases_modified += 1;
                    }
                }
            }
        }

        if modified {
            let new_arena_id = graph.arena.push_seq(&seq);
            graph.nodes[id as usize].seq_id = new_arena_id;
        }
    }

    stats
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_reverse_complement_exact() {
        let a = b"ATGCATGCATGCATGCATGCATGCATGCATG";
        let rc_a = crate::common::bitseq::reverse_complement_ascii(a);
        assert!(is_reverse_complement_exact(a, &rc_a));
        assert!(!is_reverse_complement_exact(a, a)); // Forward is not reverse complement (asymmetric)

        // Mutate one base
        let mut mutated = rc_a.clone();
        mutated[0] = if mutated[0] == b'A' { b'T' } else { b'A' };
        assert!(!is_reverse_complement_exact(a, &mutated));

        // Ambiguous base N must be rejected
        let mut with_n = rc_a.clone();
        with_n[5] = b'N';
        assert!(!is_reverse_complement_exact(a, &with_n));
    }

    #[test]
    fn test_reverse_coordinate_formula() {
        // c = 100, k = 31, r = 30
        let c = 100i64;
        let k = 31i64;
        let r = 30i64;
        let base_coord = c + k - 1 + r; // 100 + 31 - 1 + 30 = 160

        // At index i = r = 30, coordinate must be c + k - 1 = 130
        let p_at_r = base_coord - r;
        assert_eq!(p_at_r, c + k - 1);

        // At index i = r + k - 1 = 60, coordinate must be c = 100
        let p_at_end_of_seed = base_coord - (r + k - 1);
        assert_eq!(p_at_end_of_seed, c);

        // At index i = 0, coordinate is 160
        assert_eq!(base_coord - 0, 160);
    }
}
