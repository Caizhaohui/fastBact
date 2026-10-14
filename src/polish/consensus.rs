//! In-memory pileup consensus polishing.
//!
//! Replaces BWA + Pilon by anchoring reads with exact 31-mer seeds and accumulating
//! per-position nucleotide base counts, substituting reference alleles with >75% support.

use crate::common::bitseq::{ASCII_TO_2BIT, TWO_BIT_TO_ASCII};
use crate::common::kmer::nthash_canonical;
use crate::graph::cdbg::CdbgGraph;
use hashbrown::HashMap;
use rustc_hash::FxHasher;
use std::hash::BuildHasherDefault;

type FxHashMap<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;

/// Anchors reads from the read pool onto assembled contigs and performs fast in-memory
/// pileup error correction.
pub fn polish_contigs(
    graph: &mut CdbgGraph,
    reads: impl Iterator<Item = Vec<u8>>,
) {
    let k = 31;
    let n = graph.nodes.len();

    let total_active_len: usize = graph.nodes.iter().filter(|n| !n.is_deleted).map(|n| n.seq_len as usize).sum();
    let est_seeds = (total_active_len / 25).max(1024);
    let mut seed_index: FxHashMap<u64, (u32, usize)> = FxHashMap::with_capacity_and_hasher(est_seeds, Default::default());
    let mut seed_counts: FxHashMap<u64, u32> = FxHashMap::with_capacity_and_hasher(est_seeds, Default::default());

    for id in 0..n as u32 {
        if graph.nodes[id as usize].is_deleted {
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
        if node.is_deleted {
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
        let mut best_anchor: Option<(u32, usize, usize)> = None; // (contig_id, contig_pos, read_pos)

        let mut r_pos = 0;
        while r_pos + k <= read.len() {
            let h = nthash_canonical(&read[r_pos..r_pos + k]);
            if let Some(&(contig_id, contig_pos)) = seed_index.get(&h) {
                best_anchor = Some((contig_id, contig_pos, r_pos));
                break;
            }
            r_pos += 30;
        }

        if let Some((contig_id, c_pos, r_pos)) = best_anchor {
            let c_seq = &contig_seqs[contig_id as usize];
            let pileup = &mut pileups[contig_id as usize];

            let read_start_in_contig = (c_pos as i64) - (r_pos as i64);
            let anchor_kmer = &c_seq[c_pos..c_pos + k];
            let read_kmer = &read[r_pos..r_pos + k];
            let is_fwd = anchor_kmer.eq_ignore_ascii_case(read_kmer);

            let r_len = read.len();
            for (idx, &b) in read.iter().enumerate() {
                let target_pos = if is_fwd {
                    read_start_in_contig + idx as i64
                } else {
                    read_start_in_contig + (r_len - 1 - idx) as i64
                };

                if target_pos >= 0 && (target_pos as usize) < pileup.len() {
                    let code = if is_fwd {
                        ASCII_TO_2BIT[b as usize] as usize
                    } else {
                        crate::common::bitseq::COMPLEMENT_2BIT[ASCII_TO_2BIT[b as usize] as usize] as usize
                    };
                    pileup[target_pos as usize][code] =
                        pileup[target_pos as usize][code].saturating_add(1);
                }
            }
        }
    }

    // 4. Polish contig sequences based on pileup counts
    for id in 0..n as u32 {
        if graph.nodes[id as usize].is_deleted {
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
                    }
                }
            }
        }

        if modified {
            let new_arena_id = graph.arena.push_seq(&seq);
            graph.nodes[id as usize].seq_id = new_arena_id;
        }
    }
}
