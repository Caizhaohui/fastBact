//! Paired-end read mapping and empirical insert size distribution estimation.

use crate::common::bitseq::ReadBank;
use crate::common::kmer::KmerInt;
use crate::graph::cdbg::{CdbgGraph, DirUnitigId};
use hashbrown::HashMap;
use rustc_hash::FxHasher;
use std::hash::BuildHasherDefault;

type FxHashMap<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;

/// Empirical insert size distribution parameters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InsertSizeDist {
    pub median: f64,
    pub mad: f64,
    pub min_insert: usize,
    pub max_insert: usize,
}

impl Default for InsertSizeDist {
    fn default() -> Self {
        Self {
            median: 450.0,
            mad: 50.0,
            min_insert: 250,
            max_insert: 700,
        }
    }
}

/// Sparse matrix recording observed paired-end links between unitigs.
#[derive(Clone, Debug, Default)]
pub struct PairLinkMatrix {
    /// Number of read pairs linking (u1 -> u2)
    pub links: FxHashMap<(DirUnitigId, DirUnitigId), u32>,
}

impl PairLinkMatrix {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_link(&mut self, u1: DirUnitigId, u2: DirUnitigId) {
        *self.links.entry((u1, u2)).or_insert(0) += 1;
    }

    pub fn get_support(&self, u1: DirUnitigId, u2: DirUnitigId) -> u32 {
        *self.links.get(&(u1, u2)).unwrap_or(&0)
    }
}

/// Indexes unique terminal k-mers (first/last 1000 bp) of all active unitigs.
pub fn index_terminal_kmers<K: KmerInt>(
    graph: &CdbgGraph,
    k: usize,
) -> FxHashMap<u64, (DirUnitigId, u32)> {
    let n_active = graph.nodes.iter().filter(|n| !n.is_deleted).count();
    let est_kmers = (n_active * 2000).max(1024);
    let mut kmer_counts: FxHashMap<u64, u32> = FxHashMap::with_capacity_and_hasher(est_kmers, Default::default());
    let mut kmer_pos: FxHashMap<u64, (DirUnitigId, u32)> = FxHashMap::with_capacity_and_hasher(est_kmers, Default::default());

    for id in 0..graph.nodes.len() as u32 {
        let node = &graph.nodes[id as usize];
        if node.is_deleted {
            continue;
        }

        let seq = graph.get_unitig_seq(id);
        let len = seq.len();
        if len < k {
            continue;
        }

        let terminal_len = 1000.min(len - k + 1);

        // Index 5' terminal region
        for i in 0..terminal_len {
            let kmer = K::from_ascii(&seq[i..i + k]);
            let (can, is_rc) = kmer.canonical(k);
            let h = can.to_u64();
            *kmer_counts.entry(h).or_insert(0) += 1;
            kmer_pos.insert(h, (DirUnitigId::new(id, is_rc), i as u32));
        }

        // Index 3' terminal region
        if len > 1000 + k {
            let start = len - 1000;
            for i in start..=len - k {
                let kmer = K::from_ascii(&seq[i..i + k]);
                let (can, is_rc) = kmer.canonical(k);
                let h = can.to_u64();
                *kmer_counts.entry(h).or_insert(0) += 1;
                kmer_pos.insert(h, (DirUnitigId::new(id, is_rc), i as u32));
            }
        }
    }

    // Keep only uniquely occurring k-mers
    kmer_pos.retain(|h, _| *kmer_counts.get(h).unwrap_or(&0) == 1);
    kmer_pos
}

/// Maps paired-end reads to unitigs and populates `PairLinkMatrix` and `InsertSizeDist`.
pub fn map_paired_reads<K: KmerInt>(
    graph: &CdbgGraph,
    bank: &ReadBank,
    k: usize,
) -> (PairLinkMatrix, InsertSizeDist) {
    let index = index_terminal_kmers::<K>(graph, k);
    let mut matrix = PairLinkMatrix::new();
    let mut insert_sizes = Vec::new();

    let num_pairs = bank.pe_forward.len().min(bank.pe_reverse.len());

    let map_read = |r: &[u8]| -> Option<(DirUnitigId, u32)> {
        for step in [0, 20, 40, 60, 80, 100] {
            if step + k <= r.len() {
                let km = K::from_ascii(&r[step..step + k]);
                let (can, is_rc) = km.canonical(k);
                if let Some(&(uid, pos)) = index.get(&can.to_u64()) {
                    let dir = if is_rc { uid.twin() } else { uid };
                    return Some((dir, pos));
                }
            }
        }
        None
    };

    let mut r1 = Vec::with_capacity(300);
    let mut r2 = Vec::with_capacity(300);

    for i in 0..num_pairs as u32 {
        bank.pe_forward.extract_ascii_into(i, &mut r1);
        bank.pe_reverse.extract_ascii_into(i, &mut r2);

        if r1.len() < k || r2.len() < k {
            continue;
        }

        let map1 = map_read(&r1);
        let map2 = map_read(&r2);

        if let (Some((u1, pos1)), Some((u2, pos2))) = (map1, map2) {
            let u2_forward = u2.twin();

            if u1.id() == u2.id() {
                // Intra-unitig pair: estimate insert size
                let node_len = graph.nodes[u1.id() as usize].seq_len;
                if node_len >= 1500 {
                    let d = (pos2 as i64 - pos1 as i64).abs() as usize + r1.len();
                    if d >= 100 && d <= 1500 {
                        insert_sizes.push(d);
                    }
                }
            } else {
                // Inter-unitig pair: record links
                matrix.add_link(u1, u2_forward);
                matrix.add_link(u2_forward.twin(), u1.twin());
            }
        }
    }

    let dist = if insert_sizes.len() >= 20 {
        insert_sizes.sort_unstable();
        let med = insert_sizes[insert_sizes.len() / 2] as f64;
        let mut diffs: Vec<f64> = insert_sizes.iter().map(|&x| (x as f64 - med).abs()).collect();
        diffs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mad = diffs[diffs.len() / 2].max(20.0);

        let min_ins = (med - 3.5 * mad).max(50.0) as usize;
        let max_ins = (med + 3.5 * mad) as usize;

        InsertSizeDist {
            median: med,
            mad,
            min_insert: min_ins,
            max_insert: max_ins,
        }
    } else {
        InsertSizeDist::default()
    };

    (matrix, dist)
}
