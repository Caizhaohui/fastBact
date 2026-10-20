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

/// Compact, coordinate-faithful evidence of an authentic paired-end read fragment.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct ReadPairEvidence {
    /// 32-bit unique read pair index in ReadBank (enables strict deduplication)
    pub pair_id: u32,
    /// Distance in bp from R1 5' start to u1 3' exit: len(u1) - s1
    pub dist_u1_exit: u32,
    /// Distance in bp from u2 5' entry to R2 3' end: e2
    pub dist_u2_entry: u32,
    /// Read 1 physical length in bp (u16)
    pub r1_len: u16,
    /// Read 2 physical length in bp (u16)
    pub r2_len: u16,
}

/// Sparse matrix recording observed paired-end links and coordinate evidence between unitigs.
#[derive(Clone, Debug, Default)]
pub struct PairLinkMatrix {
    /// Directed unitig pair links with exact coordinate evidence
    pub links: FxHashMap<(DirUnitigId, DirUnitigId), Vec<ReadPairEvidence>>,
}

impl PairLinkMatrix {
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a synthetic link with default coordinates for testing.
    pub fn add_link(&mut self, u1: DirUnitigId, u2: DirUnitigId) {
        let pair_id = self.links.get(&(u1, u2)).map_or(0, |v| v.len() as u32);
        self.add_evidence(
            u1,
            u2,
            ReadPairEvidence {
                pair_id,
                dist_u1_exit: 120,
                dist_u2_entry: 120,
                r1_len: 150,
                r2_len: 150,
            },
        );
    }

    pub fn add_evidence(&mut self, u1: DirUnitigId, u2: DirUnitigId, ev: ReadPairEvidence) {
        self.links.entry((u1, u2)).or_default().push(ev);
    }

    pub fn get_support(&self, u1: DirUnitigId, u2: DirUnitigId) -> u32 {
        self.links.get(&(u1, u2)).map_or(0, |v| v.len() as u32)
    }

    pub fn get_evidence(&self, u1: DirUnitigId, u2: DirUnitigId) -> Option<&Vec<ReadPairEvidence>> {
        self.links.get(&(u1, u2))
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

        let total_kmers = len - k + 1;
        let mut index_kmer = |i: usize| {
            let kmer = K::from_ascii(&seq[i..i + k]);
            let (can, is_rc) = kmer.canonical(k);
            let h = can.to_u64();
            *kmer_counts.entry(h).or_insert(0) += 1;
            kmer_pos.insert(h, (DirUnitigId::new(id, is_rc), i as u32));
        };

        if total_kmers <= 2000 {
            for i in 0..total_kmers {
                index_kmer(i);
            }
        } else {
            for i in 0..1000 {
                index_kmer(i);
            }
            for i in (total_kmers - 1000)..total_kmers {
                index_kmer(i);
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

    let map_read = |r: &[u8]| -> Option<(DirUnitigId, u32, usize)> {
        for step in [0, 20, 40, 60, 80, 100] {
            if step + k <= r.len() {
                let km = K::from_ascii(&r[step..step + k]);
                let (can, is_rc) = km.canonical(k);
                if let Some(&(uid, pos)) = index.get(&can.to_u64()) {
                    let dir = if is_rc { uid.twin() } else { uid };
                    return Some((dir, pos, step));
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

        if let (Some((u1, pos1, step1)), Some((u2, pos2, step2))) = (map1, map2) {
            let u2_forward = u2.twin();

            if u1.id() == u2.id() {
                // Intra-unitig pair: estimate insert size
                let node_len = graph.nodes[u1.id() as usize].seq_len;
                if node_len >= 1500 {
                    // Paired-end Illumina reads are inward-facing (FR orientation).
                    // Subcase A: R1 forward, R2 reverse (u1.is_rc() == false, u2.is_rc() == true)
                    // Subcase B: R1 reverse, R2 forward (u1.is_rc() == true, u2.is_rc() == false)
                    // Same-strand pairs (FF / RR) are non-concordant and rejected.
                    let span: i64 = match (u1.is_rc(), u2.is_rc()) {
                        (false, true) => {
                            let s1 = pos1 as i64 - step1 as i64;
                            let e2 = pos2 as i64 + k as i64 + step2 as i64;
                            e2 - s1
                        }
                        (true, false) => {
                            let s2 = pos2 as i64 - step2 as i64;
                            let e1 = pos1 as i64 + k as i64 + step1 as i64;
                            e1 - s2
                        }
                        _ => -1,
                    };

                    let min_len = (r1.len().max(r2.len()) as i64).max(100);
                    if span >= min_len && span <= 1500 {
                        insert_sizes.push(span as usize);
                    }
                }
            } else {
                // Inter-unitig pair: record links with coordinate evidence
                let node1_len = graph.nodes[u1.id() as usize].seq_len as i64;
                let node2_len = graph.nodes[u2.id() as usize].seq_len as i64;

                let d1: i64 = if u1.is_rc() {
                    pos1 as i64 + k as i64 + step1 as i64
                } else {
                    node1_len - (pos1 as i64 - step1 as i64)
                };

                let d2: i64 = if u2_forward.is_rc() {
                    node2_len - (pos2 as i64 - step2 as i64)
                } else {
                    pos2 as i64 + k as i64 + step2 as i64
                };

                if let (Ok(d1_u32), Ok(d2_u32), Ok(r1_len_u16), Ok(r2_len_u16)) = (
                    u32::try_from(d1.max(0)),
                    u32::try_from(d2.max(0)),
                    u16::try_from(r1.len()),
                    u16::try_from(r2.len()),
                ) {
                    matrix.add_evidence(
                        u1,
                        u2_forward,
                        ReadPairEvidence {
                            pair_id: i,
                            dist_u1_exit: d1_u32,
                            dist_u2_entry: d2_u32,
                            r1_len: r1_len_u16,
                            r2_len: r2_len_u16,
                        },
                    );
                    matrix.add_evidence(
                        u2_forward.twin(),
                        u1.twin(),
                        ReadPairEvidence {
                            pair_id: i,
                            dist_u1_exit: d2_u32,
                            dist_u2_entry: d1_u32,
                            r1_len: r2_len_u16,
                            r2_len: r1_len_u16,
                        },
                    );
                }
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
