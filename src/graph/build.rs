//! Compacted De Bruijn Graph construction with solid k-mer filtering and unitig compaction.

use crate::common::bitseq::{reverse_complement_ascii, ReadBank, ASCII_TO_2BIT};
use crate::common::kmer::KmerInt;
use crate::common::sketch::CountMinSketch;
use crate::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, UnitigNode};
use hashbrown::HashMap;
use rayon::prelude::*;
use rustc_hash::FxHasher;
use std::hash::BuildHasherDefault;

type FxHashMap<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;

/// Directed representation of a k-mer: canonical k-mer + orientation flag.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct DirKmer<K: KmerInt> {
    pub kmer: K,
    pub is_rc: bool,
}

impl<K: KmerInt> DirKmer<K> {
    #[inline(always)]
    pub fn new(kmer: K, is_rc: bool) -> Self {
        Self { kmer, is_rc }
    }

    #[inline(always)]
    pub fn twin(self) -> Self {
        Self {
            kmer: self.kmer,
            is_rc: !self.is_rc,
        }
    }
}

/// Builds a Compacted De Bruijn Graph (cDBG) from the `ReadBank`.
pub fn build_cdbg<K: KmerInt>(
    bank: &ReadBank,
    k: usize,
    target_solid_cov: Option<u8>,
) -> CdbgGraph {
    build_cdbg_with_contigs::<K>(bank, &[], k, target_solid_cov)
}

/// Builds a Compacted De Bruijn Graph (cDBG) from the `ReadBank` and optional prior contigs.
pub fn build_cdbg_with_contigs<K: KmerInt>(
    bank: &ReadBank,
    contigs: &[Vec<u8>],
    k: usize,
    target_solid_cov: Option<u8>,
) -> CdbgGraph {
    let est_g = if !contigs.is_empty() {
        contigs.iter().map(|c| c.len() as u64).sum::<u64>().max(50_000)
    } else {
        (bank.total_bases() as u64 / 60).max(50_000)
    };
    let cms = CountMinSketch::for_genome_size(est_g);

    let arenas = [&bank.merged_reads, &bank.pe_forward, &bank.pe_reverse];

    // Pass 1: Populate CountMinSketch with all canonical k-mers using rolling updates in parallel
    for arena in arenas {
        (0..arena.len() as u32).into_par_iter().for_each_init(
            || Vec::with_capacity(300),
            |buf, id| {
                if arena.get_len(id) as usize >= k {
                    arena.extract_ascii_into(id, buf);
                    let seq = &buf[..];
                    let mut fwd = K::from_ascii(&seq[0..k]);
                    let mut rev = fwd.rev_comp(k);
                    let can = if rev < fwd { rev } else { fwd };
                    cms.add(can.to_u64());

                    for i in k..seq.len() {
                        let base = ASCII_TO_2BIT[seq[i] as usize];
                        let (next_fwd, next_rev) = K::roll_next(fwd, rev, base, k);
                        fwd = next_fwd;
                        rev = next_rev;
                        let can = if rev < fwd { rev } else { fwd };
                        cms.add(can.to_u64());
                    }
                }
            },
        );
    }

    // Also populate with prior contigs (with guaranteed solid weight)
    for c in contigs {
        if c.len() >= k {
            let mut fwd = K::from_ascii(&c[0..k]);
            let mut rev = fwd.rev_comp(k);
            let can = if rev < fwd { rev } else { fwd };
            for _ in 0..10 {
                cms.add(can.to_u64());
            }
            for i in k..c.len() {
                let base = ASCII_TO_2BIT[c[i] as usize];
                let (next_fwd, next_rev) = K::roll_next(fwd, rev, base, k);
                fwd = next_fwd;
                rev = next_rev;
                let can = if rev < fwd { rev } else { fwd };
                for _ in 0..10 {
                    cms.add(can.to_u64());
                }
            }
        }
    }

    // Adaptive solid threshold
    let t_solid = target_solid_cov.unwrap_or(2);

    const NUM_SHARDS: usize = 64;

    #[derive(Clone, Copy, Default)]
    struct KmerAcc {
        cov: u16,
        edge_counts: [u16; 8],
    }

    #[inline(always)]
    fn shard_idx<K: KmerInt>(kmer: K) -> usize {
        crate::common::sketch::splitmix64(kmer.to_u64()) as usize & (NUM_SHARDS - 1)
    }

    struct ThreadBuffer<K: KmerInt> {
        shards: Vec<Vec<(K, u8)>>,
    }

    impl<K: KmerInt> ThreadBuffer<K> {
        fn new() -> Self {
            let mut shards = Vec::with_capacity(NUM_SHARDS);
            for _ in 0..NUM_SHARDS {
                shards.push(Vec::with_capacity(2048));
            }
            Self { shards }
        }

        #[inline(always)]
        fn push(&mut self, s_idx: usize, kmer: K, edge: u8, mutexes: &[std::sync::Mutex<FxHashMap<K, KmerAcc>>]) {
            self.shards[s_idx].push((kmer, edge));
            if self.shards[s_idx].len() >= 2048 {
                self.flush_shard(s_idx, mutexes);
            }
        }

        fn flush_shard(&mut self, s_idx: usize, mutexes: &[std::sync::Mutex<FxHashMap<K, KmerAcc>>]) {
            let buf = &mut self.shards[s_idx];
            if buf.is_empty() {
                return;
            }
            let mut lock = mutexes[s_idx].lock().unwrap();
            for &(kmer, edge) in buf.iter() {
                let acc = lock.entry(kmer).or_default();
                acc.cov = acc.cov.saturating_add(1);
                acc.edge_counts[edge as usize] = acc.edge_counts[edge as usize].saturating_add(1);
            }
            buf.clear();
        }

        fn flush_all(&mut self, mutexes: &[std::sync::Mutex<FxHashMap<K, KmerAcc>>]) {
            for s in 0..NUM_SHARDS {
                self.flush_shard(s, mutexes);
            }
        }
    }

    let shard_mutexes: Vec<std::sync::Mutex<FxHashMap<K, KmerAcc>>> = (0..NUM_SHARDS)
        .map(|_| std::sync::Mutex::new(FxHashMap::default()))
        .collect();

    // 1. Prior contigs get heavy backbone support (guaranteed solid transitions)
    for c in contigs {
        if c.len() >= k + 1 {
            let mut u_fwd = K::from_ascii(&c[0..k]);
            let mut u_rev = u_fwd.rev_comp(k);
            let (mut u_can, mut u_rc) = if u_rev < u_fwd { (u_rev, true) } else { (u_fwd, false) };

            for i in 0..=(c.len() - k - 1) {
                let base = ASCII_TO_2BIT[c[i + k] as usize];
                let (v_fwd, v_rev) = K::roll_next(u_fwd, u_rev, base, k);
                let (v_can, v_rc) = if v_rev < v_fwd { (v_rev, true) } else { (v_fwd, false) };

                let edge_u = if !u_rc { base } else { base + 4 };
                let s_u = shard_idx(u_can);
                {
                    let mut lock = shard_mutexes[s_u].lock().unwrap();
                    let e = lock.entry(u_can).or_default();
                    e.cov = e.cov.saturating_add(100);
                    e.edge_counts[edge_u as usize] = e.edge_counts[edge_u as usize].saturating_add(1000);
                }

                let p = ASCII_TO_2BIT[c[i] as usize];
                let rc_p = 3 - p;
                let edge_v = if !v_rc { rc_p + 4 } else { rc_p };
                let s_v = shard_idx(v_can);
                {
                    let mut lock = shard_mutexes[s_v].lock().unwrap();
                    let e = lock.entry(v_can).or_default();
                    e.cov = e.cov.saturating_add(100);
                    e.edge_counts[edge_v as usize] = e.edge_counts[edge_v as usize].saturating_add(1000);
                }

                u_fwd = v_fwd;
                u_rev = v_rev;
                u_can = v_can;
                u_rc = v_rc;
            }
        }
    }

    // 2. Read transitions: chunked parallel extraction using thread buffers and 64 shards
    for arena in arenas {
        let chunk_size = 4096;
        let total_reads = arena.len() as u32;
        (0..total_reads)
            .into_par_iter()
            .chunks(chunk_size)
            .for_each(|chunk| {
                let mut buf = Vec::with_capacity(300);
                let mut tbuf = ThreadBuffer::new();
                for id in chunk {
                    if arena.get_len(id) as usize >= k + 1 {
                        arena.extract_ascii_into(id, &mut buf);
                        let seq = &buf[..];
                        let mut u_fwd = K::from_ascii(&seq[0..k]);
                        let mut u_rev = u_fwd.rev_comp(k);
                        let (mut u_can, mut u_rc) = if u_rev < u_fwd { (u_rev, true) } else { (u_fwd, false) };
                        let mut u_cnt = cms.estimate(u_can.to_u64());

                        for i in 0..=(seq.len() - k - 1) {
                            let base = ASCII_TO_2BIT[seq[i + k] as usize];
                            let (v_fwd, v_rev) = K::roll_next(u_fwd, u_rev, base, k);
                            let (v_can, v_rc) = if v_rev < v_fwd { (v_rev, true) } else { (v_fwd, false) };
                            let v_cnt = cms.estimate(v_can.to_u64());

                            if u_cnt >= t_solid && v_cnt >= t_solid {
                                let edge_u = if !u_rc { base } else { base + 4 };
                                let s_u = shard_idx(u_can);
                                tbuf.push(s_u, u_can, edge_u, &shard_mutexes);

                                let p = ASCII_TO_2BIT[seq[i] as usize];
                                let rc_p = 3 - p;
                                let edge_v = if !v_rc { rc_p + 4 } else { rc_p };
                                let s_v = shard_idx(v_can);
                                tbuf.push(s_v, v_can, edge_v, &shard_mutexes);
                            }

                            u_fwd = v_fwd;
                            u_rev = v_rev;
                            u_can = v_can;
                            u_rc = v_rc;
                            u_cnt = v_cnt;
                        }
                    }
                }
                tbuf.flush_all(&shard_mutexes);
            });
    }

    // Drop CMS immediately to release 32-256 MB RAM before graph construction
    drop(cms);

    // 3. Convert KmerAcc into compact SolidKmerNode with 8-bit extension bitmasks in parallel
    #[derive(Clone, Copy, Default, Debug)]
    struct SolidKmerNode {
        cov: u16,
        mask: u8,
        visited: u8,
    }

    let min_trans_cov = if let Some(t) = target_solid_cov { t as u16 } else { 3 };
    let mut shards: Vec<FxHashMap<K, SolidKmerNode>> = shard_mutexes
        .into_par_iter()
        .map(|mutex| {
            let map = mutex.into_inner().unwrap();
            let mut solid_shard = FxHashMap::with_capacity_and_hasher(map.len(), Default::default());
            for (kmer, acc) in map {
                let mut mask = 0u8;
                for e in 0..8 {
                    if acc.edge_counts[e] >= min_trans_cov {
                        mask |= 1 << e;
                    }
                }
                if mask != 0 {
                    solid_shard.insert(
                        kmer,
                        SolidKmerNode {
                            cov: (acc.cov / 2).max(1),
                            mask,
                            visited: 0,
                        },
                    );
                }
            }
            solid_shard
        })
        .collect();

    let mut graph = CdbgGraph::new(k);

    // 4. Compaction: Collect all boundary directed k-mers (in_deg != 1 || out_deg != 1)
    let mut boundary_nodes = Vec::new();
    for shard in &shards {
        for (&kmer, node) in shard {
            let fwd_out = (node.mask & 0x0F).count_ones() as u8;
            let fwd_in = ((node.mask >> 4) & 0x0F).count_ones() as u8;
            if fwd_in != 1 || fwd_out != 1 {
                boundary_nodes.push(DirKmer::new(kmer, false));
            }

            let rev_out = ((node.mask >> 4) & 0x0F).count_ones() as u8;
            let rev_in = (node.mask & 0x0F).count_ones() as u8;
            if rev_in != 1 || rev_out != 1 {
                boundary_nodes.push(DirKmer::new(kmer, true));
            }
        }
    }

    for u in boundary_nodes {
        let s_u = shard_idx(u.kmer);

        for b in 0..4u8 {
            let node_u = match shards[s_u].get(&u.kmer) {
                Some(&n) => n,
                None => break,
            };

            let u_out_mask = if !u.is_rc { node_u.mask & 0x0F } else { (node_u.mask >> 4) & 0x0F };
            let u_visited = if !u.is_rc { node_u.visited & 0x0F } else { (node_u.visited >> 4) & 0x0F };

            if (u_out_mask & (1 << b)) != 0 && (u_visited & (1 << b)) == 0 {
                let edge_bit_u = if !u.is_rc { b } else { b + 4 };
                shards[s_u].get_mut(&u.kmer).unwrap().visited |= 1 << edge_bit_u;

                let mut path_seq = u.kmer.to_ascii(k);
                if u.is_rc {
                    path_seq = reverse_complement_ascii(&path_seq);
                }

                let (fwd, rev) = if !u.is_rc {
                    (u.kmer, u.kmer.rev_comp(k))
                } else {
                    (u.kmer.rev_comp(k), u.kmer)
                };
                let first_b = fwd.to_ascii(k)[0];
                let rc_first = 3 - ASCII_TO_2BIT[first_b as usize];

                let (next_fwd, next_rev) = K::roll_next(fwd, rev, b, k);
                let (next_can, next_is_rc) = if next_rev < next_fwd { (next_rev, true) } else { (next_fwd, false) };
                let mut curr = DirKmer::new(next_can, next_is_rc);

                let twin_edge = if !curr.is_rc { rc_first + 4 } else { rc_first };
                if let Some(n_curr) = shards[shard_idx(curr.kmer)].get_mut(&curr.kmer) {
                    n_curr.visited |= 1 << twin_edge;
                }

                path_seq.push(match b { 0 => b'A', 1 => b'C', 2 => b'G', _ => b'T' });
                let mut total_cov = node_u.cov as u64;
                let mut cov_points = 1u64;

                loop {
                    let s_c = shard_idx(curr.kmer);
                    let n_c = match shards[s_c].get(&curr.kmer) {
                        Some(&n) => n,
                        None => break,
                    };

                    let c_out_mask = if !curr.is_rc { n_c.mask & 0x0F } else { (n_c.mask >> 4) & 0x0F };
                    let c_in_mask = if !curr.is_rc { (n_c.mask >> 4) & 0x0F } else { n_c.mask & 0x0F };

                    if c_in_mask.count_ones() == 1 && c_out_mask.count_ones() == 1 {
                        let next_b = c_out_mask.trailing_zeros() as u8;
                        let edge_bit = if !curr.is_rc { next_b } else { next_b + 4 };

                        if (n_c.visited & (1 << edge_bit)) != 0 {
                            break;
                        }
                        shards[s_c].get_mut(&curr.kmer).unwrap().visited |= 1 << edge_bit;

                        let (c_fwd, c_rev) = if !curr.is_rc {
                            (curr.kmer, curr.kmer.rev_comp(k))
                        } else {
                            (curr.kmer.rev_comp(k), curr.kmer)
                        };

                        let first_char = c_fwd.to_ascii(k)[0];
                        let rc_p = 3 - ASCII_TO_2BIT[first_char as usize];

                        let (n_fwd, n_rev) = K::roll_next(c_fwd, c_rev, next_b, k);
                        let (n_can, n_is_rc) = if n_rev < n_fwd { (n_rev, true) } else { (n_fwd, false) };
                        let next_curr = DirKmer::new(n_can, n_is_rc);

                        let twin_e = if !next_curr.is_rc { rc_p + 4 } else { rc_p };
                        if let Some(n_next) = shards[shard_idx(next_curr.kmer)].get_mut(&next_curr.kmer) {
                            n_next.visited |= 1 << twin_e;
                        }

                        path_seq.push(match next_b { 0 => b'A', 1 => b'C', 2 => b'G', _ => b'T' });
                        total_cov += n_c.cov as u64;
                        cov_points += 1;

                        curr = next_curr;
                    } else {
                        break;
                    }
                }

                let rc_seq = reverse_complement_ascii(&path_seq);
                let final_seq = if rc_seq < path_seq { rc_seq } else { path_seq };

                let seq_id = graph.arena.push_seq(&final_seq);
                let node_id = graph.nodes.len() as u32;
                let cov = (total_cov as f32 / cov_points as f32).max(1.0);
                let node = UnitigNode::new(node_id, seq_id, final_seq.len() as u32, cov);
                graph.nodes.push(node);
            }
        }
    }

    // 5. Walk isolated circular components (plasmids with all 1-in, 1-out)
    for s in 0..NUM_SHARDS {
        let kmers_in_shard: Vec<K> = shards[s].keys().copied().collect();
        for kmer in kmers_in_shard {
            let n = match shards[s].get(&kmer) {
                Some(&node) => node,
                None => continue,
            };

            let fwd_out = n.mask & 0x0F;
            let fwd_in = (n.mask >> 4) & 0x0F;
            if fwd_out.count_ones() == 1 && fwd_in.count_ones() == 1 {
                let b = fwd_out.trailing_zeros() as u8;
                if (n.visited & (1 << b)) == 0 {
                    shards[s].get_mut(&kmer).unwrap().visited |= 1 << b;

                    let start_kmer = kmer;
                    let mut path_seq = kmer.to_ascii(k);

                    let fwd = kmer;
                    let rev = kmer.rev_comp(k);
                    let (next_fwd, next_rev) = K::roll_next(fwd, rev, b, k);
                    let (next_can, next_is_rc) = if next_rev < next_fwd { (next_rev, true) } else { (next_fwd, false) };
                    let mut curr = DirKmer::new(next_can, next_is_rc);

                    let first_char = fwd.to_ascii(k)[0];
                    let rc_first = 3 - ASCII_TO_2BIT[first_char as usize];
                    let twin_edge = if !curr.is_rc { rc_first + 4 } else { rc_first };
                    if let Some(n_curr) = shards[shard_idx(curr.kmer)].get_mut(&curr.kmer) {
                        n_curr.visited |= 1 << twin_edge;
                    }

                    path_seq.push(match b { 0 => b'A', 1 => b'C', 2 => b'G', _ => b'T' });
                    let mut total_cov = n.cov as u64;
                    let mut cov_points = 1u64;

                    while curr.kmer != start_kmer {
                        let s_c = shard_idx(curr.kmer);
                        let n_c = match shards[s_c].get(&curr.kmer) {
                            Some(&nc) => nc,
                            None => break,
                        };

                        let c_out_mask = if !curr.is_rc { n_c.mask & 0x0F } else { (n_c.mask >> 4) & 0x0F };
                        let c_in_mask = if !curr.is_rc { (n_c.mask >> 4) & 0x0F } else { n_c.mask & 0x0F };

                        if c_in_mask.count_ones() == 1 && c_out_mask.count_ones() == 1 {
                            let next_b = c_out_mask.trailing_zeros() as u8;
                            let edge_bit = if !curr.is_rc { next_b } else { next_b + 4 };

                            if (n_c.visited & (1 << edge_bit)) != 0 {
                                break;
                            }
                            shards[s_c].get_mut(&curr.kmer).unwrap().visited |= 1 << edge_bit;

                            let (c_fwd, c_rev) = if !curr.is_rc { (curr.kmer, curr.kmer.rev_comp(k)) } else { (curr.kmer.rev_comp(k), curr.kmer) };
                            let first_b = c_fwd.to_ascii(k)[0];
                            let rc_p = 3 - ASCII_TO_2BIT[first_b as usize];

                            let (n_fwd, n_rev) = K::roll_next(c_fwd, c_rev, next_b, k);
                            let (n_can, n_is_rc) = if n_rev < n_fwd { (n_rev, true) } else { (n_fwd, false) };
                            let next_curr = DirKmer::new(n_can, n_is_rc);

                            let twin_e = if !next_curr.is_rc { rc_p + 4 } else { rc_p };
                            if let Some(n_next) = shards[shard_idx(next_curr.kmer)].get_mut(&next_curr.kmer) {
                                n_next.visited |= 1 << twin_e;
                            }

                            path_seq.push(match next_b { 0 => b'A', 1 => b'C', 2 => b'G', _ => b'T' });
                            total_cov += n_c.cov as u64;
                            cov_points += 1;

                            curr = next_curr;
                        } else {
                            break;
                        }
                    }

                    // Only add if it closed a genuine circular loop
                    if curr.kmer == start_kmer && path_seq.len() > k {
                        if path_seq.len() > k && path_seq[path_seq.len() - k..] == path_seq[..k] {
                            path_seq.truncate(path_seq.len() - k);
                        } else if path_seq.len() > k - 1 && path_seq[path_seq.len() - (k - 1)..] == path_seq[..(k - 1)] {
                            path_seq.truncate(path_seq.len() - (k - 1));
                        }

                        let seq_id = graph.arena.push_seq(&path_seq);
                        let node_id = graph.nodes.len() as u32;
                        let cov = (total_cov as f32 / cov_points as f32).max(1.0);
                        let mut node = UnitigNode::new(node_id, seq_id, path_seq.len() as u32, cov);
                        node.is_circular = true;
                        let self_dir = DirUnitigId::new(node_id, false);
                        node.add_edge(false, GraphEdge::new(self_dir, (k - 1) as u8, 10));
                        graph.nodes.push(node);
                    }
                }
            }
        }
    }

    // 6. Connect graph edges between unitigs by exact k and (k - 1) sequence overlap
    for &ov_len in &[k, k.saturating_sub(1)] {
        if ov_len == 0 {
            continue;
        }
        let mut prefix_map: FxHashMap<Vec<u8>, smallvec::SmallVec<[DirUnitigId; 4]>> = FxHashMap::default();

        for id in 0..graph.nodes.len() as u32 {
            let fwd_seq = graph.get_unitig_seq(id);
            if fwd_seq.len() >= ov_len {
                let fwd_prefix = fwd_seq[..ov_len].to_vec();
                prefix_map.entry(fwd_prefix).or_default().push(DirUnitigId::new(id, false));

                let rc_seq = graph.get_stranded_seq(DirUnitigId::new(id, true));
                let rc_prefix = rc_seq[..ov_len].to_vec();
                prefix_map.entry(rc_prefix).or_default().push(DirUnitigId::new(id, true));
            }
        }

        for id in 0..graph.nodes.len() as u32 {
            let fwd_seq = graph.get_unitig_seq(id);
            let rc_seq = graph.get_stranded_seq(DirUnitigId::new(id, true));

            if fwd_seq.len() >= ov_len {
                let fwd_suffix = &fwd_seq[fwd_seq.len() - ov_len..];
                if let Some(targets) = prefix_map.get(fwd_suffix) {
                    for &target in targets {
                        if target.id() != id || graph.nodes[id as usize].is_circular {
                            graph.nodes[id as usize].add_edge(false, GraphEdge::new(target, ov_len as u8, 10));
                        }
                    }
                }

                let rc_suffix = &rc_seq[rc_seq.len() - ov_len..];
                if let Some(targets) = prefix_map.get(rc_suffix) {
                    for &target in targets {
                        if target.id() != id || graph.nodes[id as usize].is_circular {
                            graph.nodes[id as usize].add_edge(true, GraphEdge::new(target, ov_len as u8, 10));
                        }
                    }
                }
            }
        }
    }

    graph
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::sketch::splitmix64;

    #[test]
    fn test_verification_gate_circular_repeat_compaction() {
        // Verification Gate for Phase 3:
        // Construct cDBG from simulated reads (with 1% random sequencing errors)
        // of a 10 kb circular sequence containing one 200 bp exact repeat:
        // Verify that singleton error k-mers are filtered out and the graph compacts into
        // exactly 3 unitigs (the repeat unitig with in-degree 2 / out-degree 2, plus 2 flanking unique unitigs).

        let k = 31;
        let repeat_len = 200;
        let flank1_len = 4800;
        let flank2_len = 4800;

        let alphabet = [b'A', b'C', b'G', b'T'];

        // Generate synthetic components using splitmix64
        let gen_seq = |seed: u64, len: usize| -> Vec<u8> {
            let mut buf = Vec::with_capacity(len);
            for i in 0..len {
                let r = splitmix64(seed + i as u64);
                buf.push(alphabet[(r % 4) as usize]);
            }
            buf
        };

        let repeat = gen_seq(0x1111_2222_3333_4444, repeat_len);
        let flank1 = gen_seq(0x5555_6666_7777_8888, flank1_len);
        let flank2 = gen_seq(0x9999_AAAA_BBBB_CCCC, flank2_len);

        // Circular genome: R -> F1 -> R -> F2 -> (circular loop back to R)
        let mut genome = Vec::new();
        genome.extend_from_slice(&repeat);
        genome.extend_from_slice(&flank1);
        genome.extend_from_slice(&repeat);
        genome.extend_from_slice(&flank2);

        let g_len = genome.len(); // 10,000 bp
        assert_eq!(g_len, 10_000);

        // Simulate 40x depth of 150bp reads with 1% random error
        let mut bank = ReadBank::new();
        let num_reads = (g_len * 40) / 150;

        let mut read_rng: u64 = 0xFEEDBEEF;
        for i in 0..num_reads {
            read_rng = splitmix64(read_rng + i as u64);
            let start = (read_rng as usize) % g_len;

            let mut read = Vec::with_capacity(150);
            for j in 0..150 {
                let base = genome[(start + j) % g_len];
                let err_roll = splitmix64(read_rng ^ (j as u64)) % 100;
                if err_roll == 0 {
                    let base_idx = match base {
                        b'A' => 0,
                        b'C' => 1,
                        b'G' => 2,
                        _ => 3,
                    };
                    let alt_offset = 1 + (splitmix64(read_rng + (j as u64) * 31) % 3) as usize;
                    let new_base = alphabet[(base_idx + alt_offset) % 4];
                    read.push(new_base);
                } else {
                    read.push(base);
                }
            }
            bank.merged_reads.push_seq(&read);
        }

        // Build cDBG with solid cutoff = 5 (filters errors while keeping 40x genome)
        let graph = build_cdbg::<u64>(&bank, k, Some(5));

        for (i, n) in graph.nodes.iter().enumerate() {
            eprintln!("node {} len {} circular={}", i, n.seq_len, n.is_circular);
        }

        // Assert that exactly 3 unitigs are produced!
        assert_eq!(
            graph.num_active_nodes(),
            3,
            "cDBG must compact into exactly 3 unitigs: repeat + 2 flanking regions. Found {}",
            graph.num_active_nodes()
        );
    }
}
