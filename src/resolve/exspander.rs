//! SPAdes-style exSPAnder paired-end repeat resolution engine.
//!
//! Traverses across repeats longer than k_max using paired-end read insert size constraints,
//! connecting unique flanking unitigs across repetitive regions.

use crate::graph::cdbg::{CdbgGraph, DirUnitigId, UnitigNode};
use crate::resolve::mapper::{InsertSizeDist, PairLinkMatrix};
use hashbrown::HashSet;

/// Configuration parameters for exSPAnder repeat resolution.
#[derive(Clone, Copy, Debug)]
pub struct ExspanderParams {
    /// Minimum raw supporting read-pair links to extend path across repeat
    pub min_support: u32,
    /// Winner-take-all confidence ratio: score(best) >= decision_threshold * score(second)
    pub decision_threshold: f64,
}

impl Default for ExspanderParams {
    fn default() -> Self {
        Self {
            min_support: 5,
            decision_threshold: 2.5,
        }
    }
}

fn extend_one_direction(
    mut curr_dir: DirUnitigId,
    graph: &CdbgGraph,
    matrix: &PairLinkMatrix,
    params: ExspanderParams,
    used_unique: &mut HashSet<u32>,
) -> Vec<DirUnitigId> {
    let mut ext = Vec::new();

    loop {
        let curr_node = &graph.nodes[curr_dir.id() as usize];
        let outgoing = curr_node.edges(curr_dir.is_rc()).clone();

        if outgoing.is_empty() {
            break;
        }

        // Case A: 1 outgoing edge
        if outgoing.len() == 1 {
            let target = outgoing[0].target;
            let target_node = &graph.nodes[target.id() as usize];
            if target_node.is_deleted || target_node.is_circular {
                break;
            }

            let next_edges = target_node.edges(target.is_rc()).clone();
            if next_edges.len() > 1 {
                // Target is a branching repeat node! Evaluate candidates exiting the repeat
                let mut best_cand: Option<DirUnitigId> = None;
                let mut best_support = 0u32;
                let mut second_support = 0u32;

                for next_e in &next_edges {
                    let cand = next_e.target;
                    let supp = matrix.get_support(curr_dir, cand);
                    if supp > best_support {
                        second_support = best_support;
                        best_support = supp;
                        best_cand = Some(cand);
                    } else if supp > second_support {
                        second_support = supp;
                    }
                }

                if let Some(cand) = best_cand {
                    let ratio = if second_support > 0 {
                        (best_support as f64) / (second_support as f64)
                    } else {
                        100.0
                    };

                    if best_support >= params.min_support && ratio >= params.decision_threshold {
                        if !used_unique.contains(&cand.id()) && !graph.nodes[cand.id() as usize].is_deleted {
                            used_unique.insert(cand.id());
                            ext.push(target);
                            ext.push(cand);
                            curr_dir = cand;
                            continue;
                        }
                    }
                }
                break;
            } else if next_edges.len() == 1 {
                // Unbranched continuation
                if target_node.seq_len >= 500 {
                    if used_unique.contains(&target.id()) {
                        break;
                    }
                    used_unique.insert(target.id());
                }
                ext.push(target);
                curr_dir = target;
                continue;
            } else {
                // Dead end at target
                ext.push(target);
                break;
            }
        }

        // Case B: multiple outgoing edges directly from curr_dir
        let mut best_cand: Option<DirUnitigId> = None;
        let mut best_support = 0u32;
        let mut second_support = 0u32;

        for e in &outgoing {
            let cand = e.target;
            let supp = matrix.get_support(curr_dir, cand);
            if supp > best_support {
                second_support = best_support;
                best_support = supp;
                best_cand = Some(cand);
            } else if supp > second_support {
                second_support = supp;
            }
        }

        if let Some(cand) = best_cand {
            let ratio = if second_support > 0 {
                (best_support as f64) / (second_support as f64)
            } else {
                100.0
            };

            if best_support >= params.min_support && ratio >= params.decision_threshold {
                if !used_unique.contains(&cand.id()) && !graph.nodes[cand.id() as usize].is_deleted {
                    if graph.nodes[cand.id() as usize].seq_len >= 500 {
                        used_unique.insert(cand.id());
                    }
                    ext.push(cand);
                    curr_dir = cand;
                    continue;
                }
            }
        }

        break;
    }

    ext
}

/// Resolves repeats in `graph` using paired-end link constraints from `PairLinkMatrix`.
pub fn resolve_repeats_exspander(
    graph: &mut CdbgGraph,
    matrix: &PairLinkMatrix,
    _dist: &InsertSizeDist,
    params: ExspanderParams,
) {
    let k = graph.k;
    let n = graph.nodes.len() as u32;

    // 1. Identify long unique seed unitigs (length >= 500 bp, non-deleted, non-circular)
    let mut seeds: Vec<u32> = (0..n)
        .filter(|&id| {
            let node = &graph.nodes[id as usize];
            !node.is_deleted && !node.is_circular && node.seq_len >= 500
        })
        .collect();

    // Sort seeds by length descending
    seeds.sort_by_key(|&id| std::cmp::Reverse(graph.nodes[id as usize].seq_len));

    let mut used_in_path: HashSet<u32> = HashSet::new();

    for seed_id in seeds {
        if used_in_path.contains(&seed_id) || graph.nodes[seed_id as usize].is_deleted {
            continue;
        }

        used_in_path.insert(seed_id);
        let seed_fwd = DirUnitigId::new(seed_id, false);

        // 1. Extend forward from 3' end
        let fwd_ext = extend_one_direction(seed_fwd, graph, matrix, params, &mut used_in_path);

        // 2. Extend backward from 5' end (which is forward from seed.twin())
        let bwd_ext = extend_one_direction(seed_fwd.twin(), graph, matrix, params, &mut used_in_path);

        // Stitch path together:
        // backward extension reversed & twinned + seed + forward extension
        let mut full_path: Vec<DirUnitigId> = Vec::new();
        for &d in bwd_ext.iter().rev() {
            full_path.push(d.twin());
        }
        full_path.push(seed_fwd);
        full_path.extend_from_slice(&fwd_ext);

        // Merge full_path if length > 1
        if full_path.len() > 1 {
            let mut assembled_seq = Vec::new();
            for (idx, &dir) in full_path.iter().enumerate() {
                let u_seq = graph.get_stranded_seq(dir);
                if idx == 0 {
                    assembled_seq.extend_from_slice(&u_seq);
                } else {
                    let overlap = (k - 1).min(u_seq.len());
                    assembled_seq.extend_from_slice(&u_seq[overlap..]);
                }
                // Mark all non-circular nodes incorporated into the path as deleted
                if !graph.nodes[dir.id() as usize].is_circular {
                    graph.nodes[dir.id() as usize].is_deleted = true;
                }
            }

            let new_seq_id = graph.arena.push_seq(&assembled_seq);
            let new_node_id = graph.nodes.len() as u32;
            let new_len = assembled_seq.len() as u32;
            let node = UnitigNode::new(new_node_id, new_seq_id, new_len, 20.0);
            graph.nodes.push(node);
        }
    }

    // Clean up any residual small non-circular unitigs
    for node in graph.nodes.iter_mut() {
        if !node.is_circular && node.seq_len < 500 {
            node.is_deleted = true;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::bitseq::{reverse_complement_ascii, ReadBank};
    use crate::common::sketch::splitmix64;
    use crate::graph::build::build_cdbg;
    use crate::resolve::mapper::map_paired_reads;

    #[test]
    fn test_verification_gate_exspander_350bp_repeat_resolution() {
        // Verification Gate for Phase 5:
        // Simulate a 30 kb chromosome containing two copies of a 350 bp IS-like repeat
        // sequenced with 2 x 150 bp reads at 450 bp insert size:
        // Verify exspander accurately traverses across both copies of the 350 bp repeat
        // without misassembly, producing a single 30 kb contig.

        let k = 31;
        let alphabet = [b'A', b'C', b'G', b'T'];
        let gen_seq = |seed: u64, len: usize| -> Vec<u8> {
            let mut buf = Vec::with_capacity(len);
            for i in 0..len {
                let r = splitmix64(seed + i as u64);
                buf.push(alphabet[(r % 4) as usize]);
            }
            buf
        };

        // 350 bp IS-like repeat
        let rep350 = gen_seq(0x1111_3500, 350);
        let u1 = gen_seq(0xAAAA_1000, 9400);
        let u2 = gen_seq(0xBBBB_2000, 9700);
        let u3 = gen_seq(0xCCCC_3000, 10200);

        // 30 kb chromosome: U1 -> Rep -> U2 -> Rep -> U3
        let mut chrom = Vec::new();
        chrom.extend_from_slice(&u1);
        chrom.extend_from_slice(&rep350);
        chrom.extend_from_slice(&u2);
        chrom.extend_from_slice(&rep350);
        chrom.extend_from_slice(&u3);

        let total_chrom_len = chrom.len(); // 9400 + 350 + 9700 + 350 + 10200 = 30,000 bp
        assert_eq!(total_chrom_len, 30_000);

        // Simulate 40x depth of 2 x 150bp paired-end reads with insert size ~450 bp
        let mut bank = ReadBank::new();
        let num_pairs = (total_chrom_len * 40) / 300;

        let insert_size = 450usize;
        let read_len = 150usize;

        let mut rng: u64 = 0x5555_6666_7777_8888;
        for i in 0..num_pairs {
            rng = splitmix64(rng + i as u64);
            if total_chrom_len <= insert_size {
                break;
            }
            let start = (rng as usize) % (total_chrom_len - insert_size);

            let r1 = &chrom[start..start + read_len];
            let r2_frag = &chrom[start + insert_size - read_len..start + insert_size];
            let r2 = reverse_complement_ascii(r2_frag);

            bank.pe_forward.push_seq(r1);
            bank.pe_reverse.push_seq(&r2);
        }

        // 1. Build cDBG graph at k=31
        let mut graph = build_cdbg::<u64>(&bank, k, Some(3));

        // 2. Map paired reads to construct link matrix
        let (matrix, dist) = map_paired_reads::<u64>(&graph, &bank, k);

        // 3. Run exSPAnder
        let params = ExspanderParams {
            min_support: 5,
            decision_threshold: 2.5,
        };
        resolve_repeats_exspander(&mut graph, &matrix, &dist, params);

        // Verify: The 30 kb chromosome assembles into a single continuous contig
        let long_contigs: Vec<&UnitigNode> = graph
            .nodes
            .iter()
            .filter(|n| !n.is_deleted && n.seq_len >= 29_000)
            .collect();

        assert_eq!(
            long_contigs.len(),
            1,
            "exSPAnder must resolve both 350 bp repeats to produce a single ~30 kb contig! Found {} contigs >= 29kb",
            long_contigs.len()
        );
        assert!(
            long_contigs[0].seq_len >= 29_500,
            "Contig length should be ~30,000 bp, found {}",
            long_contigs[0].seq_len
        );
    }
}
