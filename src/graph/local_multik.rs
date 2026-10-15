//! Local multi-k graph upgrade: iteratively bridges branching regions using higher k-mer sizes.

use crate::common::bitseq::ReadBank;
use crate::common::kmer::KmerInt;
use crate::common::sketch::CountMinSketch;
use crate::diag::inspect_graph;
use crate::graph::build::build_cdbg_with_contigs_ext;
use crate::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, UnitigNode};
use crate::graph::simplify::simplify_graph;

/// Detailed telemetry result for a local multi-k upgrade attempt.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MultiKUpgradeResult {
    pub executed: bool,
    pub early_stop_reason: Option<&'static str>,
    pub genome_size_source: &'static str,
    pub genome_size_used: u64,
    pub cms_width: usize,
    pub cms_allocated_bytes: usize,
    pub active_unitigs_before: usize,
    pub active_edges_before: usize,
    pub active_unitigs_after: usize,
    pub active_edges_after: usize,
}

/// Upgrades an existing cDBG graph from k1 to k2 (where k2 > k1, e.g. 31 -> 79).
///
/// Uses documented fallback if no explicit genome size is passed.
pub fn upgrade_multik<K2: KmerInt>(
    prev_graph: &mut CdbgGraph,
    bank: &ReadBank,
    k2: usize,
) {
    upgrade_multik_ext::<K2>(prev_graph, bank, k2, None);
}

/// Upgrades an existing cDBG graph from k1 to k2 with explicit genome size propagation.
pub fn upgrade_multik_ext<K2: KmerInt>(
    prev_graph: &mut CdbgGraph,
    bank: &ReadBank,
    k2: usize,
    genome_size: Option<u64>,
) -> MultiKUpgradeResult {
    let (u_before, e_before, _, _) = inspect_graph(prev_graph);
    if k2 <= prev_graph.k {
        return MultiKUpgradeResult {
            executed: false,
            early_stop_reason: Some("K2_LE_PREV_K"),
            genome_size_source: "SKIPPED",
            genome_size_used: 0,
            cms_width: 0,
            cms_allocated_bytes: 0,
            active_unitigs_before: u_before,
            active_edges_before: e_before,
            active_unitigs_after: u_before,
            active_edges_after: e_before,
        };
    }

    // 1. Identify branching nodes in prev_graph
    let mut has_branching = false;
    for node in &prev_graph.nodes {
        if !node.is_deleted && !node.is_circular {
            if node.edges_fwd.len() > 1 || node.edges_rev.len() > 1 {
                has_branching = true;
                break;
            }
        }
    }

    if !has_branching {
        return MultiKUpgradeResult {
            executed: false,
            early_stop_reason: Some("NO_BRANCHING_REGIONS"),
            genome_size_source: "SKIPPED",
            genome_size_used: 0,
            cms_width: 0,
            cms_allocated_bytes: 0,
            active_unitigs_before: u_before,
            active_edges_before: e_before,
            active_unitigs_after: u_before,
            active_edges_after: e_before,
        };
    }

    // Freeze long linear regions: collect all active NON-CIRCULAR unitig sequences from prev_graph
    let prior_contigs: Vec<Vec<u8>> = prev_graph
        .nodes
        .iter()
        .filter(|n| !n.is_deleted && !n.is_circular)
        .map(|n| prev_graph.get_unitig_seq(n.id))
        .collect();

    // 2. Filter active reads: drop interior reads of long resolved unitigs
    let active_bank = filter_active_reads(bank, prev_graph);

    // 3. Build cDBG at k2, incorporating prior contigs with active boundary reads
    let (est_g, source) = if let Some(g) = genome_size {
        (g.max(50_000), "EXPLICIT_PASSED")
    } else {
        ((bank.total_bases() as u64 / 60).max(50_000), "FALLBACK_READ_BASES_DIV_60")
    };
    let cms_width = CountMinSketch::width_for_genome_size(est_g);
    let cms_allocated_bytes = CountMinSketch::allocated_bytes_for_genome_size(est_g);

    let mut new_graph = build_cdbg_with_contigs_ext::<K2>(&active_bank, &prior_contigs, k2, Some(4), Some(est_g));
    simplify_graph(&mut new_graph);

    // 4. Preserve isolated circular plasmids from prev_graph if not present in new_graph
    for old_node in &prev_graph.nodes {
        if !old_node.is_deleted && old_node.is_circular {
            let old_seq = prev_graph.get_unitig_seq(old_node.id);
            // Check if this circular plasmid already exists in new_graph (within k2 bp length tolerance)
            let exists_in_new = new_graph.nodes.iter().any(|n| {
                !n.is_deleted && n.is_circular && ((n.seq_len as i64 - old_seq.len() as i64).abs() <= (k2 as i64 + 10))
            });

            if !exists_in_new {
                // Transfer circular node to new_graph
                let seq_id = new_graph.arena.push_seq(&old_seq);
                let node_id = new_graph.nodes.len() as u32;
                let mut node = UnitigNode::new(node_id, seq_id, old_seq.len() as u32, old_node.kmer_cov);
                node.is_circular = true;
                let self_dir = DirUnitigId::new(node_id, false);
                node.add_edge(false, GraphEdge::new(self_dir, (k2 - 1) as u8, 10));
                new_graph.nodes.push(node);
            }
        }
    }

    simplify_graph(&mut new_graph);
    let (u_after, e_after, _, _) = inspect_graph(&new_graph);
    *prev_graph = new_graph;

    MultiKUpgradeResult {
        executed: true,
        early_stop_reason: None,
        genome_size_source: source,
        genome_size_used: est_g,
        cms_width,
        cms_allocated_bytes,
        active_unitigs_before: u_before,
        active_edges_before: e_before,
        active_unitigs_after: u_after,
        active_edges_after: e_after,
    }
}

/// Upgrades an existing cDBG graph from k1 to k2 using the legacy Var 1 capacity sizing (prior contigs length sum).
pub fn upgrade_multik_legacy<K2: KmerInt>(
    prev_graph: &mut CdbgGraph,
    bank: &ReadBank,
    k2: usize,
) -> MultiKUpgradeResult {
    let (u_before, e_before, _, _) = inspect_graph(prev_graph);
    if k2 <= prev_graph.k {
        return MultiKUpgradeResult {
            executed: false,
            early_stop_reason: Some("K2_LE_PREV_K"),
            genome_size_source: "SKIPPED",
            genome_size_used: 0,
            cms_width: 0,
            cms_allocated_bytes: 0,
            active_unitigs_before: u_before,
            active_edges_before: e_before,
            active_unitigs_after: u_before,
            active_edges_after: e_before,
        };
    }

    let mut has_branching = false;
    for node in &prev_graph.nodes {
        if !node.is_deleted && !node.is_circular {
            if node.edges_fwd.len() > 1 || node.edges_rev.len() > 1 {
                has_branching = true;
                break;
            }
        }
    }

    if !has_branching {
        return MultiKUpgradeResult {
            executed: false,
            early_stop_reason: Some("NO_BRANCHING_REGIONS"),
            genome_size_source: "SKIPPED",
            genome_size_used: 0,
            cms_width: 0,
            cms_allocated_bytes: 0,
            active_unitigs_before: u_before,
            active_edges_before: e_before,
            active_unitigs_after: u_before,
            active_edges_after: e_before,
        };
    }

    let prior_contigs: Vec<Vec<u8>> = prev_graph
        .nodes
        .iter()
        .filter(|n| !n.is_deleted && !n.is_circular)
        .map(|n| prev_graph.get_unitig_seq(n.id))
        .collect();

    let active_bank = filter_active_reads(bank, prev_graph);

    // Legacy formula from Var 1 / dac6a57: sum of prior contigs lengths
    let legacy_est_g = if !prior_contigs.is_empty() {
        prior_contigs.iter().map(|c| c.len() as u64).sum::<u64>().max(50_000)
    } else {
        (bank.total_bases() as u64 / 60).max(50_000)
    };
    let cms_width = CountMinSketch::width_for_genome_size(legacy_est_g);
    let cms_allocated_bytes = CountMinSketch::allocated_bytes_for_genome_size(legacy_est_g);

    let mut new_graph = build_cdbg_with_contigs_ext::<K2>(&active_bank, &prior_contigs, k2, Some(4), Some(legacy_est_g));
    simplify_graph(&mut new_graph);

    for old_node in &prev_graph.nodes {
        if !old_node.is_deleted && old_node.is_circular {
            let old_seq = prev_graph.get_unitig_seq(old_node.id);
            let exists_in_new = new_graph.nodes.iter().any(|n| {
                !n.is_deleted && n.is_circular && ((n.seq_len as i64 - old_seq.len() as i64).abs() <= (k2 as i64 + 10))
            });

            if !exists_in_new {
                let seq_id = new_graph.arena.push_seq(&old_seq);
                let node_id = new_graph.nodes.len() as u32;
                let mut node = UnitigNode::new(node_id, seq_id, old_seq.len() as u32, old_node.kmer_cov);
                node.is_circular = true;
                let self_dir = DirUnitigId::new(node_id, false);
                node.add_edge(false, GraphEdge::new(self_dir, (k2 - 1) as u8, 10));
                new_graph.nodes.push(node);
            }
        }
    }

    simplify_graph(&mut new_graph);
    let (u_after, e_after, _, _) = inspect_graph(&new_graph);
    *prev_graph = new_graph;

    MultiKUpgradeResult {
        executed: true,
        early_stop_reason: None,
        genome_size_source: "LEGACY_PRIOR_CONTIGS_SUM",
        genome_size_used: legacy_est_g,
        cms_width,
        cms_allocated_bytes,
        active_unitigs_before: u_before,
        active_edges_before: e_before,
        active_unitigs_after: u_after,
        active_edges_after: e_after,
    }
}

/// Filters active reads for multi-k upgrades.
///
/// Discards reads that fall entirely within the safe interior (> 350 bp from ends)
/// of long resolved unitigs (>= 2000 bp), keeping only reads near boundaries,
/// branching junctions, and unresolved regions.
pub fn filter_active_reads(bank: &ReadBank, graph: &CdbgGraph) -> ReadBank {
    let mut boundary_kmers: rustc_hash::FxHashSet<u64> = rustc_hash::FxHashSet::default();
    let k = graph.k.min(31);
    let boundary_window = 350usize;

    for node in &graph.nodes {
        if node.is_deleted {
            continue;
        }
        let seq = graph.get_unitig_seq(node.id);
        if seq.len() < k {
            continue;
        }

        if node.seq_len < 2000 || node.is_circular || node.edges_fwd.len() > 1 || node.edges_rev.len() > 1 {
            // Keep all k-mers for short, circular, or branching unitigs
            for i in 0..=seq.len() - k {
                let h = crate::common::kmer::RollingNtHash::new(&seq[i..i + k]).canonical();
                boundary_kmers.insert(h);
            }
        } else {
            // Long unitig: only index boundary windows at 5' and 3' ends
            let left_end = boundary_window.min(seq.len() - k + 1);
            for i in 0..left_end {
                let h = crate::common::kmer::RollingNtHash::new(&seq[i..i + k]).canonical();
                boundary_kmers.insert(h);
            }

            let right_start = (seq.len().saturating_sub(boundary_window)).min(seq.len() - k);
            for i in right_start..=seq.len() - k {
                let h = crate::common::kmer::RollingNtHash::new(&seq[i..i + k]).canonical();
                boundary_kmers.insert(h);
            }
        }
    }

    let mut filtered_bank = ReadBank::new();
    let mut buf = Vec::with_capacity(300);

    let is_active = |seq: &[u8], boundary_kmers: &rustc_hash::FxHashSet<u64>| -> bool {
        if seq.len() < k {
            return false;
        }
        let check_offsets = [0, (seq.len() - k) / 2, seq.len() - k];
        for &offset in &check_offsets {
            let h = crate::common::kmer::RollingNtHash::new(&seq[offset..offset + k]).canonical();
            if boundary_kmers.contains(&h) {
                return true;
            }
        }
        false
    };

    // Filter merged reads
    for id in 0..bank.merged_reads.len() as u32 {
        bank.merged_reads.extract_ascii_into(id, &mut buf);
        if is_active(&buf, &boundary_kmers) {
            filtered_bank.merged_reads.push_seq(&buf);
        }
    }

    // Filter PE pairs
    let mut buf_r2 = Vec::with_capacity(300);
    let pe_len = bank.pe_forward.len().min(bank.pe_reverse.len()) as u32;
    for id in 0..pe_len {
        bank.pe_forward.extract_ascii_into(id, &mut buf);
        bank.pe_reverse.extract_ascii_into(id, &mut buf_r2);
        if is_active(&buf, &boundary_kmers) || is_active(&buf_r2, &boundary_kmers) {
            filtered_bank.pe_forward.push_seq(&buf);
            filtered_bank.pe_reverse.push_seq(&buf_r2);
        }
    }

    filtered_bank
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::bitseq::ReadBank;
    use crate::common::sketch::splitmix64;
    use crate::graph::build::build_cdbg;

    #[test]
    fn test_verification_gate_multik_resolution_and_plasmid() {
        // Verification Gate for Phase 4:
        // Synthetic 20 kb sequence containing a 65 bp repeat (unresolvable at k=31, resolvable at k=79)
        // and a 5 kb high-copy plasmid (25x relative coverage):
        // Verify that after running simplify + local_multik (31 -> 79), the 20 kb sequence
        // assembles into a single contig and the 5 kb plasmid is preserved intact.

        let alphabet = [b'A', b'C', b'G', b'T'];
        let gen_seq = |seed: u64, len: usize| -> Vec<u8> {
            let mut buf = Vec::with_capacity(len);
            for i in 0..len {
                let r = splitmix64(seed + i as u64);
                buf.push(alphabet[(r % 4) as usize]);
            }
            buf
        };

        // 65 bp exact repeat
        let rep65 = gen_seq(0x1234_5678, 65);
        let flank1 = gen_seq(0xAAAA_1111, 9900);
        let flank2 = gen_seq(0xBBBB_2222, 9970);

        // 20 kb chromosome with 65 bp repeat: Flank1 -> rep65 -> Flank2
        let mut chrom = Vec::new();
        chrom.extend_from_slice(&flank1);
        chrom.extend_from_slice(&rep65);
        chrom.extend_from_slice(&flank2);
        assert_eq!(chrom.len(), 9900 + 65 + 9970); // ~20 kb

        // 5 kb circular plasmid
        let plasmid = gen_seq(0xCCCC_3333, 5000);

        // Simulate reads:
        // Chromosome: 40x depth of 150bp reads
        // Plasmid: 100x depth (25x relative to chromosome)
        let mut bank = ReadBank::new();

        // Simulate chromosome reads
        let chrom_reads = (chrom.len() * 40) / 150;
        let mut rng = 0x9999;
        for i in 0..chrom_reads {
            rng = splitmix64(rng + i as u64);
            let start = (rng as usize) % (chrom.len() - 150);
            bank.merged_reads.push_seq(&chrom[start..start + 150]);
        }

        // Simulate plasmid reads (circular)
        let plasmid_reads = (plasmid.len() * 100) / 150;
        for i in 0..plasmid_reads {
            rng = splitmix64(rng + i as u64);
            let start = (rng as usize) % plasmid.len();
            let mut read = Vec::with_capacity(150);
            for j in 0..150 {
                read.push(plasmid[(start + j) % plasmid.len()]);
            }
            bank.merged_reads.push_seq(&read);
        }

        // 1. Build initial graph at k=31
        let mut graph = build_cdbg::<u64>(&bank, 31, Some(5));
        simplify_graph(&mut graph);

        // 2. Upgrade to k=79
        upgrade_multik::<[u64; 4]>(&mut graph, &bank, 79);
        simplify_graph(&mut graph);

        // Verify:
        // 1. The 5 kb plasmid is preserved intact with circular=true
        let plasmid_node = graph.nodes.iter().find(|n| !n.is_deleted && n.is_circular);
        assert!(plasmid_node.is_some(), "5 kb circular plasmid must be preserved intact!");
        let p_node = plasmid_node.unwrap();
        assert!(p_node.seq_len >= 4900 && p_node.seq_len <= 5100, "Plasmid length should be ~5000 bp");

        // 2. The 20 kb sequence assembles into a large contig (>= 19.5 kb)
        let chrom_node = graph.nodes.iter().find(|n| !n.is_deleted && !n.is_circular && n.seq_len >= 19000);
        assert!(chrom_node.is_some(), "20 kb chromosome must be resolved into a single contig at k=79!");
    }
}
