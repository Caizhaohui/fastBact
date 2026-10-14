//! Post-assembly curation and FASTA/GFA writer.
//!
//! - Trims terminal self-overlaps on circular contigs and formats FASTA headers.
//! - Filters contigs by minimum length and coverage, and discards homopolymer artifacts.
//! - Writes `contigs.fa` and GFA v1.0 `assembly_graph.gfa`.

use crate::graph::cdbg::CdbgGraph;
use anyhow::{Context, Result};
use std::fs::{create_dir_all, File};
use std::io::{BufWriter, Write};
use std::path::Path;

/// Final assembled and curated contig record.
#[derive(Clone, Debug, PartialEq)]
pub struct AssembledContig {
    pub id: usize,
    pub node_id: u32,
    pub name: String,
    pub seq: Vec<u8>,
    pub len: usize,
    pub cov: f32,
    pub is_circular: bool,
}

/// Curates contigs from the cDBG graph and returns filtered, circularity-trimmed contigs.
pub fn curate_contigs(
    graph: &CdbgGraph,
    min_len: usize,
    min_cov: f32,
) -> Vec<AssembledContig> {
    let k = graph.k;
    let mut contigs = Vec::new();

    for id in 0..graph.nodes.len() as u32 {
        let node = &graph.nodes[id as usize];
        if node.is_deleted {
            continue;
        }

        let mut seq = graph.get_unitig_seq(id);

        // 1. Terminal overlap trimming for circular contigs
        let is_circular = node.is_circular;
        if is_circular && seq.len() > k {
            if seq[seq.len() - k..] == seq[..k] {
                seq.truncate(seq.len() - k);
            } else {
                let ov = k - 1;
                if seq.len() > ov && seq[seq.len() - ov..] == seq[..ov] {
                    seq.truncate(seq.len() - ov);
                }
            }
        }

        let len = seq.len();

        // 2. Length filter
        if len < min_len {
            continue;
        }

        // 3. Coverage filter
        if node.kmer_cov < min_cov {
            continue;
        }

        // 4. Homopolymer filter: discard contigs consisting of a single repeated nucleotide
        if len > 1 && seq.iter().all(|&b| b == seq[0]) {
            continue;
        }

        contigs.push(AssembledContig {
            id: contigs.len() + 1,
            node_id: id,
            name: format!("contig_{}", contigs.len() + 1),
            seq,
            len,
            cov: node.kmer_cov,
            is_circular,
        });
    }

    // Sort contigs by length descending
    contigs.sort_by_key(|c| std::cmp::Reverse(c.len));

    // Deduplicate: remove contigs that are exact substrings of longer contigs or duplicate circular plasmids
    let mut kept_contigs: Vec<AssembledContig> = Vec::new();
    for cand in contigs {
        let cand_rc = crate::common::bitseq::reverse_complement_ascii(&cand.seq);
        let mut is_sub = false;

        for kept in &kept_contigs {
            // Check if cand is a circular duplicate
            if cand.is_circular && kept.is_circular && (kept.len as i64 - cand.len as i64).abs() <= 120 {
                is_sub = true;
                break;
            }

            if kept.len >= cand.len {
                let haystack = if kept.is_circular {
                    let mut d = kept.seq.clone();
                    d.extend_from_slice(&kept.seq);
                    d
                } else {
                    kept.seq.clone()
                };

                let probe_len = 31.min(cand.len);
                let probe = &cand.seq[..probe_len];
                let probe_rc = &cand_rc[..probe_len];

                let probe_match = haystack.windows(probe_len).any(|w| w == probe || w == probe_rc);
                if probe_match
                    && (haystack.windows(cand.len).any(|w| w == cand.seq)
                        || haystack.windows(cand.len).any(|w| w == cand_rc))
                {
                    is_sub = true;
                    break;
                }
            }
        }

        if !is_sub {
            kept_contigs.push(cand);
        }
    }

    // Re-index contig names by length rank
    for (idx, c) in kept_contigs.iter_mut().enumerate() {
        c.id = idx + 1;
        c.name = format!("contig_{}", idx + 1);
    }

    kept_contigs
}

/// Writes final contigs to `<outdir>/contigs.fa` and `<outdir>/assembly_graph.gfa`.
pub fn write_outputs<P: AsRef<Path>>(
    graph: &CdbgGraph,
    contigs: &[AssembledContig],
    outdir: P,
) -> Result<()> {
    create_dir_all(outdir.as_ref())
        .with_context(|| format!("Failed to create output directory {:?}", outdir.as_ref()))?;

    // 1. Write contigs.fa
    let fa_path = outdir.as_ref().join("contigs.fa");
    let fa_file = File::create(&fa_path)
        .with_context(|| format!("Failed to create FASTA output {:?}", fa_path))?;
    let mut fa_writer = BufWriter::new(fa_file);

    for c in contigs {
        writeln!(
            fa_writer,
            ">{} len={} cov={:.1} circular={}",
            c.name, c.len, c.cov, c.is_circular
        )?;
        // Write sequence wrapped at 80 bp
        for chunk in c.seq.chunks(80) {
            fa_writer.write_all(chunk)?;
            fa_writer.write_all(b"\n")?;
        }
    }
    fa_writer.flush()?;

    // 2. Write assembly_graph.gfa
    let gfa_path = outdir.as_ref().join("assembly_graph.gfa");
    let gfa_file = File::create(&gfa_path)
        .with_context(|| format!("Failed to create GFA output {:?}", gfa_path))?;
    let mut gfa_writer = BufWriter::new(gfa_file);

    writeln!(gfa_writer, "H\tVN:Z:1.0")?;
    for c in contigs {
        writeln!(
            gfa_writer,
            "S\t{}\t{}\tLN:i:{}\tKC:i:{}",
            c.name,
            String::from_utf8_lossy(&c.seq),
            c.len,
            c.cov as u32
        )?;
    }

    // Map active contig nodes to their contig names
    let mut node_to_name = hashbrown::HashMap::new();
    for c in contigs {
        node_to_name.insert(c.node_id, c.name.clone());
    }

    // Write links from graph
    for id in 0..graph.nodes.len() as u32 {
        let node = &graph.nodes[id as usize];
        if node.is_deleted {
            continue;
        }
        let Some(source_name) = node_to_name.get(&id) else {
            continue;
        };

        for edge in &node.edges_fwd {
            let target_id = edge.target.id();
            if (target_id as usize) < graph.nodes.len() && !graph.nodes[target_id as usize].is_deleted {
                if let Some(target_name) = node_to_name.get(&target_id) {
                    let target_strand = if edge.target.is_rc() { "-" } else { "+" };
                    writeln!(
                        gfa_writer,
                        "L\t{}\t+\t{}\t{}\t{}M",
                        source_name,
                        target_name,
                        target_strand,
                        edge.overlap_k
                    )?;
                }
            }
        }
    }
    gfa_writer.flush()?;

    Ok(())
}
