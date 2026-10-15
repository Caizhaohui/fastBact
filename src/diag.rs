//! Diagnostic instrumentation, intermediate stage dumping, and telemetry recording for fastBact.

use crate::common::bitseq::ReadBank;
use crate::common::sketch::CountMinSketch;
use crate::graph::cdbg::CdbgGraph;
use crate::polish::output::AssembledContig;
use rustc_hash::FxHashMap;
use std::fs::{create_dir_all, File};
use std::io::{BufWriter, Write};
use std::path::Path;

/// Stage-level graph statistics.
#[derive(Clone, Debug, Default)]
pub struct StageGraphTelemetry {
    pub stage_name: String,
    pub k: usize,
    pub active_unitigs: usize,
    pub active_edges: usize,
    pub isolated_cycles: usize,
    pub deleted_tips_count: usize,
    pub deleted_tips_bp: usize,
    pub popped_bubbles_count: usize,
    pub popped_bubbles_bp: usize,
    pub total_graph_bp: usize,
}

/// Ingest statistics.
#[derive(Clone, Debug, Default)]
pub struct IngestTelemetry {
    pub input_read_pairs: usize,
    pub retained_read_pairs: usize,
    pub bases_before_trim: usize,
    pub bases_after_trim: usize,
    pub merged_pairs_count: usize,
    pub merged_bases: usize,
    pub unmerged_pairs_count: usize,
    pub effective_coverage: f64,
    pub subsampling_enabled: bool,
    pub subsampling_reason: String,
    pub estimated_genome_size: u64,
    pub estimated_genome_size_note: String,
    pub estimated_coverage: f64,
    pub estimator_confidence: String,
    pub actual_input_pairs: usize,
    pub actual_retained_pairs: usize,
}

/// K-mer histogram and Count-Min Sketch agreement metrics.
#[derive(Clone, Debug, Default)]
pub struct KmerTelemetry {
    pub freq_counts: Vec<(u8, usize)>, // counts at freq 1..=20
    pub estimated_t_solid: u8,
    pub total_unique_kmers: usize,
    pub solid_kmers_retained: usize,
    pub solid_retention_rate: f64,
    pub cms_sampled_kmers: usize,
    pub cms_exact_match_count: usize,
    pub cms_overestimated_count: usize,
    pub cms_overestimation_rate: f64,
    pub cms_mean_error: f64,
}

/// Dumps all active unitigs in `graph` to a FASTA file.
pub fn dump_graph_fasta<P: AsRef<Path>>(graph: &CdbgGraph, path: P) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    let mut count = 0;
    for id in 0..graph.nodes.len() as u32 {
        let node = &graph.nodes[id as usize];
        if node.is_deleted {
            continue;
        }
        count += 1;
        let seq = graph.get_unitig_seq(id);
        writeln!(
            writer,
            ">unitig_{} len={} cov={:.2} circular={}",
            count,
            seq.len(),
            node.kmer_cov,
            node.is_circular
        )?;
        writer.write_all(&seq)?;
        writeln!(writer)?;
    }
    writer.flush()?;
    Ok(())
}

/// Dumps curated contigs to a FASTA file.
pub fn dump_contigs_fasta<P: AsRef<Path>>(contigs: &[AssembledContig], path: P) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    for c in contigs {
        writeln!(
            writer,
            ">{} len={} cov={:.2} circular={}",
            c.name, c.len, c.cov, c.is_circular
        )?;
        writer.write_all(&c.seq)?;
        writeln!(writer)?;
    }
    writer.flush()?;
    Ok(())
}

/// Inspects graph and returns (active_unitigs, active_edges, isolated_cycles, total_graph_bp).
pub fn inspect_graph(graph: &CdbgGraph) -> (usize, usize, usize, usize) {
    let mut active_unitigs = 0;
    let mut active_edges = 0;
    let mut isolated_cycles = 0;
    let mut total_graph_bp = 0;

    for id in 0..graph.nodes.len() as u32 {
        let node = &graph.nodes[id as usize];
        if node.is_deleted {
            continue;
        }
        active_unitigs += 1;
        total_graph_bp += node.seq_len as usize;
        let fwd = node.edges_fwd.len();
        let rev = node.edges_rev.len();
        active_edges += fwd + rev;
        if node.is_circular && fwd == 0 && rev == 0 {
            isolated_cycles += 1;
        }
    }
    (active_unitigs, active_edges, isolated_cycles, total_graph_bp)
}

/// Computes k-mer frequency histogram and CMS agreement on a sample of reads from `ReadBank`.
pub fn sample_kmer_telemetry(bank: &ReadBank, k: usize, genome_size: u64) -> KmerTelemetry {
    let cms = CountMinSketch::for_genome_size(genome_size.max(50_000));
    let mut sample_counts: FxHashMap<u64, u16> = FxHashMap::default();
    let sample_limit = 50_000.min(bank.merged_reads.len());

    let mut buf = Vec::with_capacity(300);
    for id in 0..sample_limit as u32 {
        if bank.merged_reads.get_len(id) as usize >= k {
            bank.merged_reads.extract_ascii_into(id, &mut buf);
            for w in buf.windows(k) {
                let mut h: u64 = 0;
                for &b in w {
                    h = h.wrapping_mul(31).wrapping_add(b as u64);
                }
                *sample_counts.entry(h).or_insert(0) += 1;
                cms.add(h);
            }
        }
    }

    if sample_counts.is_empty() {
        let pe_limit = 50_000.min(bank.pe_forward.len());
        for id in 0..pe_limit as u32 {
            if bank.pe_forward.get_len(id) as usize >= k {
                bank.pe_forward.extract_ascii_into(id, &mut buf);
                for w in buf.windows(k) {
                    let mut h: u64 = 0;
                    for &b in w {
                        h = h.wrapping_mul(31).wrapping_add(b as u64);
                    }
                    *sample_counts.entry(h).or_insert(0) += 1;
                    cms.add(h);
                }
            }
        }
    }

    let mut freq_histogram = vec![0usize; 21];
    let total_unique = sample_counts.len();
    let mut solid_retained = 0usize;

    let mut sampled_for_cms = 0usize;
    let mut cms_exact_match = 0usize;
    let mut cms_overestimated = 0usize;
    let mut total_error: f64 = 0.0;

    for (&h, &exact_cov) in &sample_counts {
        let freq_idx = (exact_cov as usize).min(20);
        freq_histogram[freq_idx] += 1;
        if exact_cov >= 2 {
            solid_retained += 1;
        }

        if sampled_for_cms < 10_000 {
            sampled_for_cms += 1;
            let cms_cov = cms.estimate(h) as u16;
            if cms_cov == exact_cov {
                cms_exact_match += 1;
            } else if cms_cov > exact_cov {
                cms_overestimated += 1;
                total_error += (cms_cov - exact_cov) as f64;
            }
        }
    }

    let mut freq_pairs = Vec::new();
    for f in 1..=20 {
        freq_pairs.push((f as u8, freq_histogram[f]));
    }

    let retention_rate = if total_unique > 0 {
        solid_retained as f64 / total_unique as f64
    } else {
        0.0
    };

    let overestimation_rate = if sampled_for_cms > 0 {
        cms_overestimated as f64 / sampled_for_cms as f64
    } else {
        0.0
    };

    let mean_error = if sampled_for_cms > 0 {
        total_error / sampled_for_cms as f64
    } else {
        0.0
    };

    KmerTelemetry {
        freq_counts: freq_pairs,
        estimated_t_solid: 2,
        total_unique_kmers: total_unique,
        solid_kmers_retained: solid_retained,
        solid_retention_rate: retention_rate,
        cms_sampled_kmers: sampled_for_cms,
        cms_exact_match_count: cms_exact_match,
        cms_overestimated_count: cms_overestimated,
        cms_overestimation_rate: overestimation_rate,
        cms_mean_error: mean_error,
    }
}

/// Writes `stage_telemetry.json` in standard JSON format.
pub fn write_stage_telemetry_json<P: AsRef<Path>>(
    path: P,
    ingest: &IngestTelemetry,
    kmer: &KmerTelemetry,
    stages: &[StageGraphTelemetry],
) -> std::io::Result<()> {
    if let Some(parent) = path.as_ref().parent() {
        create_dir_all(parent)?;
    }
    let mut writer = BufWriter::new(File::create(path)?);

    writeln!(writer, "{{")?;
    writeln!(writer, "  \"input_read_pairs\": {},", ingest.input_read_pairs)?;
    writeln!(writer, "  \"retained_read_pairs\": {},", ingest.retained_read_pairs)?;
    writeln!(writer, "  \"bases_before_trim\": {},", ingest.bases_before_trim)?;
    writeln!(writer, "  \"bases_after_trim\": {},", ingest.bases_after_trim)?;
    writeln!(writer, "  \"merged_pairs_count\": {},", ingest.merged_pairs_count)?;
    writeln!(writer, "  \"merged_bases\": {},", ingest.merged_bases)?;
    writeln!(writer, "  \"unmerged_pairs_count\": {},", ingest.unmerged_pairs_count)?;
    writeln!(writer, "  \"effective_coverage\": {:.2},", ingest.effective_coverage)?;
    writeln!(writer, "  \"subsampling_enabled\": {},", ingest.subsampling_enabled)?;
    writeln!(writer, "  \"subsampling_reason\": \"{}\",", ingest.subsampling_reason)?;
    writeln!(writer, "  \"estimated_genome_size\": {},", ingest.estimated_genome_size)?;
    writeln!(writer, "  \"estimated_genome_size_note\": \"{}\",", ingest.estimated_genome_size_note)?;
    writeln!(writer, "  \"estimated_coverage\": {:.2},", ingest.estimated_coverage)?;
    writeln!(writer, "  \"estimator_confidence\": \"{}\",", ingest.estimator_confidence)?;
    writeln!(writer, "  \"actual_input_pairs\": {},", ingest.actual_input_pairs)?;
    writeln!(writer, "  \"actual_retained_pairs\": {},", ingest.actual_retained_pairs)?;

    // kmer_histogram_summary
    writeln!(writer, "  \"kmer_histogram_summary\": {{")?;
    for (i, &(freq, count)) in kmer.freq_counts.iter().enumerate() {
        let comma = if i + 1 == kmer.freq_counts.len() { "" } else { "," };
        writeln!(writer, "    \"freq_{}\": {}{}", freq, count, comma)?;
    }
    writeln!(writer, "  }},")?;
    writeln!(writer, "  \"estimated_t_solid\": {},", kmer.estimated_t_solid)?;
    writeln!(writer, "  \"total_unique_kmers\": {},", kmer.total_unique_kmers)?;
    writeln!(writer, "  \"solid_kmers_retained\": {},", kmer.solid_kmers_retained)?;
    writeln!(writer, "  \"solid_retention_rate\": {:.4},", kmer.solid_retention_rate)?;

    // cms_vs_exact_agreement
    writeln!(writer, "  \"cms_vs_exact_agreement\": {{")?;
    writeln!(writer, "    \"sampled_kmers\": {},", kmer.cms_sampled_kmers)?;
    writeln!(writer, "    \"exact_match_count\": {},", kmer.cms_exact_match_count)?;
    writeln!(writer, "    \"overestimated_count\": {},", kmer.cms_overestimated_count)?;
    writeln!(writer, "    \"overestimation_rate\": {:.4},", kmer.cms_overestimation_rate)?;
    writeln!(writer, "    \"mean_error\": {:.4}", kmer.cms_mean_error)?;
    writeln!(writer, "  }},")?;

    // stages
    writeln!(writer, "  \"stages\": [")?;
    for (idx, st) in stages.iter().enumerate() {
        let comma = if idx + 1 == stages.len() { "" } else { "," };
        writeln!(writer, "    {{")?;
        writeln!(writer, "      \"stage_name\": \"{}\",", st.stage_name)?;
        writeln!(writer, "      \"k\": {},", st.k)?;
        writeln!(writer, "      \"active_unitigs\": {},", st.active_unitigs)?;
        writeln!(writer, "      \"active_edges\": {},", st.active_edges)?;
        writeln!(writer, "      \"isolated_cycles\": {},", st.isolated_cycles)?;
        writeln!(writer, "      \"deleted_tips_count\": {},", st.deleted_tips_count)?;
        writeln!(writer, "      \"deleted_tips_bp\": {},", st.deleted_tips_bp)?;
        writeln!(writer, "      \"popped_bubbles_count\": {},", st.popped_bubbles_count)?;
        writeln!(writer, "      \"popped_bubbles_bp\": {},", st.popped_bubbles_bp)?;
        writeln!(writer, "      \"total_graph_bp\": {}", st.total_graph_bp)?;
        writeln!(writer, "    }}{}", comma)?;
    }
    writeln!(writer, "  ]")?;
    writeln!(writer, "}}")?;
    writer.flush()?;
    Ok(())
}
