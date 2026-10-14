//! fastBact: Ultra-fast Rust bacterial short-read genome assembler
//! fusing Shovill's pre/post-processing pipeline with SPAdes' multi-k cDBG & exSPAnder.

use anyhow::Result;
use clap::Parser;
use fastbact::cli::Cli;
use fastbact::graph::build::build_cdbg;
use fastbact::graph::local_multik::upgrade_multik;
use fastbact::ingest::estimator::estimate_genome;
use fastbact::ingest::trim_merge::ingest_reads_ext;
use fastbact::polish::consensus::polish_contigs;
use fastbact::polish::output::{curate_contigs, write_outputs};
use fastbact::resolve::exspander::{resolve_repeats_exspander, ExspanderParams};
use fastbact::resolve::mapper::map_paired_reads;
use std::time::Instant;
use tracing::info;

fn read_peak_rss_mb() -> f64 {
    #[cfg(target_os = "linux")]
    {
        if let Ok(status) = std::fs::read_to_string("/proc/self/status") {
            for line in status.lines() {
                if line.starts_with("VmHWM:") {
                    let parts: Vec<&str> = line.split_whitespace().collect();
                    if parts.len() >= 2 {
                        if let Ok(kb) = parts[1].parse::<f64>() {
                            return kb / 1024.0;
                        }
                    }
                }
            }
        }
    }
    0.0
}

fn main() -> Result<()> {
    let t_start = Instant::now();
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .init();

    let args = Cli::parse();

    info!("=== fastBact v0.1.0 (Rust Bacterial Genome Assembler) ===");
    info!("R1: {:?}", args.r1);
    info!("R2: {:?}", args.r2);
    info!("Outdir: {:?}", args.outdir);
    info!("Threads: {}", args.threads);

    // Set thread pool
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .ok();
    }

    // Stage 1: Genome Estimation
    info!("[Stage 1/6] Estimating genome size and coverage from reads...");
    let t_s1 = Instant::now();
    let estimate = estimate_genome(&args.r1, &args.r2, 250_000, args.depth)?;
    info!(
        "Estimated genome size: {} bp, coverage: {:.1}x (subsample factor: {:.4})",
        estimate.genome_size, estimate.estimated_cov, estimate.subsample_factor
    );
    let subsample_factor = if args.no_subsample {
        1.0
    } else {
        estimate.subsample_factor
    };
    let kmers = if args.kmers == "auto" {
        estimate.kmers
    } else {
        args.kmers
            .split(',')
            .filter_map(|s| s.trim().parse::<usize>().ok())
            .collect()
    };
    info!(
        "Multi-k schedule: {:?} (elapsed: {:.2?})",
        kmers,
        t_s1.elapsed()
    );

    // Stage 2: Streaming Ingest (Trimming & Overlap Merging)
    info!("[Stage 2/6] Streaming FASTQ reads, adapter trimming, quality clipping & FLASH merging...");
    let t_s2 = Instant::now();
    let (bank, ingest_stats) = ingest_reads_ext(
        &args.r1,
        &args.r2,
        subsample_factor,
        30,
        20,
        !args.no_trim,
        !args.no_merge,
    )?;
    info!(
        "Ingested {} merged reads ({} bp), {} unmerged PE reads (elapsed: {:.2?})",
        bank.merged_reads.len(),
        bank.merged_reads.total_bases(),
        bank.pe_forward.len(),
        t_s2.elapsed()
    );

    let mut stage_telemetries = Vec::new();
    let k1 = kmers[0];

    // If dump_stages is active, prepare directory and initial k-mer telemetry
    let kmer_telemetry = if args.dump_stages.is_some() {
        std::fs::create_dir_all(args.dump_stages.as_ref().unwrap())?;
        fastbact::diag::sample_kmer_telemetry(&bank, k1, estimate.genome_size)
    } else {
        Default::default()
    };

    // Stage 3: Initial cDBG Build
    info!("[Stage 3/6] Building compacted De Bruijn graph at k={}...", k1);
    let t_s3 = Instant::now();
    let mut graph = build_cdbg::<u64>(&bank, k1, None);

    if let Some(ref dump_dir) = args.dump_stages {
        fastbact::diag::dump_graph_fasta(&graph, dump_dir.join("stage_01_raw_k1_cdbg.fa"))?;
        let (u, e, c, bp) = fastbact::diag::inspect_graph(&graph);
        stage_telemetries.push(fastbact::diag::StageGraphTelemetry {
            stage_name: "stage_01_raw_k1_cdbg".into(),
            k: k1,
            active_unitigs: u,
            active_edges: e,
            isolated_cycles: c,
            total_graph_bp: bp,
            ..Default::default()
        });

        // 1. Tip clipping
        let (tips, tip_bp) = if !args.no_tip_clip {
            fastbact::graph::simplify::clip_tips_tracked(&mut graph)
        } else {
            (0, 0)
        };
        fastbact::graph::simplify::recompact_graph(&mut graph);
        fastbact::diag::dump_graph_fasta(&graph, dump_dir.join("stage_02_k1_tip_clipped.fa"))?;
        let (u2, e2, c2, bp2) = fastbact::diag::inspect_graph(&graph);
        stage_telemetries.push(fastbact::diag::StageGraphTelemetry {
            stage_name: "stage_02_k1_tip_clipped".into(),
            k: k1,
            active_unitigs: u2,
            active_edges: e2,
            isolated_cycles: c2,
            deleted_tips_count: tips,
            deleted_tips_bp: tip_bp,
            total_graph_bp: bp2,
            ..Default::default()
        });

        // 2. Bubble popping and complete simplification convergence
        let mut sim_stats = fastbact::graph::simplify::SimplifyStats::default();
        if !args.no_bubble_pop {
            let (b, b_bp) = fastbact::graph::simplify::pop_bubbles_tracked(&mut graph);
            sim_stats.popped_bubbles_count += b;
            sim_stats.popped_bubbles_bp += b_bp;
        }
        fastbact::graph::simplify::recompact_graph(&mut graph);
        let rest = fastbact::graph::simplify::simplify_graph_custom(&mut graph, !args.no_tip_clip, !args.no_bubble_pop);
        sim_stats.deleted_tips_count += rest.deleted_tips_count;
        sim_stats.deleted_tips_bp += rest.deleted_tips_bp;
        sim_stats.popped_bubbles_count += rest.popped_bubbles_count;
        sim_stats.popped_bubbles_bp += rest.popped_bubbles_bp;

        fastbact::diag::dump_graph_fasta(&graph, dump_dir.join("stage_03_k1_bubble_popped.fa"))?;
        let (u3, e3, c3, bp3) = fastbact::diag::inspect_graph(&graph);
        stage_telemetries.push(fastbact::diag::StageGraphTelemetry {
            stage_name: "stage_03_k1_bubble_popped".into(),
            k: k1,
            active_unitigs: u3,
            active_edges: e3,
            isolated_cycles: c3,
            deleted_tips_count: tips + sim_stats.deleted_tips_count,
            deleted_tips_bp: tip_bp + sim_stats.deleted_tips_bp,
            popped_bubbles_count: sim_stats.popped_bubbles_count,
            popped_bubbles_bp: sim_stats.popped_bubbles_bp,
            total_graph_bp: bp3,
        });
    } else {
        fastbact::graph::simplify::simplify_graph_custom(&mut graph, !args.no_tip_clip, !args.no_bubble_pop);
    }
    info!(
        "Initial cDBG at k={}: {} active unitigs (elapsed: {:.2?})",
        k1,
        graph.num_active_nodes(),
        t_s3.elapsed()
    );

    // Stage 4: Multi-k Upgrades
    if !args.no_multik {
        info!("[Stage 4/6] Iterative local multi-k repeat bridging...");
        let t_s4 = Instant::now();
        for &next_k in &kmers[1..] {
            info!("Upgrading graph to k={}...", next_k);
            if next_k <= 63 {
                upgrade_multik::<u128>(&mut graph, &bank, next_k);
            } else {
                upgrade_multik::<[u64; 4]>(&mut graph, &bank, next_k);
            }
            if let Some(ref dump_dir) = args.dump_stages {
                fastbact::diag::dump_graph_fasta(&graph, dump_dir.join(format!("stage_04_multik_{}_raw.fa", next_k)))?;
                let (u, e, c, bp) = fastbact::diag::inspect_graph(&graph);
                stage_telemetries.push(fastbact::diag::StageGraphTelemetry {
                    stage_name: format!("stage_04_multik_{}_raw", next_k),
                    k: next_k,
                    active_unitigs: u,
                    active_edges: e,
                    isolated_cycles: c,
                    total_graph_bp: bp,
                    ..Default::default()
                });
            }
            let sim_stats = fastbact::graph::simplify::simplify_graph_custom(&mut graph, !args.no_tip_clip, !args.no_bubble_pop);
            if let Some(ref dump_dir) = args.dump_stages {
                fastbact::diag::dump_graph_fasta(&graph, dump_dir.join(format!("stage_04_multik_{}_simplified.fa", next_k)))?;
                let (u, e, c, bp) = fastbact::diag::inspect_graph(&graph);
                stage_telemetries.push(fastbact::diag::StageGraphTelemetry {
                    stage_name: format!("stage_04_multik_{}_simplified", next_k),
                    k: next_k,
                    active_unitigs: u,
                    active_edges: e,
                    isolated_cycles: c,
                    deleted_tips_count: sim_stats.deleted_tips_count,
                    deleted_tips_bp: sim_stats.deleted_tips_bp,
                    popped_bubbles_count: sim_stats.popped_bubbles_count,
                    popped_bubbles_bp: sim_stats.popped_bubbles_bp,
                    total_graph_bp: bp,
                });
            }
        }
        info!(
            "Graph after multi-k: {} active unitigs (elapsed: {:.2?})",
            graph.num_active_nodes(),
            t_s4.elapsed()
        );
    }

    // Stage 5: exSPAnder Paired-End Repeat Resolution
    if !args.no_exspander {
        info!("[Stage 5/6] Resolving repeats with exSPAnder paired-end insert size constraints...");
        let t_s5 = Instant::now();
        let (matrix, dist) = map_paired_reads::<u64>(&graph, &bank, k1);
        info!(
            "Mapped paired reads: median insert size = {:.1} +/- {:.1} bp ({} links)",
            dist.median,
            dist.mad,
            matrix.links.len()
        );
        let params = ExspanderParams::default();
        resolve_repeats_exspander(&mut graph, &matrix, &dist, params);
        let sim_stats = fastbact::graph::simplify::simplify_graph_custom(&mut graph, !args.no_tip_clip, !args.no_bubble_pop);
        info!(
            "Graph after exSPAnder: {} active unitigs (elapsed: {:.2?})",
            graph.num_active_nodes(),
            t_s5.elapsed()
        );
        if let Some(ref dump_dir) = args.dump_stages {
            fastbact::diag::dump_graph_fasta(&graph, dump_dir.join("stage_05_exspander.fa"))?;
            let (u, e, c, bp) = fastbact::diag::inspect_graph(&graph);
            stage_telemetries.push(fastbact::diag::StageGraphTelemetry {
                stage_name: "stage_05_exspander".into(),
                k: *kmers.last().unwrap_or(&k1),
                active_unitigs: u,
                active_edges: e,
                isolated_cycles: c,
                deleted_tips_count: sim_stats.deleted_tips_count,
                deleted_tips_bp: sim_stats.deleted_tips_bp,
                popped_bubbles_count: sim_stats.popped_bubbles_count,
                popped_bubbles_bp: sim_stats.popped_bubbles_bp,
                total_graph_bp: bp,
            });
        }
    }

    // Stage 6: Consensus Polishing & Final Output
    info!("[Stage 6/6] In-memory consensus polishing and FASTA/GFA curation...");
    let t_s6 = Instant::now();

    if !args.no_polish {
        let read_iter = (0..bank.merged_reads.len())
            .map(|i| bank.merged_reads.extract_ascii(i as u32))
            .chain((0..bank.pe_forward.len()).map(|i| bank.pe_forward.extract_ascii(i as u32)));
        polish_contigs(&mut graph, read_iter);
    }

    let contigs = curate_contigs(&graph, args.min_len, args.min_cov);
    write_outputs(&graph, &contigs, &args.outdir)?;

    if let Some(ref dump_dir) = args.dump_stages {
        fastbact::diag::dump_contigs_fasta(&contigs, dump_dir.join("stage_06_final_polished.fa"))?;
        let (_, e, c, bp) = fastbact::diag::inspect_graph(&graph);
        stage_telemetries.push(fastbact::diag::StageGraphTelemetry {
            stage_name: "stage_06_final_polished".into(),
            k: *kmers.last().unwrap_or(&k1),
            active_unitigs: contigs.len(),
            active_edges: e,
            isolated_cycles: c,
            total_graph_bp: bp,
            ..Default::default()
        });

        let effective_cov = if estimate.genome_size > 0 {
            ingest_stats.bases_after_trim as f64 / estimate.genome_size as f64
        } else {
            0.0
        };
        let ingest_telemetry = fastbact::diag::IngestTelemetry {
            input_read_pairs: ingest_stats.input_read_pairs,
            retained_read_pairs: ingest_stats.retained_read_pairs,
            bases_before_trim: ingest_stats.bases_before_trim,
            bases_after_trim: ingest_stats.bases_after_trim,
            merged_pairs_count: ingest_stats.merged_pairs_count,
            merged_bases: ingest_stats.merged_bases,
            unmerged_pairs_count: ingest_stats.unmerged_pairs_count,
            effective_coverage: effective_cov,
        };

        fastbact::diag::write_stage_telemetry_json(
            dump_dir.join("stage_telemetry.json"),
            &ingest_telemetry,
            &kmer_telemetry,
            &stage_telemetries,
        )?;
    }

    let total_len: usize = contigs.iter().map(|c| c.len).sum();
    let max_len = contigs.iter().map(|c| c.len).max().unwrap_or(0);
    let circular_count = contigs.iter().filter(|c| c.is_circular).count();

    // Calculate N50
    let mut sorted_lens: Vec<usize> = contigs.iter().map(|c| c.len).collect();
    sorted_lens.sort_unstable_by(|a, b| b.cmp(a));
    let half_len = total_len / 2;
    let mut acc = 0;
    let mut n50 = 0;
    for &l in &sorted_lens {
        acc += l;
        if acc >= half_len {
            n50 = l;
            break;
        }
    }

    info!(
        "Polishing and curation completed in {:.2?}. Output: {:?}",
        t_s6.elapsed(),
        args.outdir
    );
    info!("============================================================");
    info!("Assembly Summary:");
    info!("  Total Contigs: {}", contigs.len());
    info!("  Total Assembled Length: {} bp", total_len);
    info!("  Longest Contig: {} bp", max_len);
    info!("  N50: {} bp", n50);
    info!("  Circular Plasmids: {}", circular_count);
    info!("  Peak Memory (RSS): {:.1} MB", read_peak_rss_mb());
    info!("  Total Wall Time: {:.2?}", t_start.elapsed());
    info!("============================================================");

    Ok(())
}
