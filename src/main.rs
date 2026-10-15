//! fastBact: Ultra-fast Rust bacterial short-read genome assembler
//! fusing Shovill's pre/post-processing pipeline with SPAdes' multi-k cDBG & exSPAnder.

use anyhow::Result;
use clap::Parser;
use fastbact::cli::Cli;
use fastbact::graph::build::build_cdbg_ext;
use fastbact::graph::local_multik::{upgrade_multik_ext, upgrade_multik_legacy};
use fastbact::ingest::estimator::{estimate_genome, EstimatorConfidence};
use fastbact::ingest::trim_merge::ingest_reads_ext_with_seed;
use fastbact::polish::consensus::polish_contigs;
use fastbact::polish::output::{curate_contigs, write_outputs};
use fastbact::resolve::exspander::{resolve_repeats_exspander, ExspanderParams};
use fastbact::resolve::mapper::map_paired_reads;
use std::time::Instant;
use tracing::{info, warn};

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

    // B2-C1b-3B Decoupled Subsampling Policy:
    // 1) --no-subsample has highest priority -> 100% data retention.
    // 2) Subsampling is ONLY executed when explicitly requested via --enable-auto-subsample.
    //    Passing --genome-size or --depth alone provides priors/configuration, but NEVER discards reads.
    // 3) When --enable-auto-subsample is active:
    //    - If an explicit --genome-size (> 0) is provided, coverage is calculated against explicit_g.
    //    - Otherwise, coverage is calculated from HLL estimator.
    //    - If estimator confidence is LOW_COVERAGE_WARNING or UNRELIABLE, auto-subsampling
    //      is safely intercepted and downgraded to 100% data retention.
    let trigger_mult = if args.subsample_threshold_ratio > 0.0 {
        args.subsample_threshold_ratio
    } else {
        1.5
    };

    let (subsample_factor, subsampling_enabled, subsampling_reason) = if args.no_subsample {
        (1.0, false, "EXPLICIT_NO_SUBSAMPLE_FLAG".to_string())
    } else if args.enable_auto_subsample {
        if let Some(explicit_g) = args.genome_size {
            if explicit_g == 0 {
                (1.0, false, "EXPLICIT_GENOME_SIZE_ZERO_INVALID".to_string())
            } else {
                let explicit_cov = estimate.total_bases as f64 / (explicit_g as f64);
                let trigger_cov = trigger_mult * args.depth;
                if args.depth > 0.0 && explicit_cov > trigger_cov {
                    let factor = (args.depth / explicit_cov).clamp(0.01, 1.0);
                    info!(
                        "[AUTO_SUBSAMPLE] User-guided genome size cap (cov={:.1}x > trigger={:.1}x [{:.1}x]). Subsampling factor={:.4}",
                        explicit_cov, trigger_cov, trigger_mult, factor
                    );
                    (
                        factor,
                        true,
                        format!("AUTO_SUBSAMPLE_EXPLICIT_GENOME_SIZE (cov={:.1}x > trigger={:.1}x, factor={:.4})", explicit_cov, trigger_cov, factor),
                    )
                } else {
                    info!(
                        "[AUTO_SUBSAMPLE] User-guided genome size, but library coverage ({:.1}x) is within conservative trigger ({:.1}x [{:.1}x]). Subsampling bypassed.",
                        explicit_cov, trigger_cov, trigger_mult
                    );
                    (
                        1.0,
                        false,
                        format!("AUTO_SUBSAMPLE_EXPLICIT_GENOME_SIZE_WITHIN_TARGET (cov={:.1}x <= trigger={:.1}x)", explicit_cov, trigger_cov),
                    )
                }
            }
        } else {
            match estimate.confidence {
                EstimatorConfidence::High => {
                    let est_cov = estimate.estimated_cov;
                    let trigger_cov = trigger_mult * args.depth;
                    if args.depth > 0.0 && est_cov > trigger_cov {
                        let factor = (args.depth / est_cov).clamp(0.01, 1.0);
                        info!(
                            "[AUTO_SUBSAMPLE] High estimator confidence (est_cov={:.1}x > trigger={:.1}x [{:.1}x]). Auto-subsampling factor={:.4}",
                            est_cov, trigger_cov, trigger_mult, factor
                        );
                        (
                            factor,
                            true,
                            format!(
                                "AUTO_SUBSAMPLE_HIGH_CONFIDENCE (est_cov={:.1}x > trigger={:.1}x, factor={:.4})",
                                est_cov, trigger_cov, factor
                            ),
                        )
                    } else {
                        info!(
                            "[AUTO_SUBSAMPLE] High estimator confidence, but library coverage ({:.1}x) is within conservative trigger ({:.1}x [{:.1}x]). Subsampling bypassed.",
                            est_cov, trigger_cov, trigger_mult
                        );
                        (
                            1.0,
                            false,
                            format!(
                                "AUTO_SUBSAMPLE_HIGH_CONFIDENCE_WITHIN_TARGET (est_cov={:.1}x <= trigger={:.1}x)",
                                est_cov, trigger_cov
                            ),
                        )
                    }
                }
                EstimatorConfidence::LowCoverageWarning => {
                    warn!(
                        "[SAFETY_GUARD] Auto-subsampling intercepted: Estimator confidence is LOW_COVERAGE_WARNING (est_cov={:.1}x, solid_depth={:.2}x). Subsampling safely downgraded to 100% data retention.",
                        estimate.estimated_cov, estimate.sample_kmer_depth
                    );
                    (
                        1.0,
                        false,
                        "AUTO_SUBSAMPLE_INTERCEPTED_LOW_COVERAGE_WARNING".to_string(),
                    )
                }
                EstimatorConfidence::Unreliable => {
                    warn!(
                        "[SAFETY_GUARD] Auto-subsampling intercepted: Estimator confidence is UNRELIABLE (est_cov={:.1}x). Subsampling safely downgraded to 100% data retention.",
                        estimate.estimated_cov
                    );
                    (
                        1.0,
                        false,
                        "AUTO_SUBSAMPLE_INTERCEPTED_UNRELIABLE".to_string(),
                    )
                }
            }
        }
    } else if args.genome_size.is_some() {
        info!("[SUBSAMPLE] Prior genome size specified without --enable-auto-subsample. Retaining 100% reads for assembly.");
        (1.0, false, "GENOME_SIZE_PRIOR_PROVIDED_SUBSAMPLING_NOT_ENABLED".to_string())
    } else {
        // Safe Default: Subsampling is disabled unless explicitly requested via --enable-auto-subsample
        (1.0, false, "DISABLED_BY_DEFAULT_SAFE_FALLBACK".to_string())
    };

    info!("--- Input Sequencing Statistics (C1b-3 Calibrated Ingest) ---");
    info!("  exact_read_pairs: {}", estimate.exact_read_pairs);
    info!("  total_bases: {} bp", estimate.total_bases);
    info!("  median_read_len: {} bp", estimate.median_read_len);
    info!("  max_read_len: {} bp", estimate.max_read_len);
    info!("  estimated_genome_size: {} bp", estimate.genome_size);
    info!("  estimated_coverage: {:.1}x", estimate.estimated_cov);
    info!("  estimator_confidence: {}", estimate.confidence);
    info!("  solid_unique_kmers_u2: {}", estimate.u_ge_2);
    info!("  sample_solid_kmer_depth: {:.2}x", estimate.sample_kmer_depth);
    info!("--- Subsampling Policy (C1b-3 Calibrated Policy) ---");
    info!("  subsampling_enabled: {}", subsampling_enabled);
    info!("  subsampling_reason: {}", subsampling_reason);
    info!("  subsample_factor: {:.4}", subsample_factor);
    info!("-----------------------------------------------------");

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
    info!("[Stage 2/6] Streaming FASTQ reads, adapter trimming, quality clipping & FLASH merging (seed={})...", args.subsample_seed);
    let t_s2 = Instant::now();
    let (bank, ingest_stats) = ingest_reads_ext_with_seed(
        &args.r1,
        &args.r2,
        subsample_factor,
        30,
        20,
        !args.no_trim,
        !args.no_merge,
        args.subsample_seed,
    )?;
    info!(
        "Ingested {} merged reads ({} bp), {} unmerged PE reads (elapsed: {:.2?})",
        bank.merged_reads.len(),
        bank.merged_reads.total_bases(),
        bank.pe_forward.len(),
        t_s2.elapsed()
    );
    info!("--- Actual Ingestion Counts ---");
    info!("  actual_input_pairs: {}", ingest_stats.input_read_pairs);
    info!("  actual_retained_pairs: {}", ingest_stats.retained_read_pairs);
    info!("-------------------------------");

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
    let effective_g = args.genome_size.filter(|&g| g > 0).unwrap_or(estimate.genome_size);
    let g_source = if args.genome_size.filter(|&g| g > 0).is_some() {
        "CLI_USER_SPECIFIED"
    } else {
        "ESTIMATOR_HLL"
    };
    let k1_cms_width = fastbact::common::sketch::CountMinSketch::width_for_genome_size(effective_g);
    let k1_cms_bytes = fastbact::common::sketch::CountMinSketch::allocated_bytes_for_genome_size(effective_g);
    info!(
        "[Stage 3/6] Building compacted De Bruijn graph at k={} (genome sketch size: {} bp, source: {}, CMS width={}, bytes={} MB)...",
        k1, effective_g, g_source, k1_cms_width, k1_cms_bytes / (1024 * 1024)
    );
    let t_s3 = Instant::now();
    let mut graph = build_cdbg_ext::<u64>(&bank, k1, None, Some(effective_g));

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
            genome_size_source: g_source.to_string(),
            genome_size_used: effective_g,
            cms_width: k1_cms_width,
            cms_allocated_bytes: k1_cms_bytes,
            upgrade_executed: true,
            early_stop_reason: "NONE_INITIAL_BUILD".to_string(),
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
            ..Default::default()
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
            let upgrade_res = if args.legacy_multik_cms {
                info!("Upgrading graph to k={} (LEGACY Var 1 capacity path)...", next_k);
                if next_k <= 63 {
                    upgrade_multik_legacy::<u128>(&mut graph, &bank, next_k)
                } else {
                    upgrade_multik_legacy::<[u64; 4]>(&mut graph, &bank, next_k)
                }
            } else {
                info!(
                    "Upgrading graph to k={} with explicit genome size {} bp (Var 2 propagation)...",
                    next_k, effective_g
                );
                if next_k <= 63 {
                    upgrade_multik_ext::<u128>(&mut graph, &bank, next_k, Some(effective_g))
                } else {
                    upgrade_multik_ext::<[u64; 4]>(&mut graph, &bank, next_k, Some(effective_g))
                }
            };
            info!(
                "[MULTI-K TELEMETRY] k={}: executed={}, early_stop={:?}, CMS source={}, size_used={} bp, width={}, bytes={} MB, unitigs: {} -> {}",
                next_k,
                upgrade_res.executed,
                upgrade_res.early_stop_reason,
                upgrade_res.genome_size_source,
                upgrade_res.genome_size_used,
                upgrade_res.cms_width,
                upgrade_res.cms_allocated_bytes / (1024 * 1024),
                upgrade_res.active_unitigs_before,
                upgrade_res.active_unitigs_after
            );
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
                    genome_size_source: upgrade_res.genome_size_source.to_string(),
                    genome_size_used: upgrade_res.genome_size_used,
                    cms_width: upgrade_res.cms_width,
                    cms_allocated_bytes: upgrade_res.cms_allocated_bytes,
                    upgrade_executed: upgrade_res.executed,
                    early_stop_reason: upgrade_res.early_stop_reason.unwrap_or("NONE").to_string(),
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
                    ..Default::default()
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
                ..Default::default()
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
            subsampling_enabled,
            subsampling_reason: subsampling_reason.clone(),
            estimated_genome_size: estimate.genome_size,
            estimated_genome_size_note: "C1B2_HLL_SPLITMIX64_POISSON_CALIBRATED".into(),
            estimated_coverage: estimate.estimated_cov,
            estimator_confidence: estimate.confidence.to_string(),
            actual_input_pairs: ingest_stats.input_read_pairs,
            actual_retained_pairs: ingest_stats.retained_read_pairs,
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
    info!("  Subsampling Enabled: {}", subsampling_enabled);
    info!("  Subsampling Reason: {}", subsampling_reason);
    info!("  Actual Input Pairs: {}", ingest_stats.input_read_pairs);
    info!("  Actual Retained Pairs: {}", ingest_stats.retained_read_pairs);
    info!("  Peak Memory (RSS): {:.1} MB", read_peak_rss_mb());
    info!("  Total Wall Time: {:.2?}", t_start.elapsed());
    info!("============================================================");

    Ok(())
}
