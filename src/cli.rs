//! Command-line argument parser for fastBact.

use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug, Clone)]
#[command(
    name = "fastBact",
    author = "fastBact Developers",
    version = "0.1.0",
    about = "Ultra-fast bacterial short-read genome assembler fusing Shovill and SPAdes"
)]
pub struct Cli {
    /// Forward paired-end reads (FASTQ, optionally gzipped)
    #[arg(short = '1', long = "r1", required = true)]
    pub r1: PathBuf,

    /// Reverse paired-end reads (FASTQ, optionally gzipped)
    #[arg(short = '2', long = "r2", required = true)]
    pub r2: PathBuf,

    /// Output directory for assembly results
    #[arg(short = 'o', long = "outdir", default_value = "fastbact_out")]
    pub outdir: PathBuf,

    /// Number of CPU threads
    #[arg(short = 't', long = "threads", default_value_t = 8)]
    pub threads: usize,

    /// Target subsampling sequencing depth (0 to disable)
    #[arg(long = "depth", default_value_t = 120.0)]
    pub depth: f64,

    /// Minimum contig length in final FASTA output (bp)
    #[arg(long = "min-len", default_value_t = 500)]
    pub min_len: usize,

    /// Minimum k-mer coverage depth for contig retention
    #[arg(long = "min-cov", default_value_t = 2.0)]
    pub min_cov: f32,

    /// Comma-separated list of k-mers (e.g. "31,55,79"), or "auto"
    #[arg(short = 'k', long = "kmers", default_value = "auto")]
    pub kmers: String,

    /// Bypass adapter & sliding-window quality trimming
    #[arg(long, hide = true)]
    pub no_trim: bool,

    /// Bypass FLASH-style PE overlap merging (keep all reads as unmerged PE)
    #[arg(long, hide = true)]
    pub no_merge: bool,

    /// Bypass coverage subsampling (keep 100% of reads)
    #[arg(long, hide = true)]
    pub no_subsample: bool,

    /// Bypass tip clipping in simplify_graph
    #[arg(long, hide = true)]
    pub no_tip_clip: bool,

    /// Bypass bubble popping in simplify_graph
    #[arg(long, hide = true)]
    pub no_bubble_pop: bool,

    /// Bypass local_multik upgrades (run only initial k1)
    #[arg(long, hide = true)]
    pub no_multik: bool,

    /// Bypass exSPAnder PE repeat resolution
    #[arg(long, hide = true)]
    pub no_exspander: bool,

    /// Bypass consensus pileup polishing
    #[arg(long, hide = true)]
    pub no_polish: bool,

    /// Directory to dump intermediate stage FASTAs & telemetry JSON
    #[arg(long, hide = true)]
    pub dump_stages: Option<PathBuf>,
}
