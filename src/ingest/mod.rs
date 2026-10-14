//! Streaming ingest pipeline: genome size estimation, subsampling, trimming, and overlap merging.

pub mod estimator;
pub mod trim_merge;

pub use estimator::{compute_dynamic_kmers, estimate_genome, GenomeEstimate};
pub use trim_merge::{ingest_reads, merge_read_pair, trim_adapter, trim_quality, MergeResult};
