//! Paired-end mapping and exSPAnder repeat resolution.

pub mod exspander;
pub mod mapper;

pub use exspander::{
    resolve_repeats_exspander, resolve_repeats_exspander_with_stats, EdgeResolveResult,
    ExspanderParams, JunctionStats, PathStep,
};
pub use mapper::{index_terminal_kmers, map_paired_reads, InsertSizeDist, PairLinkMatrix};
