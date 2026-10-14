//! Paired-end mapping and exSPAnder repeat resolution.

pub mod exspander;
pub mod mapper;

pub use exspander::{resolve_repeats_exspander, ExspanderParams};
pub use mapper::{index_terminal_kmers, map_paired_reads, InsertSizeDist, PairLinkMatrix};
