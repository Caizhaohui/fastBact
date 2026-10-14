//! Foundational data structures: 2-bit bitseq, polymorphic k-mers, and streaming sketches.

pub mod bitseq;
pub mod kmer;
pub mod sketch;

pub use bitseq::{Base, BitSeqArena, ReadBank};
pub use kmer::{KmerInt, RollingKmer, RollingNtHash};
pub use sketch::{CountMinSketch, HyperLogLog};
