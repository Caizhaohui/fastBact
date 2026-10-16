//! Consensus polishing, curation, and FASTA/GFA output.

pub mod consensus;
pub mod output;

pub use consensus::{polish_contigs, PolishStats};
pub use output::{curate_contigs, write_outputs, AssembledContig};
