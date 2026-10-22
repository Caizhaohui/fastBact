//! fastBact: Ultra-fast Rust bacterial short-read genome assembler
//! fusing Shovill's pre/post-processing pipeline with SPAdes' multi-k cDBG & exSPAnder.

pub mod cli;
pub mod common;
pub mod diag;
pub mod graph;
pub mod ingest;
pub mod polish;
pub mod resolve;

pub use graph::ResolverLifecycleState;
pub use resolve::ExspanderError;
