//! Compacted De Bruijn Graph algorithms and data structures.

pub mod build;
pub mod cdbg;
pub mod local_multik;
pub mod simplify;

pub use build::{build_cdbg, build_cdbg_ext, build_cdbg_with_contigs, build_cdbg_with_contigs_ext, DirKmer};
pub use cdbg::{CdbgGraph, DirUnitigId, GraphEdge, ResolverLifecycleState, UnitigNode};
pub use local_multik::upgrade_multik;
pub use simplify::{clip_tips, edit_distance, pop_bubbles, recompact_graph, simplify_graph};
