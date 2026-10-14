//! Compacted De Bruijn Graph data structures and unitig container.

use crate::common::bitseq::BitSeqArena;
use smallvec::SmallVec;

/// Directed stranded node pointer (Unitig ID + Orientation).
/// Lowest bit: strand (0 = forward, 1 = reverse complement).
/// Upper 31 bits: UnitigId.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct DirUnitigId(pub u32);

impl DirUnitigId {
    #[inline(always)]
    pub fn new(id: u32, is_rc: bool) -> Self {
        Self((id << 1) | (is_rc as u32))
    }

    #[inline(always)]
    pub fn id(self) -> u32 {
        self.0 >> 1
    }

    #[inline(always)]
    pub fn is_rc(self) -> bool {
        (self.0 & 1) == 1
    }

    #[inline(always)]
    pub fn twin(self) -> Self {
        Self(self.0 ^ 1)
    }

    #[inline(always)]
    pub fn forward(self) -> Self {
        Self(self.0 & !1)
    }

    #[inline(always)]
    pub fn reverse(self) -> Self {
        Self(self.0 | 1)
    }
}

/// A directed edge connecting two unitigs with an exact k-1 overlap.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct GraphEdge {
    pub target: DirUnitigId,
    pub overlap_k: u8, // Overlap length (typically k - 1)
    pub support: u16,  // Read transition support count
}

impl GraphEdge {
    pub fn new(target: DirUnitigId, overlap_k: u8, support: u16) -> Self {
        Self {
            target,
            overlap_k,
            support,
        }
    }
}

/// A vertex representing a maximal non-branching sequence (unitig) in the cDBG.
#[derive(Clone, Debug)]
pub struct UnitigNode {
    pub id: u32,
    pub seq_id: u32,     // Sequence ID inside graph's BitSeqArena
    pub seq_len: u32,    // Length in bases (>= k)
    pub kmer_cov: f32,   // Mean k-mer coverage depth
    pub edges_fwd: SmallVec<[GraphEdge; 4]>, // Outgoing from 3' end of forward strand
    pub edges_rev: SmallVec<[GraphEdge; 4]>, // Outgoing from 3' end of reverse strand
    pub is_deleted: bool,
    pub is_circular: bool,
}

impl UnitigNode {
    pub fn new(id: u32, seq_id: u32, seq_len: u32, kmer_cov: f32) -> Self {
        Self {
            id,
            seq_id,
            seq_len,
            kmer_cov,
            edges_fwd: SmallVec::new(),
            edges_rev: SmallVec::new(),
            is_deleted: false,
            is_circular: false,
        }
    }

    #[inline(always)]
    pub fn edges(&self, is_rc: bool) -> &SmallVec<[GraphEdge; 4]> {
        if is_rc {
            &self.edges_rev
        } else {
            &self.edges_fwd
        }
    }

    #[inline(always)]
    pub fn edges_mut(&mut self, is_rc: bool) -> &mut SmallVec<[GraphEdge; 4]> {
        if is_rc {
            &mut self.edges_rev
        } else {
            &mut self.edges_fwd
        }
    }

    #[inline(always)]
    pub fn add_edge(&mut self, is_rc: bool, edge: GraphEdge) {
        let edge_list = self.edges_mut(is_rc);
        for existing in edge_list.iter_mut() {
            if existing.target == edge.target && existing.overlap_k == edge.overlap_k {
                existing.support = existing.support.saturating_add(edge.support);
                return;
            }
        }
        edge_list.push(edge);
    }
}

/// Full Compacted De Bruijn Graph containing sequence arena and flat node array.
#[derive(Clone, Debug)]
pub struct CdbgGraph {
    pub k: usize,
    pub arena: BitSeqArena,
    pub nodes: Vec<UnitigNode>,
}

impl CdbgGraph {
    pub fn new(k: usize) -> Self {
        Self {
            k,
            arena: BitSeqArena::new(),
            nodes: Vec::new(),
        }
    }

    pub fn num_active_nodes(&self) -> usize {
        self.nodes.iter().filter(|n| !n.is_deleted).count()
    }

    pub fn total_active_bases(&self) -> usize {
        self.nodes.iter().filter(|n| !n.is_deleted).map(|n| n.seq_len as usize).sum()
    }

    /// Retrieves unitig sequence as ASCII bytes in forward orientation.
    pub fn get_unitig_seq(&self, id: u32) -> Vec<u8> {
        let seq_id = self.nodes[id as usize].seq_id;
        self.arena.extract_ascii(seq_id)
    }

    /// Retrieves stranded unitig sequence (forward or reverse complement).
    pub fn get_stranded_seq(&self, target: DirUnitigId) -> Vec<u8> {
        let seq_id = self.nodes[target.id() as usize].seq_id;
        if target.is_rc() {
            self.arena.extract_rc_ascii(seq_id)
        } else {
            self.arena.extract_ascii(seq_id)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dir_unitig_id() {
        let u_fwd = DirUnitigId::new(42, false);
        let u_rev = DirUnitigId::new(42, true);

        assert_eq!(u_fwd.id(), 42);
        assert!(!u_fwd.is_rc());
        assert_eq!(u_rev.id(), 42);
        assert!(u_rev.is_rc());

        assert_eq!(u_fwd.twin(), u_rev);
        assert_eq!(u_rev.twin(), u_fwd);
    }
}
