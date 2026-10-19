use fastbact::common::bitseq::reverse_complement_ascii;
use fastbact::common::sketch::splitmix64;
use fastbact::graph::cdbg::{CdbgGraph, UnitigNode};
use fastbact::polish::consensus::polish_contigs;

fn make_deterministic_dna(seed: u64, len: usize) -> Vec<u8> {
    let alphabet = [b'A', b'C', b'G', b'T'];
    let mut buf = Vec::with_capacity(len);
    for i in 0..len {
        let r = splitmix64(seed.wrapping_add(i as u64));
        buf.push(alphabet[(r % 4) as usize]);
    }
    buf
}

fn contains_either_strand(seq: &[u8], kmer: &[u8]) -> bool {
    let rc = reverse_complement_ascii(kmer);
    seq.windows(kmer.len())
        .any(|window| window == kmer || window == rc.as_slice())
}

/// A residual node below the FASTA cutoff must not capture the first seed of a read
/// that also anchors on an output contig. Otherwise the output contig loses the pileup.
#[test]
fn residual_short_node_does_not_steal_output_contig_anchor() {
    let long = make_deterministic_dna(0xC0_02_0013, 600);
    let c_pos = 200;
    let r_pos = 60;
    let vote_pos = c_pos - 1;
    let original = long[vote_pos];
    let voted = if original == b'T' { b'A' } else { b'T' };

    let mut guard = 0u64;
    let short = loop {
        let candidate = make_deterministic_dna(0x51_05_7000 + guard, 40);
        if !contains_either_strand(&long, &candidate[..31]) {
            break candidate;
        }
        guard += 1;
        assert!(
            guard < 20,
            "could not find a short k-mer absent from the contig"
        );
    };

    let mut graph = CdbgGraph::new(31);
    let long_id = graph.arena.push_seq(&long);
    graph
        .nodes
        .push(UnitigNode::new(0, long_id, long.len() as u32, 20.0));
    let short_id = graph.arena.push_seq(&short);
    graph
        .nodes
        .push(UnitigNode::new(1, short_id, short.len() as u32, 20.0));

    let mut read = vec![b'A'; 120];
    read[..31].copy_from_slice(&short[..31]);
    read[r_pos..r_pos + 31].copy_from_slice(&long[c_pos..c_pos + 31]);
    read[r_pos - 1] = voted;
    assert_eq!(c_pos - r_pos + (r_pos - 1), vote_pos);

    let short_before = graph.get_unitig_seq(1);
    polish_contigs(&mut graph, std::iter::repeat(read).take(8), 500);

    assert_eq!(
        graph.get_unitig_seq(0)[vote_pos],
        voted,
        "output contig must receive the read pileup"
    );
    assert_eq!(
        graph.get_unitig_seq(1),
        short_before,
        "residual node below the FASTA cutoff is not a consensus target"
    );
}
