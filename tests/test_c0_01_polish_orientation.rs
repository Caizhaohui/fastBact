use fastbact::common::bitseq::reverse_complement_ascii;
use fastbact::common::sketch::splitmix64;
use fastbact::graph::cdbg::{CdbgGraph, UnitigNode};
use fastbact::polish::consensus::polish_contigs;

/// Generates a deterministic non-palindromic, non-symmetric DNA sequence of given length.
fn make_deterministic_dna(seed: u64, len: usize) -> Vec<u8> {
    let alphabet = [b'A', b'C', b'G', b'T'];
    let mut buf = Vec::with_capacity(len);
    for i in 0..len {
        let r = splitmix64(seed + i as u64);
        buf.push(alphabet[(r % 4) as usize]);
    }
    assert_ne!(buf, reverse_complement_ascii(&buf));
    buf
}

/// Helper to create a single-unitig CdbgGraph
fn create_test_graph(contig_seq: &[u8]) -> CdbgGraph {
    let mut graph = CdbgGraph::new(31);
    let seq_id = graph.arena.push_seq(contig_seq);
    graph.nodes.push(UnitigNode::new(0, seq_id, contig_seq.len() as u32, 20.0));
    graph
}

/// Helper to generate a reverse read anchored at a specific contig seed (c_pos)
/// with a specific read seed position (r_pos).
fn make_reverse_read(contig: &[u8], c_pos: usize, r_pos: usize, read_len: usize) -> Vec<u8> {
    let k = 31;
    let max_contig_pos = c_pos + k - 1 + r_pos;
    let min_contig_pos = max_contig_pos - (read_len - 1);
    assert!(max_contig_pos < contig.len(), "Contig overflow: max_pos={}", max_contig_pos);
    assert!(min_contig_pos <= c_pos, "Seed must be inside read window: min_pos={}", min_contig_pos);

    let contig_window = &contig[min_contig_pos..=max_contig_pos];
    assert_eq!(contig_window.len(), read_len);
    let read = reverse_complement_ascii(contig_window);

    let contig_seed = &contig[c_pos..c_pos + k];
    let read_seed = &read[r_pos..r_pos + k];
    assert_eq!(*read_seed, *reverse_complement_ascii(contig_seed));

    read
}

/// Helper to generate a forward read anchored at a specific contig seed (c_pos)
/// with a specific read seed position (r_pos).
fn make_forward_read(contig: &[u8], c_pos: usize, r_pos: usize, read_len: usize) -> Vec<u8> {
    let k = 31;
    assert!(c_pos >= r_pos, "Contig underflow");
    let contig_start = c_pos - r_pos;
    let contig_end = contig_start + read_len;
    assert!(contig_end <= contig.len(), "Contig overflow");

    let read = contig[contig_start..contig_end].to_vec();
    let contig_seed = &contig[c_pos..c_pos + k];
    let read_seed = &read[r_pos..r_pos + k];
    assert_eq!(contig_seed, read_seed);

    read
}

fn run_matrix_case(len: usize, offset: usize, mode: &str) {
    let contig_seq = make_deterministic_dna(0x9876_5432, 500);
    let mut graph = create_test_graph(&contig_seq);
    let c_pos = 200;

    let fwd = make_forward_read(&contig_seq, c_pos, offset, len);
    let rev = make_reverse_read(&contig_seq, c_pos, offset, len);

    let reads = match mode {
        "fwd" => vec![fwd; 10],
        "rev" => vec![rev; 10],
        "mixed" => {
            let mut r = vec![fwd; 5];
            r.extend(vec![rev; 5]);
            r
        }
        _ => unreachable!(),
    };

    let stats = polish_contigs(&mut graph, reads.into_iter(), 0);
    assert_eq!(
        graph.get_unitig_seq(0),
        contig_seq,
        "Failed preservation for len={}, offset={}, mode={}",
        len,
        offset,
        mode
    );
    assert_eq!(stats.bases_modified, 0);
    match mode {
        "fwd" => {
            assert!(stats.fwd_mapped_reads >= 10);
            assert_eq!(stats.rev_mapped_reads, 0);
        }
        "rev" => {
            assert_eq!(stats.fwd_mapped_reads, 0);
            assert!(stats.rev_mapped_reads >= 10);
        }
        "mixed" => {
            assert!(stats.fwd_mapped_reads >= 5);
            assert!(stats.rev_mapped_reads >= 5);
        }
        _ => unreachable!(),
    }
}

// 1. Minimal failure reproduction test (passes on fixed code, fails on old code)
#[test]
fn test_reverse_coordinate_repro_minimal_failure() {
    let contig_seq = make_deterministic_dna(0xABCD_1234, 400);
    let mut graph = create_test_graph(&contig_seq);

    let rev_read = make_reverse_read(&contig_seq, 150, 30, 75);
    let reads = vec![rev_read; 10];

    let stats = polish_contigs(&mut graph, reads.into_iter(), 0);
    let polished_seq = graph.get_unitig_seq(0);

    assert_eq!(
        polished_seq, contig_seq,
        "Polishing with perfect reverse reads corrupted the contig sequence due to shifted coordinates!"
    );
    assert_eq!(stats.bases_modified, 0);
    assert_eq!(stats.rev_mapped_reads, 10);
}

// 2. SNV correction test with reverse reads
#[test]
fn test_reverse_snv_correction() {
    let truth_seq = make_deterministic_dna(0x1234_5678, 400);

    let mut contig_seq = truth_seq.clone();
    let err_pos = 160;
    let true_base = truth_seq[err_pos];
    let mutated_base = match true_base {
        b'A' => b'C',
        _ => b'A',
    };
    contig_seq[err_pos] = mutated_base;

    let mut graph = create_test_graph(&contig_seq);

    // c_pos = 125, r_pos = 30, len = 100 covers pos [86..185], seed [125..156] untouched
    let rev_read = make_reverse_read(&truth_seq, 125, 30, 100);
    let reads = vec![rev_read; 10];

    let stats = polish_contigs(&mut graph, reads.into_iter(), 0);
    let polished_seq = graph.get_unitig_seq(0);

    assert_eq!(
        polished_seq[err_pos], true_base,
        "Reverse read pileup failed to correct SNV at pos {} from {} to {}",
        err_pos, mutated_base as char, true_base as char
    );
    assert_eq!(polished_seq, truth_seq, "Entire polished contig must match truth");
    assert_eq!(stats.bases_modified, 1);
    assert_eq!(stats.rev_mapped_reads, 10);
}

// 3. Matrix tests: Read length (75, 100, 150) x Seed offset (0, 30, 60) x Orientation (fwd, rev, mixed)
#[test]
fn test_matrix_len075_offset00_fwd() { run_matrix_case(75, 0, "fwd"); }
#[test]
fn test_matrix_len075_offset00_rev() { run_matrix_case(75, 0, "rev"); }
#[test]
fn test_matrix_len075_offset00_mixed() { run_matrix_case(75, 0, "mixed"); }

#[test]
fn test_matrix_len075_offset30_fwd() { run_matrix_case(75, 30, "fwd"); }
#[test]
fn test_matrix_len075_offset30_rev() { run_matrix_case(75, 30, "rev"); }
#[test]
fn test_matrix_len075_offset30_mixed() { run_matrix_case(75, 30, "mixed"); }

#[test]
fn test_matrix_len100_offset00_fwd() { run_matrix_case(100, 0, "fwd"); }
#[test]
fn test_matrix_len100_offset00_rev() { run_matrix_case(100, 0, "rev"); }
#[test]
fn test_matrix_len100_offset00_mixed() { run_matrix_case(100, 0, "mixed"); }

#[test]
fn test_matrix_len100_offset30_fwd() { run_matrix_case(100, 30, "fwd"); }
#[test]
fn test_matrix_len100_offset30_rev() { run_matrix_case(100, 30, "rev"); }
#[test]
fn test_matrix_len100_offset30_mixed() { run_matrix_case(100, 30, "mixed"); }

#[test]
fn test_matrix_len100_offset60_fwd() { run_matrix_case(100, 60, "fwd"); }
#[test]
fn test_matrix_len100_offset60_rev() { run_matrix_case(100, 60, "rev"); }
#[test]
fn test_matrix_len100_offset60_mixed() { run_matrix_case(100, 60, "mixed"); }

#[test]
fn test_matrix_len150_offset00_fwd() { run_matrix_case(150, 0, "fwd"); }
#[test]
fn test_matrix_len150_offset00_rev() { run_matrix_case(150, 0, "rev"); }
#[test]
fn test_matrix_len150_offset00_mixed() { run_matrix_case(150, 0, "mixed"); }

#[test]
fn test_matrix_len150_offset30_fwd() { run_matrix_case(150, 30, "fwd"); }
#[test]
fn test_matrix_len150_offset30_rev() { run_matrix_case(150, 30, "rev"); }
#[test]
fn test_matrix_len150_offset30_mixed() { run_matrix_case(150, 30, "mixed"); }

#[test]
fn test_matrix_len150_offset60_fwd() { run_matrix_case(150, 60, "fwd"); }
#[test]
fn test_matrix_len150_offset60_rev() { run_matrix_case(150, 60, "rev"); }
#[test]
fn test_matrix_len150_offset60_mixed() { run_matrix_case(150, 60, "mixed"); }

// 4. Sequencing errors tolerance test
#[test]
fn test_reverse_reads_with_small_sequencing_errors() {
    let contig_seq = make_deterministic_dna(0x5555_AAAA, 400);
    let mut graph = create_test_graph(&contig_seq);

    let base_read = make_reverse_read(&contig_seq, 150, 30, 100);
    let mut reads = Vec::new();
    for i in 0..10 {
        let mut r = base_read.clone();
        let mutate_idx = if i % 2 == 0 { 5 } else { 85 };
        r[mutate_idx] = match r[mutate_idx] {
            b'A' => b'C',
            _ => b'A',
        };
        reads.push(r);
    }

    let stats = polish_contigs(&mut graph, reads.into_iter(), 0);
    assert_eq!(
        graph.get_unitig_seq(0),
        contig_seq,
        "Consensus threshold (>75%) should reject sporadic sequencing errors"
    );
    assert_eq!(stats.bases_modified, 0);
    assert_eq!(stats.rev_mapped_reads, 10);
}

// 5. Repeated seed rejection test
#[test]
fn test_repeated_seed_and_unanchored_reads() {
    let mut contig_seq = make_deterministic_dna(0x1111_2222, 500);
    let seed50 = contig_seq[50..81].to_vec();
    contig_seq[250..281].copy_from_slice(&seed50);

    let mut graph = create_test_graph(&contig_seq);

    let rev_read = make_reverse_read(&contig_seq, 50, 30, 75);
    let reads = vec![rev_read; 10];

    let stats = polish_contigs(&mut graph, reads.into_iter(), 0);
    assert_eq!(graph.get_unitig_seq(0), contig_seq);
    assert_eq!(stats.bases_modified, 0);
    assert_eq!(stats.fwd_mapped_reads, 0);
    assert_eq!(stats.rev_mapped_reads, 0);
}

// 6. Unknown orientation / mismatch rejection test
#[test]
fn test_unknown_orientation_rejected() {
    let contig_seq = make_deterministic_dna(0x3333_4444, 400);
    let mut graph = create_test_graph(&contig_seq);

    // Read with unrelated sequence that does not match contig
    let bogus_read = make_deterministic_dna(0x7777_8888, 100);
    let reads = vec![bogus_read; 10];

    let stats = polish_contigs(&mut graph, reads.into_iter(), 0);
    assert_eq!(graph.get_unitig_seq(0), contig_seq);
    assert_eq!(stats.bases_modified, 0);
    assert_eq!(stats.fwd_mapped_reads, 0);
    assert_eq!(stats.rev_mapped_reads, 0);
}
