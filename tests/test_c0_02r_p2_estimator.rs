use fastbact::common::bitseq::{reverse_complement_ascii, ReadBank};
use fastbact::common::sketch::splitmix64;
use fastbact::graph::cdbg::{CdbgGraph, UnitigNode};
use fastbact::resolve::mapper::{index_terminal_kmers, map_paired_reads};

fn make_deterministic_dna(seed: u64, len: usize) -> Vec<u8> {
    let alphabet = [b'A', b'C', b'G', b'T'];
    (0..len)
        .map(|i| alphabet[(splitmix64(seed.wrapping_add(i as u64)) % 4) as usize])
        .collect()
}

fn setup_single_unitig_graph(len: usize) -> (CdbgGraph, Vec<u8>) {
    let mut graph = CdbgGraph::new(31);
    let seq = make_deterministic_dna(0xCAFE_BABE_0001, len);
    let seq_id = graph.arena.push_seq(&seq);
    graph
        .nodes
        .push(UnitigNode::new(0, seq_id, seq.len() as u32, 30.0));
    (graph, seq)
}

/// Helper to build a pair of reads from a fragment on the unitig forward sequence.
/// frag_start: start coordinate on unitig
/// frag_len: physical fragment span
/// r1_len, r2_len: read lengths
/// r1_offset: offset of non-mutated seed in r1 (0, 20, 40, etc.)
/// r2_offset: offset of non-mutated seed in r2
/// If r1_fwd_r2_rev is true: R1 is forward strand (frag_start..frag_start+r1_len),
///   R2 is reverse strand (RC of frag_end-r2_len..frag_end).
/// If false: R1 is reverse strand, R2 is forward strand.
fn make_fragment_pair(
    unitig_seq: &[u8],
    frag_start: usize,
    frag_len: usize,
    r1_len: usize,
    r2_len: usize,
    r1_seed_step: usize,
    r2_seed_step: usize,
    r1_fwd_r2_rev: bool,
) -> (Vec<u8>, Vec<u8>) {
    let frag_end = frag_start + frag_len;
    assert!(frag_end <= unitig_seq.len());

    let (mut r1, mut r2) = if r1_fwd_r2_rev {
        let r1 = unitig_seq[frag_start..frag_start + r1_len].to_vec();
        let r2 = reverse_complement_ascii(&unitig_seq[frag_end - r2_len..frag_end]);
        (r1, r2)
    } else {
        let r1 = reverse_complement_ascii(&unitig_seq[frag_end - r1_len..frag_end]);
        let r2 = unitig_seq[frag_start..frag_start + r2_len].to_vec();
        (r1, r2)
    };

    // To force map_read to match at r1_seed_step, mutate the k-mers at steps < r1_seed_step
    // Note: step in [0, 20, 40, 60, 80, 100].
    // If r1_seed_step > 0, we can alter earlier step k-mers so they don't match the index.
    for step in [0, 20, 40, 60, 80, 100] {
        if step < r1_seed_step && step + 31 <= r1.len() {
            // mutate base at step to 'N' or mismatch
            r1[step] = if r1[step] == b'A' { b'T' } else { b'A' };
        }
        if step < r2_seed_step && step + 31 <= r2.len() {
            r2[step] = if r2[step] == b'A' { b'T' } else { b'A' };
        }
    }

    (r1, r2)
}

#[test]
fn test_insert_size_estimation_seed_offset_zero() {
    let (graph, unitig) = setup_single_unitig_graph(2000);
    let mut bank = ReadBank::new();

    // 25 identical pairs with exact true fragment length = 450, r1_len=150, r2_len=150, seed_step=0
    for i in 0..25 {
        let start = 100 + i * 2;
        let (r1, r2) = make_fragment_pair(&unitig, start, 450, 150, 150, 0, 0, true);
        bank.pe_forward.push_seq(&r1);
        bank.pe_reverse.push_seq(&r2);
    }

    let (_matrix, dist) = map_paired_reads::<u64>(&graph, &bank, 31);
    // True span is 450.
    // Old implementation computes: |pos2 - pos1| + r1.len()
    // Here pos1 = start. pos2 = frag_end - 31 (since r2 is RC and seed_step=0, r2[0..31] matches unitig[frag_end-31..frag_end]).
    // pos2 - pos1 = (start + 450 - 31) - start = 419.
    // d = 419 + 150 = 569! (an error of +119 bp!)
    assert_eq!(
        dist.median.round() as usize,
        450,
        "Median insert size must match true fragment span 450, got {}",
        dist.median
    );
}

#[test]
fn test_insert_size_estimation_different_seed_offsets() {
    let (graph, unitig) = setup_single_unitig_graph(2000);
    let mut bank = ReadBank::new();

    // 25 pairs with true span = 500, but R1 matches at seed_step=40, R2 matches at seed_step=20
    for i in 0..25 {
        let start = 100 + i * 2;
        let (r1, r2) = make_fragment_pair(&unitig, start, 500, 150, 150, 40, 20, true);
        bank.pe_forward.push_seq(&r1);
        bank.pe_reverse.push_seq(&r2);
    }

    let (_matrix, dist) = map_paired_reads::<u64>(&graph, &bank, 31);
    assert_eq!(
        dist.median.round() as usize,
        500,
        "Median insert size must be robust to seed offset, expected 500, got {}",
        dist.median
    );
}

#[test]
fn test_insert_size_reverse_complement_pair() {
    let (graph, unitig) = setup_single_unitig_graph(2000);
    let mut bank = ReadBank::new();

    // 25 pairs with R1 on reverse strand, R2 on forward strand (r1_fwd_r2_rev = false)
    for i in 0..25 {
        let start = 100 + i * 2;
        let (r1, r2) = make_fragment_pair(&unitig, start, 480, 150, 150, 0, 0, false);
        bank.pe_forward.push_seq(&r1);
        bank.pe_reverse.push_seq(&r2);
    }

    let (_matrix, dist) = map_paired_reads::<u64>(&graph, &bank, 31);
    assert_eq!(
        dist.median.round() as usize,
        480,
        "Reverse-orientation pair must correctly yield true fragment span 480, got {}",
        dist.median
    );
}

#[test]
fn test_outward_and_same_strand_pairs_rejected() {
    let (graph, unitig) = setup_single_unitig_graph(2000);
    let mut bank = ReadBank::new();

    // 25 pairs on the SAME strand (both forward)
    for i in 0..25 {
        let start = 100 + i * 2;
        let r1 = unitig[start..start + 150].to_vec();
        let r2 = unitig[start + 200..start + 350].to_vec();
        bank.pe_forward.push_seq(&r1);
        bank.pe_reverse.push_seq(&r2);
    }

    let (_matrix, dist) = map_paired_reads::<u64>(&graph, &bank, 31);
    // Because same-strand pairs are not concordant inward pairs, they must be rejected.
    // Fewer than 20 qualifying pairs -> default fallback!
    assert_eq!(
        dist,
        fastbact::resolve::mapper::InsertSizeDist::default(),
        "Same-strand pairs must not contribute to insert size estimation"
    );
}

#[test]
fn test_different_read_lengths() {
    let (graph, unitig) = setup_single_unitig_graph(2000);
    let mut bank = ReadBank::new();

    // 25 pairs with R1=150 bp, R2=100 bp, true span = 400
    for i in 0..25 {
        let start = 100 + i * 2;
        let (r1, r2) = make_fragment_pair(&unitig, start, 400, 150, 100, 0, 0, true);
        bank.pe_forward.push_seq(&r1);
        bank.pe_reverse.push_seq(&r2);
    }

    let (_matrix, dist) = map_paired_reads::<u64>(&graph, &bank, 31);
    assert_eq!(
        dist.median.round() as usize,
        400,
        "Different read lengths must accurately yield span 400, got {}",
        dist.median
    );
}

#[test]
fn test_terminal_kmer_index_overlap_fix_f1_p2_002() {
    let k = 31;

    // Test cases: (unitig_length, expected_indexed_kmers)
    // For len <= 2000 + k - 1 (total_kmers <= 2000): all total_kmers must be indexed.
    // For total_kmers > 2000: exactly 2000 terminal k-mers (1000 from 5', 1000 from 3') must be indexed.
    let test_cases = [
        (100, 70),       // Short unitig: 100 - 31 + 1 = 70
        (1031, 1001),    // Boundary: 1000 + 31 = 1031 -> 1001 kmers
        (1032, 1002),    // Catastrophic overlap boundary in unpatched code!
        (1500, 1470),    // Mid-range unitig: 1500 - 31 + 1 = 1470 kmers
        (2000, 1970),    // 2000 - 31 + 1 = 1970 kmers
        (2500, 2000),    // Long unitig: 1000 from 5' + 1000 from 3' = 2000 kmers
    ];

    for &(len, expected_kmers) in &test_cases {
        let (graph, _seq) = setup_single_unitig_graph(len);
        let index = index_terminal_kmers::<u64>(&graph, k);
        assert_eq!(
            index.len(),
            expected_kmers,
            "Failed for unitig len {}: expected {} unique indexed k-mers, got {}",
            len,
            expected_kmers,
            index.len()
        );
    }
}
