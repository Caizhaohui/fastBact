use fastbact::common::sketch::splitmix64;
use fastbact::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, UnitigNode};
use fastbact::resolve::exspander::{
    resolve_repeats_exspander_with_stats, score_candidate_path, ExspanderParams,
};
use fastbact::resolve::mapper::{InsertSizeDist, PairLinkMatrix, ReadPairEvidence};

fn dna(seed: u64, len: usize) -> Vec<u8> {
    (0..len)
        .map(|i| b"ACGT"[(splitmix64(seed + i as u64) % 4) as usize])
        .collect()
}

fn add_node(graph: &mut CdbgGraph, seq: &[u8]) -> DirUnitigId {
    let id = graph.nodes.len() as u32;
    let seq_id = graph.arena.push_seq(seq);
    graph
        .nodes
        .push(UnitigNode::new(id, seq_id, seq.len() as u32, 25.0));
    DirUnitigId::new(id, false)
}

fn add_edge(graph: &mut CdbgGraph, from: DirUnitigId, to: DirUnitigId, overlap: u8) {
    graph.nodes[from.id() as usize].add_edge(from.is_rc(), GraphEdge::new(to, overlap, 10));
    graph.nodes[to.id() as usize].add_edge(to.twin().is_rc(), GraphEdge::new(from.twin(), overlap, 10));
}

fn link_both_strands(matrix: &mut PairLinkMatrix, from: DirUnitigId, to: DirUnitigId, n: usize) {
    for _ in 0..n {
        matrix.add_link(from, to);
        matrix.add_link(to.twin(), from.twin());
    }
}

/// Two-copy repeat R: A1,A2 -> R -> B1,B2. Only the A1/B1 pair has read-pair support.
fn two_copy_repeat(mid_len: usize) -> (CdbgGraph, PairLinkMatrix, Vec<u8>) {
    let ov_a = dna(0x100, 30);
    let ov_b = dna(0x200, 30);
    let mid = dna(0x300, mid_len);

    let pa1 = dna(0x400, 670);
    let mut a1 = pa1.clone();
    a1.extend_from_slice(&ov_a);
    let mut a2 = dna(0x500, 620);
    a2.extend_from_slice(&ov_a);
    let mut r = ov_a.clone();
    r.extend_from_slice(&mid);
    r.extend_from_slice(&ov_b);
    let mut b1 = ov_b.clone();
    b1.extend_from_slice(&dna(0x600, 570));
    let mut b2 = ov_b.clone();
    b2.extend_from_slice(&dna(0x700, 520));

    let mut expected = pa1;
    expected.extend_from_slice(&ov_a);
    expected.extend_from_slice(&mid);
    expected.extend_from_slice(&b1);

    let mut graph = CdbgGraph::new(31);
    let na1 = add_node(&mut graph, &a1);
    let na2 = add_node(&mut graph, &a2);
    let nr = add_node(&mut graph, &r);
    let nb1 = add_node(&mut graph, &b1);
    let nb2 = add_node(&mut graph, &b2);
    add_edge(&mut graph, na1, nr, 30);
    add_edge(&mut graph, na2, nr, 30);
    add_edge(&mut graph, nr, nb1, 30);
    add_edge(&mut graph, nr, nb2, 30);

    let mut matrix = PairLinkMatrix::new();
    link_both_strands(&mut matrix, na1, nb1, 10);
    (graph, matrix, expected)
}

#[test]
fn link_support_is_exact_per_directed_pair() {
    let u0 = DirUnitigId::new(0, false);
    let u1 = DirUnitigId::new(1, false);
    let mut matrix = PairLinkMatrix::new();
    for _ in 0..10 {
        matrix.add_link(u0, u1);
    }
    assert_eq!(matrix.get_support(u0, u1), 10);
    assert_eq!(matrix.get_support(u0, u1.twin()), 0);
    assert_eq!(matrix.get_support(u0.twin(), u1), 0);
    assert_eq!(matrix.get_support(u0.twin(), u1.twin()), 0);
}

#[test]
fn repeat_wider_than_insert_size_is_not_resolved_by_pair_links() {
    let (mut graph, matrix, _) = two_copy_repeat(1440);
    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    );
    assert_eq!(stats.valid_k_joins + stats.valid_k_minus_1_joins + stats.valid_other_joins, 0);
    assert_eq!(graph.nodes.len(), 5, "no merged node may be created");
    assert!(graph.nodes.iter().all(|n| !n.is_deleted));
}

#[test]
fn repeat_within_insert_size_is_resolved_with_exact_sequence() {
    let (mut graph, matrix, expected) = two_copy_repeat(240);
    resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    );
    let merged: Vec<_> = graph
        .nodes
        .iter()
        .filter(|n| !n.is_deleted && n.seq_len as usize == expected.len())
        .collect();
    assert_eq!(merged.len(), 1);
    assert_eq!(graph.get_unitig_seq(merged[0].id), expected);
}

#[test]
fn test_short_repeat_cycle_terminates_without_looping_f1_p2_003() {
    let mut graph = CdbgGraph::new(31);
    let ov = dna(0x999, 30);

    // S is a long seed unitig (600 bp >= 500 bp)
    let mut s_seq = dna(0x101, 570);
    s_seq.extend_from_slice(&ov);
    let s_id = add_node(&mut graph, &s_seq);

    // R1 is a short repeat (100 bp < 500 bp)
    let mut r1_seq = ov.clone();
    r1_seq.extend_from_slice(&dna(0x202, 40));
    r1_seq.extend_from_slice(&ov);
    let r1_id = add_node(&mut graph, &r1_seq);

    // R2 is a short repeat (100 bp < 500 bp)
    let mut r2_seq = ov.clone();
    r2_seq.extend_from_slice(&dna(0x303, 40));
    r2_seq.extend_from_slice(&ov);
    let r2_id = add_node(&mut graph, &r2_seq);

    // Edges forming S -> R1 -> R2 -> R1 (unbranched cycle of short nodes)
    add_edge(&mut graph, s_id, r1_id, 30);
    add_edge(&mut graph, r1_id, r2_id, 30);
    add_edge(&mut graph, r2_id, r1_id, 30);

    let matrix = PairLinkMatrix::new();
    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    );

    // The repeat cycle must terminate deterministically without infinite loop!
    // Path validation may reject splicing or terminate cleanly.
    // Graph must remain in a valid state.
    assert!(
        stats.valid_k_joins + stats.valid_k_minus_1_joins + stats.valid_other_joins <= 2,
        "Short repeat cycle must not produce infinite joins"
    );
}

/// Case 8: 64-step prefix must not orphan frontier R65 (F3-P1-002).
///
/// Chain S(700 bp) -> R1(100 bp) -> ... -> R65(100 bp), overlap 31 bp, reciprocal
/// edges, empty pair matrix. The unbranched walk stops when `ext.len()` reaches
/// MAX_EXTENSION_STEPS (64), so the committed prefix is S+R1..=R64 (spliced
/// length 5116 bp = 700 + 64*69). R65 is not on that path and is not inserted
/// into `rejected_nodes`. Its only neighbor R64 is deleted, and cleanup removes
/// the isolated 100 bp node. The 31 bp overlap with R64 is already in the prefix;
/// the remaining 69 bp of R65 is lost.
///
/// Oracle: R65's 69 bp tail is present in an active unitig, and the sum of active
/// unitig lengths is at least the overlap-collapsed genome 700 + 65*69 = 5185 bp.
/// Keeping R65 beside the prefix, merging the whole chain, or merging nothing all
/// satisfy this. Prefix commit plus orphan deletion does not.
///
/// ON UNPATCHED PRODUCTION CODE: This test FAILS (RED). Surviving bases are 5116
/// and the R65 tail is absent.
#[test]
fn test_case_8_bounded_prefix_frontier_sequence_retained_f3_p1_002() {
    let k = 31usize;
    let n_repeat = 65usize;
    let mut graph = CdbgGraph::new(k);

    let junctions: Vec<Vec<u8>> = (0..=n_repeat).map(|i| dna(0x6100_0000 + i as u64, k)).collect();
    for pair in junctions.windows(2) {
        assert_ne!(pair[0], pair[1], "junction k-mers must be distinct");
    }

    let mut s_seq = dna(0x5300_0001, 700 - k);
    s_seq.extend_from_slice(&junctions[0]);
    assert_eq!(s_seq.len(), 700);
    let mut prev = add_node(&mut graph, &s_seq);

    let mut repeat_seqs = Vec::with_capacity(n_repeat);
    for i in 1..=n_repeat {
        let mut seq = junctions[i - 1].clone();
        seq.extend(dna(0x7200_0000 + i as u64, 100 - 2 * k));
        seq.extend_from_slice(&junctions[i]);
        assert_eq!(seq.len(), 100);
        let id = add_node(&mut graph, &seq);
        add_edge(&mut graph, prev, id, k as u8);
        prev = id;
        repeat_seqs.push(seq);
    }

    assert_eq!(graph.nodes.len(), 1 + n_repeat);
    assert!(graph.nodes[1..].iter().all(|node| node.seq_len == 100));
    assert_eq!(&s_seq[s_seq.len() - k..], &repeat_seqs[0][..k]);
    for pair in repeat_seqs.windows(2) {
        assert_eq!(&pair[0][pair[0].len() - k..], &pair[1][..k]);
    }

    let r65_tail = repeat_seqs[n_repeat - 1][k..].to_vec();
    assert_eq!(r65_tail.len(), 69);
    let expected_unique = 700 + n_repeat * (100 - k);
    assert_eq!(expected_unique, 5185);

    resolve_repeats_exspander_with_stats(
        &mut graph,
        &PairLinkMatrix::new(),
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    );

    let tail_retained = graph.nodes.iter().filter(|node| !node.is_deleted).any(|node| {
        let seq = graph.get_unitig_seq(node.id);
        seq.windows(r65_tail.len()).any(|window| window == r65_tail.as_slice())
    });
    let surviving: usize = graph
        .nodes
        .iter()
        .filter(|node| !node.is_deleted)
        .map(|node| node.seq_len as usize)
        .sum();

    assert!(
        tail_retained && surviving >= expected_unique,
        "REGRESSION FAILURE (F3-P1-002): frontier R65 unique tail retained = {}, \
         surviving active bases = {}, authentic genome lower bound = {} (700 + 65*69). \
         Unpatched code commits the 64-step prefix (5116 bp) and cleanup deletes isolated R65, \
         losing 69 bp.",
        tail_retained,
        surviving,
        expected_unique
    );
}

/// Case 9: Duplicate `pair_id` rows count once in support and weight (F3-P2-003).
///
/// Nine `ReadPairEvidence` rows are stored, interleaved:
///   pair_id 1, four copies, span 200 + 300 + 200 = 700 bp (Z = 0, weight 1)
///   pair_id 2, two copies, same span
///   pair_id 3, three copies, span 4300 bp, outside [350, 1050]
/// `score_candidate_path` must count the two in-range pair ids once each.
///
/// ON UNPATCHED PRODUCTION CODE: This test FAILS (RED) with support 6 and weight 6.
#[test]
fn test_case_9_duplicate_pair_id_counted_once_f3_p2_003() {
    let dist = InsertSizeDist {
        median: 700.0,
        mad: 100.0,
        min_insert: 350,
        max_insert: 1050,
    };
    let u1 = DirUnitigId::new(1, false);
    let u2 = DirUnitigId::new(2, false);
    let path_separation = 300i64;

    let rows = [
        (1u32, 200u32, 200u32),
        (2, 200, 200),
        (3, 2_000, 2_000),
        (1, 200, 200),
        (1, 200, 200),
        (2, 200, 200),
        (3, 2_000, 2_000),
        (1, 200, 200),
        (3, 2_000, 2_000),
    ];

    let mut matrix = PairLinkMatrix::new();
    for (pair_id, dist_u1_exit, dist_u2_entry) in rows {
        matrix.add_evidence(
            u1,
            u2,
            ReadPairEvidence {
                pair_id,
                dist_u1_exit,
                dist_u2_entry,
                r1_len: 150,
                r2_len: 150,
            },
        );
    }

    let stored = matrix.get_evidence(u1, u2).expect("duplicate rows stored");
    assert_eq!(
        stored.len(),
        rows.len(),
        "precondition: PairLinkMatrix keeps every inserted row"
    );

    let (support, weight) = score_candidate_path(&matrix, &dist, u1, u2, path_separation);
    assert!(
        support == 2 && (weight - 2.0).abs() < 1e-9,
        "REGRESSION FAILURE (F3-P2-003): each pair_id must contribute once. \
         expected support 2 and weight 2.0, got support {} and weight {}. \
         Unpatched scoring counts all 6 in-range copies of pair_id 1 and 2.",
        support,
        weight
    );
}
