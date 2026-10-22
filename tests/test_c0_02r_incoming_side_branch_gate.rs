//! F7-A / F7-B gates for incoming side-branch preservation and overlap-adjusted
//! search distance.
//!
//! F5-P1-001: committing a path deletes every constituent node. An incoming
//! side branch X → R then loses its only active neighbor. If X is shorter than
//! 500 bp and was never inserted into `rejected_nodes`, Stage 5 cleanup deletes
//! it. `protect_unmerged_frontier` only shields outgoing neighbors of the final
//! extension node at the 64-step cap.
//!
//! F5-P2-002: `find_candidate_paths` and `depth_limit_hides_continuation` reject
//! a hop when `sep_u + len(u) > max_separation` before subtracting `overlap_k`.
//! A hop with `sep_u + len(u) - overlap <= max_separation` is still admissible.
//!
//! Cases A–E and both F7-B tests encode the corrected contract. On commit
//! `0e513a3` they fail. Cases F and G lock behavior that already holds:
//! validation failure must not mutate topology or reject an external neighbor,
//! and an ordinary unbranched path must still splice exactly.

use fastbact::common::bitseq::reverse_complement_ascii;
use fastbact::common::sketch::splitmix64;
use fastbact::graph::cdbg::{CdbgGraph, DirUnitigId, GraphEdge, UnitigNode};
use fastbact::polish::output::curate_contigs;
use fastbact::resolve::exspander::{
    resolve_edge_between, resolve_repeats_exspander_with_stats, EdgeResolveResult, ExspanderParams,
};
use fastbact::resolve::mapper::{InsertSizeDist, PairLinkMatrix, ReadPairEvidence};

fn dna(seed: u64, len: usize) -> Vec<u8> {
    // Mix seed and index before splitmix64. `splitmix64(seed + i) % 4` is only a
    // shift of a nearby seed, so a "unique" tail would be a substring of another node.
    (0..len)
        .map(|i| {
            let mixed = seed
                .wrapping_mul(0x9E37_79B9_7F4A_7C15)
                .wrapping_add((i as u64).wrapping_mul(0xBF58_476D_1CE4_E5B9));
            b"ACGT"[(splitmix64(mixed) % 4) as usize]
        })
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
    graph.nodes[to.id() as usize]
        .add_edge(to.twin().is_rc(), GraphEdge::new(from.twin(), overlap, 10));
}

fn add_forward_only(graph: &mut CdbgGraph, from: DirUnitigId, to: DirUnitigId, overlap: u8) {
    graph.nodes[from.id() as usize].add_edge(from.is_rc(), GraphEdge::new(to, overlap, 10));
}

fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    !needle.is_empty()
        && haystack
            .windows(needle.len())
            .any(|window| window == needle)
}

fn active_summary(graph: &CdbgGraph) -> String {
    graph
        .nodes
        .iter()
        .filter(|node| !node.is_deleted)
        .map(|node| format!("{}:{}", node.id, node.seq_len))
        .collect::<Vec<_>>()
        .join(",")
}

fn resolve(graph: &mut CdbgGraph, matrix: &PairLinkMatrix, dist: &InsertSizeDist) {
    resolve_repeats_exspander_with_stats(graph, matrix, dist, ExspanderParams::default()).expect("repeat resolution failed");
}

fn assert_valid_edge(graph: &CdbgGraph, from: DirUnitigId, to: DirUnitigId, overlap: u8) {
    match resolve_edge_between(graph, from, to) {
        EdgeResolveResult::Valid(edge) => {
            assert_eq!(edge.target, to);
            assert_eq!(edge.overlap_k, overlap);
        }
        other => panic!("expected a sequence-valid edge {from:?} -> {to:?}, got {other:?}"),
    }
}

fn splice(nodes: &[&[u8]], overlaps: &[usize]) -> Vec<u8> {
    assert_eq!(nodes.len(), overlaps.len() + 1);
    let mut assembled = nodes[0].to_vec();
    for (node, &overlap) in nodes.iter().skip(1).zip(overlaps) {
        assert!(overlap > 0 && overlap < node.len());
        assembled.extend_from_slice(&node[overlap..]);
    }
    assembled
}

fn assert_not_homopolymer(seq: &[u8]) {
    assert!(seq.len() > 1 && seq.iter().any(|base| *base != seq[0]));
}

fn structural_fingerprint(
    graph: &CdbgGraph,
) -> Vec<(
    u32,
    u32,
    bool,
    bool,
    Vec<GraphEdge>,
    Vec<GraphEdge>,
    Vec<u8>,
)> {
    graph
        .nodes
        .iter()
        .map(|node| {
            (
                node.id,
                node.seq_len,
                node.is_deleted,
                node.is_circular,
                node.edges_fwd.to_vec(),
                node.edges_rev.to_vec(),
                graph.get_unitig_seq(node.id),
            )
        })
        .collect()
}

fn side_branch_survives(
    graph: &CdbgGraph,
    x_id: u32,
    unique: &[u8],
    spliced: &[u8],
) -> (bool, bool, bool, bool) {
    let spliced_ok = graph
        .nodes
        .iter()
        .filter(|node| !node.is_deleted)
        .any(|node| graph.get_unitig_seq(node.id) == spliced);
    let x_alive = graph
        .nodes
        .get(x_id as usize)
        .is_some_and(|node| !node.is_deleted);
    let in_graph = graph
        .nodes
        .iter()
        .filter(|node| !node.is_deleted)
        .any(|node| contains(&graph.get_unitig_seq(node.id), unique));
    let raw = curate_contigs(graph, 0, 2.0);
    let in_raw = raw.iter().any(|contig| contains(&contig.seq, unique));
    (spliced_ok, x_alive, in_graph, in_raw)
}

/// Case A: S(700) → R(169), with incoming X(100) → R.
///
/// S → R is an unbranched terminal join and commits. X's 69 bp prefix is absent
/// from S and R. Cleanup must not delete X, and raw curation (`min_len = 0`)
/// must still emit that 69 bp.
///
/// ON UNPATCHED CODE THIS FAILS: X is orphaned when R is deleted and Stage 5
/// removes it.
#[test]
fn test_case_a_classic_incoming_side_branch_preserves_unique_sequence() {
    let mut graph = CdbgGraph::new(31);
    let junction = dna(0xA001, 31);
    let mut s_seq = dna(0xA002, 669);
    s_seq.extend_from_slice(&junction);
    let mut r_seq = junction.clone();
    r_seq.extend(dna(0xA003, 138));
    let unique = dna(0xA004, 69);
    let mut x_seq = unique.clone();
    x_seq.extend_from_slice(&junction);

    assert_eq!(s_seq.len(), 700);
    assert_eq!(r_seq.len(), 169);
    assert_eq!(x_seq.len(), 100);
    assert_eq!(&s_seq[s_seq.len() - 31..], &r_seq[..31]);
    assert_eq!(&x_seq[x_seq.len() - 31..], &r_seq[..31]);
    assert!(!contains(&s_seq, &unique));
    assert!(!contains(&r_seq, &unique));
    assert_not_homopolymer(&x_seq);

    let s_id = add_node(&mut graph, &s_seq);
    let r_id = add_node(&mut graph, &r_seq);
    let x_id = add_node(&mut graph, &x_seq);
    add_edge(&mut graph, s_id, r_id, 31);
    add_edge(&mut graph, x_id, r_id, 31);
    assert_valid_edge(&graph, s_id, r_id, 31);
    assert_valid_edge(&graph, x_id, r_id, 31);
    assert!(graph.nodes[r_id.id() as usize].edges_fwd.is_empty());

    let spliced = splice(&[&s_seq, &r_seq], &[31]);
    assert_eq!(spliced.len(), 838);
    assert!(!contains(&spliced, &unique));

    resolve(
        &mut graph,
        &PairLinkMatrix::new(),
        &InsertSizeDist::default(),
    );

    let (spliced_ok, x_alive, in_graph, in_raw) =
        side_branch_survives(&graph, x_id.id(), &unique, &spliced);
    assert!(
        spliced_ok && x_alive && in_graph && in_raw,
        "F5-P1-001 case A RED: S→R spliced={spliced_ok}, X alive={x_alive}, \
         unique 69 bp in graph={in_graph}, in raw curation={in_raw}, active=[{}]. \
         Unpatched cleanup deletes the incoming 100 bp branch after R is consumed.",
        active_summary(&graph)
    );
}

/// Case B: S(700) → R1(100) → R2(100), with X(100) → R1.
///
/// The side branch enters an intermediate consumed node, not the path sink.
#[test]
fn test_case_b_intermediate_incoming_side_branch_preserved() {
    let mut graph = CdbgGraph::new(31);
    let j0 = dna(0xB001, 31);
    let j1 = dna(0xB002, 31);
    assert_ne!(j0, j1);

    let mut s_seq = dna(0xB003, 669);
    s_seq.extend_from_slice(&j0);
    let mut r1_seq = j0.clone();
    r1_seq.extend(dna(0xB004, 38));
    r1_seq.extend_from_slice(&j1);
    let mut r2_seq = j1.clone();
    r2_seq.extend(dna(0xB005, 69));
    let unique = dna(0xB006, 69);
    let mut x_seq = unique.clone();
    x_seq.extend_from_slice(&j0);

    assert_eq!(
        (s_seq.len(), r1_seq.len(), r2_seq.len(), x_seq.len()),
        (700, 100, 100, 100)
    );
    assert!(!contains(&s_seq, &unique));
    assert!(!contains(&r1_seq, &unique));
    assert!(!contains(&r2_seq, &unique));

    let s_id = add_node(&mut graph, &s_seq);
    let r1_id = add_node(&mut graph, &r1_seq);
    let r2_id = add_node(&mut graph, &r2_seq);
    let x_id = add_node(&mut graph, &x_seq);
    add_edge(&mut graph, s_id, r1_id, 31);
    add_edge(&mut graph, r1_id, r2_id, 31);
    add_edge(&mut graph, x_id, r1_id, 31);

    let spliced = splice(&[&s_seq, &r1_seq, &r2_seq], &[31, 31]);
    assert_eq!(spliced.len(), 838);
    assert!(!contains(&spliced, &unique));

    resolve(
        &mut graph,
        &PairLinkMatrix::new(),
        &InsertSizeDist::default(),
    );

    let (spliced_ok, x_alive, in_graph, in_raw) =
        side_branch_survives(&graph, x_id.id(), &unique, &spliced);
    assert!(
        spliced_ok && x_alive && in_graph && in_raw,
        "F5-P1-001 case B RED: path spliced={spliced_ok}, X alive={x_alive}, \
         unique 69 bp in graph={in_graph}, in raw={in_raw}, active=[{}]",
        active_summary(&graph)
    );
}

/// Case C: X(100) has forward edges into two different consumed nodes.
///
/// Both R1 and R2 start with the same 31-mer, so one suffix of X is a valid
/// incoming edge to each of them. Consuming the path deletes both targets.
#[test]
fn test_case_c_multi_connection_side_branch_preserved() {
    let mut graph = CdbgGraph::new(31);
    let junction = dna(0xC001, 31);
    let mut s_seq = dna(0xC002, 669);
    s_seq.extend_from_slice(&junction);
    let mut r1_seq = junction.clone();
    r1_seq.extend(dna(0xC003, 38));
    r1_seq.extend_from_slice(&junction);
    let mut r2_seq = junction.clone();
    r2_seq.extend(dna(0xC004, 69));
    let unique = dna(0xC005, 69);
    let mut x_seq = unique.clone();
    x_seq.extend_from_slice(&junction);

    assert_eq!(r1_seq.len(), 100);
    assert_eq!(x_seq.len(), 100);
    assert_eq!(&r1_seq[..31], &r1_seq[r1_seq.len() - 31..]);
    assert!(!contains(&s_seq, &unique));
    assert!(!contains(&r1_seq, &unique));
    assert!(!contains(&r2_seq, &unique));

    let s_id = add_node(&mut graph, &s_seq);
    let r1_id = add_node(&mut graph, &r1_seq);
    let r2_id = add_node(&mut graph, &r2_seq);
    let x_id = add_node(&mut graph, &x_seq);
    add_edge(&mut graph, s_id, r1_id, 31);
    add_edge(&mut graph, r1_id, r2_id, 31);
    add_edge(&mut graph, x_id, r1_id, 31);
    add_edge(&mut graph, x_id, r2_id, 31);
    assert_eq!(graph.nodes[x_id.id() as usize].edges_fwd.len(), 2);
    assert_valid_edge(&graph, x_id, r1_id, 31);
    assert_valid_edge(&graph, x_id, r2_id, 31);

    let spliced = splice(&[&s_seq, &r1_seq, &r2_seq], &[31, 31]);
    assert!(!contains(&spliced, &unique));

    resolve(
        &mut graph,
        &PairLinkMatrix::new(),
        &InsertSizeDist::default(),
    );

    let (spliced_ok, x_alive, in_graph, in_raw) =
        side_branch_survives(&graph, x_id.id(), &unique, &spliced);
    assert!(
        spliced_ok && x_alive && in_graph && in_raw,
        "F5-P1-001 case C RED: path spliced={spliced_ok}, multi-connected X alive={x_alive}, \
         unique 69 bp in graph={in_graph}, in raw={in_raw}, active=[{}]",
        active_summary(&graph)
    );
}

/// Case D: the incoming branch is stored on the opposite strand.
///
/// The node sequence is the reverse complement of the biological branch
/// `unique 69 bp + R prefix`. The edge leaves `X.twin()` and enters R.
/// The 69 bp unique to the stored forward sequence must survive curation.
#[test]
fn test_case_d_reverse_orientation_side_branch_preserved() {
    let mut graph = CdbgGraph::new(31);
    let junction = dna(0xD001, 31);
    let mut s_seq = dna(0xD002, 669);
    s_seq.extend_from_slice(&junction);
    let mut r_seq = junction.clone();
    r_seq.extend(dna(0xD003, 138));
    assert_eq!(r_seq.len(), 169);

    let biological_unique = dna(0xD004, 69);
    let mut biological = biological_unique.clone();
    biological.extend_from_slice(&junction);
    let x_seq = reverse_complement_ascii(&biological);
    let stored_unique = x_seq[31..].to_vec();
    assert_eq!(stored_unique.len(), 69);
    assert_eq!(stored_unique, reverse_complement_ascii(&biological_unique));
    assert!(!contains(&s_seq, &stored_unique));
    assert!(!contains(&r_seq, &stored_unique));
    assert_not_homopolymer(&x_seq);

    let s_id = add_node(&mut graph, &s_seq);
    let r_id = add_node(&mut graph, &r_seq);
    let x_id = add_node(&mut graph, &x_seq);
    add_edge(&mut graph, s_id, r_id, 31);
    add_edge(&mut graph, x_id.twin(), r_id, 31);
    assert!(graph.nodes[x_id.id() as usize].edges_fwd.is_empty());
    assert_eq!(graph.nodes[x_id.id() as usize].edges_rev.len(), 1);
    assert_valid_edge(&graph, x_id.twin(), r_id, 31);
    assert_valid_edge(&graph, r_id.twin(), x_id, 31);

    let spliced = splice(&[&s_seq, &r_seq], &[31]);
    assert!(!contains(&spliced, &stored_unique));

    resolve(
        &mut graph,
        &PairLinkMatrix::new(),
        &InsertSizeDist::default(),
    );

    let (spliced_ok, x_alive, in_graph, in_raw) =
        side_branch_survives(&graph, x_id.id(), &stored_unique, &spliced);
    assert!(
        spliced_ok && x_alive && in_graph && in_raw,
        "F5-P1-001 case D RED: path spliced={spliced_ok}, RC side branch alive={x_alive}, \
         stored unique 69 bp in graph={in_graph}, in raw={in_raw}, active=[{}]",
        active_summary(&graph)
    );
}

/// Case E: 64-step cap commits S → R1 → … → R64 and protects outgoing R65.
/// Incoming X → R64 is outside that outgoing-only shield and must also survive.
///
/// ON UNPATCHED CODE THIS FAILS: R65 remains, X is deleted.
#[test]
fn test_case_e_cap_preserves_outgoing_frontier_and_incoming_branch() {
    let k = 31usize;
    let n_repeat = 65usize;
    let mut graph = CdbgGraph::new(31);
    let junctions: Vec<Vec<u8>> = (0..=n_repeat)
        .map(|i| dna(0xE100_0000 + i as u64, k))
        .collect();

    let mut s_seq = dna(0xE200_0001, 700 - k);
    s_seq.extend_from_slice(&junctions[0]);
    let mut prev = add_node(&mut graph, &s_seq);
    let mut repeat_seqs = Vec::with_capacity(n_repeat);
    for i in 1..=n_repeat {
        let mut seq = junctions[i - 1].clone();
        seq.extend(dna(0xE300_0000 + i as u64, 100 - 2 * k));
        seq.extend_from_slice(&junctions[i]);
        assert_eq!(seq.len(), 100);
        let id = add_node(&mut graph, &seq);
        add_edge(&mut graph, prev, id, k as u8);
        prev = id;
        repeat_seqs.push(seq);
    }

    // R64 is node 64 and repeat_seqs[63]. X enters its prefix.
    let r64_prefix = &repeat_seqs[63][..k];
    let unique = dna(0xE400_0001, 69);
    let mut x_seq = unique.clone();
    x_seq.extend_from_slice(r64_prefix);
    assert_eq!(x_seq.len(), 100);
    assert!(!contains(&s_seq, &unique));
    assert!(repeat_seqs.iter().all(|seq| !contains(seq, &unique)));

    let x_id = add_node(&mut graph, &x_seq);
    let r64 = DirUnitigId::new(64, false);
    let r65 = DirUnitigId::new(65, false);
    add_edge(&mut graph, x_id, r64, k as u8);
    assert_eq!(graph.nodes[64].edges_fwd.len(), 1);
    assert_eq!(graph.nodes[64].edges_fwd[0].target, r65);
    assert_valid_edge(&graph, x_id, r64, k as u8);

    let r65_unique = repeat_seqs[64][k..].to_vec();
    assert_eq!(r65_unique.len(), 69);

    resolve(
        &mut graph,
        &PairLinkMatrix::new(),
        &InsertSizeDist::default(),
    );

    let prefix_committed = graph
        .nodes
        .iter()
        .any(|node| !node.is_deleted && node.seq_len == 5116);
    let r65_alive = !graph.nodes[65].is_deleted;
    let x_alive = !graph.nodes[x_id.id() as usize].is_deleted;
    let r65_seq_ok = r65_alive && graph.get_unitig_seq(65) == repeat_seqs[64];
    let raw = curate_contigs(&graph, 0, 2.0);
    let x_in_raw = raw.iter().any(|contig| contains(&contig.seq, &unique));
    let r65_in_raw = raw.iter().any(|contig| contains(&contig.seq, &r65_unique));
    let x_in_graph = x_alive && contains(&graph.get_unitig_seq(x_id.id()), &unique);

    assert!(
        prefix_committed
            && r65_alive
            && r65_seq_ok
            && r65_in_raw
            && x_alive
            && x_in_graph
            && x_in_raw,
        "F5-P1-001 case E RED: prefix5116={prefix_committed}, R65 alive={r65_alive}, \
         R65 sequence={r65_seq_ok}, R65 in raw={r65_in_raw}, incoming X alive={x_alive}, \
         X unique in graph={x_in_graph}, X unique in raw={x_in_raw}, active=[{}]. \
         Outgoing R65 is protected at the cap; incoming X→R64 is not.",
        active_summary(&graph)
    );
}

/// Case F: a candidate that fails reciprocal validation must not change topology,
/// and must not put an external neighbor into `rejected_nodes`.
///
/// The second component is a later, shorter seed whose only supported branch
/// exit is that external neighbor. `find_candidate_paths` skips rejected node
/// ids, so a polluted neighbor would make the legal branch unresolvable.
#[test]
fn test_case_f_validation_failure_is_atomic_for_external_neighbors() {
    // Scenario 1: broken reciprocal, external neighbor, nothing else to commit.
    {
        let mut graph = CdbgGraph::new(31);
        let junction = dna(0xF001, 31);
        let mut s_seq = dna(0xF002, 769);
        s_seq.extend_from_slice(&junction);
        let mut r_seq = junction.clone();
        r_seq.extend(dna(0xF003, 138));
        let mut x_seq = dna(0xF004, 69);
        x_seq.extend_from_slice(&junction);
        assert_eq!((s_seq.len(), r_seq.len(), x_seq.len()), (800, 169, 100));

        let s_id = add_node(&mut graph, &s_seq);
        let r_id = add_node(&mut graph, &r_seq);
        let x_id = add_node(&mut graph, &x_seq);
        add_forward_only(&mut graph, s_id, r_id, 31);
        add_edge(&mut graph, x_id, r_id, 31);
        assert_valid_edge(&graph, s_id, r_id, 31);
        assert!(matches!(
            resolve_edge_between(&graph, r_id.twin(), s_id.twin()),
            EdgeResolveResult::NoEdge
        ));
        assert_valid_edge(&graph, x_id, r_id, 31);

        let before = structural_fingerprint(&graph);
        let stats = resolve_repeats_exspander_with_stats(
            &mut graph,
            &PairLinkMatrix::new(),
            &InsertSizeDist::default(),
            ExspanderParams::default(),
        ).expect("repeat resolution failed");
        assert_eq!(
            stats.valid_k_joins + stats.valid_k_minus_1_joins + stats.valid_other_joins,
            0
        );
        assert!(stats.rejected_invalid_joins >= 1);
        assert_eq!(
            structural_fingerprint(&graph),
            before,
            "failed reciprocal validation must leave nodes, edges, sequences, and deleted flags untouched"
        );
    }

    // Scenario 2: the same failure happens first, then a shorter seed resolves
    // through the external neighbor. That branch is reachable only if X was not
    // inserted into rejected_nodes.
    {
        let mut graph = CdbgGraph::new(31);
        let j_fail = dna(0xF101, 31);
        let j_qx = dna(0xF102, 31);
        let j_sq = dna(0xF103, 31);
        assert_ne!(j_fail, j_qx);
        assert_ne!(j_qx, j_sq);

        let mut s1 = dna(0xF111, 869);
        s1.extend_from_slice(&j_fail);
        let mut r1 = j_fail.clone();
        r1.extend(dna(0xF112, 138));
        let x_middle = dna(0xF113, 38);
        let mut x_seq = j_qx.clone();
        x_seq.extend_from_slice(&x_middle);
        x_seq.extend_from_slice(&j_fail);
        let mut s2 = dna(0xF114, 669);
        s2.extend_from_slice(&j_sq);
        let q_middle = dna(0xF115, 38);
        let mut q_seq = j_sq.clone();
        q_seq.extend_from_slice(&q_middle);
        q_seq.extend_from_slice(&j_qx);
        let mut d_seq = j_qx.clone();
        d_seq.extend(dna(0xF116, 69));
        // Second exit from X. Without it, the unbranched walker would append R1
        // after choosing X and consume the failed path's node. R1 is already in
        // rejected_nodes, but that set is only consulted by candidate search.
        let mut z_seq = j_fail.clone();
        z_seq.extend(dna(0xF117, 69));

        assert_eq!(s1.len(), 900);
        assert_eq!(r1.len(), 169);
        assert_eq!(x_seq.len(), 100);
        assert_eq!(s2.len(), 700);
        assert_eq!(q_seq.len(), 100);
        assert_eq!(d_seq.len(), 100);
        assert_eq!(z_seq.len(), 100);

        let s1_id = add_node(&mut graph, &s1);
        let r1_id = add_node(&mut graph, &r1);
        let x_id = add_node(&mut graph, &x_seq);
        let s2_id = add_node(&mut graph, &s2);
        let q_id = add_node(&mut graph, &q_seq);
        let d_id = add_node(&mut graph, &d_seq);
        let z_id = add_node(&mut graph, &z_seq);

        add_forward_only(&mut graph, s1_id, r1_id, 31);
        add_edge(&mut graph, x_id, r1_id, 31);
        add_edge(&mut graph, x_id, z_id, 31);
        add_edge(&mut graph, s2_id, q_id, 31);
        add_edge(&mut graph, q_id, x_id, 31);
        add_edge(&mut graph, q_id, d_id, 31);
        assert_eq!(graph.nodes[x_id.id() as usize].edges_fwd.len(), 2);
        assert!(matches!(
            resolve_edge_between(&graph, r1_id.twin(), s1_id.twin()),
            EdgeResolveResult::NoEdge
        ));
        assert_valid_edge(&graph, s2_id, q_id, 31);
        assert_valid_edge(&graph, q_id, x_id, 31);

        let mut matrix = PairLinkMatrix::new();
        for _ in 0..8 {
            matrix.add_link(s2_id, x_id);
        }
        // sep(X) = len(Q) - ov(S2→Q) - ov(Q→X) = 100 - 31 - 31 = 38.
        // add_link distances are 120 + 120, so span = 278.
        let dist = InsertSizeDist {
            median: 278.0,
            mad: 40.0,
            min_insert: 100,
            max_insert: 700,
        };

        let failed_edge = graph.nodes[s1_id.id() as usize].edges_fwd.clone();
        resolve(&mut graph, &matrix, &dist);

        assert!(
            !graph.nodes[s1_id.id() as usize].is_deleted,
            "failed seed must stay active"
        );
        assert!(
            !graph.nodes[r1_id.id() as usize].is_deleted,
            "failed path node must stay active"
        );
        assert_eq!(graph.nodes[s1_id.id() as usize].edges_fwd, failed_edge);
        assert_eq!(graph.get_unitig_seq(s1_id.id()), s1);
        assert_eq!(graph.get_unitig_seq(r1_id.id()), r1);

        let expected = splice(&[&s2, &q_seq, &x_seq], &[31, 31]);
        assert_eq!(expected.len(), 838);
        assert!(contains(&expected, &x_middle));
        let resolved = graph
            .nodes
            .iter()
            .filter(|node| !node.is_deleted)
            .any(|node| graph.get_unitig_seq(node.id) == expected);
        assert!(
            resolved,
            "external neighbor X must remain searchable after the failed path; \
             rejecting it would hide the only supported exit. active=[{}]",
            active_summary(&graph)
        );
        let forbidden = splice(&[&s1, &r1], &[31]);
        assert!(graph
            .nodes
            .iter()
            .all(|node| graph.get_unitig_seq(node.id) != forbidden));
    }
}

/// Case G: an ordinary unbranched path with no side branch still splices exactly.
#[test]
fn test_case_g_unbranched_path_splices_completely() {
    let mut graph = CdbgGraph::new(31);
    let j0 = dna(0x6701, 31);
    let j1 = dna(0x6702, 31);
    let mut s_seq = dna(0x6703, 669);
    s_seq.extend_from_slice(&j0);
    let mut r_seq = j0.clone();
    r_seq.extend(dna(0x6704, 38));
    r_seq.extend_from_slice(&j1);
    let mut t_seq = j1.clone();
    t_seq.extend(dna(0x6705, 519));
    assert_eq!((s_seq.len(), r_seq.len(), t_seq.len()), (700, 100, 550));

    let s_id = add_node(&mut graph, &s_seq);
    let r_id = add_node(&mut graph, &r_seq);
    let t_id = add_node(&mut graph, &t_seq);
    add_edge(&mut graph, s_id, r_id, 31);
    add_edge(&mut graph, r_id, t_id, 31);
    assert_valid_edge(&graph, s_id, r_id, 31);
    assert_valid_edge(&graph, r_id, t_id, 31);

    let expected = splice(&[&s_seq, &r_seq, &t_seq], &[31, 31]);
    assert_eq!(expected.len(), 1288);

    resolve(
        &mut graph,
        &PairLinkMatrix::new(),
        &InsertSizeDist::default(),
    );

    let active: Vec<_> = graph.nodes.iter().filter(|node| !node.is_deleted).collect();
    assert_eq!(
        active.len(),
        1,
        "unbranched S→R→T must collapse to one contig, active=[{}]",
        active_summary(&graph)
    );
    assert_eq!(graph.get_unitig_seq(active[0].id), expected);
    assert!(graph.nodes[s_id.id() as usize].is_deleted);
    assert!(graph.nodes[r_id.id() as usize].is_deleted);
    assert!(graph.nodes[t_id.id() as usize].is_deleted);
}

fn add_centered_pairs(
    matrix: &mut PairLinkMatrix,
    from: DirUnitigId,
    to: DirUnitigId,
    n: u32,
    dist1: u32,
    dist2: u32,
) {
    for pair_id in 0..n {
        matrix.add_evidence(
            from,
            to,
            ReadPairEvidence {
                pair_id,
                dist_u1_exit: dist1,
                dist_u2_entry: dist2,
                r1_len: 150,
                r2_len: 150,
            },
        );
    }
}

/// F7-B / line 282: the supported flank sits one hop past a node whose raw exit
/// exceeds the budget, but the overlap-adjusted exit does not.
///
/// Audit counterexample embedded in the first search hop:
///   max_separation = 600
///   sep_u          = 410
///   len(u)         = 200
///   overlap        = 127
///   sep_u_exit     = 610
///   sep_next       = 483
///   dist1 = dist2  = 50, span = 583
///
/// R branches, so the walker must use candidate search. The decoy has no
/// pair support. The admissible continuation is U → B.
///
/// ON UNPATCHED CODE THIS FAILS: the pre-check drops U → B before subtracting
/// 127 bp, the search returns no candidate, and A+R+U+B is never spliced.
#[test]
fn test_f7b_overlap_adjusted_transition_is_explored() {
    const MAX_SEP: i64 = 600;
    let mut graph = CdbgGraph::new(31);
    let j_ar = dna(0x71_0001, 31);
    let j_ru = dna(0x71_0002, 31);
    let j_ub = dna(0x71_0003, 127);
    assert_ne!(j_ar, j_ru);
    assert_ne!(j_ru, j_ub);
    assert_ne!(j_ar, j_ub);

    let mut a_seq = dna(0x71_0010, 669);
    a_seq.extend_from_slice(&j_ar);
    let mut r_seq = j_ar.clone();
    r_seq.extend(dna(0x71_0011, 410));
    r_seq.extend_from_slice(&j_ru);
    let mut u_seq = j_ru.clone();
    u_seq.extend(dna(0x71_0012, 42));
    u_seq.extend_from_slice(&j_ub);
    let b_tail = dna(0x71_0013, 293);
    let mut b_seq = j_ub.clone();
    b_seq.extend_from_slice(&b_tail);
    let mut d_seq = j_ru.clone();
    d_seq.extend(dna(0x71_0014, 49));

    assert_eq!(a_seq.len(), 700);
    assert_eq!(r_seq.len(), 472);
    assert_eq!(u_seq.len(), 200);
    assert_eq!(b_seq.len(), 420);
    assert_eq!(d_seq.len(), 80);
    assert!(
        b_seq.len() < 500,
        "B must not become a second seed and walk the path backwards"
    );

    let sep_u = r_seq.len() as i64 - 31 - 31;
    let sep_exit = sep_u + u_seq.len() as i64;
    let sep_next = sep_exit - j_ub.len() as i64;
    assert_eq!(sep_u, 410);
    assert_eq!(sep_exit, 610);
    assert_eq!(sep_next, 483);
    assert!(sep_exit > MAX_SEP);
    assert!(sep_next <= MAX_SEP);
    let span = 50 + sep_next + 50;
    assert_eq!(span, 583);

    let a_id = add_node(&mut graph, &a_seq);
    let r_id = add_node(&mut graph, &r_seq);
    let u_id = add_node(&mut graph, &u_seq);
    let b_id = add_node(&mut graph, &b_seq);
    let d_id = add_node(&mut graph, &d_seq);
    add_edge(&mut graph, a_id, r_id, 31);
    add_edge(&mut graph, r_id, u_id, 31);
    add_edge(&mut graph, r_id, d_id, 31);
    add_edge(&mut graph, u_id, b_id, 127);
    assert_eq!(graph.nodes[r_id.id() as usize].edges_fwd.len(), 2);
    assert!(graph.nodes[b_id.id() as usize].edges_fwd.is_empty());
    assert_valid_edge(&graph, u_id, b_id, 127);
    assert_valid_edge(&graph, b_id.twin(), u_id.twin(), 127);

    let expected = splice(&[&a_seq, &r_seq, &u_seq, &b_seq], &[31, 31, 127]);
    assert_eq!(expected.len(), 1603);
    assert!(contains(&expected, &b_tail));
    let partial = splice(&[&a_seq, &r_seq, &u_seq], &[31, 31]);
    assert!(!contains(&partial, &b_tail));

    let mut matrix = PairLinkMatrix::new();
    add_centered_pairs(&mut matrix, a_id, b_id, 8, 50, 50);
    let dist = InsertSizeDist {
        median: 583.0,
        mad: 50.0,
        min_insert: 250,
        max_insert: 600,
    };

    resolve(&mut graph, &matrix, &dist);

    let found = graph
        .nodes
        .iter()
        .filter(|node| !node.is_deleted)
        .any(|node| graph.get_unitig_seq(node.id) == expected);
    assert!(
        found,
        "F5-P2-002 RED: transition sep_u+len(u)=610 > 600 but sep_next=483 <= 600 \
         was not explored, so A→R→U→B (1603 bp, including B's unique tail) is absent. \
         active=[{}]",
        active_summary(&graph)
    );
}

/// F7-B / line 327: at candidate depth 4, C4 still has a sequence-valid neighbor
/// whose overlap-adjusted separation fits the budget.
///
/// Intermediates C1..C3 satisfy `sep + len <= 600` even before overlap
/// subtraction, so both the patched and unpatched searches reach C4. Only C4
/// is in the audit window (`410 + 200 = 610`, `610 - 127 = 483`). The direct
/// shortcut R → B has pair support. A complete search would commit it. A search
/// that sees the hidden continuation must refuse that shortcut and leave the
/// graph unchanged.
///
/// ON UNPATCHED CODE THIS FAILS: the depth helper returns false before
/// subtracting overlap, the search is reported complete, and the shortcut is
/// committed.
#[test]
fn test_f7b_depth_limit_overlap_adjusted_continuation_stays_hidden() {
    const MAX_SEP: i64 = 600;
    let mut graph = CdbgGraph::new(31);

    let j_ar = dna(0x72_0001, 31);
    let j_rb = dna(0x72_0002, 31);
    let j_12 = dna(0x72_0003, 31);
    let j_23 = dna(0x72_0004, 31);
    let j_34 = dna(0x72_0005, 31);
    let j_45 = dna(0x72_0006, 127);
    for (i, left) in [&j_ar, &j_rb, &j_12, &j_23, &j_34, &j_45]
        .iter()
        .enumerate()
    {
        for right in [&j_ar, &j_rb, &j_12, &j_23, &j_34, &j_45]
            .iter()
            .skip(i + 1)
        {
            assert_ne!(left, right);
            assert_ne!(*left, &reverse_complement_ascii(right));
        }
    }

    let mut a_seq = dna(0x72_0010, 669);
    a_seq.extend_from_slice(&j_ar);
    let mut r_seq = j_ar.clone();
    r_seq.extend(dna(0x72_0011, 38));
    r_seq.extend_from_slice(&j_rb);
    let mut b_seq = j_rb.clone();
    b_seq.extend(dna(0x72_0012, 369));
    let mut c1 = j_rb.clone();
    c1.extend(dna(0x72_0021, 133));
    c1.extend_from_slice(&j_12);
    let mut c2 = j_12.clone();
    c2.extend(dna(0x72_0022, 58));
    c2.extend_from_slice(&j_23);
    let mut c3 = j_23.clone();
    c3.extend(dna(0x72_0023, 88));
    c3.extend_from_slice(&j_34);
    let mut c4 = j_34.clone();
    c4.extend(dna(0x72_0024, 42));
    c4.extend_from_slice(&j_45);
    let mut c5 = j_45.clone();
    c5.extend(dna(0x72_0025, 33));

    assert_eq!(a_seq.len(), 700);
    assert_eq!(r_seq.len(), 100);
    assert_eq!(b_seq.len(), 400);
    assert_eq!(
        (c1.len(), c2.len(), c3.len(), c4.len(), c5.len()),
        (195, 120, 150, 200, 160)
    );
    assert!(graph_nodes_except_seed_are_short(&[
        r_seq.len(),
        b_seq.len(),
        c1.len(),
        c2.len(),
        c3.len(),
        c4.len(),
        c5.len(),
    ]));

    let sep_c1 = (r_seq.len() - 31 - 31) as i64;
    assert_eq!(sep_c1, 38);
    let sep_c2 = assert_enqueued(sep_c1, c1.len() as i64, 31);
    assert_eq!(sep_c2, 202);
    let sep_c3 = assert_enqueued(sep_c2, c2.len() as i64, 31);
    assert_eq!(sep_c3, 291);
    let sep_c4 = assert_enqueued(sep_c3, c3.len() as i64, 31);
    assert_eq!(sep_c4, 410);
    let sep_exit = sep_c4 + c4.len() as i64;
    let sep_next = sep_exit - 127;
    assert_eq!(sep_exit, 610);
    assert_eq!(sep_next, 483);
    assert!(sep_exit > MAX_SEP && sep_next <= MAX_SEP);

    let a_id = add_node(&mut graph, &a_seq);
    let r_id = add_node(&mut graph, &r_seq);
    let b_id = add_node(&mut graph, &b_seq);
    let c1_id = add_node(&mut graph, &c1);
    let c2_id = add_node(&mut graph, &c2);
    let c3_id = add_node(&mut graph, &c3);
    let c4_id = add_node(&mut graph, &c4);
    let c5_id = add_node(&mut graph, &c5);
    add_edge(&mut graph, a_id, r_id, 31);
    add_edge(&mut graph, r_id, b_id, 31);
    add_edge(&mut graph, r_id, c1_id, 31);
    add_edge(&mut graph, c1_id, c2_id, 31);
    add_edge(&mut graph, c2_id, c3_id, 31);
    add_edge(&mut graph, c3_id, c4_id, 31);
    add_edge(&mut graph, c4_id, c5_id, 127);
    assert_eq!(graph.nodes[r_id.id() as usize].edges_fwd.len(), 2);
    assert!(graph.nodes[b_id.id() as usize].edges_fwd.is_empty());
    assert_valid_edge(&graph, c4_id, c5_id, 127);
    assert_valid_edge(&graph, c5_id.twin(), c4_id.twin(), 127);

    let mut matrix = PairLinkMatrix::new();
    // sep(B) = 38. dist 120 + 120 ⇒ span 278, centered on the supplied median.
    add_centered_pairs(&mut matrix, a_id, b_id, 10, 120, 120);
    let dist = InsertSizeDist {
        median: 278.0,
        mad: 40.0,
        min_insert: 100,
        max_insert: 600,
    };
    let before_len = graph.nodes.len();
    let before = structural_fingerprint(&graph);

    resolve(&mut graph, &matrix, &dist);

    assert!(
        graph.nodes.len() == before_len && structural_fingerprint(&graph) == before,
        "F5-P2-002 RED: depth-4 node C4 hides C5 only after overlap subtraction \
         (610 > 600, 483 <= 600). Unpatched search treats that as complete and commits \
         shortcut R→B. nodes {} -> {}, active=[{}]",
        before_len,
        graph.nodes.len(),
        active_summary(&graph)
    );
}

fn assert_enqueued(sep_u: i64, len_u: i64, overlap: i64) -> i64 {
    let sep_exit = sep_u + len_u;
    assert!(
        sep_exit <= 600,
        "fixture error: intermediate exit {sep_exit} must stay inside the budget before overlap subtraction"
    );
    sep_exit - overlap
}

fn graph_nodes_except_seed_are_short(lengths: &[usize]) -> bool {
    lengths.iter().all(|len| *len < 500)
}
