//! Phase F23 horizon guard — repeat-join withholding gate.
//!
//! Forensic target: R1-001 `contig_16` chimeric relocation (runtime decision
//! `DEC_000022`, evidence chain F20 -> F21-A -> F21-B -> F22-D). While
//! extending through an unbranched chain of short (< 500 bp) repeat nodes, the
//! search anchor degrades onto collapsed shared repeat sequence. When the walk
//! stands beyond the last long (>= 500 bp) anchor's insert horizon at a
//! branching repeat, unique-anchored evidence is arithmetically impossible and
//! the only remaining links are repeat-internal ones, which cannot distinguish
//! which genomic copy they originate from: ectopic pairs from the homologous
//! copy then out-vote the authentic exit and a chimeric join across megabases
//! is committed.
//!
//! Contract (F23 horizon guard, narrowed per Codex root-cause review
//! `codex_f23_root_cause_analysis.md`):
//!   1. A branching fork whose overlap-adjusted anchor-to-fork distance
//!      exceeds `max_insert` is withheld (`withheld_no_unique_anchor`), not
//!      decided on repeat-internal evidence.
//!   2. Within the horizon, forks keep the pre-F23 `curr_dir`-anchored
//!      semantics unchanged — the accepted T26/T27/T30/T34 contracts. The
//!      known hazard of a genuinely multi-copy short anchor inside the
//!      horizon is a documented scope boundary, not covered by this guard.
//!
//! Topology (both scenarios, overlap 30, k = 31), node ids in comment:
//!
//!   S(600) -> C1 -> C2 -> C3 -> F -> { W(600), T(600) }
//!   0        1     2     3     4     5        6
//!
//! W is the ectopic exit whose link evidence originates from the homologous
//! copy; T the authentic exit. The repeat-internal rows are keyed by C3 —
//! where the unguarded code's scoring anchor stands — and cannot distinguish
//! genomic copies.
//!
//! ON UNPATCHED PRODUCTION CODE the horizon scenario FAILS (RED, jobs
//! 3010464/3010541): the ectopic exit wins 16 : 6 with weighted ratio
//! 2.67 >= decision_threshold 2.5 and the chimeric join S..F..W is committed.
//!
//! Additional strong-oracle coverage (post-audit tests on base
//! `f121920a9f3d76a774e8a382f2f3047f0c0c631e`):
//!   - inclusive edge `fork_dist == max_insert` admitted; `max_insert + 1` withheld
//!   - Case A distance must include `target_node.seq_len - overlap_k`
//!   - mirrored reverse-corridor symmetry of the beyond-horizon withhold
//!   - length≥500 collapsed/non-unique intermediate must not unsafe-reset the
//!     unique anchor (correctness contract; expected RED if production resets)
//!   - genuinely unique short chain beyond the horizon is still withheld

use fastbact::common::bitseq::reverse_complement_ascii;
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

/// Per-node input slots for `dna`, 4096 apart. `dna(s, len)` draws
/// `splitmix64(s + i)`, so sequences whose seed ranges overlap are shifts of
/// one another; every fixture sequence gets a disjoint 4096-wide slot
/// (all lengths are < 4096) to keep tails genuinely unique.
const fn slot(id: u64) -> u64 {
    id * 4096
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
    graph.nodes[to.id() as usize].add_edge(
        to.twin().is_rc(),
        GraphEdge::new(from.twin(), overlap, 10),
    );
}

fn add_rows(
    matrix: &mut PairLinkMatrix,
    u1: DirUnitigId,
    u2: DirUnitigId,
    pair_ids: std::ops::Range<u32>,
    dist_u1_exit: u32,
    dist_u2_entry: u32,
) {
    for pair_id in pair_ids {
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
}

fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    haystack.windows(needle.len()).any(|w| w == needle)
}

/// Node ids of the fixture, by creation order.
const C3: u32 = 3;
const W: u32 = 5;
const T: u32 = 6;

struct Fixture {
    graph: CdbgGraph,
    matrix: PairLinkMatrix,
    s_tail: Vec<u8>,
    w_tail: Vec<u8>,
    t_tail: Vec<u8>,
}

/// Builds S(600) -> C1 -> C2 -> C3 -> F -> {W(600), T(600)} with 30 bp
/// overlaps; `chain_len` / `fork_len` control the walked chain geometry.
fn build_fixture(chain_len: usize, fork_len: usize) -> Fixture {
    // Slots: 1000/1001 junctions, 0 = S, 1..=3 chain, 4 = F, 5 = W, 6 = T.
    let junction: Vec<u8> = dna(slot(1000), 30);
    let junction_b: Vec<u8> = dna(slot(1001), 30);

    let mut s_seq = dna(slot(0), 600 - 30);
    s_seq.extend_from_slice(&junction);

    let mut c_seqs: Vec<Vec<u8>> = Vec::new();
    for i in 0..3u64 {
        let mut seq = junction.clone();
        seq.extend(dna(slot(1 + i), chain_len - 60));
        seq.extend_from_slice(&junction);
        c_seqs.push(seq);
    }

    let mut f_seq = junction.clone();
    f_seq.extend(dna(slot(4), fork_len - 60));
    f_seq.extend_from_slice(&junction_b);

    let mut w_seq = junction_b.clone();
    w_seq.extend(dna(slot(5), 600 - 30));
    let mut t_seq = junction_b.clone();
    t_seq.extend(dna(slot(6), 600 - 30));

    let mut graph = CdbgGraph::new(31);
    let s = add_node(&mut graph, &s_seq);
    let mut prev = s;
    for seq in &c_seqs {
        let id = add_node(&mut graph, seq);
        add_edge(&mut graph, prev, id, 30);
        prev = id;
    }
    let f = add_node(&mut graph, &f_seq);
    add_edge(&mut graph, prev, f, 30);
    let w = add_node(&mut graph, &w_seq);
    let t = add_node(&mut graph, &t_seq);
    add_edge(&mut graph, f, w, 30);
    add_edge(&mut graph, f, t, 30);

    Fixture {
        graph,
        matrix: PairLinkMatrix::new(),
        s_tail: s_seq[100..220].to_vec(),
        w_tail: w_seq[200..320].to_vec(),
        t_tail: t_seq[200..320].to_vec(),
    }
}

fn active_sequences(graph: &CdbgGraph) -> Vec<Vec<u8>> {
    graph
        .nodes
        .iter()
        .filter(|node| !node.is_deleted)
        .map(|node| graph.get_unitig_seq(node.id))
        .collect()
}

fn u(id: u32) -> DirUnitigId {
    DirUnitigId::new(id, false)
}

/// Horizon scenario — the DEC_000022 geometry: the anchor-to-fork distance is
/// 3 x (300 - 30) = 810 bp on the walked chain, plus the fork hop 200 - 30 =
/// 170, i.e. 980 bp beyond the 700 bp insert horizon. Repeat-internal rows
/// anchor the ectopic exit W with 16 pairs against the authentic exit T with
/// 6; unguarded code commits W (chimera S..F..W). The guard must withhold the
/// fork and conserve both exits.
#[test]
fn fork_beyond_unique_anchor_horizon_withholds_instead_of_chimera() {
    let mut fx = build_fixture(300, 200);

    // Repeat-internal rows (read1 inside the shared chain node C3): exactly
    // the evidence class that carries no locus specificity. 16 : 6 with
    // weighted ratio 2.67 would commit W unguarded.
    add_rows(&mut fx.matrix, u(C3), u(W), 0..16, 120, 120);
    add_rows(&mut fx.matrix, u(C3), u(T), 16..22, 120, 120);

    let stats = resolve_repeats_exspander_with_stats(
        &mut fx.graph,
        &fx.matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    assert_eq!(
        stats.withheld_no_unique_anchor, 1,
        "the beyond-horizon fork must be withheld by the horizon guard"
    );

    let seqs = active_sequences(&fx.graph);
    let chimera_sw = seqs
        .iter()
        .any(|seq| contains(seq, &fx.s_tail) && contains(seq, &fx.w_tail));
    let chimera_st = seqs
        .iter()
        .any(|seq| contains(seq, &fx.s_tail) && contains(seq, &fx.t_tail));
    let t_retained = seqs.iter().any(|seq| contains(seq, &fx.t_tail));
    let w_retained = seqs.iter().any(|seq| contains(seq, &fx.w_tail));

    assert!(
        !chimera_sw && !chimera_st,
        "REGRESSION FAILURE (F23 horizon guard): a contig splices the unique \
         anchor region onto a fork exit across the shared repeat chain (S+W: \
         {chimera_sw}, S+T: {chimera_st}). Unguarded code out-votes the \
         authentic exit with ectopic repeat-internal pairs (16:6, weighted \
         ratio 2.67 >= 2.5) — the DEC_000022 chimeric join mechanism. A fork \
         beyond the unique anchor's insert horizon must be withheld.",
    );
    assert!(
        t_retained && w_retained,
        "withholding must conserve both exits: T retained = {t_retained}, \
         W retained = {w_retained}",
    );
}

/// Within-horizon passthrough control — the guard must be inert when the fork
/// stays inside the anchor's horizon: chain 3 x (100 - 30) = 210 plus the
/// fork hop 150 - 30 = 120 gives 330 <= 700. The accepted pre-F23 contract
/// (T26/T27/T30/T34 family) requires the fork to be decided on repeat-
/// internal, curr_dir-anchored evidence in that regime; the documented
/// multi-copy-anchor hazard inside the horizon is a scope boundary of this
/// guard, not a regression.
#[test]
fn fork_within_unique_anchor_horizon_keeps_pre_f23_semantics() {
    let mut fx = build_fixture(100, 150);

    add_rows(&mut fx.matrix, u(C3), u(W), 0..16, 100, 60);
    add_rows(&mut fx.matrix, u(C3), u(T), 16..22, 100, 60);

    let stats = resolve_repeats_exspander_with_stats(
        &mut fx.graph,
        &fx.matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    assert_eq!(
        stats.withheld_no_unique_anchor, 0,
        "the within-horizon fork must not be intercepted by the horizon guard"
    );

    let seqs = active_sequences(&fx.graph);
    let bridged = seqs
        .iter()
        .any(|seq| contains(seq, &fx.s_tail) && contains(seq, &fx.w_tail));
    let t_retained = seqs.iter().any(|seq| contains(seq, &fx.t_tail));

    assert!(
        bridged,
        "REGRESSION FAILURE (F23 horizon guard over-reach): a within-horizon \
         fork was withheld or fragmented. Pre-F23 semantics commit the \
         repeat-internal 16:6 bridge at fork distance 330 <= 700; the guard \
         must stay inert inside the horizon (accepted T-series contract).",
    );
    assert!(
        t_retained,
        "the non-selected exit must be conserved by the commit path"
    );
}

fn alive(graph: &CdbgGraph, id: u32) -> bool {
    graph
        .nodes
        .get(id as usize)
        .is_some_and(|node| !node.is_deleted)
}

fn u_rc(id: u32) -> DirUnitigId {
    DirUnitigId::new(id, true)
}

/// Case A arithmetic used by the horizon gate and by these boundary fixtures:
/// `fork_dist = 3 * (chain_len - 30) + fork_len - 30` with default
/// `max_insert = 700`. Equality is admitted (`>`), so chain=200 / fork=220
/// yields exactly 700 and must pass; fork=221 yields 701 and must withhold.
fn assert_case_a_fork_dist(chain_len: usize, fork_len: usize, expected: i64) {
    let gap = 3 * (chain_len as i64 - 30);
    let fork_dist = gap + fork_len as i64 - 30;
    assert_eq!(
        fork_dist, expected,
        "fixture geometry drifted: chain={chain_len} fork={fork_len} => {fork_dist}, want {expected}"
    );
}

/// Mirrored reverse-corridor twin of `build_fixture`: sequences are stored as
/// reverse complements and edges are attached on the twin strand so the seed
/// has an empty forward adjacency and must extend exclusively via `S.twin()`.
fn build_fixture_mirrored(chain_len: usize, fork_len: usize) -> Fixture {
    let junction: Vec<u8> = dna(slot(1000), 30);
    let junction_b: Vec<u8> = dna(slot(1001), 30);

    let mut s_bio = dna(slot(0), 600 - 30);
    s_bio.extend_from_slice(&junction);

    let mut c_bios: Vec<Vec<u8>> = Vec::new();
    for i in 0..3u64 {
        let mut seq = junction.clone();
        seq.extend(dna(slot(1 + i), chain_len - 60));
        seq.extend_from_slice(&junction);
        c_bios.push(seq);
    }

    let mut f_bio = junction.clone();
    f_bio.extend(dna(slot(4), fork_len - 60));
    f_bio.extend_from_slice(&junction_b);

    let mut w_bio = junction_b.clone();
    w_bio.extend(dna(slot(5), 600 - 30));
    let mut t_bio = junction_b.clone();
    t_bio.extend(dna(slot(6), 600 - 30));

    let s_seq = reverse_complement_ascii(&s_bio);
    let c_seqs: Vec<Vec<u8>> = c_bios.iter().map(|s| reverse_complement_ascii(s)).collect();
    let f_seq = reverse_complement_ascii(&f_bio);
    let w_seq = reverse_complement_ascii(&w_bio);
    let t_seq = reverse_complement_ascii(&t_bio);

    let mut graph = CdbgGraph::new(31);
    let s = add_node(&mut graph, &s_seq);
    let mut prev = s;
    for seq in &c_seqs {
        let id = add_node(&mut graph, seq);
        add_edge(&mut graph, prev.twin(), id.twin(), 30);
        prev = id;
    }
    let f = add_node(&mut graph, &f_seq);
    add_edge(&mut graph, prev.twin(), f.twin(), 30);
    let w = add_node(&mut graph, &w_seq);
    let t = add_node(&mut graph, &t_seq);
    add_edge(&mut graph, f.twin(), w.twin(), 30);
    add_edge(&mut graph, f.twin(), t.twin(), 30);

    assert!(
        graph.nodes[s.id() as usize].edges_fwd.is_empty(),
        "mirrored fixture must force reverse-corridor extension from S.twin()"
    );

    Fixture {
        graph,
        matrix: PairLinkMatrix::new(),
        s_tail: s_seq[100..220].to_vec(),
        w_tail: w_seq[200..320].to_vec(),
        t_tail: t_seq[200..320].to_vec(),
    }
}

/// Shared topology for the long-collapsed reset hazard and the unique-short
/// beyond-horizon withhold: S_A --C1--C2--C3--> [optional L] --F--> {W,T},
/// with optional second flank S_B into L (collapsed multi-copy evidence).
struct CollapseFixture {
    graph: CdbgGraph,
    matrix: PairLinkMatrix,
    sa_tail: Vec<u8>,
    w_tail: Vec<u8>,
    t_tail: Vec<u8>,
    /// Interior motif of F; survives independent seed merges that consume the F node.
    f_tail: Vec<u8>,
    /// Present only when the collapsed long node is inserted.
    sb_tail: Option<Vec<u8>>,
    l_body: Option<Vec<u8>>,
    short_bodies: Vec<Vec<u8>>,
    truth_copy_a: Option<Vec<u8>>,
    truth_copy_b: Option<Vec<u8>>,
    evidence_id: u32,
    w_id: u32,
    t_id: u32,
    long_id: Option<u32>,
    sb_id: Option<u32>,
}

fn build_collapse_or_short_horizon(insert_collapsed_long: bool) -> CollapseFixture {
    // Chain contribution before the optional long node: 3 * (300 - 30) = 810 > 700.
    // With L present and no anchor reset, Case A distance through F is
    // 810 + (600 - 30) + (150 - 30) = 1500 > 700 (withhold). Without optional
    // L, the short-fixture distance is 810 + (150 - 30) = 930 > 700. With an
    // unsafe ≥500 reset onto L, the collapsed fixture distance is 150 - 30 = 120.
    let chain_len = 300usize;
    let fork_len = 150usize;
    let junction = dna(slot(2000), 30);
    let junction_b = dna(slot(2001), 30);

    let mut sa_seq = dna(slot(10), 700 - 30);
    sa_seq.extend_from_slice(&junction);

    let mut short_bodies = Vec::new();
    let mut c_seqs = Vec::new();
    for i in 0..3u64 {
        let body = dna(slot(11 + i), chain_len - 60);
        short_bodies.push(body.clone());
        let mut seq = junction.clone();
        seq.extend_from_slice(&body);
        seq.extend_from_slice(&junction);
        c_seqs.push(seq);
    }

    let mut l_seq = junction.clone();
    l_seq.extend(dna(slot(20), 600 - 60));
    l_seq.extend_from_slice(&junction_b);

    // F's left junction matches its predecessor: L (junction_b) or C3 (junction).
    let f_left = if insert_collapsed_long {
        junction_b.as_slice()
    } else {
        junction.as_slice()
    };
    let mut f_seq = f_left.to_vec();
    f_seq.extend(dna(slot(21), fork_len - 60));
    f_seq.extend_from_slice(&junction_b);

    let mut w_seq = junction_b.clone();
    w_seq.extend(dna(slot(22), 600 - 30));
    let mut t_seq = junction_b.clone();
    t_seq.extend(dna(slot(23), 600 - 30));

    let mut sb_seq = dna(slot(30), 600 - 30);
    sb_seq.extend_from_slice(&junction);

    let (truth_copy_a, truth_copy_b) = if insert_collapsed_long {
        let mut copy_a = sa_seq.clone();
        for seq in &c_seqs {
            copy_a.extend_from_slice(&seq[30..]);
        }
        copy_a.extend_from_slice(&l_seq[30..]);
        copy_a.extend_from_slice(&f_seq[30..]);
        copy_a.extend_from_slice(&t_seq[30..]);

        let mut copy_b = sb_seq.clone();
        copy_b.extend_from_slice(&l_seq[30..]);
        copy_b.extend_from_slice(&f_seq[30..]);
        copy_b.extend_from_slice(&w_seq[30..]);
        (Some(copy_a), Some(copy_b))
    } else {
        (None, None)
    };

    let mut graph = CdbgGraph::new(31);
    let sa = add_node(&mut graph, &sa_seq);
    let mut prev = sa;
    for seq in &c_seqs {
        let id = add_node(&mut graph, seq);
        add_edge(&mut graph, prev, id, 30);
        prev = id;
    }

    let (evidence_id, long_id, sb_id) = if insert_collapsed_long {
        let l = add_node(&mut graph, &l_seq);
        add_edge(&mut graph, prev, l, 30);
        let sb = add_node(&mut graph, &sb_seq);
        // Second unique flank into the same long node: pooled collapsed-repeat
        // topology (two genomic copies share L).
        add_edge(&mut graph, sb, l, 30);
        prev = l;
        (l.id(), Some(l.id()), Some(sb.id()))
    } else {
        (prev.id(), None, None)
    };

    let f = add_node(&mut graph, &f_seq);
    add_edge(&mut graph, prev, f, 30);
    let w = add_node(&mut graph, &w_seq);
    let t = add_node(&mut graph, &t_seq);
    add_edge(&mut graph, f, w, 30);
    add_edge(&mut graph, f, t, 30);

    // F is only fork_len bp; take an interior window of the unique body.
    let f_tail = f_seq[40..100].to_vec();
    let (sb_tail, l_body) = if insert_collapsed_long {
        (
            Some(sb_seq[100..220].to_vec()),
            Some(l_seq[100..220].to_vec()),
        )
    } else {
        (None, None)
    };

    CollapseFixture {
        graph,
        matrix: PairLinkMatrix::new(),
        sa_tail: sa_seq[100..220].to_vec(),
        w_tail: w_seq[200..320].to_vec(),
        t_tail: t_seq[200..320].to_vec(),
        f_tail,
        sb_tail,
        l_body,
        short_bodies,
        truth_copy_a,
        truth_copy_b,
        evidence_id,
        w_id: w.id(),
        t_id: t.id(),
        long_id,
        sb_id,
    }
}

/// Inclusive Case A edge: `fork_dist == max_insert` (700) must remain admitted.
/// Geometry: gap 3*(200-30)=510, fork contribution 220-30=190, total 700.
#[test]
fn fork_dist_equal_max_insert_is_allowed() {
    assert_case_a_fork_dist(200, 220, 700);
    let mut fx = build_fixture(200, 220);

    add_rows(&mut fx.matrix, u(C3), u(W), 0..16, 100, 60);
    add_rows(&mut fx.matrix, u(C3), u(T), 16..22, 100, 60);

    let stats = resolve_repeats_exspander_with_stats(
        &mut fx.graph,
        &fx.matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    assert_eq!(
        stats.withheld_no_unique_anchor, 0,
        "fork_dist == max_insert must use inclusive `>` and must not withhold"
    );

    let seqs = active_sequences(&fx.graph);
    let bridged = seqs
        .iter()
        .any(|seq| contains(seq, &fx.s_tail) && contains(seq, &fx.w_tail));
    let t_retained = alive(&fx.graph, T) && seqs.iter().any(|seq| contains(seq, &fx.t_tail));
    let w_consumed = !alive(&fx.graph, W);

    assert!(
        bridged && w_consumed,
        "inclusive-edge fork must keep pre-F23 commit semantics (S+W bridge, W consumed); \
         bridged={bridged}, W alive={}",
        alive(&fx.graph, W)
    );
    assert!(
        t_retained,
        "non-selected exit T must remain an observable preserved candidate"
    );
}

/// Strict exterior: `fork_dist == max_insert + 1` (701) must withhold.
/// Geometry: gap 510 + fork contribution 221-30=191 = 701.
#[test]
fn fork_dist_max_insert_plus_one_is_withheld() {
    assert_case_a_fork_dist(200, 221, 701);
    let mut fx = build_fixture(200, 221);

    add_rows(&mut fx.matrix, u(C3), u(W), 0..16, 100, 60);
    add_rows(&mut fx.matrix, u(C3), u(T), 16..22, 100, 60);

    let stats = resolve_repeats_exspander_with_stats(
        &mut fx.graph,
        &fx.matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    assert_eq!(
        stats.withheld_no_unique_anchor, 1,
        "fork_dist == max_insert + 1 must be withheld"
    );

    let seqs = active_sequences(&fx.graph);
    let chimera_sw = seqs
        .iter()
        .any(|seq| contains(seq, &fx.s_tail) && contains(seq, &fx.w_tail));
    let chimera_st = seqs
        .iter()
        .any(|seq| contains(seq, &fx.s_tail) && contains(seq, &fx.t_tail));
    let t_retained = seqs.iter().any(|seq| contains(seq, &fx.t_tail));
    let w_retained = seqs.iter().any(|seq| contains(seq, &fx.w_tail));

    assert!(
        !chimera_sw && !chimera_st,
        "max_insert+1 withhold must not splice S onto either exit (S+W={chimera_sw}, S+T={chimera_st})"
    );
    assert!(
        t_retained && w_retained,
        "withholding must conserve W and T sequence in completed contigs (nodes may be \
         consumed into independent valid paths): T retained = {t_retained}, \
         W retained = {w_retained}",
    );
}

/// Case A must add `target_node.seq_len - overlap_k`. Fixture: gap
/// 3*(260-30)=690 <= 700, but full Case A distance 690+(100-30)=760 > 700.
/// Omitting the target contribution would admit the fork and commit the
/// 16:6 chimera; the correct formula withholds.
#[test]
fn case_a_fork_dist_includes_target_seq_len_minus_overlap() {
    let chain_len = 260usize;
    let fork_len = 100usize;
    let gap = 3 * (chain_len as i64 - 30);
    let target_contrib = fork_len as i64 - 30;
    let fork_dist = gap + target_contrib;
    assert_eq!(gap, 690);
    assert_eq!(target_contrib, 70);
    assert_eq!(fork_dist, 760);
    assert!(
        gap <= 700 && fork_dist > 700,
        "fixture must distinguish gap-only ({gap}) from full Case A ({fork_dist})"
    );

    let mut fx = build_fixture(chain_len, fork_len);
    add_rows(&mut fx.matrix, u(C3), u(W), 0..16, 120, 120);
    add_rows(&mut fx.matrix, u(C3), u(T), 16..22, 120, 120);

    let stats = resolve_repeats_exspander_with_stats(
        &mut fx.graph,
        &fx.matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    assert_eq!(
        stats.withheld_no_unique_anchor, 1,
        "Case A must count target_node.seq_len - overlap_k; omitting it yields \
         gap={gap} <= 700 and would wrongly admit this fork"
    );

    let seqs = active_sequences(&fx.graph);
    let chimera_sw = seqs
        .iter()
        .any(|seq| contains(seq, &fx.s_tail) && contains(seq, &fx.w_tail));
    let t_retained = seqs.iter().any(|seq| contains(seq, &fx.t_tail));
    let w_retained = seqs.iter().any(|seq| contains(seq, &fx.w_tail));
    assert!(
        !chimera_sw,
        "REGRESSION: Case A missing target contribution would admit gap={gap} and \
         commit the ectopic S+W chimera under 16:6 repeat-internal evidence"
    );
    assert!(
        t_retained && w_retained,
        "Case A withhold must conserve W and T sequence in completed contigs \
         (independent seed path consumption is allowed): T retained = {t_retained}, \
         W retained = {w_retained}",
    );
}

/// Reverse-corridor mirror of the beyond-horizon withhold (seed extends only
/// via `S.twin()`; evidence keyed on twin orientations).
#[test]
fn mirrored_reverse_corridor_beyond_horizon_withholds() {
    let mut fx = build_fixture_mirrored(300, 200);

    add_rows(&mut fx.matrix, u_rc(C3), u_rc(W), 0..16, 120, 120);
    add_rows(&mut fx.matrix, u_rc(C3), u_rc(T), 16..22, 120, 120);

    let stats = resolve_repeats_exspander_with_stats(
        &mut fx.graph,
        &fx.matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    assert_eq!(
        stats.withheld_no_unique_anchor, 1,
        "mirrored reverse corridor must withhold the beyond-horizon fork"
    );

    let seqs = active_sequences(&fx.graph);
    let chimera_sw = seqs
        .iter()
        .any(|seq| contains(seq, &fx.s_tail) && contains(seq, &fx.w_tail));
    let chimera_st = seqs
        .iter()
        .any(|seq| contains(seq, &fx.s_tail) && contains(seq, &fx.t_tail));
    let t_retained = seqs.iter().any(|seq| contains(seq, &fx.t_tail));
    let w_retained = seqs.iter().any(|seq| contains(seq, &fx.w_tail));

    assert!(
        !chimera_sw && !chimera_st,
        "reverse-corridor withhold must not splice S onto W/T (S+W={chimera_sw}, S+T={chimera_st})"
    );
    assert!(
        t_retained && w_retained,
        "mirrored withhold must conserve both exits in completed contigs: \
         T retained = {t_retained}, W retained = {w_retained}",
    );
}

/// Correctness contract: a ≥500 bp node that represents a collapsed /
/// non-unique repeat must not reset `unique_anchor`. Pooled topology gives L
/// two distinct long flanks (S_A and S_B); after walking a short chain already
/// past S_A's insert horizon, resetting onto L would re-open a within-horizon
/// decision and accept the ectopic 16:6 join S_A..L..W.
///
/// Synthetic truth assignment for this fixture: copy A follows
/// S_A -> C1 -> C2 -> C3 -> L -> F -> T; copy B follows S_B -> L -> F -> W.
/// Pair IDs 0..16 are assigned to copy B's L->W evidence and IDs 16..22 to
/// copy A's L->T evidence. PairLinkMatrix pools both groups at shared node L,
/// so neither group is locus-specific evidence for the S_A fork. This is an
/// abstract evidence fixture, not a claim that raw reads were generated.
///
/// EXPECTED TO FAIL on production base `f121920a...`: current bookkeeping
/// treats every seq_len ≥ 500 node as unique and resets the gap, which is the
/// unsafe bypass this oracle refuses to encode as desired behavior.
#[test]
fn collapsed_long_node_must_not_reset_anchor_and_false_accept() {
    let mut fx = build_collapse_or_short_horizon(true);
    let long_id = fx.long_id.expect("collapsed long node");
    assert!(
        fx.sb_id.is_some(),
        "collapsed fixture must include the second unique flank S_B"
    );
    let sb_tail = fx.sb_tail.as_ref().expect("S_B motif");
    let l_body = fx.l_body.as_ref().expect("L motif");

    // Fixture preconditions: L is a ≥500 bp pooled repeat with two distinct
    // incoming flanks (C3 from the S_A chain and S_B). On the forward
    // orientation those predecessors appear as reverse-strand out-edges.
    let l_node = &fx.graph.nodes[long_id as usize];
    assert!(
        l_node.seq_len >= 500,
        "collapsed long node must be ≥500 bp (got {})",
        l_node.seq_len
    );
    let incoming_flank_ids: Vec<u32> = l_node.edges_rev.iter().map(|e| e.target.id()).collect();
    assert_eq!(
        incoming_flank_ids.len(),
        2,
        "collapsed L must have exactly two distinct incoming flank edges; got {incoming_flank_ids:?}"
    );
    assert_ne!(
        incoming_flank_ids[0], incoming_flank_ids[1],
        "incoming flank predecessors must be distinct nodes"
    );

    // Pooled repeat-internal evidence keyed at the collapsed long node — the
    // scoring curr_dir after an unsafe reset — favoring ectopic exit W.
    add_rows(
        &mut fx.matrix,
        u(fx.evidence_id),
        u(fx.w_id),
        0..16,
        120,
        120,
    );
    add_rows(
        &mut fx.matrix,
        u(fx.evidence_id),
        u(fx.t_id),
        16..22,
        120,
        120,
    );

    let truth_copy_a = fx.truth_copy_a.as_ref().expect("copy A reference path");
    let truth_copy_b = fx.truth_copy_b.as_ref().expect("copy B reference path");
    assert!(
        contains(truth_copy_a, &fx.sa_tail)
            && contains(truth_copy_a, &fx.t_tail)
            && !contains(truth_copy_a, &fx.w_tail),
        "copy A truth must contain S_A and T, but not W"
    );
    assert!(
        contains(truth_copy_b, sb_tail)
            && contains(truth_copy_b, &fx.w_tail)
            && !contains(truth_copy_b, &fx.sa_tail)
            && !contains(truth_copy_b, &fx.t_tail),
        "copy B truth must contain S_B and W, but not S_A or T"
    );

    let dist = InsertSizeDist::default();
    let params = ExspanderParams::default();
    let candidate_separation = 150 - 2 * 30;
    let w_score = score_candidate_path(
        &fx.matrix,
        &dist,
        u(fx.evidence_id),
        u(fx.w_id),
        candidate_separation,
    );
    let t_score = score_candidate_path(
        &fx.matrix,
        &dist,
        u(fx.evidence_id),
        u(fx.t_id),
        candidate_separation,
    );
    assert_eq!(w_score.0, 16, "all 16 pooled L→W PE rows must be in-window");
    assert_eq!(t_score.0, 6, "all 6 pooled L→T PE rows must be in-window");
    assert!(
        w_score.1 / t_score.1 >= params.decision_threshold,
        "pooled L evidence should favor W by the configured decision threshold"
    );

    let stats = resolve_repeats_exspander_with_stats(
        &mut fx.graph,
        &fx.matrix,
        &dist,
        params,
    )
    .expect("repeat resolution failed");

    let seqs = active_sequences(&fx.graph);
    let chimera_sa_w = seqs
        .iter()
        .any(|seq| contains(seq, &fx.sa_tail) && contains(seq, &fx.w_tail));
    let sa_retained = seqs.iter().any(|seq| contains(seq, &fx.sa_tail));
    let w_retained = seqs.iter().any(|seq| contains(seq, &fx.w_tail));
    let t_retained = seqs.iter().any(|seq| contains(seq, &fx.t_tail));
    let f_retained = seqs.iter().any(|seq| contains(seq, &fx.f_tail));
    let sb_retained = seqs.iter().any(|seq| contains(seq, sb_tail));
    let l_retained = seqs.iter().any(|seq| contains(seq, l_body));
    let withholding_ok = stats.withheld_no_unique_anchor >= 1;
    let chimera_ok = !chimera_sa_w;
    let pe_summary = |target| {
        fx.matrix
            .get_evidence(u(fx.evidence_id), u(target))
            .map(|rows| {
                let min_u1 = rows.iter().map(|row| row.dist_u1_exit).min();
                let max_u1 = rows.iter().map(|row| row.dist_u1_exit).max();
                let min_u2 = rows.iter().map(|row| row.dist_u2_entry).min();
                let max_u2 = rows.iter().map(|row| row.dist_u2_entry).max();
                format!(
                    "count={},dist_u1_exit_min_max={min_u1:?}:{max_u1:?},dist_u2_entry_min_max={min_u2:?}:{max_u2:?}",
                    rows.len()
                )
            })
            .unwrap_or_else(|| "count=0,no_rows".to_owned())
    };

    eprintln!(
        "[AUD01_DIAG] withheld_no_unique_anchor={}",
        stats.withheld_no_unique_anchor
    );
    eprintln!("[AUD01_DIAG] withholding_ok={withholding_ok}");
    eprintln!("[AUD01_DIAG] chimera_sa_w={chimera_sa_w}");
    eprintln!("[AUD01_DIAG] chimera_ok={chimera_ok}");
    eprintln!(
        "[AUD01_DIAG] motif_retained S_A={sa_retained} L={l_retained} W={w_retained} T={t_retained} F={f_retained} S_B={sb_retained}"
    );
    eprintln!(
        "[AUD01_DIAG] fixture_truth copy_A=S_A-C1-C2-C3-L-F-T copy_B=S_B-L-F-W pe_pair_ids_0_16=L-W_from_copy_B pe_pair_ids_16_22=L-T_from_copy_A"
    );
    eprintln!(
        "[AUD01_DIAG] truth_reference copy_A_len={} copy_A_has_S_A={} copy_A_has_T={} copy_A_has_W={} copy_B_len={} copy_B_has_S_B={} copy_B_has_W={} copy_B_has_S_A={} copy_B_has_T={}",
        truth_copy_a.len(),
        contains(truth_copy_a, &fx.sa_tail),
        contains(truth_copy_a, &fx.t_tail),
        contains(truth_copy_a, &fx.w_tail),
        truth_copy_b.len(),
        contains(truth_copy_b, sb_tail),
        contains(truth_copy_b, &fx.w_tail),
        contains(truth_copy_b, &fx.sa_tail),
        contains(truth_copy_b, &fx.t_tail)
    );
    eprintln!(
        "[AUD01_DIAG] fixture_geometry max_insert={} overlap=30 expected_no_reset_fork_dist=1500 expected_after_L_reset_fork_dist=120",
        InsertSizeDist::default().max_insert
    );
    eprintln!(
        "[AUD01_DIAG] candidate_pe separation={} span={} W_count={} T_count={} weight_ratio={} decision_threshold={}",
        candidate_separation,
        120 + candidate_separation + 120,
        w_score.0,
        t_score.0,
        w_score.1 / t_score.1,
        params.decision_threshold
    );
    eprintln!(
        "[AUD01_DIAG] pe_L_to_W {}",
        pe_summary(fx.w_id)
    );
    eprintln!(
        "[AUD01_DIAG] pe_L_to_T {}",
        pe_summary(fx.t_id)
    );
    eprintln!("[AUD01_DIAG] candidate_trace=CANDIDATE_TRACE_UNAVAILABLE");
    for (index, seq) in seqs.iter().enumerate() {
        eprintln!(
            "[AUD01_DIAG] contig index={index} length={} has_S_A={} has_L={} has_W={} has_T={} has_F={} has_S_B={} has_S_A_and_W={}",
            seq.len(),
            contains(seq, &fx.sa_tail),
            contains(seq, l_body),
            contains(seq, &fx.w_tail),
            contains(seq, &fx.t_tail),
            contains(seq, &fx.f_tail),
            contains(seq, sb_tail),
            contains(seq, &fx.sa_tail) && contains(seq, &fx.w_tail)
        );
    }

    assert!(
        stats.withheld_no_unique_anchor >= 1,
        "CORRECTNESS CONTRACT (expected FAIL on base f121920a): collapsed long \
         node L must not reset the unique anchor; the S_A-anchored walk is \
         already beyond max_insert before L, so the fork must withhold \
         (withheld_no_unique_anchor>=1, got {}). Unsafe ≥500 reset yields a \
         false within-horizon decision on pooled L-keyed 16:6 evidence.",
        stats.withheld_no_unique_anchor
    );
    assert!(
        !chimera_sa_w,
        "CORRECTNESS CONTRACT (expected FAIL on base f121920a): must not commit \
         chimeric S_A..L..W (chimera={chimera_sa_w}). Length≥500 is not genomic \
         uniqueness; L is pooled across S_A and S_B. W may still be consumed into \
         an independent F+W contig — that is not the forbidden S_A..W join."
    );
    assert!(
        sa_retained && w_retained && t_retained && f_retained && sb_retained && l_retained,
        "fixture motifs must remain represented after resolution (independent \
         seed merges may consume nodes): S_A={sa_retained} W={w_retained} \
         T={t_retained} F={f_retained} S_B={sb_retained} L={l_retained}"
    );
}

/// Genuinely unique short-node chain (disjoint dna slots / distinct bodies)
/// beyond the insert horizon must still be withheld — the guard is
/// conservative on length heuristics and does not promote short unique nodes
/// into fresh anchors. Shares topology with the collapsed-long counterexample
/// minus the ≥500 intermediate.
#[test]
fn unique_short_chain_beyond_horizon_conservatively_withheld() {
    let mut fx = build_collapse_or_short_horizon(false);
    assert!(fx.long_id.is_none());

    // Independent diagnosability: short interiors are pairwise distinct, so
    // any withhold here is not explained by sequence collapse among C1..C3.
    for i in 0..fx.short_bodies.len() {
        for j in (i + 1)..fx.short_bodies.len() {
            assert_ne!(
                fx.short_bodies[i], fx.short_bodies[j],
                "short chain bodies C{} and C{} must be unique for this oracle",
                i + 1,
                j + 1
            );
        }
    }

    add_rows(
        &mut fx.matrix,
        u(fx.evidence_id),
        u(fx.w_id),
        0..16,
        120,
        120,
    );
    add_rows(
        &mut fx.matrix,
        u(fx.evidence_id),
        u(fx.t_id),
        16..22,
        120,
        120,
    );

    let stats = resolve_repeats_exspander_with_stats(
        &mut fx.graph,
        &fx.matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    assert_eq!(
        stats.withheld_no_unique_anchor, 1,
        "unique short chain with Case A distance 810+(150-30)=930 > 700 must withhold"
    );

    let seqs = active_sequences(&fx.graph);
    let chimera_sa_w = seqs
        .iter()
        .any(|seq| contains(seq, &fx.sa_tail) && contains(seq, &fx.w_tail));
    let chimera_sa_t = seqs
        .iter()
        .any(|seq| contains(seq, &fx.sa_tail) && contains(seq, &fx.t_tail));
    let sa_retained = seqs.iter().any(|seq| contains(seq, &fx.sa_tail));
    let w_retained = seqs.iter().any(|seq| contains(seq, &fx.w_tail));
    let t_retained = seqs.iter().any(|seq| contains(seq, &fx.t_tail));
    let f_retained = seqs.iter().any(|seq| contains(seq, &fx.f_tail));

    assert!(
        !chimera_sa_w && !chimera_sa_t,
        "conservative withhold must not splice S_A onto W/T despite unique shorts \
         (S_A+W={chimera_sa_w}, S_A+T={chimera_sa_t})"
    );
    assert!(
        sa_retained && w_retained && t_retained && f_retained,
        "unique-short beyond-horizon withhold must conserve S_A, F, and both exits \
         as observable sequence (nodes may be consumed into independent valid \
         paths): S_A={sa_retained} F={f_retained} W={w_retained} T={t_retained}"
    );
}

fn f23_same_contig(seqs: &[Vec<u8>], left: &[u8], right: &[u8]) -> bool {
    let left_rc = reverse_complement_ascii(left);
    let right_rc = reverse_complement_ascii(right);
    seqs.iter().any(|seq| {
        (contains(seq, left) && contains(seq, right))
            || (contains(seq, &left_rc) && contains(seq, &right_rc))
    })
}

fn f23_add_pair_rows(
    matrix: &mut PairLinkMatrix,
    from: DirUnitigId,
    to: DirUnitigId,
    pair_ids: std::ops::Range<u32>,
    from_exit: u32,
    to_entry: u32,
    read_len: u16,
) {
    for pair_id in pair_ids {
        matrix.add_evidence(
            from,
            to,
            ReadPairEvidence {
                pair_id,
                dist_u1_exit: from_exit,
                dist_u2_entry: to_entry,
                r1_len: read_len,
                r2_len: read_len,
            },
        );
        matrix.add_evidence(
            to.twin(),
            from.twin(),
            ReadPairEvidence {
                pair_id,
                dist_u1_exit: to_entry,
                dist_u2_entry: from_exit,
                r1_len: read_len,
                r2_len: read_len,
            },
        );
    }
}

/// F2: declared truth is A=S_A-C1-C2-C3-L-F-T and B=S_B-L-F-W.
/// The F1 graph and geometry are unchanged, but pooled L→T has 16 rows and
/// L→W has 6. Both groups remain non-attributable to the S_A walk. The
/// source-key separation is 90 and the synthetic span is 120+90+120=330;
/// the last attributable flank is 1500 bp away, beyond max_insert=700.
/// The old implementation is expected to join S_A to T and fail this oracle.
#[test]
fn f23_opus_f2_flipped_pooled_weights_still_withhold() {
    let mut fx = build_collapse_or_short_horizon(true);
    let l = u(fx.evidence_id);
    assert!(fx.truth_copy_a.is_some() && fx.truth_copy_b.is_some());
    f23_add_pair_rows(&mut fx.matrix, l, u(fx.t_id), 0..16, 120, 120, 100);
    f23_add_pair_rows(&mut fx.matrix, l, u(fx.w_id), 16..22, 120, 120, 100);

    let stats = resolve_repeats_exspander_with_stats(
        &mut fx.graph,
        &fx.matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");
    let seqs = active_sequences(&fx.graph);
    assert!(
        !f23_same_contig(&seqs, &fx.sa_tail, &fx.t_tail)
            && !f23_same_contig(&seqs, &fx.sa_tail, &fx.w_tail),
        "flipped pooled evidence must not join S_A to either exit across F"
    );
    assert!(seqs.iter().any(|seq| contains(seq, &fx.sa_tail)));
    assert!(stats.withheld_no_unique_anchor >= 1);
}

fn f23_seed_node(seed: u64, len: usize, junction: &[u8]) -> Vec<u8> {
    let mut seq = junction.to_vec();
    seq.extend(dna(slot(seed), len - 2 * junction.len()));
    seq.extend_from_slice(junction);
    seq
}

fn f23_add_bio_node(graph: &mut CdbgGraph, seq: &[u8], mirror: bool) -> DirUnitigId {
    let stored = if mirror {
        reverse_complement_ascii(seq)
    } else {
        seq.to_vec()
    };
    let id = add_node(graph, &stored);
    if mirror {
        id.twin()
    } else {
        id
    }
}

/// F3: L is the sole eligible long seed. Declared copies are
/// S_A→L→F→T and S_B→L→F→W. Forward L→W/T and reverse
/// L.twin→S_A.twin/S_B.twin evidence are pooled across those copies.
/// F has length 150, overlap 30, so the forward separation is 90 and
/// 150+90+150=390; reverse direct separation is −30 and span 200−30+200=370.
/// No source-specific anchor exists at L. Both forks must abstain and retain
/// S_A and S_B. Current code is expected to select pooled links from L.
fn f23_check_shared_long_seed(mirror: bool) {
    let junction = dna(slot(3000), 30);
    let sa_seq = f23_seed_node(100, 400, &junction);
    let sb_seq = f23_seed_node(101, 400, &junction);
    let l_seq = f23_seed_node(102, 600, &junction);
    let f_seq = f23_seed_node(103, 150, &junction);
    let t_seq = f23_seed_node(104, 400, &junction);
    let w_seq = f23_seed_node(105, 400, &junction);
    let mut graph = CdbgGraph::new(31);
    let sa = f23_add_bio_node(&mut graph, &sa_seq, mirror);
    let sb = f23_add_bio_node(&mut graph, &sb_seq, mirror);
    let l = f23_add_bio_node(&mut graph, &l_seq, mirror);
    let f = f23_add_bio_node(&mut graph, &f_seq, mirror);
    let t = f23_add_bio_node(&mut graph, &t_seq, mirror);
    let w = f23_add_bio_node(&mut graph, &w_seq, mirror);
    for (from, to) in [(sa, l), (sb, l), (l, f), (f, t), (f, w)] {
        add_edge(&mut graph, from, to, 30);
    }
    assert_eq!(
        graph
            .nodes
            .iter()
            .filter(|node| node.seq_len >= 500)
            .count(),
        1
    );
    assert_eq!(
        graph.nodes[l.id() as usize].edges(l.twin().is_rc()).len(),
        2
    );
    let mut matrix = PairLinkMatrix::new();
    f23_add_pair_rows(&mut matrix, l, w, 0..16, 150, 150, 150);
    f23_add_pair_rows(&mut matrix, l, t, 16..22, 150, 150, 150);
    f23_add_pair_rows(&mut matrix, l.twin(), sa.twin(), 22..38, 200, 200, 150);
    f23_add_pair_rows(&mut matrix, l.twin(), sb.twin(), 38..44, 200, 200, 150);

    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");
    let seqs = active_sequences(&graph);
    let l_motif = &l_seq[100..220];
    for flank in [
        &sa_seq[100..220],
        &sb_seq[100..220],
        &t_seq[100..220],
        &w_seq[100..220],
    ] {
        assert!(
            !f23_same_contig(&seqs, l_motif, flank),
            "shared L seed must not join a fork flank (mirror={mirror})"
        );
    }
    assert!(
        !graph.nodes[sa.id() as usize].is_deleted,
        "S_A must remain unconsumed"
    );
    assert!(
        !graph.nodes[sb.id() as usize].is_deleted,
        "S_B must remain unconsumed"
    );
    assert!(stats.withheld_no_unique_anchor >= 1);
}

#[test]
fn f23_opus_f3_shared_long_seed_forward() {
    f23_check_shared_long_seed(false);
}

#[test]
fn f23_opus_f3_shared_long_seed_mirrored() {
    f23_check_shared_long_seed(true);
}

/// F4: Shared short predecessor bypass.
/// Topology: S_A -> C1 -> C2 -> C3 -> X -> L -> F -> {T, W}
///           S_B -> X
/// where X is under 500 bp (120 bp).
/// L is a 600 bp long node whose immediate predecessor is X (in-degree 1 to L).
/// Because X has in-degree 2 from C3 and S_B, L inherits shared context.
/// Total distance from S_A exceeds max_insert (700 bp).
/// Pooled L evidence favors W (16:6).
/// Expected RED on unpatched code: production checks only immediate predecessors of L,
/// resets the anchor at L, and falsely accepts the chimera S_A..W.
#[test]
fn f23_opus_f4_shared_short_predecessor_bypass() {
    let junction = dna(slot(4000), 30);
    let sa_seq = f23_seed_node(200, 600, &junction);
    let c1_seq = f23_seed_node(201, 300, &junction);
    let c2_seq = f23_seed_node(202, 300, &junction);
    let c3_seq = f23_seed_node(203, 300, &junction);
    let sb_seq = f23_seed_node(204, 600, &junction);
    let x_seq = f23_seed_node(205, 120, &junction);
    let l_seq = f23_seed_node(206, 600, &junction);
    let f_seq = f23_seed_node(207, 150, &junction);
    let t_seq = f23_seed_node(208, 400, &junction);
    let w_seq = f23_seed_node(209, 400, &junction);

    let mut graph = CdbgGraph::new(31);
    let sa = add_node(&mut graph, &sa_seq);
    let c1 = add_node(&mut graph, &c1_seq);
    let c2 = add_node(&mut graph, &c2_seq);
    let c3 = add_node(&mut graph, &c3_seq);
    let sb = add_node(&mut graph, &sb_seq);
    let x = add_node(&mut graph, &x_seq);
    let l = add_node(&mut graph, &l_seq);
    let f = add_node(&mut graph, &f_seq);
    let t = add_node(&mut graph, &t_seq);
    let w = add_node(&mut graph, &w_seq);

    for (from, to) in [
        (sa, c1), (c1, c2), (c2, c3), (c3, x),
        (sb, x),
        (x, l),
        (l, f),
        (f, t), (f, w),
    ] {
        add_edge(&mut graph, from, to, 30);
    }

    assert_eq!(graph.nodes[x.id() as usize].edges(x.twin().is_rc()).len(), 2, "X must have in-degree 2");
    assert_eq!(graph.nodes[l.id() as usize].edges(l.twin().is_rc()).len(), 1, "L has immediate in-degree 1");

    let mut matrix = PairLinkMatrix::new();
    f23_add_pair_rows(&mut matrix, l, w, 0..16, 120, 120, 100);
    f23_add_pair_rows(&mut matrix, l, t, 16..22, 120, 120, 100);

    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    let seqs = active_sequences(&graph);
    let sa_tail = &sa_seq[100..220];
    let w_tail = &w_seq[100..220];
    let t_tail = &t_seq[100..220];

    assert!(
        !f23_same_contig(&seqs, sa_tail, w_tail) && !f23_same_contig(&seqs, sa_tail, t_tail),
        "short predecessor X must propagate shared context to L; fork must withhold"
    );
    assert!(stats.withheld_no_unique_anchor >= 1, "fork must be withheld");
}

/// F5: Legal resolution through a short repeat (Positive Control).
/// Topology: S_A(600) -> C1(80) -> L(200) -> F(100) -> {T(400), W(400)}
///           S_B(600) -> L(200)
/// Total distance from S_A through C1, L, F to exit is (80-30)+(200-30)+(100-30) = 50+170+70 = 290 bp <= max_insert (700 bp).
/// S_A is an attributable anchor within the insert horizon!
/// Attributable evidence: S_A -> T (20 pairs, span 460 in window [250, 700]).
/// Misleading pooled evidence: L -> W (16 pairs, span 300 in window [250, 700]), L -> T (6 pairs, span 300).
/// Expected: S_A..T is resolved and committed; W is rejected.
#[test]
fn f23_opus_f5_legal_resolution_through_short_repeat() {
    let junction = dna(slot(5000), 30);
    let sa_seq = f23_seed_node(300, 600, &junction);
    let c1_seq = f23_seed_node(301, 80, &junction);
    let sb_seq = f23_seed_node(302, 600, &junction);
    let l_seq = f23_seed_node(303, 200, &junction);
    let f_seq = f23_seed_node(304, 100, &junction);
    let t_seq = f23_seed_node(305, 400, &junction);
    let w_seq = f23_seed_node(306, 400, &junction);

    let mut graph = CdbgGraph::new(31);
    let sa = add_node(&mut graph, &sa_seq);
    let c1 = add_node(&mut graph, &c1_seq);
    let sb = add_node(&mut graph, &sb_seq);
    let l = add_node(&mut graph, &l_seq);
    let f = add_node(&mut graph, &f_seq);
    let t = add_node(&mut graph, &t_seq);
    let w = add_node(&mut graph, &w_seq);

    for (from, to) in [
        (sa, c1), (c1, l), (sb, l), (l, f), (f, t), (f, w),
    ] {
        add_edge(&mut graph, from, to, 30);
    }

    let mut matrix = PairLinkMatrix::new();
    // Attributable evidence from S_A to authentic exit T (span = 100 + 260 + 100 = 460):
    f23_add_pair_rows(&mut matrix, sa, t, 0..20, 100, 100, 100);
    // Misleading repeat-internal pooled evidence at L favoring W (span = 130 + 40 + 130 = 300):
    f23_add_pair_rows(&mut matrix, l, w, 20..36, 130, 130, 100);
    f23_add_pair_rows(&mut matrix, l, t, 36..42, 130, 130, 100);

    let _stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    let seqs = active_sequences(&graph);
    let sa_tail = &sa_seq[100..220];
    let t_tail = &t_seq[100..220];
    let w_tail = &w_seq[100..220];

    assert!(
        f23_same_contig(&seqs, sa_tail, t_tail),
        "in-horizon attributable anchor S_A must legally resolve S_A..T"
    );
    assert!(
        !f23_same_contig(&seqs, sa_tail, w_tail),
        "misleading pooled evidence at L must not commit S_A..W"
    );
}

/// F6(a): Continuity control: Unique long node with no convergence resets anchor.
/// Topology: S(600) -> U(600) -> F(150) -> {T(400), W(400)}.
/// U has in-degree 1 (no convergence). U is unique and >= 500 bp.
/// Evidence keyed at U favors T: U -> T (20 pairs), U -> W (2 pairs).
/// Expected: S..U..F..T commits cleanly.
#[test]
fn f23_opus_f6a_unique_long_node_continuity() {
    let junction = dna(slot(6000), 30);
    let s_seq = f23_seed_node(400, 600, &junction);
    let u_seq = f23_seed_node(401, 600, &junction);
    let f_seq = f23_seed_node(402, 150, &junction);
    let t_seq = f23_seed_node(403, 400, &junction);
    let w_seq = f23_seed_node(404, 400, &junction);

    let mut graph = CdbgGraph::new(31);
    let s = add_node(&mut graph, &s_seq);
    let u_node = add_node(&mut graph, &u_seq);
    let f = add_node(&mut graph, &f_seq);
    let t = add_node(&mut graph, &t_seq);
    let w = add_node(&mut graph, &w_seq);

    for (from, to) in [(s, u_node), (u_node, f), (f, t), (f, w)] {
        add_edge(&mut graph, from, to, 30);
    }

    let mut matrix = PairLinkMatrix::new();
    f23_add_pair_rows(&mut matrix, u_node, t, 0..20, 100, 100, 100);
    f23_add_pair_rows(&mut matrix, u_node, w, 20..22, 100, 100, 100);

    let _stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    let seqs = active_sequences(&graph);
    let s_tail = &s_seq[100..220];
    let t_tail = &t_seq[100..220];
    assert!(
        f23_same_contig(&seqs, s_tail, t_tail),
        "unique long node U must reset anchor and resolve S..U..F..T"
    );
}

/// F6(b): Continuity control: SNP bubble diverging and reconverging onto L within insert horizon.
/// Topology: S(600) -> B1(80) -> L(200) -> F(100) -> {T(400), W(400)}
///           S(600) -> B2(80) -> L(200)
/// Bubble branches B1 and B2 diverge from single locus S and reconverge at L.
/// Bubble span = 80 <= max_insert (700). Reconvergence check certifies same locus.
/// Evidence from S favors T: S -> T (20 pairs, span 460 in [250, 700]), S -> W (0 pairs).
/// Expected: bubble does not mark shared context; S..T resolves.
#[test]
fn f23_opus_f6b_snp_bubble_continuity() {
    let junction = dna(slot(7000), 30);
    let s_seq = f23_seed_node(500, 600, &junction);
    let b1_seq = f23_seed_node(501, 80, &junction);
    let b2_seq = f23_seed_node(502, 80, &junction);
    let l_seq = f23_seed_node(503, 200, &junction);
    let f_seq = f23_seed_node(504, 100, &junction);
    let t_seq = f23_seed_node(505, 400, &junction);
    let w_seq = f23_seed_node(506, 400, &junction);

    let mut graph = CdbgGraph::new(31);
    let s = add_node(&mut graph, &s_seq);
    let b1 = add_node(&mut graph, &b1_seq);
    let b2 = add_node(&mut graph, &b2_seq);
    let l = add_node(&mut graph, &l_seq);
    let f = add_node(&mut graph, &f_seq);
    let t = add_node(&mut graph, &t_seq);
    let w = add_node(&mut graph, &w_seq);

    for (from, to) in [
        (s, b1), (s, b2), (b1, l), (b2, l), (l, f), (f, t), (f, w),
    ] {
        add_edge(&mut graph, from, to, 30);
    }

    let mut matrix = PairLinkMatrix::new();
    // Primary allele support bridging S across bubble arm B1 into L (span = 150 - 30 + 150 = 270 in [250, 700]):
    f23_add_pair_rows(&mut matrix, s, b1, 0..20, 150, 150, 100);
    // Attributable evidence from S to authentic exit T (span = 100 + 260 + 100 = 460 in [250, 700]):
    f23_add_pair_rows(&mut matrix, s, t, 20..40, 100, 100, 100);
    f23_add_pair_rows(&mut matrix, l, t, 40..55, 130, 130, 100);
    f23_add_pair_rows(&mut matrix, l, w, 55..58, 130, 130, 100);

    let _stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    let seqs = active_sequences(&graph);
    let s_tail = &s_seq[100..220];
    let t_tail = &t_seq[100..220];
    assert!(
        f23_same_contig(&seqs, s_tail, t_tail),
        "bubble reconvergence within horizon must be exempted from shared context; S..T must resolve"
    );
}

/// F7: Case B boundary distances (`anchor_gap > max_insert`).
#[test]
fn f23_opus_f7_case_b_boundary_distances() {
    let gap_admit_exact = (263 - 30) + (263 - 30) + (264 - 30);
    assert_eq!(gap_admit_exact, 700);
    assert!(gap_admit_exact <= 700, "gap 700 <= 700 must be admitted under Case B");

    let gap_withhold = (263 - 30) + (264 - 30) + (264 - 30);
    assert_eq!(gap_withhold, 701);
    assert!(gap_withhold > 700, "gap 701 > 700 must be withheld under Case B");
}

/// F8: Reverse-complement mirror of F4.
#[test]
fn f23_opus_f8_reverse_mirror_f4() {
    let junction = dna(slot(8000), 30);
    let sa_seq = f23_seed_node(600, 600, &junction);
    let c1_seq = f23_seed_node(601, 300, &junction);
    let c2_seq = f23_seed_node(602, 300, &junction);
    let c3_seq = f23_seed_node(603, 300, &junction);
    let sb_seq = f23_seed_node(604, 600, &junction);
    let x_seq = f23_seed_node(605, 120, &junction);
    let l_seq = f23_seed_node(606, 600, &junction);
    let f_seq = f23_seed_node(607, 150, &junction);
    let t_seq = f23_seed_node(608, 400, &junction);
    let w_seq = f23_seed_node(609, 400, &junction);

    let mut graph = CdbgGraph::new(31);
    let sa = f23_add_bio_node(&mut graph, &sa_seq, true);
    let c1 = f23_add_bio_node(&mut graph, &c1_seq, true);
    let c2 = f23_add_bio_node(&mut graph, &c2_seq, true);
    let c3 = f23_add_bio_node(&mut graph, &c3_seq, true);
    let sb = f23_add_bio_node(&mut graph, &sb_seq, true);
    let x = f23_add_bio_node(&mut graph, &x_seq, true);
    let l = f23_add_bio_node(&mut graph, &l_seq, true);
    let f = f23_add_bio_node(&mut graph, &f_seq, true);
    let t = f23_add_bio_node(&mut graph, &t_seq, true);
    let w = f23_add_bio_node(&mut graph, &w_seq, true);

    for (from, to) in [
        (sa, c1), (c1, c2), (c2, c3), (c3, x),
        (sb, x),
        (x, l),
        (l, f),
        (f, t), (f, w),
    ] {
        add_edge(&mut graph, from, to, 30);
    }

    let mut matrix = PairLinkMatrix::new();
    f23_add_pair_rows(&mut matrix, l, w, 0..16, 120, 120, 100);
    f23_add_pair_rows(&mut matrix, l, t, 16..22, 120, 120, 100);

    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    let seqs = active_sequences(&graph);
    let sa_tail = &sa_seq[100..220];
    let w_tail = &w_seq[100..220];
    let t_tail = &t_seq[100..220];

    assert!(
        !f23_same_contig(&seqs, sa_tail, w_tail) && !f23_same_contig(&seqs, sa_tail, t_tail),
        "mirrored F4 must also withhold"
    );
    assert!(stats.withheld_no_unique_anchor >= 1);
}

/// F9(a): Merge remap: After an earlier merge absorbs copy B's context, the convergence at L remains visible.
#[test]
fn f23_opus_f9a_merge_remap_retains_convergence() {
    let junction = dna(slot(9000), 30);
    let sa_seq = f23_seed_node(700, 600, &junction);
    let c1_seq = f23_seed_node(701, 300, &junction);
    let c2_seq = f23_seed_node(702, 300, &junction);
    let c3_seq = f23_seed_node(703, 300, &junction);
    let sb_merged_seq = f23_seed_node(704, 700, &junction);
    let l_seq = f23_seed_node(705, 600, &junction);
    let f_seq = f23_seed_node(706, 150, &junction);
    let t_seq = f23_seed_node(707, 400, &junction);
    let w_seq = f23_seed_node(708, 400, &junction);

    let mut graph = CdbgGraph::new(31);
    let sa = add_node(&mut graph, &sa_seq);
    let c1 = add_node(&mut graph, &c1_seq);
    let c2 = add_node(&mut graph, &c2_seq);
    let c3 = add_node(&mut graph, &c3_seq);
    let sb_m = add_node(&mut graph, &sb_merged_seq);
    let l = add_node(&mut graph, &l_seq);
    let f = add_node(&mut graph, &f_seq);
    let t = add_node(&mut graph, &t_seq);
    let w = add_node(&mut graph, &w_seq);

    for (from, to) in [
        (sa, c1), (c1, c2), (c2, c3), (c3, l),
        (sb_m, l),
        (l, f),
        (f, t), (f, w),
    ] {
        add_edge(&mut graph, from, to, 30);
    }

    let mut matrix = PairLinkMatrix::new();
    f23_add_pair_rows(&mut matrix, l, w, 0..16, 120, 120, 100);
    f23_add_pair_rows(&mut matrix, l, t, 16..22, 120, 120, 100);

    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    let seqs = active_sequences(&graph);
    let sa_tail = &sa_seq[100..220];
    let w_tail = &w_seq[100..220];
    let t_tail = &t_seq[100..220];
    assert!(
        !f23_same_contig(&seqs, sa_tail, w_tail) && !f23_same_contig(&seqs, sa_tail, t_tail),
        "merged predecessor must preserve in-degree convergence at L; fork must withhold"
    );
    assert!(stats.withheld_no_unique_anchor >= 1);
}

/// F9(b): Qualified cycle re-entry (Astra Architectural Ruling).
/// Re-entry into a multi-context collapsed cycle must not restore attribution or permit a fork decision
/// from pooled evidence when no admissible attributable support exists.
/// Topology: S_A -> C1..C3 -> L, S_B -> L
/// Cycle: L(200) -> R(100) -> L(200) (L < 500 bp so repeat cycle exploration is exercised)
/// Fork: L -> F(150) -> {T(400), W(400)}
/// S_A is beyond insert horizon (> 1500 bp).
/// Evidence at L: L -> W (16 pairs), L -> T (6 pairs).
/// Expected RED on baseline: cycle re-entry at multi-context L must withhold, not commit chimera S_A..W.
#[test]
fn f23_opus_f9b_qualified_cycle_reentry_withholds() {
    let junction = dna(slot(10000), 30);
    let sa_seq = f23_seed_node(800, 600, &junction);
    let c1_seq = f23_seed_node(801, 300, &junction);
    let c2_seq = f23_seed_node(802, 300, &junction);
    let c3_seq = f23_seed_node(803, 300, &junction);
    let sb_seq = f23_seed_node(804, 600, &junction);
    let l_seq = f23_seed_node(805, 200, &junction);
    let r_seq = f23_seed_node(806, 100, &junction);
    let f_seq = f23_seed_node(807, 150, &junction);
    let t_seq = f23_seed_node(808, 400, &junction);
    let w_seq = f23_seed_node(809, 400, &junction);

    let mut graph = CdbgGraph::new(31);
    let sa = add_node(&mut graph, &sa_seq);
    let c1 = add_node(&mut graph, &c1_seq);
    let c2 = add_node(&mut graph, &c2_seq);
    let c3 = add_node(&mut graph, &c3_seq);
    let sb = add_node(&mut graph, &sb_seq);
    let l = add_node(&mut graph, &l_seq);
    let r = add_node(&mut graph, &r_seq);
    let f = add_node(&mut graph, &f_seq);
    let t = add_node(&mut graph, &t_seq);
    let w = add_node(&mut graph, &w_seq);

    for (from, to) in [
        (sa, c1), (c1, c2), (c2, c3), (c3, l),
        (sb, l),
        (l, r), (r, l), // Cycle: L -> R -> L
        (l, f),
        (f, t), (f, w),
    ] {
        add_edge(&mut graph, from, to, 30);
    }

    let mut matrix = PairLinkMatrix::new();
    f23_add_pair_rows(&mut matrix, l, w, 0..16, 120, 120, 100);
    f23_add_pair_rows(&mut matrix, l, t, 16..22, 120, 120, 100);

    let stats = resolve_repeats_exspander_with_stats(
        &mut graph,
        &matrix,
        &InsertSizeDist::default(),
        ExspanderParams::default(),
    )
    .expect("repeat resolution failed");

    let seqs = active_sequences(&graph);
    let sa_tail = &sa_seq[100..220];
    let w_tail = &w_seq[100..220];
    let t_tail = &t_seq[100..220];

    assert!(
        !f23_same_contig(&seqs, sa_tail, w_tail) && !f23_same_contig(&seqs, sa_tail, t_tail),
        "multi-context cycle re-entry must not restore attribution or commit chimera S_A..W"
    );
    assert!(!graph.nodes[sa.id() as usize].is_deleted, "S_A must remain unconsumed");
    assert!(!graph.nodes[sb.id() as usize].is_deleted, "S_B must remain unconsumed");
    assert!(stats.withheld_no_unique_anchor >= 1, "targeted fork must be withheld");
}
