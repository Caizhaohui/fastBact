# Phase F9-D2 Architecture & State Semantics Review Report

## 1. Executive Summary & Design Verdict

- Target Commit: `37a7c0b`. HEAD at session start was `3484972`, not `5a85f67`. The commits after `37a7c0b` are docs-only according to their messages and the F8 report. I did not re-check this with `git diff` in this session.
- Model: `claude-opus-5-5` (Claude Code CLI, High reasoning effort)
- Provider: OmniaKey Claude Gateway. This comes from the request; I can't see the provider from inside the session.
- Design Verdict: ARCHITECTURE_CHANGE_REQUIRED. The change needed is small (Option C, about 5 lines of logic plus a rename). A full state machine (Option B) is not needed.

Findings:

- F8-P1-001 is real. I traced it by hand through the source (§5). The commit-time neighbor scan (`exspander.rs:997`) writes into the same set that the BFS checks when it dequeues a node (`exspander.rs:246`).
- The two-set split fixes that one write site, but the same failure still happens through older writers that predate F7:
  - Decision-failure writes (`:477–487`, `:560–562`, `:612`, `:619`, `:683`).
  - Validation-failure writes (`:934–940`, `:1023–1026`).
  - The frontier-protection write (`:363`), which is a cleanup write filed as a rejection.

  §5.6 gives a new counterexample (Test J). On `37a7c0b`, and also after the split, it commits a weaker endpoint while the search reports `Complete`.
- Root cause: the BFS silently drops nodes based on rejections recorded earlier in the run, possibly by other seeds. Removing a candidate is not a conservative move. It can create an accept that the full candidate set would not have produced (Lemma, Q-D).
- Recommendation (Option C):
  - Remove `rejected_nodes` from the BFS filter.
  - Rename the set to `cleanup_protected_nodes` and keep every current writer.
  - Stage 5 cleanup output stays exactly as it is in F7, so Cases A–E are unaffected.
  - Search eligibility then depends only on graph state and the current seed's state.
- How I verified: I read and searched the source only. I did not compile anything or run any tests in this session. Every trace below was done by hand from the code.

## 2. Formal State Semantics Matrix

The target contract below is Option C. Where F7 behaves differently, I say so.

### 2a. Ownership and lifecycle

| State | Representation | Owner | Writers | Readers | Lifecycle |
|---|---|---|---|---|---|
| SEARCH_REJECTED | F7: `rejected_nodes: HashSet<u32>` (physical id), shared across the whole resolver run. Target: no stored state; it is just the `break` at the failed fork. | `extend_one_direction`, plus the validation branch in `resolve_repeats_exspander_with_stats` | Edge-level: `:412`, `:416`, `:450`, `:454`. Decision failure: `:477/479`, `:485/487`, `:560/562`, `:612`, `:619`, `:683`. Path failure: `:934–940`, `:1023/1026` | F7: `:246` (BFS), `:1051` (Stage 5). Target: no search reader; the ids go to CLEANUP_PROTECTED. | F7: only grows, from the first write until return. Target: ends at the fork decision. |
| CLEANUP_PROTECTED | `cleanup_protected_nodes: HashSet<u32>` | `resolve_repeats_exspander_with_stats` | Commit scan `:977–999` (Ok branch, before deletion). Under C, also every SEARCH_REJECTED writer and `:363`. | `:1051` only | Only grows; dropped at return; never stored on the graph. |
| CONSUMED | `consumed_nodes` (`:977–980`, local to the commit). Long nodes also go into `used_in_path` (`:1017`). | Commit block `:958–1018` | Ok branch only | `:991` (scan exclusion). Via `used_in_path`/`local_used`: `:252`, `:523`, `:655`, `:849` | Lasts an instant; becomes DELETED at `:1005`. |
| DELETED | `UnitigNode.is_deleted` | Graph. Inside the resolver: commit `:1005`, Stage 5 `:1052` | `:1005`, `:1052` (upstream simplify runs outside this function) | `:246`, `:430`, `:522`, `:654`, `:729`, `:849`, `:994`, `:1042–1044` | Terminal; survives past return. |
| VISITED_IN_SEARCH | `visited: SmallVec<[u32; 8]>` per queued state | `find_candidate_paths` | `:229–230` (curr, fork), `:256` (on dequeue) | `:291` (enqueue filter), `:334` (depth helper) | One per BFS branch (cloned at `:300`); dropped when the call returns. |
| USED_UNIQUE (related) | `local_used` per seed, `used_in_path` across seeds | Seed loop | `:534`, `:540`, `:574`, `:662`; merged only on Ok (`:1017`) | `:252`, `:523`, `:655`, `:849` | `local_used` is discarded on failure, so this state is atomic. |

### 2b. Transitions and effects

| State | Allowed transitions | Invalid transitions | BFS eligibility | Stage 5 cleanup | Sequence preservation |
|---|---|---|---|---|---|
| SEARCH_REJECTED | Can become CONSUMED by a later path. F7 already allows this through the unbranched walker `:565–583`, the fork target `:527` and the dead-end push `:586`, none of which check rejection. Can coexist with CLEANUP_PROTECTED. | Must not change any later fork's candidate set. F7 breaks this. | F7: silently excluded (the defect). Target: no effect. | Shields the node. This is redundant under F7 (see Q-B). | No guarantee. |
| CLEANUP_PROTECTED | Can become CONSUMED and then DELETED by a later commit; its membership is then stale and has no effect. | Must not be read by BFS, scoring, the decision rule or the validator. Must not be removed before Stage 5. Under A, must not be written on the failure path. | None | Shields the node | Survives the resolver. Downstream curation can still drop it (§8). |
| CONSUMED | Becomes DELETED (non-circular nodes only). In practice circular nodes never enter a path: `:246`, `:430`, and seeds are non-circular. | A long node in two paths (blocked at `:732`, `:252`, `:523`). Un-consuming. | Excluded through DELETED or USED_UNIQUE | Not applicable | Contained in the validated assembled contig |
| DELETED | None (terminal) | Un-deleting. Stage 5 deleting a CLEANUP_PROTECTED node. | Skipped silently (residual R-2) | Not applicable | Kept only if the node was CONSUMED first |
| VISITED_IN_SEARCH | None outside its own path | Leaking across BFS calls or seeds | Blocks re-entering the same physical id (both strands) within one path | None | None |

## 3. Detailed Answers to Questions Q-A through Q-G

**Q-A: Does `cleanup_protected_nodes` stay out of BFS?**
Under Option A it does, with two conditions:
1. `find_candidate_paths` and `extend_one_direction` must never receive the new set. Enforce this through their signatures.
2. `protect_unmerged_frontier` (`:363`) must stop writing to `rejected_nodes`.

`:363` is a cleanup write: its doc comment at `:347–350` says so. Under F7 it also hides R65 from every later seed's BFS. If Grok moves it into `cleanup_protected_nodes`, the write happens before validation (it is called at `:385`, during exploration), and that breaks Option A's atomicity. The fix is to delete the `:385` call. It is redundant: at the step cap, R64 is the last step in `full_path`, and the commit scan reads `R64.edges_fwd` and protects R65. If validation fails, R64 stays active and R65 stays connected. By my trace, Case E and the N=66 sweep stay GREEN either way. Under Option C the question goes away, because nothing gates the BFS.

**Q-B: Should `rejected_nodes` keep its cleanup role?**
It can keep it. The role is harmless and redundant now that F7 adds the commit scan. Every exploration-time write targets a node N that was reached through an edge stored on an active node M (`curr.edges(...)` or `target.edges(...)`). At Stage 5 there are two cases:
- M is still active: N is in `connected` (`:1043–1046`).
- M was consumed by some commit: that commit's scan read M's own adjacency list and protected N (`:981–998`).

Nothing else deletes nodes inside the resolver before Stage 5. The Err-branch writes add `local_used` (seeds and nodes ≥500 bp, which Stage 5 never deletes, `:1051`) and path nodes (each reached by an edge from the step before it). So the old writes do not change the Stage 5 result. Keep them as cleanup-only, so the cleanup diff against F7 is zero. Failed candidates should not be singled out for deletion.

**Q-C: Can a node be in both states?**
Yes. Example: X is rejected at S1's failed fork (`:683`) and is also an external neighbor of S0's commit (`:997`).
- Under A, the two readers are separate. The search excludes X because of `rejected_nodes`. Cleanup shields X because it is in `protected ∪ rejected`. Nothing conflicts, but the search exclusion is still the hazard from Q-D.
- Under C there is one set, so the question doesn't arise.

**Q-D: Is cross-seed coupling valid?**
No. None of the writers records a reason that justifies a ban on the whole node across seeds:

| Writer | What actually failed | Is a node-wide ban across seeds justified? |
|---|---|---|
| `:412/416/450/454` | One edge (source→target) | No. The BFS already enqueues only `Valid` edges (`:225`, `:294`), so the bad edge is never walked. The ban only blocks other, valid edges into the same node. |
| `:477–562`, `:612–683` | Pair evidence from one anchor (`curr_dir`) | No. Another anchor's evidence is independent. |
| `:934–940`, `:1023–1026` | One junction in one specific path | No. `validate_and_splice_path` re-checks every future path. Concrete harm: `:1022` and `:934` put the failed seed and its claimed long nodes into the set. A later seed whose real partner across a repeat is that failed seed can then never score it. |

Lemma (exclusion is not conservative). Let C be the full candidate set and C′ ⊂ C the set left after exclusion. The decision rule (`:511–518`, `:643–650`) looks only at the top two candidates. There are two ways exclusion goes wrong:
- Removing the true winner promotes a weaker candidate, which can then pass the single-candidate test (this is F8-P1-001).
- Removing the runner-up turns a failed ratio test into a single-candidate accept.

So C′ can produce accepts that C would not. Only returning `Truncated` (reject) is a safe response to uncertainty. Silently skipping a node never is.

**Q-E: What happens when a protected node is later consumed?**
It becomes CONSUMED and then DELETED at `:1005`. Stage 5 only ever sets `is_deleted = true` (`:1052`), so a deleted node is never brought back. Its leftover entry in the protected set has no effect. No contradiction. The later commit's scan protects that node's own neighbors in turn.

**Q-F: Can noise get over-protected?**
Yes, within limits. The commit scan protects every active neighbor regardless of coverage, including error tips that used to be removed as orphans. These orphans also survive the simplify pass that runs after the resolver (`main.rs:436`): `clip_tips_tracked` builds predecessors only from non-deleted nodes (`simplify.rs:165–176`), so `has_strong_alternative` stays false (`:183–192`). In practice, only the `curate_contigs` length and coverage filters remove them. The retained set is at most Σ deg(consumed) and is dropped at return. I accept this cost: losing real sequence is worse than keeping some noise.

**Q-G: Is failure atomic?**
- The graph, the arena and `used_in_path` are atomic. The validator takes `&CdbgGraph` (`:705`), and mutation happens only at `:1002–1017`.
- `stats` changes, which is intended.
- `rejected_nodes` is not atomic. Exploration writes happen before validation, and the Err branch adds more. Under F7 and A these writes feed into later BFS runs, so a failed attempt changes future decisions. Under C they reach only Stage 5, where they have no effect (Q-B). Under A, `cleanup_protected_nodes` stays atomic only if `:363` is not moved into it.

## 4. Formal Correctness Invariants (Invariants 1–7)

Definition. Let Elig(u) = ¬deleted(u) ∧ ¬circular(u) ∧ u.id ≠ curr.id ∧ u.id ∉ ext ∧ ¬(len(u) ≥ 500 ∧ u.id ∈ used_unique).

| # | Invariant | F7 | A | C |
|---|---|---|---|---|
| 1 | Search eligibility. A dequeue skip happens if and only if ¬Elig(u). Eligibility depends only on graph state and the current seed's state. | FAIL (`:246`) | FAIL (`:246` still reads cross-seed history) | PASS |
| 2 | Completeness. If the search returns `Complete`, every path from the fork with depth ≤ 4, within the separation budget, through Elig nodes, ending at an endpoint with in-range evidence appears among the candidates. | FAIL (F8 case) | FAIL (Test J) | PASS, relative to Inv 1 |
| 3 | Boundary preservation. Any N ∉ P that is adjacent to some v ∈ P and active at commit time is still active at return, unless a later commit consumed it. | PASS | PASS (Stage 5 reads the union) | PASS |
| 4 | Atomicity. (a) Err leaves the graph, arena and `used_in_path` unchanged. (b) No state written during a failed attempt changes a later decision. | (a) PASS, (b) FAIL | (a) PASS, (b) FAIL | PASS |
| 5 | Twin consistency. All membership is keyed by `id()` (`cdbg.rs:19`). The scan reads both `edges_fwd` and `edges_rev`, so an incoming X→R, stored as `R.rc→X.rc`, is found. | PASS | PASS | PASS |
| 6 | Decision independence. A fork's decision is a function of (graph, matrix, dist, params, curr, used_unique, ext) only. | FAIL | PASS | PASS |
| 7 | Cleanup precision. A short non-circular u is deleted at Stage 5 if and only if ¬connected(u) ∧ u ∉ Protected. Protected ⊆ V; it is dropped at return; isolated nodes that no path touched are deleted. | PASS | PASS | PASS |

Proof sketch for Inv 2 under C. The BFS enumerates paths, not nodes: there is no global visited set, and `visited` is per branch (`:300`). The only dequeue skips are `:246–254` with the rejection term removed, which is exactly ¬Elig. Every limit sets `Truncated`: `:240–242`, `:270–272` and `:304–306`. Note that Inv 2 is stated relative to Inv 1, so the silent DELETED skip is outside it. That remaining risk is R-2 in §8.

## 5. Counterexample Forensic Audit & Production Reachability

**5.1 Arithmetic.** Both endpoints have separation `100 − 31 − 31 = 38` and span `50 + 38 + 50 = 138` = median. So every pair has Z = 0 and weight 1, and the ratio is exactly 20/5 = 4.0 ≥ 2.5, whatever the MAD. D on its own scores W/C = 1.0 ≥ 0.05 and count 5 ≥ `min_support` 5.

**5.2 Trace on `37a7c0b`:**
1. S1 (900 bp) is processed first.
2. S1 → R is a dead end (`:584–593`). The path validates and commits.
3. `R.edges_rev` holds `X.twin`, so `:997` inserts X.
4. S2 → Q. Q has two continuations, so the BFS runs (`:460–474`).
5. X is dequeued and skipped at `:246`. D is scored.
6. The search returns `Complete` with D as the only candidate, which is accepted.
7. `S2→Q→D` commits.

With L246 ignoring the cleanup set, X wins 4.0 to 1. Confirmed.

**5.3 Topology.** Junction a has inputs {S1, X} and output {R}. Junction b has input {Q} and outputs {X, D}. The edges are twin-complete and consistent with a dBG. One problem: S2 (out-degree 1) and Q (in-degree 1) are not maximally compacted, so `recompact_graph` would merge them before exSPAnder runs (`main.rs:334`, `:397`). The merged variant still reproduces the bug. The fork moves to S2 (Case B, `:597–608`), and both separations become −31, which are still equal.

**5.4 Physical realism.** The 50 bp reads and 138 bp fragments are internally consistent (d₂ = 50 ≤ len(X) = 100), but they don't match fastBact's usual 2×150 input:
- Short fragments are probably FLASH-merged in Stage 2 and produce no PE links. I did not read the merge thresholds.
- A 100 bp X cannot hold a 150 bp R2. For a realistic fixture, use X and D around 300 bp (still below 500, so they don't become seeds) and keep d₁ ≤ len(u₁), d₂ ≤ len(u₂).

The same realism issue exists in Case F (`add_link` gives d₂ = 120 > len(X) = 100, `mapper.rs:67–68`). `score_candidate_path` doesn't check this, so the code paths are still exercised correctly.

**5.5 Reachability preconditions and unstated assumptions:**
- S1 is processed before S2. Seeds are sorted by length (`:843`).
- X is adjacent to a node S1 consumes, and S1 does not consume X itself.
- X is within S2's BFS limits (depth ≤ 4, separation ≤ `max_insert`).
- D has at least `min_support` pairs. This sits exactly on the boundary: with 4 pairs the candidate list is empty and the result is under-assembly, not misassembly. Even with no competitor, F7 still makes the result worse than F5B, because S2 stays unresolved (Test B).
- D survives tip clipping before the resolver runs. With k = 31, len(D) = 100 ≥ 2k = 62, so D is not a tip candidate (`simplify.rs:148`). At larger k1 values D becomes a candidate and survives only if cov(D) > 2 and cov(Q) < 3·cov(D) (`:187`).
- The skip at `:246` ignores node length, so Codex's "short neighbor" understates how far the bug reaches.

Biologically, this is a short shared segment with a minority alternative branch (strain heterogeneity, an IS variant, or mis-mapped pairs at repeat boundaries). A 20%-minority branch is plausible.

**5.6 New counterexample through older writers (Test J).** This one stays RED on F7 and after Option A:
- Nodes: S1′ (900) → {X, Y}; S2′ (700) → {M, D}; M (100) → {X, Y}. Junction j2 has inputs {S1′, M} and outputs {X, Y}, consistent with a dBG.
- S1′ has only 3 pairs to Y (below `min_support`), so its fork decision fails and `:683` inserts X and Y.
- S2′'s BFS goes M → X (separation 38) and to D (separation −31). X has 20 pairs and D has 5, both centered.
- On F7 and A, X is skipped, D is the only candidate, and `S2′→D` commits (misassembly).
- On C, X beats D 4.0 to 1 and `S2′→M→X` commits.

## 6. Architectural Decision Evaluation (Option A vs B vs C)

| Criterion | A (split, plus removing `:385`) | B (state enum) | C (cleanup-only set, no BFS gate) |
|---|---|---|---|
| Correctness | Inv 1, 2, 4 FAIL (Test J, Test G2) | Depends on its transition rules. An enum alone doesn't fix Inv 2, and Q-C shows the states overlap, so it would need bitflags. | All 7 PASS |
| LoC | about 10 | 150–300 | about 4 of logic plus a mechanical rename |
| API | none (all private) | internal refactor | none (`find_candidate_paths` is private) |
| Runtime / memory | +1 HashSet | O(n) bytes | Same memory. The BFS may process more states, capped at 256 per call (`:178`). |
| Verification | Low | High | Medium. Search behavior changes in multi-seed graphs, so it needs the full suite and the integration genome. |
| Rollback | Trivial | Moderate | Trivial (one revert) |

I considered a variant where reaching a rejected node returns `Truncated` instead of skipping it, and rejected it: it breaks Case F scenario 2, where the rejected R1 sits under X in S2's search.

I checked Option C against the existing tests I read, by trace:
- Case F scenario 2 is unchanged: R1 has no support from X, so it yields no candidate either way.
- g2 and g6 are unchanged: they use unbranched walks, or the rejected node is never consulted.
- Cases A–E are unchanged: the cleanup diff is zero.

I did not read `graph_atomicity`, `edge_gate`, `directed_edge_identity`, `junction_overlap`, `pe_evidence`, `integration_test` or the in-module tests. Test I covers them.

Option C sketch:

```rust
// find_candidate_paths: remove the `rejected_nodes` parameter and its two call sites (:472, :606).
// Eligibility now depends only on graph state and current-seed state (Invariant 1).
if u_node.is_deleted || u_node.is_circular {
    continue;
}
```

Also rename `rejected_nodes` → `cleanup_protected_nodes` everywhere and update Contract 2 at `:1037`. Keep every writer.

Recommendation: Option C, in two commits:
1. Tests only: A, B, J and G2. All should be RED on `37a7c0b`.
2. The fix.

If an existing test fails under C, first check whether that test encodes the old exclusion behavior. Fall back to A only if the owner signs off on Test J as a known pre-existing P1.

## 7. RED → GREEN Regression Test Contract for Grok

Every test must assert its own preconditions (scores, ratios, W/C, edge validity) before it runs the resolver. Use explicit `add_evidence` with distinct `pair_id`s, and keep d₁ ≤ len(u₁) and d₂ ≤ len(u₂).

| Test | Fixture | Oracle | On `37a7c0b` | After A | After C |
|---|---|---|---|---|---|
| A1 | Codex topology. X and D are 300 bp. 20 vs 5 pairs at span 450 (d₁ = 200, d₂ = 212, default dist). | Contig `S2+Q+X` is present, `S2+Q+D` is absent, D is still active | RED | GREEN | GREEN |
| A2 | Same with S2+Q recompacted (Case B fork, separation −31) | Same | RED | GREEN | GREEN |
| B | A1 with D at 0 pairs. Variant B2: the supported endpoint is one hop past the protected X. | The S2 contig through X exists | RED | GREEN | GREEN |
| C | The existing Cases A–G, plus a check that X stays alive when a later seed doesn't consume it | Existing assertions unchanged | GREEN | GREEN | GREEN |
| D | Fork failures (support, tie, truncation) and Ambiguous/Invalid edges | No committed path contains a rejected edge or decision. A later seed with the same weak evidence also fails. | GREEN | GREEN | GREEN |
| E | X is rejected at S1's fork and also protected by S0's commit | X survives Stage 5. A: X is ineligible (document this). C: X is eligible. | n/a | GREEN | GREEN |
| F | A1 mirrored: X stored as RC (Case D style), and S2's fork reached through the backward extension | Same as A | RED | GREEN | GREEN |
| G1 | Validation failure fingerprint (Case F, g1) | Graph and `used_in_path` unchanged | GREEN | GREEN | GREEN |
| G2 | The failed path S1→X shares X with S2's fork (Q→M→X vs D). Differential: G vs G minus S1. | S2's committed contigs are identical in both runs | RED | RED | GREEN |
| H | Three seeds: S1 protects X, S2 consumes X, S3 commits nearby | Conservation: every short node is either consumed (its unique bases appear in exactly one contig) or still active. No deleted node with unique bases is missing from the output. | GREEN | GREEN | GREEN |
| I | Full `cargo test --release`; record commands, exit codes and logs | 138 reported plus new tests, 0 failed; diff touches only `exspander.rs` and new test files | n/a | GREEN | GREEN |
| J | §5.6 (fork decision failure → later seed) | `S2′+M+X` is present, `S2′+D` is absent. Also run with D at 6–8 pairs so the result doesn't depend on the boundary. | RED | RED | GREEN |

As written in the brief, Test D ("genuine failed candidates remain excluded from BFS") conflicts with Option C. I restated it as "a rejected decision is never committed". Literal exclusion is a property only Option A has, and the Lemma in Q-D shows it isn't a safety property.

Carried over from F8: equality and one-over-budget checks (`sep == max`, `max + 1`) for both `:297` and `:339` are still outstanding.

## 8. Potential Secondary Risks & Residual Vulnerabilities

- **R-1 (P1, existed before F7):** the decision-failure coupling in Test J and the validation-failure coupling in G2. Option C fixes these; Option A does not.
- **R-2 (P2, existed before F7):** if an earlier seed already consumed a supported endpoint, the BFS skips it silently at `:246` or `:252`. The search can then accept a weaker alternative and still report `Complete`. This comes from the walker taking short repeats regardless of in-degree (`:584–593`). Suggested fix for F10: when a skipped node has pair support from `curr`, return `Truncated`.
- **R-3 (P2, existed before F7):** commits don't rewire edges to the assembled node, so a protected neighbor keeps its sequence but loses its link (noted by Codex).
- **R-4:** protected orphans survive the post-resolver simplify pass (Q-F), and dangling edges to deleted nodes inflate the out-degree counts in `clip_tips_tracked`. I haven't checked whether `write_outputs` puts these nodes in the GFA.
- **R-5:** because the visited set uses physical ids, hairpin and inverted-repeat paths can't be explored. This is conservative and existed before F7.
- **R-6 (Option C cost):** more BFS states means more `MAX_DFS_STATES` truncations, which can reduce contiguity. Measure it on the integration dataset.
- **R-7:** fixture realism (§5.4), in both Case F and the Codex parameterization.
- **R-8:** the commit scan can't see malformed, non-twin incoming edges (noted by Codex; existed before F7).
- **Unverified:** no part of this report was compiled or executed. The 138-test count comes from the F7 report and has not been reproduced.
