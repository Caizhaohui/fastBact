# Phase F9-D2 Round 2 Adversarial Design Audit Report

Scope and limits: The Round 1 report (`phase_f9/opus55_round1_architecture.md`) is not in this worktree. `phase_f9/` contains only the prompt file. I therefore judged Option C, Test J and Test G2 from the summary in your prompt. I did inspect the code: `src/resolve/exspander.rs`, `src/resolve/mapper.rs`, `src/graph/cdbg.rs`, `src/graph/simplify.rs`, `src/main.rs`, `tests/test_c0_02r_incoming_side_branch_gate.rs`, and the Codex F8 report. I did not build anything or run tests; every counterexample below comes from tracing the source. Your prompt promises 10 vectors but lists 8, so I added vectors 9 and 10 myself.

## 1. Executive Summary & Adversarial Verdict

- Target Evaluation: Round 1 Option C vs Option A
- Model: `claude-opus-5-5` (Claude Code CLI, High reasoning effort)
- Verdict: **REFINED_OPTION_REQUIRED**

Option C is right about the part it changes. Removing `rejected_nodes` from BFS is safe for termination, keeps reverse-complement symmetry, and is strictly better than Option A. Option A still lets failed fork decisions, truncated searches and validation errors ban nodes for later seeds (lines 477–487, 560–562, 612–619, 682–683, 1022–1027).

Option C is not enough, though. Round 1 defines BFS eligibility as "physical validity (`!is_deleted && !is_circular`)". That definition carries the same flaw as F8-P1-001. A node consumed by an earlier seed is tombstoned (lines 1002–1006), but its edges remain and its paired-end evidence is still in the matrix. BFS silently skips it at line 246 (`is_deleted`) and line 252 (`used_unique`) and still returns `Complete`. That hides a competitor, and the lenient single-candidate rule can then commit a weaker alternative (lines 513–515, 645–647). This is Risk R-2, it has the same failure class as F8-P1-001, and neither A nor C fixes it. Counterexample CE-3 below is the Codex 6-node topology with different evidence. It fails on HEAD, on Option A and on Option C.

Recommended Option D has three parts:
1. Option C.
2. Nodes consumed by earlier seeds are traversed as non-committable "shadow" candidates. They can block a fork decision but can never be spliced.
3. Make the separation structural by removing the `rejected_nodes` parameter from `find_candidate_paths`, so the separation does not rely on naming.

## 2. Attack Vector Analysis (Vectors 1–10)

### V1. BFS explosion and cycles under C: SAFE, with one unmeasured cost (P3)

- Hard bounds hold whatever is in the rejected set. `states_visited` is counted before any filter (239–243), giving at most 256 dequeues, at most 32 candidates, and depth at most 4 (281). The queue is iterative, so there is no recursion and no stack exhaustion. `visited` holds at most 6 IDs (curr, fork, plus 4 steps), so it never spills out of `SmallVec<[u32; 8]>`.
- Cycles are blocked by physical ID (228–230, 256, 291). That covers `A→R→A`, a self-loop `R→R`, and a hairpin `R→R.rc`, because `twin()` keeps `id()` the same (cdbg.rs:19–31). `depth_limit_hides_continuation` uses the same `visited` (334).
- Forward and reverse traversal are the same code: `edges(u.is_rc())` (282), the twin-complete graph, and mapper evidence written in both orientations (mapper.rs:231–252).
- Real cost: nodes that were rejected used to be dequeued and dropped. Now they are expanded and enqueue children, so the 256-state cap is hit more often. That produces more `Truncated` results and therefore more fragmentation, never wrong joins. Nothing counts this today (see V10).

### V2. Candidate competition and ratio distortion: conservative by construction (P2, pre-existing)

- Adding candidates to one decision can do only two things. It can make a candidate with higher weight the best one, which is the intended fix. Or it can turn the lenient single-candidate rule into the stricter ratio test, which can only move pass to fail. The one fail-to-pass transition happens when a newly visible candidate beats a previously ambiguous pair by at least 2.5×. That is correct, because each seed's paired-end evidence is independent of why some other seed's fork failed.
- Stray evidence is limited at the source. The mapper only indexes k-mers that are unique within the graph (mapper.rs:134). A stray candidate with 1–2 pairs, weight 2 or less, blocks only when the best candidate's weight is about 5 or less. That is right at `min_support = 5` (exspander.rs:24), so blocking there is statistically correct.
- Pre-existing defect that C exposes a little more: candidates are compared by endpoint, not by fork exit. If `Q→X` has 20 pairs and `Q→X→Y` has 15, the ratio is 1.33 and the decision fails, even though both point to exit X. The result is fragmentation, not misassembly. Measure it first (V10) and fix it later.

### V3. Cleanup over-protection: unchanged by C, bounded downstream (P3)

- Option C uses the same writers as F7 HEAD, so how much is protected does not change. Codex A5 already accepted that.
- Adversarial detail: a protected neighbor whose only edges point at tombstones also survives the tip clipping after exSPAnder (main.rs:436). `clip_tips_tracked` builds predecessors from live nodes only (simplify.rs:166–176), so `has_strong_alternative` stays false. These nodes reach `curate_contigs`, where `min_len`, `min_cov` and substring/RC deduplication limit them (Codex A7; I did not re-read `output.rs`).
- Tips with strong-coverage neighbors were already clipped before Stage 5 (main.rs:334, 397). What survives is sequence that upstream simplification did not class as a tip. I could not verify the effect on R1-001/R1-003; I ran no benchmarks.

### V4. Validation rollback: C improves atomicity

- On `Err` (1019–1029), Option C writes only to the cleanup set. The graph is untouched (the validator borrows immutably), and under C the search state is untouched as well. That means Invariant 4 holds for search eligibility, which Option A does not achieve.
- Protecting the nodes of a failed path is what atomicity requires: deleting orphaned path nodes after a failed commit would itself be a mutation. Inserting `local_used`, which includes `used_in_path` (856), is a no-op for nodes already deleted.

### V5. Tombstone re-visitation (R-2): CONFIRMED, P1 latent, pre-existing, not fixed by A or C

- No deleted node can ever be spliced: the validator rejects it (729) and so does `all_steps_valid` (520–524, 652–656).
- The hole is in competition. Line 246 `continue`s on `is_deleted` and leaves `search_status` unchanged, so a supported route through a consumed node vanishes and the search reports `Complete`. Line 252 does the same for unique nodes already used.
- Why this matters in production: nodes shorter than 500 bp may be reused within one path (523, 655), but every commit deletes all its constituents (1002–1006). Collapsed short repeats therefore become single-use across seeds. These repeats are also where the mapper puts evidence from both copies, because uniqueness is checked across the graph, not the genome (mapper.rs:89–136). The first seed consumes the repeat. The second seed's BFS then cannot see it, and any alternative exit with at least 5 in-range pairs and average weight of at least 0.05 commits without a ratio test.
- How often this happens on real isolates is unknown. It needs a decoy exit with at least 5 consistent pairs. It has to be measured (V10), not argued away.

### V6. Twin symmetry: preserved, and C removes one asymmetry

- Every exclusion set is keyed by physical ID, which is strand-agnostic. Evidence is mirrored, so `(u1,u2)` and `(u2.twin(),u1.twin())` hold swapped distances (mapper.rs:231–252).
- On HEAD, the forward extension writes `rejected_nodes` before the backward extension runs (860–863). That means flipping which strand of a seed is stored can change the outcome. Under C, the only coupling left between directions is `local_used` (unique-node reuse), which is order-dependent but semantically required.

### V7. Test J and the Codex 6-node case: genuine, reachable

- In the Codex topology, Q has out-degree 2 and R has in-degree 2, so X cannot be merged away. `S2→Q` may fuse during recompaction, which just renames the case. I did not read `recompact_graph`.
- At production k, the dead ends R and D (169 bp and 100 bp) can fall under 2k and become tip candidates. They are clipped only when the predecessor's coverage is at least 3× theirs or their own coverage is 2 or less (simplify.rs:148–192), so with comparable coverage they survive.
- The current fixtures call the resolver directly at k=31 with no simplification passes. A pipeline-level fixture would show reachability directly.

### V8. A vs C vs D: D

Covered in section 4.

### V9 (added). Malformed reciprocal edges under C (P3)

On HEAD, when `bwd_failed` triggers (931–941), the seed and all its extension nodes are banned. Under C, a later seed can route through them. BFS does not check reciprocal edges, but the validator does (790–798), so the later seed's whole path is rejected instead of committing a shorter alternative. This is conservative and can only happen in non-canonical graphs.

### V10 (added). Observability gap (P2, process)

`JunctionStats` (31–43) has no counters for truncated searches or for consumed-node conflicts. As a result, neither the cost of C (V1, V2) nor the exposure to R-2 (V5) can be measured on R1-001 or R1-003.

## 3. Stress-Testing Counterexamples

All use the Codex topology: `S1(900)→R(169)`, `X(100)→R`, `S2(700)→Q(100)→{X, D(100)}`, overlap 31, twin-complete, coverage 25. S1 commits `S1+R` first, then S2 searches the fork at Q with `init_sep = 69`, giving sep(X) = sep(D) = 38 and sep(R) = 38 + 100 − 31 = 107.

| Case | S2 evidence | HEAD | Option A | Option C | Option D |
|---|---|---|---|---|---|
| CE-1 (Codex F8) | X: 20, D: 5 | commits D ✗ | X ✓ | X ✓ | X ✓ |
| CE-2 (cross-seed fork failure)¹ | X: 20, D: 5 | commits D ✗ | commits D ✗ | X ✓ | X ✓ |
| CE-3 (tombstone R-2)² | R: 20, X: 0, D: 5 | commits D ✗ | commits D ✗ | **commits D ✗** | declines ✓ |
| CE-4 (same-exit shadow)³ | X: 20, R: 18, D: 2 | no commit | X ✓ | X ✓ | X ✓ |

1. CE-2: S1 instead reaches the fork `P→{X,Y}` with 6 vs 5 pairs. The ratio is 1.2, so the decision fails and lines 560–562 write P, X and Y. S1's full path has length 1, so nothing is spliced. X now also has the in-edge `Q→X`.
2. CE-3 trace under C: X is dequeued but has no support. It expands to R, which is skipped silently at line 246. D is the only candidate, the search returns `Complete`, and D passes the single-candidate rule. Centered evidence: R uses d1=50, d2=50, so span 207. D uses d1=119, d2=50, so span 207. Median 207, MAD 40, insert range 100–400.
3. CE-4 is a guard against D blocking too much. R is a shadow under the same exit (X) as the winner, so it must not block.

Termination check: a directed cycle `A→R→A` plus a hairpin `R→R.rc` terminates under C and D, because of the physical-ID `visited` check (291).

## 4. Final Architectural Recommendation

Adopt Option D. It has the following property: for any single fork decision with the same prior state, D commits only what C would commit, minus the decisions that consumed evidence contradicts. Shadow candidates are filtered out before ranking and are also stopped by `all_steps_valid`, so they can never be spliced.

**D1: decouple structurally (this is Option C).** Rename `rejected_nodes` to `cleanup_protected_nodes` in `extend_one_direction`, `protect_unmerged_frontier` and `resolve_repeats_exspander_with_stats`, at Stage 5 (1051), and in the Contract 2 comment (1037–1038). Remove the parameter from `find_candidate_paths`, so the compiler stops the set from ever being read during search again.

**D2: track consumed nodes per resolver call.** Use an explicit set rather than `is_deleted`. Nodes pruned upstream by tip clipping or bubble popping keep being skipped, and only nodes spliced by this resolver become shadows.

```rust
// resolve_repeats_exspander_with_stats
let mut consumed_nodes: HashSet<u32> = HashSet::new();
// ... in Ok(assembled_seq), replacing the per-commit set at 977–980:
consumed_nodes.extend(full_path.iter().map(|step| step.dir.id()));
```

Pass `&consumed_nodes` through `extend_one_direction` into `find_candidate_paths`, in the slot that `rejected_nodes` frees.

**D3: shadow traversal in `find_candidate_paths`.**

```rust
struct CandidatePath {
    // existing fields ...
    /// Route enters a node spliced by an earlier seed or a unique node this
    /// seed already uses. Never committed; its evidence still competes.
    through_consumed: bool,
}

// queue entries gain a bool; initial pushes use `false`
while let Some((u, sep_u, steps, mut visited, through_consumed)) = queue.pop_front() {
    // ... state cap unchanged ...
    let u_node = &graph.nodes[u.id() as usize];
    if u_node.is_circular || (u_node.is_deleted && !consumed_nodes.contains(&u.id())) {
        continue;
    }
    if u.id() == curr_dir.id() || ext_nodes.contains(&u.id()) {
        continue;
    }
    let through_consumed = through_consumed
        || u_node.is_deleted
        || (u_node.seq_len >= 500 && used_unique.contains(&u.id()));
    visited.push(u.id());
    // scoring unchanged; set `through_consumed` on the pushed CandidatePath
    // expansion unchanged (tombstone edges remain); children inherit the flag
}
```

**D4: shadow gate at both decision sites** (around 484–518 and 617–650). Put this after the `Truncated` check and before the empty check:

```rust
let (shadow_paths, mut candidate_paths): (Vec<_>, Vec<_>) =
    candidate_paths.into_iter().partition(|c| c.through_consumed);
// existing empty-check, sort, best/second unchanged ...
let is_decision_valid = !shadow_blocks_decision(best, &shadow_paths, params.decision_threshold)
    && /* existing expression */;
```

```rust
/// True when evidence for a consumed route through a different fork exit is
/// too strong for `best` to win by `threshold`.
fn shadow_blocks_decision(best: &CandidatePath, shadows: &[CandidatePath], threshold: f64) -> bool {
    shadows.iter().any(|s| {
        s.steps[0].0 != best.steps[0].0 && best.total_weight < threshold * s.total_weight
    })
}
```

**D5: counters.** Add `truncated_searches` and `consumed_conflicts` to `JunctionStats`. Check first whether any tests build `JunctionStats` as a struct literal, because adding fields would break them. If they do, use `tracing::debug!` instead.

Deferred: grouping candidates by fork exit (V2). Only do this once the D5 counters show the fragmentation cost.

Known trade-off: under D, the depth-cap continuation check (304) also fires on shadow branches. That is deliberately conservative, and its cost should be measured with D5.

## 5. Verification & Testing Requirements for Grok

1. **T1 (CE-1):** RED on HEAD, GREEN on D. Assert the exact `S2+Q+X` splice, and that no contig contains D's unique tail joined to S2.
2. **T2 (CE-2):** RED on HEAD and on Option A, GREEN on D. This is the gate that rules out Option A.
3. **T3 (CE-3):** RED on HEAD, A and C; GREEN on D. Assert no splice at Q, that S2, Q, X and D are still active, and that D's unique sequence is present in the raw curation output.
4. **T4 (CE-4):** RED on HEAD, GREEN on D. Guards against D blocking too much.
5. **T5:** a bubble node popped upstream sits next to the fork and must not become a shadow. It needs a deleted node that is not in `consumed_nodes`, and the result must be identical to C.
6. **T6:** Cases A–G in `test_c0_02r_incoming_side_branch_gate.rs` stay GREEN with unchanged assertions.
7. **T7 (atomicity):** extend Case F scenario 2 so the later seed's only supported exit is a node of the failed path (R1). RED on HEAD, GREEN on D. Also keep asserting that the graph fingerprint is unchanged.
8. **T8 (twin metamorphic):** for T1–T4, build the mirror fixture: reverse-complement every stored sequence and mirror the edges and evidence. The output contig set must be identical up to reverse complement.
9. **T9 (termination):** a directed cycle `A→R→A`, a hairpin `R→R.rc`, and a dense 8×8 fork with every node pre-protected. Each must finish without panicking, and every fork that reports `Truncated` must leave the graph fingerprint unchanged.
10. **T10:** budget boundary tests at `sep_next == max_separation` and `max + 1`, for both BFS and the depth helper. This is still owed from F8 transition item 4.
11. **T11:** a full rebuild plus the complete release suite (all 138 tests), with commands, exit codes and output saved. Then run R1-001 and R1-003 against HEAD and report contig count, N50, misassemblies, and the D5 counters. A rise in `consumed_conflicts` measures real-data exposure to R-2 directly.

I haven't saved this report anywhere, since I only had read access in this session. If you want it next to the Round 1 prompt, the natural location is `docs/development/evidence/C0_02R/phase_f9/`.
