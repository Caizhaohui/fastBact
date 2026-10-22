# fastBact — Phase F18 Independent Codex Audit Report

**Audit ID**: AUDIT-20261022-013  
**Auditor**: Codex CLI (Model: `gpt-6.1-sol`, Reasoning Effort: `medium`, Sandbox: `read-only`)  
**Target Commit**: `c37ef6e03475230dc7f8596e0fc6788ba4fc6355`  
**Baseline Commit**: `4d08b380d573b568e67498ef71cff632f8a0ef67`  
**Target Branch**: `feature/c0-02r-phase-f17-grok`  
**Target Worktree**: `fastBact-phase-f17-grok-worktree`

## 1. Executive Summary

**Overall verdict: CONDITIONAL PASS. C0-02R remains NOT ACCEPTED.**

No confirmed P0 or P1 production correctness blocker was found. Under the production graph and scoring contracts, the evidence floor β is sound, repeat exploration conservatively accounts for unfinished work, and both candidate decision sites reject `SearchStatus::Truncated` before applying acceptance thresholds.

Five P2 findings prevent an unconditional PASS:

1. The repaired oracle still accepts inconsistent reciprocal edges.
2. The oracle does not fully enforce long-node uniqueness at search initialization.
3. T30 withholds at an earlier S-anchored search, bypassing its intended downstream decision.
4. T37 does not exercise β pruning or its advertised two-sided boundary.
5. T39 never approaches the queue-capacity boundary.

Independent verification performed:

| Verification | Result |
|---|---|
| HEAD and branch | Match requested target |
| Production diff | `exspander.rs`: +385 / −275 |
| State-semantics release binary | 42 passed |
| Oracle release binary | 6 passed |
| Complete available release-test inventory | 181 tests: 175 passed; 6 blocked by read-only temporary-file creation |
| Stored S0-001/S0-003 output hashes | Match recorded baseline hashes |
| Stored contig sequence identity | Independently confirmed exact reverse-complement reference substrings |

These were **existing release binaries**, not a fresh compilation. No repository files were modified. The pre-existing untracked `phase_f18` directory was preserved.

The supplied F16 audit and Opus architecture-review files are absent from this target checkout. Their conclusions were not used as correctness evidence.

## 2. Mathematical Soundness Proof of Evidence Floor β(curr)

### Evidence domain

[Source: separation_beta](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f17-grok-worktree/src/resolve/exspander.rs:225)

The implementation computes:

\[
E_c=\bigcup_{(c,t)\in\text{matrix.links}}\text{matrix.links}[(c,t)]
\]

\[
\beta(c)=
\begin{cases}
\min_{e\in E_c}(e.d_1+e.d_2),&E_c\ne\varnothing\\
0,&E_c=\varnothing.
\end{cases}
\]

It examines every evidence row whose **first directed key equals `curr_dir`**. There is no filtering by:

- endpoint reachability or deletion;
- support threshold;
- insert window;
- read length;
- pair-ID uniqueness;
- shadow status.

Consequently, no evidence row available to `score_candidate_path(matrix, dist, curr_dir, endpoint, sep)` is omitted from the β domain.

Unreachable or duplicate evidence can lower β and weaken pruning. It cannot inflate β above the offset of an evidence row subsequently used to score a candidate.

There is no β cache. The matrix is immutably borrowed throughout the search. Edge-target deduplication does not alter evidence.

### Coordinate consistency

`score_candidate_path`, at `exspander.rs:161`, computes:

\[
\operatorname{span}_e(v)=d_1(e)+s(v)+d_2(e).
\]

The mapper records orientation-adjusted endpoint distances at `mapper.rs:213–250`. With coordinates expressed in the directed traversal frame:

\[
d_1=|L_1|-p_1,\qquad d_2=p_2+r_2.
\]

Here \(p_2\) is the inferred start of the complete second-read alignment in that frame. The mapper reconstructs the relevant outer endpoint using the mapped k-mer position and read offset; `pos2` alone is not \(p_2\).

β and scoring consume the **same stored distances**. Neither reconstructs coordinates independently from `r1_len` or `r2_len`. Thus, the lower-bound proof does not depend on an additional read-length approximation.

Twin evidence reverses the directed endpoints and swaps the two distances. This preserves their sum. Different directed anchors can nevertheless have different β values because they have different evidence sets.

Consumed-node traversal uses the original graph nodes and the same separation arithmetic. Shadow status introduces no coordinate adjustment.

### Separation invariant

The implementation requires an initialization distinction.

For the first edge from the fork:

\[
s(v)=\operatorname{init\_sep}-\operatorname{overlap}(fork,v).
\]

For subsequent edges:

\[
s(v)=s(u)+\operatorname{len}(u)-\operatorname{overlap}(u,v).
\]

For Case B, `init_sep = 0`, so an adjacent endpoint can have separation \(-\operatorname{overlap}\). Separation is therefore **not monotonically non-decreasing relative to the initial zero value**.

After initialization, it is non-decreasing:

\[
s(v)-s(u)=\operatorname{len}(u)-\operatorname{overlap}(u,v)\ge0.
\]

`resolve_edge_between`, at `exspander.rs:101`, rejects zero overlaps and overlaps exceeding either actual stranded sequence length. Production constructors and sequence-merging code maintain `seq_len` consistently with the stored sequence.

Accordingly, the monotonicity claim applies to admitted transitions after initialization on valid production graphs. A public, manually malformed node whose `seq_len` understates its arena sequence could violate this premise; that is not a valid production graph.

### Lower-bound proof

Let \(x\) be the node reached by the hop being considered for pruning, with post-hop separation \(s_x\). Let \(v\) be a descendant legal candidate endpoint, and let \(e\) be any evidence row that can score it.

By the evidence-domain definition:

\[
d_1(e)+d_2(e)\ge\beta(c).
\]

By non-negative subsequent separation increments:

\[
s(v)\ge s_x.
\]

Therefore:

\[
\operatorname{span}_e(v)
=s(v)+d_1(e)+d_2(e)
\ge s_x+\beta(c).
\]

Hence:

\[
s_x+\beta(c)>\operatorname{max\_insert}
\Longrightarrow
\operatorname{span}_e(v)>\operatorname{max\_insert}.
\]

No descendant evidence row can score in-window.

**The implemented β prune is sound.** The correct bound uses the **post-hop** separation `sep_v`, as the code does, rather than the predecessor separation before initialization.

### Other `SafePrune` cases

`classify_hop` also returns `SafePrune` for:

- raw separation exceeding `max_insert`;
- reuse of a long physical node already in the current path/prefix.

The raw bound is sound because evidence offsets are non-negative. Long-node reuse is structural rejection: no completion through that reuse is legal under the splicer’s uniqueness contract. Its safety does not arise from β.

Thus, “every `SafePrune` is an evidence-floor proof” would be an inaccurate description of the actual enum semantics.

### Boundaries and integer safety

| Case | Source conclusion |
|---|---|
| `sep_v + beta == max_insert` | Admitted by the β test |
| `sep_v + beta == max_insert + 1` | Pruned on an eligible repeat re-entry |
| No anchored evidence | β = 0; no candidate can score from that anchor |
| Reads shorter than mapping k | Mapper skips them |
| Synthetic short-read evidence | No fixed 31-bp assumption in β or scoring |
| Single-base nodes | Only sequence-valid overlaps can pass; a valid full-length overlap gives zero subsequent increment |
| Distance addition | `u32` operands widen before addition; bounded search depth keeps sums far below `i64::MAX` |

Production insert distributions are also far below the `usize`-to-`i64` conversion boundary. Arbitrary externally constructed distributions exceeding `i64::MAX` are outside this proof.

### Counterexample search outcome

I traced short self-loops, alternating cycles, twin reuse, consumed-node paths, zero-increment transitions, and equality boundaries. No valid production counterexample was found where β removes an in-window strong competitor, the search returns `Complete`, and a weak decoy commits.

## 3. classify_hop Single Authority & Call Graph Tracing

[Source: classify_hop](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f17-grok-worktree/src/resolve/exspander.rs:250)

The call graph is:

```text
find_candidate_paths
 ├─ classify_hop: initial fork edges
 ├─ classify_hop: subsequent edges
 ├─ pruned_repeat_hides_supported_endpoint
 │   ├─ classify_hop: probe entry
 │   ├─ classify_hop: subsequent probe edges
 │   └─ depth_limit_hides_continuation
 │       └─ classify_hop
 └─ depth_limit_hides_continuation
     └─ classify_hop
```

All three walkers delegate edge-extension admission to this function. No surviving walker-local `ext_nodes` or `visited` exclusion silently removes short-node routes.

The main walker’s physical-anchor check changes endpoint treatment: it prevents candidate emission, scores self evidence, and continues traversal.

| Class | Main BFS | Repeat probe | Depth check |
|---|---|---|---|
| `Live` | Enqueue and score when dequeued | Enqueue and score | Unfinished continuation |
| `Probe` | Invoke repeat probe | Enqueue with bounded depth | Unfinished continuation |
| `Unresolved` | Set `Truncated` | Return `true` | Return `true` |
| `SafePrune` | Discard | Discard | Discard |
| `Inadmissible` | Discard | Discard | Discard |

Qualifications:

- Long-node reuse is labeled `SafePrune`, not `Inadmissible`.
- `Inadmissible` includes unresolved-by-contract edges such as tied edge selection, invalid sequence overlaps, and upstream-deleted nodes.
- The main depth check deliberately supplies β = 0; the probe depth check supplies the actual β. They use one classifier with different conservative policies.
- The initial caller ignores `Probe`, but `step_depth == 0` makes that variant unreachable there.

Shadow inheritance is monotonic:

\[
shadow_v=shadow_u\lor consumed(v)\lor previouslyUsedUnique(v).
\]

Both decision sites partition shadows before choosing a committable candidate. Their step checks reject deleted nodes and previously used long nodes, and final validation rejects deleted constituents again.

No unsafe walker `continue` was identified. Probe depth exhaustion continues only after checking for unfinished admissible successors.

## 4. SearchStatus Complete / Truncated Correctness Ledger

| Event | Handling | Assessment |
|---|---|---|
| Reachable circular node within admitted budget | `Unresolved` propagates outward | Sound |
| Main state count exceeds 256 | `Truncated` | Sound |
| Probe state count exceeds 256 | Returns `true`; outer search truncates | Sound |
| Candidate count reaches 32 | `Truncated` | Conservative |
| Queue cannot accept a state | `Truncated`, or probe returns `true` | Sound |
| Depth cap with `Live`, `Probe`, or `Unresolved` successor | Search truncates | Sound |
| Probe finds supported endpoint | Returns `true`; outer search truncates | Conservative |
| Revisited physical anchor has in-window support | No candidate emitted; search truncates | Sound |
| All remaining branches rejected structurally or by a proven bound | `Complete` permitted | Sound |

A repeat probe returns `false` only after exhausting its accounted search without supported endpoints or unresolved continuation. Budget exhaustion does not return `false`.

The outer BFS cannot ignore probe incompleteness: a `true` result sets `Truncated` and exits the BFS.

Both decision sites check status first, at `exspander.rs:823` and `:963`. The single-candidate average-weight threshold of 0.05 cannot bypass truncation.

Pair-ID deduplication remains candidate-local and occurs after insert-window validation. An out-of-range duplicate does not prevent a later in-range row from contributing.

The examined historical defect mechanisms are structurally addressed:

- **F12-P1-001:** subsequent short repeat transitions remain explorable inside the probe.
- **F14-P1-001:** short nodes in the accumulated extension are not categorically excluded.
- **F16-P1-001:** short current-anchor re-entry is traversable; it is not an unconditional dequeue rejection.

This establishes completeness relative to the current selected-edge, legal-path, and scoring contracts. It does not establish biological completeness beyond the evidence represented in `PairLinkMatrix`.

## 5. Self-Anchor Handling & T27 / T33 / T34 / T35 Traces

`find_candidate_paths`, at `exspander.rs:391`, excludes both orientations of the current physical node from candidate emission. It still explores their continuations.

The pre-search self-evidence guard at `:728` abstains when `(curr_dir, curr_dir)` scores at separation zero. That branch does **not literally return `SearchStatus::Truncated`**; it stops extension before invoking the search. Its effect is conservative withholding.

### T27: legitimate prefix reuse

All relevant short nodes have length 100; overlap is 31.

First decision:

\[
s(C)=100-31-31=38,\qquad span(S,C)=80+38+82=200.
\]

After selecting C, extension is `[P,C]`.

At the C-anchored Case A search through P:

\[
\beta(C)=80+82=162.
\]

The P→X exit has:

\[
s(X)=38,\qquad span(C,X)=200.
\]

The alternative P→C reaches C at separation 38 without self evidence. Its next P re-entry would have:

\[
s(P)=38+100-31=107,
\]

\[
107+162=269>250.
\]

That continuation is safely pruned. The search can return `Complete` and commit **S·P·C·P·X**.

The existing release test passed.

### T33 and T34: current-node self-loop

At the C-anchored decision:

| Route | Separation | Span |
|---|---:|---:|
| C→D | −31 | 200 |
| C→X | −31 | 119, outside window |
| C→C→X | 38 | 188 |

The first self-loop is `Live` because it is an initial hop. C is traversed without becoming an endpoint candidate.

The 20-pair X candidate has weight approximately 18.43, versus D’s weight 5, exceeding the 2.5 ratio.

Further repeat exploration does not force false truncation: with β = 150, the later C re-entry at separation 107 satisfies:

\[
107+150=257>250.
\]

T33 therefore withholds D; T34 confirms the assembled sequence **S·P·C·C·X**. Both passed.

### T35: alternating cycle

At C:

\[
s(A)=-31,\quad s(C_{\rm revisit})=38,\quad s(X)=107.
\]

The competitor span is:

\[
75+107+78=260.
\]

The repeat probe reaches this supported endpoint, returns `true`, and forces outer truncation. D remains uncommitted.

The existing release test passed.

## 6. Queue Budgeting & Resource Bounds (Q_MAX = 1024)

The main queue checks capacity before its initial and subsequent pushes, at `exspander.rs:363` and `:455`.

The repeat probe checks capacity before subsequent pushes at `:608`. Its initial push is into a newly empty queue, so it cannot violate the limit.

Therefore:

- Neither queue can exceed 1024 entries.
- Refusing an enqueue conservatively marks unfinished work.
- Main and probe state limits do not silently discard branches as complete.
- Depth limits terminate zero-increment cycles.
- The outer extension loop also has a finite step budget.

The queue bound is not a constant bound on total search cost. `distinct_outgoing` materializes distinct targets using linear membership checks, giving quadratic deduplication work in high-degree adjacency lists. Sequence extraction and repeated probes add graph-dependent cost. No unbounded cyclic execution was found.

**T39 does not test the capacity boundary.** Its 50 terminal targets produce a queue of at most 50 states. With no evidence, it naturally makes zero joins even if queue-capacity handling is absent.

## 7. Reference Oracle & Differential Verification Audit

[Source: independent_resolve_edge](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f17-grok-worktree/tests/test_c0_02r_repeat_completeness_gate.rs:103)

The oracle independently implements sequence-overlap checks and candidate scoring without calling the corresponding production functions.

Verified improvements:

- Zero overlap is rejected.
- Forward-edge selection chooses a unique highest-support sequence-valid edge.
- Ties reject.
- Self-anchor evidence affects completeness without emitting a self candidate.
- Scoring independently deduplicates in-window pair IDs.

Two validity defects remain.

**Reciprocal selection is incomplete.** At line 141, reverse candidates are filtered to `chosen_ov` before support comparison. A higher-support valid reverse edge with another overlap is invisible. Reverse support is also never compared with forward support.

**Search initialization omits prefix uniqueness.** Initial paths contain only `edge.target`; the already-used fork is absent. A long fork self-loop can therefore be admitted and scored despite representing a repeated physical long node.

Other oracle limitations include:

- no consumed-node/prefix context;
- no upstream-deletion filtering;
- circular successors treated as continuation before edge and budget validation;
- iteration over edge records rather than distinct directed targets.

The six differential tests genuinely invoke both the oracle and production resolver. Their comparisons are principally endpoint/deletion assertions, however; they do not compare production candidate inventories, weights, and `SearchStatus` directly.

They can detect invalid pruning that causes a tested decoy commit or loses a required positive-control commit. They do not establish general mutation sensitivity. Invalid pruning of an untested route, or pruning masked by another truncation cause, can survive.

## 8. Reverse-Complement Mirroring & Fixture Sensitivity Audit

### T22 and T24

The repair is substantive:

- `push_bio(..., true)` stores reverse-complemented sequences.
- `link_bio(..., true)` connects twinned directed nodes with reciprocal edges.
- `walk(..., true)` anchors evidence in the matching directed frame.
- Both tests assert that S has no forward outgoing edges.

These fixtures genuinely exercise backward extension from `S.twin()`. Both passed.

### T28 and T36

Both explicitly construct:

```text
C → P.twin() → Y
```

C’s suffix matches `rc(j_b)`, and Y’s prefix matches `rc(j_a)`. These are sequence-valid opposite-orientation routes.

The source establishes actual twin reachability; both tests passed. Their final assertions alone do not provide traversal telemetry.

### T29

Evidence is correctly anchored at S2, the active Case A anchor:

\[
span(S2,D)=80+38+82=200,
\]

\[
span(S2,X)=40+107+53=200.
\]

The control preserves topology and removes X evidence. D then commits. This is a meaningful evidence-sensitivity control.

### T30

The intended downstream decision is not reached.

Its first search is anchored at S, through branching P1:

\[
s(P2)=80-31-31=18,
\]

\[
span(S,P2)=60+18+60=138<180.
\]

Following P2→C→P1→P2 gives separations:

\[
18,\ 67,\ 116,\ 165.
\]

The later P2 scores:

\[
60+165+60=285\in[180,320].
\]

The repeat probe therefore truncates the **initial S-anchored search**. The P2-anchored evidence does not control the asserted withholding.

The negative control removes C and X and makes the graph linear. Linear extension can consume D without exercising a competing-evidence decision.

T30 therefore does not verify the requested decision-anchored sensitivity repair.

## 9. Code Scope & Production Regression Audit

Only `src/resolve/exspander.rs` changes in production. Mapper, graph-building, and polishing files are unchanged.

The substantive production changes concern admission classification, repeat/depth exploration, queue handling, and self evidence. Changes around splicing, reciprocal path reconstruction, lifecycle code, and final cleanup are formatting-only.

The F11-B contracts remain intact:

- checked lifecycle entry;
- `Result<JunctionStats, ExspanderError>`;
- guard transition to `Completed` on success;
- `Failed` on unwinding while running;
- validation before constituent deletion and assembled-node insertion.

The guard is not a general rollback transaction. Its unchanged comment explicitly states that partial edits survive panic/OOM.

### Test inventory and execution

Source inventory confirms **181 tests**, including 42 state-semantics tests versus 33 at baseline: nine additions.

Existing release artifacts were executed across 16 suites:

- **175 passed**;
- **6 failed at temporary-file creation because the filesystem is read-only**.

The six environment-blocked tests were:

- two frontier-preservation boundary sweeps;
- three ingest-estimator FASTX tests;
- the end-to-end synthetic assembly test.

No fresh build, compiler diagnostics, or complete writable-environment release run was performed. Existing binaries cannot independently establish a fresh source-to-binary build.

### S0-001 and S0-003

The available release assembler hash matches the execution report:

```text
0e708ab19ba43bd7a1af12340e56a1cff41d6ca6af593b861ce61a3d45c6e5ea
```

Stored output hashes independently match the baseline-recorded values:

```text
S0-001  bba42d4d86a010bff26ef11d2022f20004fa8237bf4dd655c7f273f61a1b5dae
S0-003  587f72eb5b9ffb9b77de1aa2fdf2b2a6110196ca710c4fe600780e0431b7ba93
```

Direct sequence comparison established:

| Fixture | Reference length | Contig length | Exact reverse-complement match |
|---|---:|---:|---|
| S0-001 | 5000 | 4957 | Reference offset 5 |
| S0-003 | 8000 | 7969 | Reference offset 13 |

This confirms **100% identity across each stored contig**, including zero substitutions or indels. It does not mean 100% reference coverage.

Fresh assemblies and historical alignment telemetry were not reproduced.

## 10. Subsystem Verdicts Ledger (Subsystems A through K)

| Subsystem | Verdict | Basis |
|---|---|---|
| A: β Mathematical Soundness | PASS | Complete evidence domain; valid post-hop lower-bound proof |
| B: Single Authority Consistency | PASS | All three walkers delegate admission; documented policy distinctions |
| C: Complete / Truncated Soundness | PASS | No confirmed false-`Complete` path; exhaustion propagates |
| D: Self-Anchor & Cycle Handling | PASS | Intermediate traversal preserved; endpoints excluded |
| E: Repeat Admissibility & T27 | PASS | Source trace and passing positive control |
| F: Shadow Candidate Safety | PASS | Monotonic shadow inheritance; commit exclusion |
| G: Queue & Complexity Safety | CONDITIONAL PASS | Source bound sound; T39 misses the boundary |
| H: Reference Oracle Validity | CONDITIONAL PASS | Independent primitives, but reciprocal and initialization defects |
| I: RC / Twin Reachability | PASS | Genuine mirrored and twin graph routes |
| J: Regression Integrity | CONDITIONAL PASS | 175 passes; six environment-blocked tests; stored artifacts verified |
| K: Real-Data Benchmark Readiness | CONDITIONAL PASS | Selective evaluation reasonable; acceptance evidence incomplete |

## 11. Findings Inventory (P0 / P1 / P2 with Root Cause, File/Line, Counterexample, and Evidence Status)

**P0: none confirmed. P1: none confirmed.**

### F18-P2-001 — Oracle accepts inconsistent reciprocal selection

**File/line:** `tests/test_c0_02r_repeat_completeness_gate.rs:141`, `:153`.

**Root cause:** Reverse candidates are restricted to the forward-selected overlap before resolving reverse ambiguity; reciprocal support equality is omitted.

**Counterexample:** Sequence-valid homopolymer junctions with:

```text
Forward: overlap 31/support 10; overlap 30/support 9
Reverse: overlap 31/support 10; overlap 30/support 20
```

The oracle returns overlap 31. Production reverse selection chooses overlap 30, so final path validation rejects the junction.

Even a sole reverse overlap-31 edge with support 1 is accepted by the oracle against forward support 10, while production validation rejects it.

**Impact:** The oracle can certify physically uncommittable paths.

**Evidence status:** Source-confirmed counterexample; no new executable fixture compiled.

### F18-P2-002 — Oracle initialization permits long-fork reuse

**File/line:** `tests/test_c0_02r_repeat_completeness_gate.rs:230`, `:235`, `:300`.

**Root cause:** Initial visited paths omit the already-used fork.

**Counterexample:** S→F, F→F, F→D, with F length 600 and overlap 31. With `init_sep = 569`, the oracle admits the initial F self-loop at separation 538 and can score D after it at separation 1107.

The represented assembly path contains F twice. Production rejects that reuse because F is already in `visited`; the splicer also prohibits it.

**Impact:** Reference candidate inventories can contain structurally illegal paths.

**Evidence status:** Source-confirmed counterexample; no new executable fixture compiled.

### F18-P2-003 — T30 bypasses its intended evidence decision

**File/line:** `tests/test_c0_02r_state_semantics_gate.rs:3084`, `:3086`, `:3098`.

**Root cause:** First-decision evidence has span 138, below the 180 minimum. A later cycle revisit scores at 285 and truncates the initial search. The control removes topology instead of isolating evidence.

**Counterexample/sensitivity:** Removing `(P2,X)` competitor evidence leaves the S-anchored truncation mechanism unchanged.

**Impact:** T30 does not establish downstream competitor-evidence sensitivity; this part of F16-P2-003 remains unverified.

**Evidence status:** Exact source/arithmetic trace; existing test passes without establishing the intended property.

### F18-P2-004 — T37 does not exercise β boundaries

**File/line:** `tests/test_c0_02r_state_semantics_gate.rs:3475`.

**Root cause:** The fixture is linear S→P→D. No candidate search or repeat re-entry invokes the β prune.

**Counterexample/sensitivity:** Changing the β comparison from `>` to `>=` would not affect this test’s execution path. No `max_insert + 1` case exists.

**Impact:** The advertised β boundary regression gate is ineffective.

**Evidence status:** Source-confirmed reachability analysis; existing test passes.

### F18-P2-005 — T39 does not exercise queue overflow

**File/line:** `tests/test_c0_02r_state_semantics_gate.rs:3559`.

**Root cause:** Only 50 terminal branches are created; capacity is 1024. No evidence supports a commit.

**Counterexample/sensitivity:** Removing queue-capacity enforcement leaves this fixture’s maximum queue size and zero-join result unchanged.

**Impact:** Queue saturation and its required truncation behavior lack the advertised regression coverage.

**Evidence status:** Source-confirmed state bound; existing test passes.

## 12. Overall Verdict & Final Recommendation

**CONDITIONAL PASS for the audited production implementation. C0-02R remains NOT ACCEPTED.**

The evidence-floor proof holds under valid production graph and scoring contracts. No unsafe β pruning, surviving silent short-competitor exclusion, truncation bypass, or false-`Complete` production route was confirmed.

Selective real-data benchmarking is reasonable as further evaluation. Before unconditional acceptance:

1. Repair reciprocal selection and initial uniqueness in the oracle.
2. Rebuild T30 around a demonstrably reached decision with identical-topology evidence controls.
3. Replace T37 and T39 with fixtures that actually cross their stated boundaries.
4. Perform a fresh build and complete release regression in a writable environment, including fresh S0-001/S0-003 assemblies.

The present evidence supports continued evaluation, not a claim that every F17 verification mandate has been satisfied.