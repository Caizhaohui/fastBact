# F24_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT

Date: 2026-10-29  
Project: fastBact  
Phase: F24 — Provenance-Safe Contiguity Recovery  
Decision class: V3.2 §4, L3 architectural specification  
Disposition: **SPECIFIED; IMPLEMENTATION AND VALIDATION NOT PERFORMED**  
Worktree: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f24-contiguity-worktree`  
Branch: `feature/f24-provenance-safe-contiguity`  
Frozen baseline: `f23-accepted-2026-10-29^{commit}` = `68b19f06d22537cf199e8fe124ea6b62bb0bb50a`

## 0. Authority, scope, and decision

This document is the F24 technical implementation and acceptance contract requested by the user. MUST, MUST NOT, and SHALL are normative. It does not grant scientific acceptance, merge, push, or release permission. The user remains FINAL_AUTHORIZATION_AUTHORITY. No production code or frozen test oracle is changed by issuing it.

**Chosen design: separate graph-preserving traversal, locus attribution, and physical evidence reach.** Replace the overloaded SI-2 early return with two capabilities: structural corridor traversal and evidence-anchored fork resolution. Introduce directed entry/exit metadata and certified absorption barriers. Preserve intrinsic pooled provenance and monotonic contamination on every existing walk. Use linear conservative propagation instead of claiming weighted Dijkstra is linear.

**Explicit qualification of the mandate:** a long node can isolate its far exit from a specified upstream influence. Length cannot prove genomic uniqueness. A 10 kb collapsed repeat can carry pooled evidence at its last 100 bp just as a 600 bp repeat can. Consequently:
- forward absorption clears only an entry-shadow obligation, and only at an independently attributable exit;
- it never clears intrinsic shared-copy provenance;
- SCC membership alone establishes neither uniqueness nor a new genomic copy;
- deterministic traversal does not authorize a later fork or an evidence-anchor reset.

The unconditional statement “length >= max_insert implies a 100% locus-specific exit” is rejected. The conditional distance statement in §3.2 is adopted. Without an independent attribution witness, ambiguous long nodes remain unable to decide forks. This is a correctness requirement, even if the numerical recovery gates are not reached.

The specification authorizes a conservative recovery implementation. It does **not** assert that the present graph and pair matrix contain enough information to recover every withheld locus, or that recovery targets are achievable without additional scientific information.

### 0.1 Baseline and status precedence

The annotated tag object resolves to `0a74fc496d511526753c365b02cdea0ecb9754ae`; its peeled commit, not the tag-object hash, is the baseline above. HEAD in the designated worktree matched that commit during preparation.

Root AGENTS.md retains `F23_ALGORITHM_CONTRACT_BLOCKED` and NOT ACCEPTED text. The dated acceptance block in HANDOFF and ACTIVE_TASK records user-authorized F23 completion, and the current user explicitly freezes that accepted baseline. Record this stale-status conflict in HANDOFF; preserve historical text. F24 does not reopen F23 acceptance or declare F24 accepted.

### 0.2 Evidence inventory and strength

| Evidence | Observation | Strength / limitation |
|---|---|---|
| `src/resolve/exspander.rs:355-390` at baseline | Distinct external origins are aggregated per SCC and copied to all member orientations | Inspected source; comments elsewhere describe different root-lineage semantics |
| `src/resolve/exspander.rs:394-477` | Reverse BCI absorbs at long single-context predecessors; forward BCI has no corresponding absorption | Inspected source |
| `src/resolve/exspander.rs:545-574,1585-1600` | Seed gate requires length >=500, context count 1, no downstream BCI; rejection precedes structural traversal | Inspected source; backward BCI is NOT itself a seed-gate predicate |
| `src/resolve/exspander.rs:688-778,1339-1405` | PE span, PB filtering, pair deduplication and bounded search | Inspected source |
| `tests/test_c0_02r_unique_anchor_gate.rs` | Shared-long, short-predecessor bypass, attributable short bridge, bubble, remap, cycle controls | Frozen behavioral obligations, not inferred from old architecture prose |
| `tests/test_c0_02r_state_semantics_gate.rs` | T26–T39, lifecycle, shadow, re-entry and truncation contracts | Frozen behavioral obligations |
| HANDOFF / ACTIVE_TASK dated 2026-10-29 | F23 safety results, 118 targeted and 194 workspace tests, accepted NGA50 trade-off | Recorded historical evidence; not rerun here |
| User-provided F24 telemetry | R1-001: 689/716 SI-2; R1-003: 874/876; 495 seeds with fail_ctx=438, fail_bci_ds=443, fail_bci=438 | Supplied observations; raw run manifests/logs must be attached and hashed before numerical acceptance |

Historical rulings consulted: `~/.local/state/fastbact-f23-post-audit-tests/validation/F23_ASTRA_ARCHITECTURAL_RULING.md` and `F23_FINAL_PROVENANCE_ARCHITECTURE.md`. They contain superseded status and inconsistent formula prose. The accepted source/test coordinate semantics govern; do not reintroduce the old reversed PB inequality or add read length twice to span.

## 1. Root cause analysis

### 1.1 A seed capability was conflated with a traversal capability

For orientation d, baseline SI-2 effectively requires:

```text
seed_ok(d) = len(d) >= 500
             AND ctx_count(d) == 1
             AND NOT bci_downstream(d)
```

Failure returns an empty extension before examining outgoing structure. This makes “cannot originate a PE fork decision” equivalent to “cannot preserve any deterministic graph adjacency.” That implication is unnecessary for a verified 1-out/1-in corridor.

The supplied abort fractions are 96.23% and 99.77%. They establish where control flow stops, not which rejected seeds are biologically recoverable. The three seed-failure counts overlap and must not be added; fail_bci is diagnostic of ARG/contamination, not an additional conjunct of the implemented SI-2 predicate. Report physical seeds separately from the two oriented extension attempts.

### 1.2 Entry and exit were conflated

Let an oriented node u have sequence length L(u), entry coordinate 0 and exit coordinate L(u). An edge u→v has validated overlap o(u,v).

The existing forward update is:

```text
entry_distance(v) = entry_distance(u) + L(u) - o(u,v)
```

Thus v is flagged according to its entry, independently of L(v). A 10 kb v entered 100 bp after a confluence is marked just like a 100 bp v. SI-2 then rejects the 10 kb node's far exit despite the very different local geometry.

Entry contamination, exit exposure, and shared-copy identity are different facts. Correcting the first two cannot erase the third.

### 1.3 SCC reachability is too coarse as a physical locality model

Tarjan correctly computes reachability classes. The failure is interpreting the union of all SCC entry origins as a local, bounded contamination fact at every port. A bacterial chromosome cycle can contain hundreds of long corridors; one shared junction then affects distant ports under that interpretation.

Conversely, if the semantic definition is unrestricted reachability of source identities, all vertices of an SCC really do receive the same union. Replacing Tarjan with a worklist without changing that semantic definition cannot improve attribution. F24 therefore separates:
1. structural cycle membership;
2. inherited locus attribution;
3. port-local physical shadow.

A cycle return must not mint a new lineage. Two genuinely independent entries must not be collapsed into one merely because their paths return to the same cycle.

### 1.4 The linear-time claim needs correction

The baseline uses two BinaryHeap Dijkstra passes, repeated `find_branch_origin` walks, and linear searches in origin lists. Its unconditional O(V+E) comment is not a valid bound. A distance-limited weighted search is not ordinary BFS. F24 must supply an actual linear analysis, with conservative loss of precision where necessary (§6).

## 2. Observable recoverability

All decisions use the live oriented graph, retained provenance, accepted edge geometry, and actual pair rows. Reference alignments are acceptance oracles, never resolver inputs.

### 2.1 Recoverable classes

**R-S: structural corridor.** Each accepted edge u→v has exactly one live structural successor of u and one live structural predecessor of v; a unique validated directed edge identity; no circular/self/re-entry exception; and all lifecycle, ownership and conservation checks. Contraction creates no new graph walk or branch pairing. It needs no PE winner. It may combine sequence that remains unattributable as an evidence source.

**R-A: attributed continuation.** The deciding source has a valid directed single-lineage certificate, its source interval is admissible, and the complete candidate satisfies unchanged span, support, rivalry, shadow and search-completeness gates.

**R-I: physically isolated, attributed exit.** R-A's attribution requirement is independently established, and every entry-only shadow reaching the proposed evidence interval is excluded by the physical test in §3.2. A fresh walk may originate there; an already contaminated walk does not acquire a retroactive bridge through its upstream repeat.

**R-C: cycle corridor or legal cycle bridge.** A structural corridor inside a cycle can be traversed up to the existing re-entry boundary. A fork/copy-count decision still requires R-A. T34 remains a legal evidence-based bridge; F9b remains a pooled-evidence withholding case.

### 2.2 Non-recoverable classes

Withhold fork decisions for unknown or shared attribution; insufficient/pooled/out-of-horizon evidence; ambiguity or incomplete candidate enumeration; consumed divergent shadows; invalid overlap; unsupported repeat multiplicity; and composite intervals whose attribution cannot be represented safely.

A long shared node, a sole valid edge remaining after discarding a structural competitor, and a source-free cycle without an attribution witness are not new unique anchors. A deterministic corridor ending at an unresolved fork can be retained only under graph-conservation and frozen rollback rules, without selecting the fork exit.

“Recoverable” is a capability established by a witness, not a promise that every graph fragment will merge. Resource caps, protected frontiers and later validation may still prevent a commit.

## 3. Chosen architectural design and core invariants

### 3.1 Directed identity and attribution certificates

Use d=(physical_node_id, orientation), with explicit entry and exit ports. Query the current PathStep orientation at every step; do not hold a global boolean strand constant while following edges that flip orientation. Reverse-complement tests transform both graph and evidence.

Attribution is a lattice:

```text
Unknown
Single(lineage_id, witness_id)
Shared
```

Unknown is not Single. Two unequal Single identities join to Shared; Shared is absorbing under ordinary propagation. Identical identities do not multiply on re-entry. Pending/unvisited analysis state must be separate from finalized Unknown so that an unresolved cyclic predecessor is not silently ignored.

A witness is one of:
- a graph-relative rooted context, propagated without meeting a distinct or unknown context;
- an explicitly verified same-locus bubble equivalence, preserving its incoming lineage;
- retained constituent provenance from a previously validated transformation, with exact orientation/interval projection.

These witnesses retain F23's graph-relative assumptions; none proves the absence of an unobserved genomic copy. A new certificate cannot be manufactured from length, coverage, mapping to a single collapsed graph node, the current SCC ID, a winning PE ratio, or deleting a competing entry. No reference-derived certificate may enter production.

For rootless single-entry cycle fixtures, an external attributable entry grounds the cycle; internal self-loops contribute no new identity. For a wholly rootless cycle without retained provenance, attribution remains Unknown. Structural traversal is still available.

**No new uniqueness oracle is silently assumed.** If present metadata cannot supply a certificate, the implementation must return Unknown/Shared. Additional certificate producers require a separately specified, reviewed correctness contract. A test-injected “unique” bit is not production evidence.

### 3.2 Locus-length physical isolation

Physical isolation is relative to a specified upstream uncertainty boundary, not to all possible repeat copies.

For an influence interval ending at coordinate b of U, a candidate observation at coordinate x>b is unreachable from that interval under horizon H only if a conservative lower bound on the physical separation is **strictly greater than H**. Equality is admissible under the existing inclusive span predicate. Every path or internal entry that could shorten this separation must be considered; otherwise the result is Unknown.

For the simplified point-entry model b=0, a mate observation strictly beyond a length-L exit has distance >L. Therefore L>=H excludes such a crossing. This says nothing about a source read already inside a collapsed copy of U. It also says nothing about a downstream read whose mapped end lies in a shared overlap before that exit.

The operational fork oracle remains:

```text
span(row, candidate) =
    row.dist_u1_exit + signed_path_separation + row.dist_u2_entry

admissible_span = min_insert <= span <= max_insert
```

Direct edges can have negative separation because overlap is subtracted. Do not clamp that signed span calculation to zero or add r2_len again. Validate read coordinates against node lengths.

For a conservative, constant-size exit certificate:
- retain the greatest possibly contaminated coordinate b in the walk orientation;
- an evidence-source start p must satisfy p>b+H to be outside that entry influence;
- require its complete read interval to fit the certified source interval;
- use the actual candidate span to validate the pair;
- if a less conservative bound is desired, it requires a new proven coordinate predicate, not raw L>=H.

This deliberately conservative p>b+H rule may miss safe rows. It cannot admit an entry-exposed row. A 10 kb unique node can have a usable terminal interval; a 700 bp node with overlap or a nonzero contaminated prefix need not have one. Boundary positives at L=H must explicitly place the observation beyond the exit and prove the exact span exclusion; they must not assume the entire node is clean.

### 3.3 Forward anchor absorption

Separate **intrinsic sharedness** from **entry shadow**. Intrinsic sharedness traverses long nodes unchanged. Forward absorption applies only to entry-shadow propagation at an exit with both:
1. independently established Single attribution for the evidence interval;
2. a valid physical isolation certificate.

Mark entry exposure even when exit propagation is absorbed. Never replace both with one false bci_downstream flag. At an unclassified or shared exit, propagate shadow conservatively.

Backward Q8 retains its accepted role: a long attributable predecessor may serve as an upstream anchor despite a future confluence. Forward symmetry is semantic, not a copy of “len>=500 ⇒ stop”: the forward side must cross the node and account for overlaps and affected intervals. Backward and forward predicates share coordinate helpers and orientation tests, not necessarily identical length thresholds.

### 3.4 Cycle-safe attribution

SCC IDs are used for termination, unknown-component detection, and diagnostics. They are not evidence-source identities. Propagate actual capped lineage states at oriented ports. An SCC-wide count may remain as a diagnostic summary but MUST NOT alone prohibit R-S traversal or erase a certified local interval.

Every independently sourced ingress is processed at its actual port. Re-entry propagates the same state and adds actual path distance. Genuine inherited Shared state is not removed by leaving and re-entering a node. In the absence of a certified barrier/equivalence, all mutually reachable ports may still converge to Shared; this is the correct conservative result.

A claim that local propagation automatically cleans all long nodes of a multi-entry SCC is forbidden. If the accepted observations do not distinguish two collapsed copies from one circular locus, F24 withholds evidence-based resolution there. Safe corridor contraction remains possible.

### 3.5 Deterministic traversal before SI-2

Treat structural continuation and PE resolution as separate modes.

For R-S, count **all live structural distinct directed neighbors** before evidence, ownership or overlap filtering. Parallel overlap edges to one neighbor still require unambiguous edge resolution. A bad competing edge does not disappear to create a degree-one proof.

Require out_degree(u)=1 and in_degree(v)=1 in the same directed graph. A 1-out edge entering a confluence is not R-S. A self-loop, physical-node re-entry (including opposite orientation), circular target, self-support stop, or already consumed target is not handled by this shortcut. Existing cycle search handles applicable exceptions.

R-S can append without an evidence anchor. It does not reset attribution, contamination, horizon or repeat multiplicity. At a fork, an unanchored walk stops. A shared seed may compact a common corridor, but cannot select either flank. No already contaminated anchored walk may use R-S to avoid PAG rollback.

Structural degree certificates use the immutable analysis snapshot plus current liveness/ownership checks. Deletion during resolution cannot make a formerly branching port newly degree-one for purposes of this pass.

### 3.6 Invariant register

| ID | Mandatory invariant |
|---|---|
| F24-I01 | Length proves neither Single attribution nor repeat copy number |
| F24-I02 | Entry shadow, exit isolation, inherited sharedness and walk contamination are separate states |
| F24-I03 | Unknown never grants anchor or PE-source eligibility |
| F24-I04 | Shared evidence has zero deciding weight regardless of weight magnitude or reversal |
| F24-I05 | Walk contamination is monotone; physical absorption never clears it |
| F24-I06 | Each deciding row is attributed in its actual source orientation and interval |
| F24-I07 | Span, deduplication, exclusive rival scoring, consumed shadows and truncation remain unchanged |
| F24-I08 | R-S preserves graph walks and cannot choose a branch or copy count |
| F24-I09 | A rejected path leaves no ownership, graph or pair-matrix mutation; conservative protected-frontier additions remain required |
| F24-I10 | Remap preserves conservative provenance; lost competitors do not create uniqueness |
| F24-I11 | Static analysis is O(V+E); runtime enumeration retains all frozen bounds |
| F24-I12 | Safety and completeness gates cannot be traded for N50/NGA50 improvements |

## 4. Rejected alternatives and counterexamples

| Rejected rule | Counterexample / reason |
|---|---|
| Delete SI-2 | F3 shared L seeds pooled downstream and reverse-flank decisions |
| Any len>=H node is unique | Lengthen F4's pooled L to 10 kb; its terminal rows still mix copies A and B |
| Copy backward >=500 absorption into forward propagation | F4 has X(shared)→L(600), with H=700; both attribution and physical test fail |
| Only inspect immediate in-degree | X is shared but L has one predecessor; inherited ambiguity survives |
| All SCC members are clean / one lineage | F9b's two independent entries then gain pooled deciding evidence |
| All SCC members are permanently blocked | T34 and structural cycle corridors are unnecessarily withheld |
| Ignore all cycle-return edges | An external lineage arriving through a return path disappears |
| A worklist alone fixes SCC poisoning | Unrestricted reachability has the same SCC fixed point |
| Retain one valid outgoing edge after pruning | Invalid/ambiguous competitor or consumed shadow can be the true alternative |
| Treat every 1-out edge as safe | A→R and B→R makes R a confluence; arbitrary consumption loses locus context |
| Reset a contaminated walk at an isolated exit | Attaches an unresolved upstream copy to a new downstream decision; F1/F4 attack |
| Raise max_insert or lower support/ratio | Converts pooled or physically impossible evidence into apparent support |
| Use whole-composite evidence after a clean suffix appears | F9a and unique/shared/unique composites launder pooled intervals |
| Ordinary BFS for exact weighted distance | Different overlap/length costs violate shortest-path assumptions |
| Accept N50 increase alone | A longer chimera increases N50 while reducing aligned correctness |

## 5. Algorithm and data-structure specifications

### 5.1 ProvenanceTable schema and lifecycle

Keep a side table keyed by physical node ID, with two directed records. Each directed record must represent:

```text
lineage: Unknown | Single(id, witness) | Shared
entry_shadow: bool
exit_state: Exposed | Isolated(certificate) | Unknown
trusted_source_interval: optional [lo, hi) with lineage
backward_contamination: bool
structural_in_degree, structural_out_degree
scc_id: diagnostic only
generation: graph-analysis generation
```

Use explicit validity tags; absent interval is not “whole node trusted.” A physical intrinsic-shared flag, where established, dominates both orientation records. Directional incoming ambiguity need not be identical in both orientations. Retain old counters only as derived diagnostics with documented meanings.

The interval is a conservative subset of attributable positions, never the convex hull across a pooled gap. At most one retained interval per orientation is permitted for F24. If multiple disjoint intervals arise, retain the interval adjacent to that orientation's exit if certified, otherwise retain no interval; never join them through unknown sequence. Existing PB filtering remains an additional restriction for frozen composite cases.

Analysis is immutable within one resolver pass. New composites receive projected metadata atomically with their graph and pair remap. They do not trigger full graph reanalysis or a second resolver invocation. Metadata length equals graph.nodes length after every successful commit. Invalid generation or missing metadata fails closed.

### 5.2 Linear static analysis

1. Build live oriented CSR adjacency, reciprocal-edge validation metadata, and structural degree counts once. Preserve self-loops for traversal; identify them separately for lineage identity.
2. Compute SCCs iteratively in O(V+E). Build maximal nonbranching corridors once, assigning every edge to at most one corridor; do not repeat predecessor walks per edge.
3. Initialize grounded lineage witnesses. Same-source closed bubble arms propagate the same identity; any unknown/independent ingress prevents an unsupported equivalence. The simple F6b bubble is covered by shared root identity, without an all-pairs reachability test.
4. Propagate capped lineage labels using a monotone worklist: unset→one identity→Shared, at most two upward changes. Afterwards mark ungrounded source SCCs Unknown and propagate unknown-obligation bits; Single is usable only when that bit is absent. This second monotone boolean pass prevents treating an unresolved predecessor as an empty set. Never replace a Single ID by a different Single ID.
5. Propagate backward and forward shadow obligations as conservative boolean reachability, stopping only at the corresponding certified absorption ports. Each port changes once per pass. This intentionally overapproximates distance-limited contamination.
6. Refine only within maximal nonbranching corridors by prefix/suffix scans over validated overlap-adjusted coordinates. No repeated path search. At a junction, or where no sound finite entry bound is available, retain Exposed/Unknown rather than guessing a distance. Independent internal contamination invalidates any earlier absorption certificate.
7. Freeze metadata and diagnostic counts. Check RC-equivalent fixtures and graph-generation agreement.

This specifies a strict linear conservative implementation, not exact shortest distances everywhere. The earlier bounded Dijkstra results may be used as a test oracle on small fixtures, never as an unaccounted production pass. If the conservative replacement loses a frozen positive control, repair the certificate/analysis within these bounds; do not waive that control.

A fresh source interval downstream of intrinsic Shared requires an independent retained witness; neither this worklist nor its physical shadow refinement produces one from nothing. This is an explicit data limitation, not an implementation TODO to be filled with a length heuristic.

### 5.3 seed_anchor_eligible

Make eligibility directed and interval-aware. The boolean is true iff:
- node is live, ordinary and length >=500;
- attribution is Single with valid generation/witness;
- a nonempty admissible source interval exists;
- relevant entry shadow is either absent or excluded for that interval by a physical certificate;
- inherited intrinsic Shared/PB restrictions do not invalidate the interval.

The function must not decide whether structural traversal is possible. For a rejected anchor, initialize `anchor=None`, `mode=StructuralOnly`; do not automatically return before R-S inspection.

Backward BCI continues to govern anchor reset/contamination under accepted Q8 semantics; it is not silently added as a new blanket seed rejection. Seed permission and subsequent fork permission are distinct.

### 5.4 extend_one_direction

Maintain anchor ID/orientation and certified interval, overlap-adjusted anchor_gap, anchor_steps, monotone contaminated flag, structural-only mode, path-local visited physical IDs, and a tentative ownership delta.

At each iteration:
1. Apply existing extension cap, self-support stop, liveness and lifecycle checks.
2. Inspect full structural neighborhood and validate directed overlap.
3. If R-S holds, tentatively append the step; accumulate geometry and provenance; reset the anchor only under unchanged clean-walk ARG plus the new certificate checks.
4. If R-S does not hold, use the existing resolution path only with a usable evidence anchor/source. A fork must not be “simplified” by deleting rejected candidates.
5. In contaminated regions, exclude head-keyed pooled evidence; score from the attributable anchor. Else preserve accepted single-context source semantics, including T34. Apply PB/source interval filtering before support counts, weights, beta bounds and rival comparison.
6. No anchor or no sufficient admissible support: stop and record a reason. Apply existing PAG rollback for tainted anchored walks. A structural-only prefix may survive only if its contraction independently passes graph conservation.
7. Commit only through the existing path verification/atomic mutation mechanism. Protect unresolved outgoing and incoming frontiers. Release exactly this attempt's tentative ownership entries on rollback.

Do not concatenate forward and reverse structural-only results across an ambiguous seed unless the combined contraction itself meets the directed degree/conservation proof. A fresh certified suffix seed is independent of an older ambiguous prefix; the two must not be spliced merely because they touch.

### 5.5 Evidence and merge remap

Source interval filtering precedes ALL support, weight, deduplication, beta-minimum and shadow calculations. The same admissibility predicate must serve primary candidates and consumed shadows. Pair IDs count once per candidate; re-entry does not multiply evidence. Preserve exclusive-support rival scoring and the single-candidate quality bar.

For a forward source read start p in a length-L node, p=L-dist_u1_exit. Baseline PB boundary b admits starts p<b, equivalently dist_u1_exit>L-b. For newly certified intervals, require the full read footprint inside the interval. Reverse intervals transform as [lo,hi)→[L-hi,L-lo), with evidence coordinates transformed consistently.

For a merge placement offset a, interval endpoints shift by a after orientation transformation. a is the overlap-adjusted path placement, not a sum of raw lengths. Classify overlap positions conservatively if either constituent marks them unknown/shared. Existing earliest pooled-boundary semantics remain conservative; multiple alternating zones must never collapse into a falsely clean whole node.

Contraction of a shared corridor preserves every external predecessor/successor and the fact that the composite is shared. Retain raw nodes/evidence in the audit ledger; no historical artifact is overwritten.

### 5.6 Telemetry

Record mutually exclusive terminal reason per oriented attempt and separate nonexclusive diagnostics:
- anchor_rejected_shared, anchor_rejected_unknown, anchor_rejected_entry_shadow;
- structural_steps_committed, structurally_recovered_bp;
- exit_absorption_certificates, certificate_source_kind;
- fork_withheld_no_anchor, horizon, pooled_empty, insufficient_score, truncated;
- attempted/committed fork source interval and lineage witness;
- SCC size/entry count as diagnostics;
- original and new SI-2 summaries with explicit denominator.

A structural-only walk that reaches no fork is not an SI-2 fork abort. Preserve legacy withheld counters where frozen tests require them, but do not equate their changed totals to biological recovery. Trace output must be noninterfering.

## 6. Complexity, memory and bounded traversal

Let V be physical nodes, E directed stored edges, P pair rows, B sequence bytes and H=max_insert. Oriented ports multiply V by a constant.

| Component | Bound |
|---|---|
| CSR, SCC, corridor decomposition | O(V+E) time and memory |
| Capped lineage + unknown propagation | O(V+E), constant state changes per port |
| Conservative shadow propagation and corridor scans | O(V+E), independent of numeric H |
| New static metadata | O(V); target <=128 bytes/physical node including both orientations, measured with actual layout |
| Transient static storage | O(V+E); no per-seed graph arrays or all-source bitsets |
| Sequence-overlap validation | Charge sequence work separately, O(B) for cached scan of required sequence data, not falsely O(V+E) for arbitrary strings |
| Pair filtering/remap | Charged to rows inspected/remapped; not included in graph-static bound |
| Extension state | Existing bounded path storage and tentative delta; no global walk-history bitsets |

Preserve MAX_EXTENSION_STEPS=64, MAX_CANDIDATE_DEPTH=4, MAX_CANDIDATE_PATHS=32, MAX_DFS_STATES=256 and MAX_QUEUE_CAPACITY=1024. Caps remain conservative truncation, never “best route so far is complete.” A candidate can still require inspecting its evidence rows and structural neighbors; do not claim O(1) total work independent of degree or P.

Exact weighted multi-source Dijkstra costs O((V+E) log V); a bucket implementation costs O(V+E+H) under suitable integer weights. Neither satisfies the requested unconditional O(V+E) bound for arbitrary H. They are rejected as the mandatory static implementation. Test-only exhaustive distance oracles are permitted.

Performance acceptance: static edge relaxations/visits must scale linearly on doubled chains, bubbles, high-indegree junctions and large SCCs; no count proportional to V×E. Whole-assembly median runtime and peak RSS must each be <=1.20 times the frozen baseline on identical resources. Also preserve the recorded integration-test limit (<3.5 s under its frozen harness). Failure is diagnosed, not hidden by more CPUs or changed read subsets.

## 7. Deterministic test obligations

### 7.1 Frozen gates

Preserve all assertions, truth labels, edge geometry, support thresholds and outcomes in the frozen baseline suite, including:
- F1 collapsed-long hijack, F2 flipped pooled weights, F3 shared seed in both orientations;
- F4 short-predecessor bypass, F5 attributable in-horizon repeat bridge;
- F6a long-node continuity, F6b SNP-bubble continuity;
- F7 boundary geometry, F8 reverse F4, F9a merge-remap, qualified F9b;
- T26–T39, especially T26/T27/T30/T34 positive commits, T33/T35/T36 decoys, T37 beta boundaries, T38 self-anchor stop and T39 queue truncation;
- remaining baseline directed-edge, conservation, cleanup, lifecycle, frontier, consumed-shadow and integration tests, including T40/T41.

F9b means pooled multi-context re-entry must not restore attribution, not “all cycles abort.” Require a new instrumented control to demonstrate actual repeated-node exploration; the existing F9b final-output assertion alone does not prove that runtime event. F7 currently includes arithmetic-only checks; add end-to-end boundary controls without replacing it.

### 7.2 New F24 matrix

Each case must assert exact sequence/motif outcomes, graph/evidence conservation, deciding-source eligibility and trace reachability. “No chimera” alone is insufficient for a positive control.

| ID | Construction | Required result |
|---|---|---|
| P01 | Independently attributed 10 kb node with entry-only shadow; terminal source rows beyond certified exclusion boundary | Eligible fresh forward anchor; authentic supported fork commits |
| N01 | Same geometry as P01 but F4 inherited Shared, including 10 kb L and reversed/amplified pooled weights | No new anchor or fork commit; pooled deciding weight zero |
| P02/N02 | Isolation bounds H-1, H, H+1; lengths 499/500, H and 10 kb; nonzero overlap and contaminated prefix | Exact span/interval oracle decides; equality never automatically excluded |
| P03 | Unanchored deterministic corridor ending before a fork, valid 1-out/1-in edges | Corridor retained where conservation permits; neither fork exit selected |
| N03 | 1-out into a 2-in confluence; reverse case; one invalid structural competitor | No structural shortcut through ambiguous junction |
| P04 | Single-entry cyclic fixture and same-lineage return | T34 behavior and legal corridor traversal preserved; no new lineage on return |
| N04 | Multi-entry cycle, including entry through inherited-shared predecessor chain | Shared stays shared; no pooled fork decision |
| P05 | Large SCC containing a corridor with an independently retained local certificate | Certificate provenance retained; SCC summary alone cannot reject structural corridor |
| N05 | Same SCC with certificate removed / unsupported “unique” label | Unknown/Shared, no fresh evidence anchor; structural capability tested separately |
| P06 | True reverse-complement mirrors of P01/P03/P04, with read coordinates transformed | Canonically equal sequence decisions and reasons |
| N06 | Opposite-orientation physical-node reuse, hairpin and internal opposite-strand ingress | No cross-strand certificate leakage or unbounded traversal |
| P07/N07 | Unique/shared/unique composite, repeated merge, reverse constituent, interval boundary reads | Only certified interval rows score; pooled gap never becomes clean |
| N08 | Tainted A walk reaches a separately certified U | A cannot reset/clear taint; independent U seed cannot retroactively join A |
| N09 | Parallel overlaps tied; invalid edge; consumed divergent shadow | Frozen rejection/competition semantics preserved |
| N10 | Zero-progress edges/cycles and every existing depth/queue/step cap | Termination, correct Truncated status, protected frontiers |
| P08/N11 | Allowed R-S prefix followed by verification failure | No partial ownership/graph mutation; another legal seed remains usable |
| P09 | F5 external in-horizon support against arbitrarily large shared internal weights | Authentic bridge only; attributable support controls decision |
| N12 | Repeated pair IDs across loop, candidates and composite rows | Deduplication and exclusive rival support unchanged |
| P10/N13 | Whole-source vs partial interval anchor; invalid/stale witness/generation | Correct row filtering or fail closed in every scoring/beta path |

P01/P05 require a real supported witness producer, not an unverified fixture flag used only to satisfy tests. If only a unit-level synthetic certificate is available, label that unit coverage and leave the production integration gate unsatisfied. Do not invent a production uniqueness oracle to meet the target.

### 7.3 Verification protocol

Before editing, inventory/hashes of all frozen tests and source. Add new fixtures separately (suggested target: `tests/test_c0_02r_f24_contiguity_gate.rs`). Recovery-positive tests that expose baseline withholding must demonstrably fail on 68b19f0; safety negatives should pass there. Record the expected failure mechanism, not merely an exit code.

Run, on Slurm qfree compute nodes only:
1. new F24 target and the unique-anchor/state-semantics targets;
2. all targeted F23 gates and full release workspace tests;
3. trace-enabled/default equivalence and integration harness;
4. linearity/resource controls;
5. real-data validation in §8.

Use unpiped command logs or explicit pipefail and capture Cargo plus Slurm exit status. Record full source/test/binary hashes, command, job ID, node, resources, compiler and all failures. Historical 118/194 counts are inventory references, not a new fixed total after adding tests. Sol independently reviews algorithm and oracle integrity; Luna performs deterministic verification. Neither model judgment substitutes for a command result.

## 8. Real-data acceptance contract

### 8.1 Comparators and reproducibility

Evaluate R1-001 (K-12 MG1655) and R1-003 (BW25113) separately. Rebuild/replay the frozen F23 commit as comparator A and F24 as B with identical input hashes, parameters, read ordering, threads, k schedule, stage cutoffs and QUAST configuration. Use the accepted job 3052717 artifacts where manifests establish equivalence; otherwise rerun both.

Attach raw telemetry behind the supplied seed counts. Extract exact integer N50 and NGA50 from QUAST report.tsv. Do not label 58.9/62.3 kb as N50: the supplied values are NGA50, and historical F20 values differ. No comparison may silently mix those runs.

Use QUAST v5.2.0 and the accepted reference/configuration, including ambiguity handling and contig cutoffs. Preserve Stage 5 and final FASTA, reports, alignments, traces, input/source/binary checksums and Slurm accounting. Three identical-resource paired runs establish deterministic canonical sequence/metric equality and median runtime; use separate output directories.

### 8.2 Numerical recovery gates

The following are F24 engineering acceptance targets issued by this specification, not claims of achieved performance or retroactive changes to F23 acceptance.

For EACH dataset, both N50 and NGA50 must be >=1.50× its exact frozen F23 value (integer ceiling). In addition, NGA50 must meet the dataset floor:

| Dataset | Supplied F23 NGA50 (approx.) | Supplied historical NGA50 (approx.) | F24 minimum NGA50 | Aspirational NGA50 |
|---|---:|---:|---:|---:|
| R1-001 | 27.0 kb | 58.9 kb | 44,200 bp | >=58,900 bp |
| R1-003 | 17.5 kb | 62.3 kb | 46,800 bp | >=62,300 bp |

The floors recover approximately 75% of the supplied historical aligned contiguity. The exact 1.50× F23 gate also applies if baseline replay differs from rounded figures. Aspirational values do not substitute for mandatory gates. N50 has a relative gate because exact frozen N50 was not supplied; no number is fabricated.

Require no decrease in genome fraction or aligned reference coverage, no increase in unaligned sequence attributable to new joins, and no loss of previously represented unique loci. Check Stage 5 separately so output filtering or polishing cannot masquerade as resolver recovery. SI-2 reduction is explanatory telemetry, not a pass threshold.

### 8.3 Zero-new-error safety gates

All of the following are mandatory, with no cross-dataset averaging:
- R1-003 extensive misassemblies remain **0**.
- R1-001 extensive misassemblies do not exceed the replayed frozen F23 count. Its known pre-existing events are matched by reference breakpoint coordinates and alignment context, not contig names.
- **Zero new validated extensive or local misassembly events on either dataset**, even if another old event disappears and total count stays constant.
- No increase in reported local-misassembly count, mismatch rate or indel rate relative to frozen F23 under identical evaluation.
- No newly introduced repeat-crossing junction supported only by pooled evidence.
- No reappearance of AUD-01: specifically audit the historical cross-locus connection around reference positions 3,666,175 and 1,570,641, its reverse complement, and nearby equivalent breakpoint placements. Sequence/reference-coordinate evidence governs, not `contig_16` naming.
- Re-examine the previously documented R1-001 ambiguity-filter artifact using both candidates' alignments. It is not a general waiver for new events.

A new QUAST event, including an apparent ambiguity artifact, prevents PASS until event-level evidence reconciles it under the frozen evaluation policy. Do not change ambiguity thresholds, references, error labels or input truth to obtain acceptance.

### 8.4 Attribution audit

For every newly committed fork, retain its source witness, permitted interval, candidate paths, admissible pair IDs/spans, pooled rows excluded, rival scores and completeness status. Independently recompute a deterministic sample plus every join near a misassembly/long repeat. Structural-only joins carry a degree/edge/conservation certificate instead of a fabricated PE decision.

Require that newly recovered bases are traceable to R-S, R-A or R-I. A falling withheld counter without recovered correct sequence is not recovery. An unexplained fork has no acceptance credit.

### 8.5 Result states

- Safety failure: F24_REJECTED_SAFETY; restore no scientific claim, preserve all evidence.
- Safety passes but recovery or resource gates fail: F24_SAFE_BUT_TARGET_NOT_MET.
- Missing comparator, raw metrics, witness producer or required test evidence: F24_VALIDATION_INCOMPLETE.
- All deterministic, safety, recovery and resource gates pass with independent review: F24_TECHNICALLY_VALIDATED.
- Scientific acceptance/integration remains the user's reserved decision.

No target may be lowered automatically after seeing results.

## 9. Implementation order, limitations and handoff

Sequential implementation order:
1. Freeze baseline/test manifests; add failing recovery and passing adversarial controls.
2. Separate structural traversal from anchor eligibility, preserving transactional commit and rollback.
3. Add directed interval/port metadata, linear attribution/shadow passes, and certified absorption.
4. Audit all scoring, beta, shadow and merge-remap consumers.
5. Sol review; Luna deterministic/resource verification; then real-data validation.
6. Publish evidence and the actual result state for user acceptance.

One writer per isolated worktree. No Rust compilation, tests, assembly or QUAST on login nodes. No main merge, push, release, paid-provider escalation or expanded compute allocation is authorized by this document.

Known limitations:
- Observable topology and a collapsed pair matrix can be identical for one circular locus and several genomic copies. No algorithm can distinguish those hidden truths from these inputs alone.
- Global lineage union around an SCC remains shared without valid barriers. F24 does not solve that information deficit by relabeling nodes.
- Strict linear conservative shadow analysis may retain false withholding. Exact weighted distances are a different complexity contract.
- Single-interval metadata may discard valid evidence in complex composites; it must not admit pooled evidence.
- Physical isolation relies on the accepted bounded-insert model and valid coordinates. It is not a proof that real libraries contain no outliers.
- The required NGA50 recovery is an empirical objective. If it needs unsupported attribution, retain safety and report target failure; escalate with the minimal counterexample and missing evidence.
- This architecture has been checked against source and frozen test definitions, not executed. No new benchmark, runtime or scientific PASS is asserted.

**Final technical ruling:** implement only provenance-preserving recovery. Permit deterministic corridor traversal without a seed evidence anchor; absorb forward entry shadows at physically isolated, independently attributable exits; keep lineage and orientation explicit; retain pooled-copy and contaminated-walk restrictions. Accept F24 only after all unchanged safety oracles and the separate real-data recovery gates pass.

