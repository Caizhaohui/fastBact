# Phase F8 Independent Codex Re-Audit Report (AUDIT-20261021-005)

## 1. Executive Summary & Overall Verdict

- **Target Commit:** `37a7c0b`; verified HEAD: `8ceb3e5`.
- **Branch:** `feature/c0-02r-phase-f7-grok`.
- **Auditor Model:** Codex CLI. Requested designation: `gpt-6.1-sol`, medium reasoning effort; exact runtime model and effort are not independently exposed to this session.
- **Overall Verdict: FAIL.**

**One P1 regression was identified:** external-neighbor preservation inserts active nodes into `rejected_nodes`, which also excludes nodes from subsequent candidate searches. A successful splice can therefore hide a stronger candidate from a later seed and cause a weaker candidate to be committed while search returns `Complete`.

The original incoming-neighbor deletion defect is locally repaired. Both overlap-budget calculations are corrected, and the strengthened N=66 assertions are present without weakening the RED tests. Those improvements do not resolve the newly introduced search exclusion.

**Verification scope:** strictly the incremental changes from `0e513a3` to `37a7c0b`, plus their immediate callers and consumers. This session executed **37 existing release-binary tests: 37 passed, 0 failed**. The environment is read-only, so no fresh compilation, historical checkout execution, new counterexample compilation, or complete write-producing regression run was performed. The counterexample below is established by source tracing, not a newly executed harness.

## 2. Subsystem Verdicts

| Subsystem | Verdict | Basis |
|---|---|---|
| A. External Neighbor Preservation | **CONDITIONAL PASS** | Complete immediate-neighbor registration for canonical twin-complete graphs; registration has harmful search side effects. |
| B. Graph Cleanup Correctness | **PASS** | Registered active neighbors survive this resolver’s orphan cleanup. |
| C. Failure Atomicity | **PASS** | New registration executes only after successful validation. |
| D. Search Budget Correctness | **FAIL** | Arithmetic is repaired, but F7 protection can exclude admissible candidates while returning `Complete`. |
| E. Sequence Preservation | **CONDITIONAL PASS** | Tested branch sequences survive direct raw curation; `min_len = 0` is not an unconditional preservation mode. |
| F. Regression Integrity | **CONDITIONAL PASS** | No tests removed or weakened; historical execution evidence and full-suite reproduction remain limited. |

## 3. Detailed Audit Findings & Answers to A1–A8

### A1 — Does the mechanism protect all active external neighbors?

**Yes, within the stated canonical graph invariants and this resolver’s cleanup pass.**

The new code builds a set of consumed physical node IDs, scans both adjacency lists of every path constituent, excludes consumed and already deleted targets, and registers the remainder before deletion. See [exspander.rs:977](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f7-grok-worktree/src/resolve/exspander.rs:977).

For an incoming edge `X → R`, its canonical twin `R.rc → X.rc` makes `X` visible in one of `R`’s two adjacency lists. The same argument applies to intermediate path nodes and either orientation.

Cleanup excludes registered nodes from orphan deletion at [exspander.rs:1051](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f7-grok-worktree/src/resolve/exspander.rs:1051).

This preserves node sequence; it does **not** reconnect the neighbor to the assembled replacement node.

### A2 — Can incoming branches still be missed or deleted?

No missed immediate incoming neighbor was identified for an in-bounds, twin-complete graph.

Important limits:

- A malformed incoming-only edge without its twin is invisible to this scan. That violates the canonical graph assumption and is not an F7-introduced defect.
- Protection prevents orphan cleanup, not legitimate incorporation into a later splice.
- Protection is local to this invocation. It is not a persistent graph attribute or a guarantee against every downstream filtering operation.

### A3 — Can orientation protect the wrong physical node?

**No orientation error was identified in the new registration.**

`DirUnitigId::id()` removes the strand bit, while `twin()` changes only that bit. Forward and reverse representations share the same physical node and deletion flag. Scanning both adjacency lists and storing physical IDs is therefore appropriate. See [cdbg.rs:18](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f7-grok-worktree/src/graph/cdbg.rs:18).

The regression described below concerns search eligibility, not strand conversion.

### A4 — Does validation failure leave spurious external-neighbor protection?

**The new block does not.**

It is inside the `Ok(assembled_seq)` branch following `validate_and_splice_path`. An `Err` never reaches it.

The existing error branch intentionally adds attempted path nodes and `local_used` IDs to `rejected_nodes`. Exploration can also register rejected candidates before validation. Consequently, rollback does **not** mean that `rejected_nodes` remains identical to its pre-attempt state. F7 introduces no additional external-neighbor registration on validation failure.

The validator borrows the graph immutably; this failure path does not push assembled sequence, delete constituents, or commit `used_in_path`.

### A5 — Is there harmful over-protection?

**Yes: protection also becomes search exclusion.**

The insertion at line 997 affects the later BFS filter:

```rust
if u_node.is_deleted || u_node.is_circular || rejected_nodes.contains(&u.id()) {
    continue;
}
```

See [exspander.rs:245](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f7-grok-worktree/src/resolve/exspander.rs:245).

This conflates two meanings:

1. Retain this node’s sequence during cleanup.
2. Exclude this node from candidate evaluation.

A neighbor requiring the first treatment does not necessarily warrant the second. **F8-P1-001** below supplies a six-node counterexample.

The new policy also retains active artifacts adjacent to committed paths without examining their coverage or edge validity. That is a conservative consequence of “protect every active neighbor,” not independently sufficient evidence of another correctness defect. These nodes are not permanently immortal: the set is invocation-local, and downstream curation still filters them.

### A6 — Are shared neighbors deduplicated with bounded storage?

**Yes.**

- `consumed_nodes` and `rejected_nodes` are hash sets.
- Multiple edges or consumed nodes referencing the same physical neighbor produce one retained set entry.
- Each `neighbor_ids` vector is temporary and proportional to that constituent’s stored degree.
- The protection set is bounded by graph node identities and is dropped when the resolver returns.

Repeated scanning and temporary degree-sized allocations exist, but no unbounded retention independent of graph size was identified.

### A7 — Can a surviving graph node be omitted or altered by raw curation?

**Yes. `min_len = 0` disables only the length threshold.**

[output.rs:42](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f7-grok-worktree/src/polish/output.rs:42) still performs:

- Terminal-overlap trimming and circularity inference.
- Minimum coverage filtering.
- Homopolymer filtering.
- Substring/reverse-complement deduplication.
- Circular-contig deduplication.

A concrete example is Case A with `X.kmer_cov = 1.0`: F7 retains `X` in the graph, but `curate_contigs(&graph, 0, 2.0)` omits it.

These are pre-existing output behaviors, **not new F7 findings**. The new fixtures use coverage 25 and non-homopolymer unique sequences, so their direct curation checks meaningfully avoid those filters.

The incoming-side-branch tests check returned contig sequences; they do not themselves exercise physical FASTA writing.

### A8 — Are graph invariants maintained on success and failure?

**The new code preserves the existing graph mutation model, but violates search-state semantics.**

On success:

- Validation precedes the new registration.
- Registration does not modify sequence or adjacency.
- Existing commit logic marks constituents deleted and appends assembled sequence.

On validation failure:

- The new registration is bypassed.
- Graph mutation and `used_in_path` commit are bypassed.

Successful splicing still leaves stored edges referencing tombstoned constituents and creates an assembled node without boundary rewiring. That is existing behavior, so an assertion that every active node has only active neighbors would be incorrect.

The F7-specific invariant failure is that a cleanup-protected active node becomes unavailable to subsequent BFS candidate evaluation.

## 4. Search Budget Consistency Analysis

### Overlap arithmetic

Both changed locations now evaluate:

```text
sep_next = sep_u + len(u) - overlap
```

before applying the upper bound. This repairs F5-P2-002.

For the fixture:

```text
410 + 200 − 127 = 483 ≤ 600
```

The transition is correctly admitted despite the unadjusted exit being 610.

### Equality boundary

Both locations use `<= max_separation`. Therefore, **`sep_next == max_separation` is admitted**, including as a hidden continuation at the depth cap.

The new F7-B fixtures use `483 < 600`; they do not dynamically test exact equality or `max + 1`.

### Signed arithmetic

`seq_len: u32` and `overlap_k: u8` convert losslessly to `i64`. Negative separation is intentional for overlapping endpoints; unsigned subtraction is not used.

A valid selected overlap is bounded by the source sequence length, making subsequent transition increments nonnegative. With candidate depth four and current call-site initialization, no production-reachable `i64` overflow was identified.

`max_insert` remains an existing `usize as i64` conversion. Arbitrary public inputs above `i64::MAX` are not represented correctly, but that conversion predates F7; no new reachable overflow defect was established.

### BFS versus depth-limit helper

The two paths share these transition rules:

| Rule | BFS expansion | Depth helper |
|---|---|---|
| Source orientation | `edges(u.is_rc())` | Same |
| Target deduplication | Directed target identity | Same |
| Revisit exclusion | Physical node ID | Same |
| Edge selection | `resolve_edge_between` | Same |
| Overlap-adjusted distance | Yes | Yes |
| Inclusive upper bound | Yes | Yes |

They do **not** have identical complete eligibility filters. BFS additionally rejects deleted, circular, rejected, already extended, and unavailable unique nodes when dequeuing. The helper does not receive all those sets.

This can conservatively report truncation for a continuation BFS would later discard. It does not, by itself, falsely certify completeness.

Directed cycles and twin revisits are bounded using physical node IDs. Multiple overlap edges to one directed target are resolved using the same sequence/support selection function in both paths. These restrictions are pre-existing.

### Truncation safety

The existing safeguards remain:

- More than 256 processed states → `Truncated`.
- Reaching 32 candidates → `Truncated`.
- An admissible unvisited continuation beyond depth four → `Truncated`.

However, **F7-P1-001 breaks the broader completeness guarantee**: a newly cleanup-protected node is discarded at dequeue without changing search status. A supported candidate can therefore be omitted while returning `Complete`.

Both forward and backward extensions invoke the same search implementation. No F7 arithmetic divergence between directions was identified.

## 5. Regression Suite Audit & Assertion Rigor

### Commit integrity

Verified:

- `0e513a3 → 4ed4f8c`: only the two regression-test files change.
- `4ed4f8c → 37a7c0b`: only `src/resolve/exspander.rs` changes.
- `37a7c0b → 8ceb3e5`: evidence/documentation additions.

**No pre-existing tests were deleted, skipped, or weakened by F7.** The repaired commit preserves the RED-test assertions.

### Historical RED evidence

The committed [Step 1 record](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f7-grok-worktree/docs/development/evidence/C0_02R/phase_f7/grok_step1_raw.jsonl) reports:

- Cases A–E: failed on unpatched production code.
- Both F7-B tests: failed.
- Cases F and G: passed.
- Strengthened N=66 gate: **already passed**.

Thus, **F5-P2-003 was an assertion-strengthening change, not a demonstrated production RED-to-GREEN repair**. R65 and R66 already survived that baseline fixture.

The files named “raw” contain agent result summaries, not individual command/output transcripts. They support historical attribution but do not independently prove the reported process results. Source tracing is consistent with the seven claimed failures.

### Actual case mapping

The implemented letters differ from the requested mapping:

| Implemented case | Actual coverage |
|---|---|
| A | Terminal incoming `X → R` |
| B | Incoming branch into intermediate `R1` |
| C | One external node connected to two consumed nodes |
| D | Reverse-oriented branch `X.rc → R` |
| E | 64-step prefix, outgoing R65, incoming `X → R64` |
| F | Validation failure, graph fingerprint preservation, later search availability |
| G | Exact ordinary path splice |

Raw-contig nucleotide checks are included across A–E. Reverse orientation is covered, although the exact requested `X → R2.rev` intermediate fixture is not present.

### Assertion strength

The tests are substantive:

- A–D require the expected assembled sequence, branch survival, unique DNA in active graph sequences, and unique DNA in curated output.
- C exercises shared-neighbor registration.
- F compares graph fingerprints and checks that a later supported path remains resolvable after failure.
- G requires one active contig with exact assembled sequence and deleted constituents.
- F7-B ordinary search requires the exact 1603-bp assembly.
- F7-B depth safety requires an unchanged graph fingerprint.

The strengthened N=66 gate strictly checks:

- `short_survivors == 2`.
- R65/R66 identities and exact sequences.
- Their retained nucleotide tails in raw contigs.
- Exact reference-substring placement and full reference union coverage.

Some tail assertions repeat implications of stronger sequence-equality assertions, but the gate as a whole is not tautological.

### Execution evidence

This session executed existing release binaries:

| Suite | Passed |
|---|---:|
| Incoming side branches | 9 |
| Graph atomicity | 6 |
| Cleanup | 6 |
| P1 regression | 7 |
| Directed edge identity | 8 |
| Frontier node-count audit only | 1 |
| **Total** | **37** |

`git diff --check 0e513a3 37a7c0b` also passed.

The [F7 execution report](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f7-grok-worktree/docs/development/evidence/C0_02R/phase_f7/phase_f7_grok_execution_report.md) claims **138 passed, 0 failed, 0 ignored**. That full result was **not independently reproduced** here. Prebuilt binaries also lack a fresh source-to-build attestation from this session.

## 6. Secondary Defects Ledger (P0/P1/P2)

### F8-P1-001 — Cleanup protection suppresses valid competing candidates

- **Severity:** P1 correctness.
- **File:line:** [exspander.rs:997](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f7-grok-worktree/src/resolve/exspander.rs:997); downstream exclusion at line 246.
- **Root cause:** Successful-commit neighbor protection writes into the same set used to reject BFS candidates.
- **Production reachability:** Multiple length-ordered seeds sharing a short external neighbor. No malformed edges, missing twins, cycles, or invalid overlaps are required.

**Minimal counterexample**

Use six active, non-circular nodes, all with coverage 25:

```text
S1(900) ──→ R(169)
              ↑
            X(100)
              ↑
S2(700) ──→ Q(100) ──→ D(100)
```

Edges are:

```text
S1 → R
X  → R
S2 → Q
Q  → X
Q  → D
```

Every edge has overlap 31 and its authentic twin.

A sequence-valid construction uses distinct 31-mers `a`, `b`, `c`:

```text
S1 = random669?  [use 869 bp filler] + a     = 900 bp
R  = a + 138 bp filler                     = 169 bp
X  = b + 38 bp filler + a                  = 100 bp
S2 = 669 bp filler + c                     = 700 bp
Q  = c + 38 bp filler + b                  = 100 bp
D  = b + 69 bp filler                      = 100 bp
```

For paired evidence from S2:

- X: 20 distinct supporting pairs.
- D: 5 distinct supporting pairs.
- Both have endpoint distances 50 and 50, with physically compatible 50-bp reads.
- Both path separations are `100 − 31 − 31 = 38`.
- Both spans are `50 + 38 + 50 = 138`.
- Use median 138, MAD 40, and insert bounds 50–278.

The execution follows directly from the code:

1. S1 is processed first because it is longer.
2. `S1 → R` validates and commits.
3. The new neighbor scan inserts X into `rejected_nodes`.
4. S2 searches the fork at Q.
5. BFS dequeues X and skips it solely because it is protected.
6. D remains the only scored candidate.
7. Search returns `Complete`; D’s five centered pairs satisfy the single-candidate acceptance rule.
8. F7 commits `S2 → Q → D`.

With F5B production code, X is not newly registered by the first commit. Both candidates are scored; X wins with weight ratio `20 / 5 = 4`, exceeding the default 2.5 threshold.

**Impact:** F7 changes a supported winner into an excluded candidate and commits the weaker alternative.

**Existing test coverage:** Case F checks availability after a **failed** first attempt. It does not test availability after a **successful** first commit, which is where the new block executes.

**Minimal repair direction:** Keep successful-commit neighbor retention in a separate cleanup-protection set. Consult that set during orphan cleanup without using it to exclude BFS candidates. Add the successful-first-commit counterexample as a regression gate.

**Additional findings:** No additional F7-introduced P0 or independent P2 defect was established. Output filtering and tombstone-edge behavior described above are existing limitations, not newly attributed defects.

## 7. Final Recommendations & Transition Gate

**Keep C0-02R NOT ACCEPTED. Phase F9 should repair F8-P1-001 before acceptance.**

Required transition evidence:

1. Separate successful-commit cleanup retention from candidate rejection.
2. Demonstrate RED then GREEN for the six-node successful-first-commit counterexample.
3. Preserve all current F7 assertions and rerun the full release suite from rebuilt sources.
4. Add exact-budget equality and one-base-over-budget checks for both search expansion and depth-cap detection.
5. Save actual commands, exit codes, and test output alongside commit identities.

The incoming-neighbor enumeration and overlap-budget repairs can remain. Acceptance is blocked by the new coupling between sequence retention and candidate eligibility.