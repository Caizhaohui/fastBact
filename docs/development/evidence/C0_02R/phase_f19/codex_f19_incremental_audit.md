## 1. Executive Summary & Audit Baseline

| Field | Value |
|---|---|
| Audit ID | **AUDIT-20261022-014** |
| Auditor | Codex CLI |
| Model / Reasoning Effort | `gpt-6.1-sol` / `medium` — requested audit metadata |
| Sandbox | `read-only` |
| Overall Verdict | **FAIL** |
| C0-02R status | **NOT ACCEPTED** |

**Four findings are resolved. F18-P2-005 remains partially unresolved because T39 reaches queue saturation but cannot detect removal of the queue-capacity guard.** The independent 256-state limit produces the same asserted result.

No new P0 or P1 defect was confirmed. Production code is unchanged.

### Commit identity discrepancy

Both full commit hashes supplied in the request are unavailable in this repository. The short IDs resolve to commits whose subjects match the intended audit:

| Role | Actual repository commit |
|---|---|
| Baseline | `39f11ca90a5de13a5b3d5572d29c191c40d9c03d` |
| Audited HEAD | `a9930471a73612159a83694d3d661f060042efc8` |
| Branch | `feature/c0-02r-phase-f19-codex` |

This report audits **these actual commits**, using `git diff 39f11ca..HEAD`. It cannot certify the unavailable full hashes.

### Evidence and execution limits

The audit used direct diff inspection, Rust source tracing, graph invariants, separation arithmetic, and test sensitivity analysis.

Existing test executables were run successfully:

| Suite | Result |
|---|---:|
| Repeat completeness gate | **8 passed** |
| State semantics gate | **42 passed** |

These were prebuilt binaries. No fresh compilation or executable mutation test was performed in the read-only sandbox; their correspondence to the audited source was not independently rebuilt. The conclusions below rest on source analysis, with binary execution as supporting evidence.

---

## 2. F18-P2-001: Oracle Reciprocal Selection & Tie Rejection

**Disposition: RESOLVED.**

Source: `tests/test_c0_02r_repeat_completeness_gate.rs:101–163`; controls begin at line 760.

The repaired algorithm:

1. Finds all sequence-valid forward edges to the requested target.
2. Requires exactly one forward maximum-support edge.
3. Inspects all reverse edges from `to.twin()` targeting `from.twin()`.
4. Validates each reverse overlap independently, without prefiltering by `chosen_ov`.
5. Requires exactly one reverse maximum-support edge.
6. Requires its overlap and support to equal the forward winner.

Writing \(F\) and \(R\) for the sequence-valid forward and reverse candidate sets, acceptance now requires:

\[
|\operatorname{argmax}_{e\in F}\operatorname{support}(e)|=1,
\qquad
|\operatorname{argmax}_{e\in R}\operatorname{support}(e)|=1,
\]

and equality of both selected overlap and support.

This agrees with the production path validator’s reciprocal-selection requirement at `src/resolve/exspander.rs:1153`.

### Control coverage

| Case | Fixture | Expected |
|---|---|---|
| Symmetric positive | Forward and reverse `(31,10)` | `Some(31)` |
| Asymmetric support | Forward `(31,10)`, reverse `(31,5)` | `None` |
| Forward tie | Forward `(31,10)` and `(30,10)` | `None` |
| Stronger competing reverse overlap | Reverse `(31,10)` and `(30,20)` | `None` |
| Reverse tie | Reverse `(31,10)` and `(30,10)` | `None` |

The homopolymer sequences make both overlap lengths sequence-valid. Thus the negative cases exercise selection and reciprocity rather than failing incidentally on sequence mismatch.

All five mandated controls are present and meaningful.

---

## 3. F18-P2-002: Oracle Initialization Long-Fork Reuse

**Disposition: RESOLVED.**

Source: `tests/test_c0_02r_repeat_completeness_gate.rs:236`, `:285`, `:311`; new differential test begins at line 785.

Initial admission now rejects a target when:

\[
L_{\text{target}}\ge500
\quad\land\quad
(\text{target.id}=\text{curr.id}\ \lor\ \text{target.id}=\text{fork.id}).
\]

Normal expansion and depth-bound continuation checks also prohibit long targets already represented by `curr_dir`, `fork_node`, or the candidate path.

These comparisons use physical node IDs, so changing strand does not evade the uniqueness restriction.

### Differential fixture

The graph contains:

```text
S → F
    ↺
    → D
```

For \(L_F=600\):

- The initial `F → F` hop is rejected.
- No oracle candidate path contains `F`.
- The direct `[D]` candidate remains.
- The oracle reports completeness.
- Production commits the valid direct route, consuming D.

For \(L_F=100\):

- The oracle admits the unrolled `[F,D]` candidate.
- Production withholds D because repeat exploration exposes supported alternatives.

The API distinction matters: `reference_oracle_search` returns an `OracleReport`, not `Ambiguous` or `SafePrune`. The long **self-loop route** is rejected; the entire search is not rejected. Retaining the legal direct D route is correct.

---

## 4. F18-P2-003: T30 Decision Reachability & Sensitivity Control

**Disposition: RESOLVED.**

Source: `tests/test_c0_02r_state_semantics_gate.rs:3066`.

### Junction validity

`P1` ends with `j_b`, and X now begins with `j_b`. The `P1 → X` edge therefore has a valid 31-base overlap.

C now has length 129, giving:

\[
L_C-31=98.
\]

### Decision 1 reaches P2

The initial S-anchored search through P1 reaches P2 at:

\[
\operatorname{sep}(P2)=(80-31)-31=18.
\]

Its evidence span is:

\[
110+18+110=238\in[180,320].
\]

The potentially troublesome return to P1 occurs at separation:

\[
18+(80-31)+(129-31)=165.
\]

For this decision, \(\beta(S)=220\), so:

\[
165+220=385>320.
\]

The repeat return is safely pruned, avoiding the earlier S-anchored false truncation. The test also explicitly asserts that P1 was consumed.

### Decision 2 has a real competitor

With P2 as anchor:

| Endpoint | Separation | Evidence span |
|---|---:|---:|
| D | \(-31\) | \(140-31+140=249\) |
| X via C, P1 | \(-31+98+49=116\) | \(65+116+65=246\) |

Both spans lie in `[180,320]`. X has 20 supporting pairs versus D’s 5. Their Gaussian weights are near their maxima around median 250; X provides substantial competing evidence.

### Sensitivity control

The control uses `graph.clone()`, preserving all six directed transitions and their reciprocal edges. It retains the S→P2 and P2→D evidence and omits only P2→X evidence.

Assertions establish:

- Main fixture: P1 consumed; D withheld.
- Identical-topology control: D consumed.

This removes the previous confounding change in topology and establishes downstream evidence sensitivity.

---

## 5. F18-P2-004: T37 β Boundary Sensitivity

**Disposition: RESOLVED.**

Source: `tests/test_c0_02r_state_semantics_gate.rs:3492`; production classification at `src/resolve/exspander.rs:288`.

The new topology is:

```text
S → D
└ → P ↺
```

For both tested lengths, P is short. On `P → P`, P is already visited and `step_depth >= 1`, satisfying the production repeat-reentry predicate.

The evidence floor is:

\[
\beta(S)=\min(110+130,\;80+82)=162.
\]

Direct S→P evidence is outside the insert window:

\[
80-31+82=131<180.
\]

D remains a valid candidate:

\[
110-31+130=209\in[180,250].
\]

For the repeated P hop:

\[
\operatorname{sep}_v=L_P-62,
\qquad
\operatorname{sep}_v+\beta=L_P+100.
\]

| \(L_P\) | \(\operatorname{sep}_v+\beta\) | Classification | Consequence |
|---:|---:|---|---|
| 150 | 250 | `Probe` | Repeated P scores at span 250; search truncates; D withheld |
| 151 | 251 | `SafePrune` | Repeat route cannot score; search completes; D commits |

The raw separation bound does not distinguish these cases: both 88 and 89 are below 250. The distinguishing condition is specifically the β comparison.

Changing the production comparison from `>` to `>=` would prune the equality case and permit D to commit, contradicting the test assertion. This is a valid sharp-boundary regression test.

---

## 6. F18-P2-005: T39 Queue Capacity Limit Saturation

**Disposition: PARTIALLY RESOLVED; remains open at P2.**

Source: `tests/test_c0_02r_state_semantics_gate.rs:3587`.

### Saturation is now reachable

The fixture creates exactly N outgoing targets: D plus \(N-1\) terminal branches.

All initial hops have separation \(-31\), valid overlaps, and admissible targets. Production fills the queue before processing it.

For \(N=1025\), the attempt to enqueue target 1025 encounters:

```rust
if queue.len() >= MAX_QUEUE_CAPACITY {
    search_status = SearchStatus::Truncated;
    break;
}
```

The capacity is 1024. Therefore the queue-capacity truncation branch is reached.

D has authentic support:

\[
160-31+171=300,
\]

equal to the insert median. The test asserts D remains unconsumed and total joins equal zero.

For \(N=50\), the queue and state limits are not reached; D is the sole supported candidate and commits.

### Remaining defect: another limit masks failure of the queue guard

Production also defines:

```rust
const MAX_DFS_STATES: usize = 256;
```

and truncates when `states_visited > MAX_DFS_STATES` at `src/resolve/exspander.rs:381`.

Consider removing only the initial queue-capacity guard:

1. All 1025 terminal targets enter the queue.
2. D is processed first and becomes a valid candidate.
3. Processing state 257 triggers the independent state limit.
4. Search returns `Truncated`.
5. The caller withholds D and records zero joins.

**Every current T39 assertion still passes.** The 50-branch control also remains unchanged.

This is a source-derived mutation counterexample, not an executed mutation experiment.

T39 now covers the queue-saturation branch, correcting the former reachability defect. It still does not verify that queue-capacity enforcement causes the observed withholding. The advertised regression gate remains insensitive to removal of that enforcement.

To close the finding, add a focused check that distinguishes queue-capacity truncation from state-budget truncation, such as an internal search test that observes the truncation cause or queue admission bound. The production public contract need not change.

---

## 7. Regression & Contract Integrity Assessment

### Production immutability

Confirmed:

```text
git diff 39f11ca..HEAD -- src/ | wc -l
0
```

Both commits have the same `src/` tree object:

```text
30967bf4940844aeb6b744732bcad6d2a72daab4
```

Only the two scoped test files differ. Tracked `src/` and `tests/` working-tree contents match HEAD.

Consequently, this incremental change does not weaken production contracts concerning reciprocal edges, long-node uniqueness, repeat admission, separation bounds, search truncation, or commit validation.

### Regression evidence

- Existing oracle binary: **8/8 passed**.
- Existing state-semantics binary: **42/42 passed**.
- `git diff --check 39f11ca..HEAD`: passed.
- No new P0 or P1 defect confirmed.
- No additional incremental P2 defect confirmed beyond the incomplete T39 repair.

Fresh build, compiler diagnostics, broader regression execution, and compiled mutation testing remain unverified in this read-only audit.

The pre-existing untracked `docs/development/evidence/C0_02R/phase_f19/` directory was left untouched.

---

## 8. Remaining Findings Ledger

| Severity | Finding | Status | Assessment |
|---|---|---|---|
| P0 | None confirmed | — | Within audited scope |
| P1 | None confirmed | — | Within audited scope |
| P2 | F18-P2-001 | **Closed** | Reciprocal selection and all five controls repaired |
| P2 | F18-P2-002 | **Closed** | Long-fork reuse rejected; short-loop control retained |
| P2 | F18-P2-003 | **Closed** | T30 reaches the intended decision and isolates competitor evidence |
| P2 | F18-P2-004 | **Closed** | T37 distinguishes the exact β boundary |
| P2 | F18-P2-005 | **Open, partially repaired** | Queue saturation reached, but state-limit truncation masks removal of the queue guard |

**Remaining open findings: P0 = 0, P1 = 0, P2 = 1.**

---

## 9. Final Recommendation for Stage C0-02R

**FAIL the claim that all five F18 P2 verification findings are fully repaired.** Accept F19-A through F19-D as resolving their respective findings; require a queue-specific sensitivity assertion for F19-E.

Reconcile the full commit hashes before recording this report as evidence for an exact revision. After T39 is strengthened, rerun both suites from freshly compiled audited sources and demonstrate that disabling the queue-capacity guard fails the focused check.

**C0-02R remains NOT ACCEPTED.**