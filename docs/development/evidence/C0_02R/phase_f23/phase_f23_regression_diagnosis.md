# Phase F23 — Release-Suite Regression Diagnosis (Pre-Repair Record)

**Worktree**: `fastBact-phase-f23-unique-anchor-worktree` (branch `feature/c0-02r-phase-f23-unique-anchor`, base `2ae4acc`, production src tree `30967bf4...`)
**Patch**: uncommitted `src/resolve/exspander.rs` diff, archived at `phase_f23_option_a_uncommitted.diff` (SHA-256 `87f3514b...`)
**Release suite job**: `3010914` (vnode1) — 38/42 state-semantics pass, **4 FAILED**; all other suites green; integration e2e green (1.54s).
**Status**: `PHASE_F23_REGRESSION_BLOCKED` per coordinator directive — no GREEN acceptance, no merge.

## Common root cause (working hypothesis, Codex to verify)

The F23 patch applies the unique-anchor substitution **unconditionally** at every
fork reached after walking through any short (<500 bp) node: candidate exits are
scored exclusively with pairs keyed by the last long (≥500 bp) anchor
(`score_from = unique_anchor`, `sep + anchor_gap`), both in
`find_candidate_paths` and in `pruned_repeat_hides_supported_endpoint`.

The accepted F15/F16/F18-era semantics contract instead requires that a fork
**within** the unique anchor's insert horizon be scored anchored at the **current
walk node** (`curr_dir`), even when that node is short. T30 encodes this
explicitly ("Decision 2 from P2: evaluated anchored at P2!"); T26/T27/T34 are
positive controls whose evidence rows are keyed only by the short current node.
The DEC_000022 hijack that F23 must eliminate arises specifically when the walk
has traversed past the anchor's insert horizon (anchor→fork overlap-adjusted
distance 936 bp > max_insert 559 bp), where anchored evidence is
arithmetically impossible and repeat-internal links are the only — locus-
ambiguous — evidence (F22 counterfactual §3.1 proves the case by horizon
arithmetic alone; Rules 1/2 length-gating was never exercised by a counterfactual
within the horizon).

## Per-failure runtime reconstruction (deterministic fixtures; scoring model
## verified against execution by probe job 3010539)

All four: `InsertSizeDist` windows per test; `ExspanderParams::default()`
(min_support 5, decision_threshold 2.5); every failure has
`search_status = Complete`; every rejection under the patch is
`WITHHELD_NO_UNIQUE_ANCHOR` (candidate set empties under anchor-keyed scoring).

### T26 `test_t26_prefix_reuse_exceeding_max_budget_resolves_cleanly` (assert :2865)
- Topology S(700)→P(100)→C(100); C→P(reuse)→X(300); C→D(300). Window [150,250].
- PE rows: (S,C)×20 span 200 ✓; (C,D)×10 span 200 ✓; (C,X)×20 span 251 ✗(>250).
- Fork 1 @P, curr=S(anchor inert): CAND_01=[C] 20 pairs w=1.0; CAND_02=[X] 0 → commit [C]. Anchor state: S, gap=138.
- Fork 2 @C (Case B, C→{P,D}), curr=C, gap=138:
  - Pre-patch: [D] sep −31 → (C,D) 10 in-window → single supported → D commits (deleted ✓).
  - Patched: [D]→(S,D,sep 107) 0 rows; [P,X]→(S,X,sep 176) 0 rows → empty → withhold.
- Assert expects D deleted; observed D alive → FAIL.

### T27 `test_t27_legitimate_prefix_reuse_bridge_commits_correctly` (assert :2910)
- Topology S(700)→P(100)→C(100); C→P(reuse)→X(300). Window [150,250].
- PE rows: (S,C)×20 span 200 ✓; (C,X)×25 span 200 ✓.
- Fork 1 @P, curr=S: commit [C]; gap=138, anchor=S.
- Fork 2 Case A @target=P, curr=C (P→{C(revisit→Probe),X}):
  - Pre-patch: [X] sep 38 → (C,X) 25 in-window w=1 → commit → contig S-P-C-P-X, X consumed ✓.
  - Patched: (S,X,sep 176) 0 rows; probe reentry check scores (S,C) 0 → no candidates → withhold.
- Assert expects X consumed; observed withheld → FAIL.

### T30 `test_t30_multiple_short_prefix_nodes_in_cycle` (assert :3132, control leg)
- Topology S(700)→P1(80)→P2(80); P2→{D(300),C(129)}; C→P1(reuse); P1→X(300). Window [180,320].
- PE main: (S,P2)×20 span 238; (P2,D)×5 span 249; (P2,X)×20 span 246. Control: (P2,X) removed.
- Fork 1 @P1, curr=S: [P2] 20 ✓ vs [X] 0 → commit [P2]; gap=98, anchor=S.
- Fork 2 @P2 (Case B): contract comment: "evaluated anchored at P2!".
  - Pre-patch main: [D] 5 w≈4.99 vs [C,P1,X] 20 w≈20 → X route commits; D alive ✓ (assert passes).
  - Pre-patch control: [D] 5 sole supported → D commits (deleted) ✓.
  - Patched (both): (S,D,67)/(S,C,67)/(S,X,214) all 0 rows → withhold → control D never commits.
- Control assert expects D deleted; observed D alive → FAIL.

### T34 `test_t34_short_current_node_reentry_positive_control_commits` (assert :3361)
- Topology S(700)→P(100)→{C(100),X(300)}; C→{C(self),D(300),X}. Window [150,250].
- PE rows: (S,C)×20 span 200; (C,D)×5 span 200; (C,X)×20 span 188.
- Fork 1 @P, curr=S: commit [C]; gap=138.
- Fork 2 @C (Case B, self-reentry): pre-patch commits the C→C→X bridge at 3.69× over D → contig S-P-C-C-X ✓.
  Patched: (S,·) rows absent → withhold → contig never created → FAIL.

## Candidate minimal repair (for Codex review)

1. Revert the unconditional `score_from`/`anchor_offset` substitution in
   `find_candidate_paths` (both scoring sites) and in
   `pruned_repeat_hides_supported_endpoint` — within-horizon forks return to
   pre-F23 `curr_dir`-anchored scoring bit-for-bit.
2. Keep the anchor bookkeeping (`unique_anchor`, `anchor_gap`) and
   `JunctionStats::withheld_no_unique_anchor`.
3. Add a horizon gate at both fork sites (Case A and Case B): withhold with
   `withheld_no_unique_anchor += 1` iff `unique_anchor != curr_dir` and the
   overlap-adjusted anchor→fork distance exceeds `dist.max_insert`
   (Case A: `anchor_gap + target.len − ov`; Case B: `anchor_gap`).
   This fires exactly where the F22 counterfactual proved DEC_000022
   (936 > 559) and nowhere in the accepted T-series (fork_dist ≤ 138 ≤ 250).
4. β(curr) / `classify_hop` / decision thresholds untouched.
5. Revise this session's own gate test (`test_c0_02r_unique_anchor_gate.rs`,
   not a historical fixture): horizon-exceeded scenario (chain 3×300, gap 810 >
   700) must withhold; add a within-horizon passthrough control asserting the
   T34-style bridge still commits. No historical test is modified, weakened,
   or skipped.

Known scope boundary to document: within-horizon forks anchored on a genuinely
multi-copy short node remain decided on repeat-internal evidence (pre-existing
baseline behavior; copy-number awareness was Option B/D territory, rejected in
the F22 strategy matrix).
