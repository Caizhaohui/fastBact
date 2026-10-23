ROOT_CAUSE_CONFIRMED

**Q1. Yes—the diagnosis identifies the common cause.**

The patch changes the scoring source whenever `unique_anchor != curr_dir`, without checking distance: [exspander.rs:353](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f23-unique-anchor-worktree/src/resolve/exspander.rs:353). It applies that substitution to candidate scoring, self-return handling, and the repeat probe at lines 412–420 and 584–585.

That contradicts the historical fixtures’ current-node evidence contract. T30 states it explicitly at [state_semantics_gate.rs:3106](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f23-unique-anchor-worktree/tests/test_c0_02r_state_semantics_gate.rs:3106); T26/T27/T34 supply the decisive evidence under C or P2, not S.

One diagnosis detail needs correction: **T27’s second fork distance is 207 bp, not ≤138 bp**. Its `anchor_gap` is 138, but Case A adds `100 − 31`. This remains within 250 bp and does not change the root cause.

**Q2. Yes—the repair restores these four contracts without weakening their assertions.**

This requires restoring **all** substituted scoring and self-return behavior, including the probe, and placing the gate before candidate search.

| Test | Repaired control flow |
|---|---|
| **T26** | First decision remains S-anchored and commits P,C. At C, Case B distance is `2 × (100−31) = 138 ≤250`; gate passes. C→D has 10 pairs at span 200. C→P→X has span 251 and contributes nothing; the unchanged β-bound safely prunes the unsupported reentry. Search completes and D commits. |
| **T27** | First decision commits P,C. Second decision is Case A through P: `138 +100−31 =207 ≤250`. Restored C-keyed scoring gives X 25 pairs at span 200. The C-return has no self evidence, and further reentry is safely bounded. Search completes; S-P-C-P-X commits, satisfying both consumption and exact-sequence assertions. |
| **T30** | Initial S-anchored decision commits P1,P2. At P2, Case B distance is `2 ×(80−31)=98 ≤320`. Main fixture: D has about 4.997 weight; C-P1-X has about 19.819, ratio ≈3.97, so the X route wins and D survives. Control: removing P2→X evidence leaves D as the sole supported candidate, so D commits. Both assertions retain their intended sensitivity. |
| **T34** | Initial decision commits P,C. At C, `138 ≤250`; gate passes. Direct X is out of window, but C→C→X has span 188 and weight ≈18.43 against D’s 5: ratio ≈3.69. Restored self-return/probe semantics preserve completeness, and S-P-C-C-X commits while D remains protected. |

These are static control-flow conclusions, not a claim that I executed the repaired implementation.

**Q3. Yes, given the supplied 936 bp geometry; independent trace verification is unavailable.**

The F22 matrix and TSV are absent at the supplied paths and were not found in the nearby F22 diagnostics worktree. I therefore cannot independently verify DEC_000022’s row or reconstruct its component lengths.

The formulas agree with the bookkeeping:

- **Case B:** `fork_dist = anchor_gap =936`; `936−559=377 bp` beyond the horizon.
- **Case A:** `fork_dist = anchor_gap + target_len − incoming_overlap`. If 936 measures distance to the branching target’s end, this expression must equal 936, and the same rejection follows.
- If 936 instead measures distance to Case A’s current node, the target’s overlap-adjusted contribution must additionally be included.

Thus the supplied mechanism is intercepted when the walk still carries u241 and `unique_anchor != curr_dir`. **The 2.819 winning ratio never reaches the decision stage.** Component-level verification of the Case A sum remains outstanding because the trace row is missing.

**Q4. No historical regression fixture was identified, with limits on audit completeness.**

The named families behave as follows:

- **`two_copy_repeat`:** both lengths are examined from a long flank, so `unique_anchor == curr_dir` at the initial Case A decision. The 240 bp interior produces a 300 bp repeat and 270 bp fork distance. The 1440 bp interior remains rejected by existing candidate separation logic. Neither newly trips the gate. See [pe_evidence.rs:36](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f23-unique-anchor-worktree/tests/test_c0_02r_pe_evidence.rs:36).
- **65-node chain and frontier sweeps:** accumulated gaps exceed the horizon, but these are unbranched chains. A fork-only gate does not fire.
- **Depth-chain tests:** the depth cap concerns candidate exploration from the first seed-anchored fork; explored candidate length is not committed `anchor_gap`. The gate must not be inserted inside that exploration.
- **Repeat-completeness suite:** the tandem/depth/separation fixtures start their decision at S→Q with S still current; direct forks also start on S. Existing completeness and β behavior therefore remains decisive.
- **Other inspected hand-built suites:** edge identity, atomicity, cleanup, incoming-side-branch, and regression fixtures supplied no identified new horizon rejection. In the cleanup reuse control, reaching long B resets the anchor before the next fork.
- **Exspander unit module:** the lifecycle test never reaches extension. The 350 bp repeat test’s intended geometry fits its approximately 450 bp library; however, its graph and distribution are generated at runtime, so source inspection alone does not certify every intermediate fork.

There is also a request inconsistency: **T26–T39 are in `state_semantics_gate.rs`, not `pe_evidence.rs`**. The explicit 2600–3400 bound excludes the full T35–T39 bodies. Consequently, I cannot honestly certify an exhaustive all-suite audit under that bound. No repaired-suite run was performed.

**Q5. No objection, provided the narrowed contract is stated explicitly.**

Changing only the new session test preserves every historical assertion.

The corrected scenario should distinguish:

- gap to C3: `3 ×(300−30)=810 bp`;
- Case A distance to F, retaining F=200: `810+200−30=980 bp >700`.

Keep the no-cross-fork-join and exit-retention assertions, and assert the withholding counter. Add the T34-style positive control.

However, replacing the current within-horizon unique-anchor test **does relinquish a new F23 guarantee**: that unique-anchor evidence overrides repeat-internal evidence within the horizon. It weakens no historical contract, but the resulting implementation must be documented as a **horizon guard**, rather than full enforcement of the currently stated Rules 1–2.

**Q6. Residual risks**

- **Boundary over-withholding — acceptable-with-documentation.** Fork distance measures distance to the fork’s end; candidate separation subtracts the outgoing overlap. Therefore `fork_dist > max_insert` alone is not universally proof that every possible anchor-keyed pair is out of window. The gate is a conservative policy and may sacrifice legitimate resolution near the boundary. Preserve strict `>` so equality passes.
- **Within-horizon multi-copy hazard — acceptable-with-documentation for this bounded repair.** Restoring current-node scoring preserves the pre-existing possibility of a repeat-internal winner choosing the wrong locus. The current new scenario 2 demonstrates the distinction. This would be **blocking for a claim of complete Option A enforcement**, but not for the explicitly narrowed beyond-horizon repair.
- **Length as a uniqueness proxy — acceptable-with-documentation.** The bookkeeping records the last ≥500 bp node; the reviewed evidence does not establish that every such node is biologically unique.

**VERDICT: ENDORSE_REPAIR**

This endorses the repair design, with the scope and verification limits above; it is not a release-green certification.

Minimal patch sketch:

```rust
// Case A: inside next_continuations.len() > 1, before search.
let fork_dist =
    anchor_gap + i64::from(target_node.seq_len)
        - i64::from(edge_to_target.overlap_k);

if unique_anchor != curr_dir && fork_dist > dist.max_insert as i64 {
    stats.withheld_no_unique_anchor += 1;
    cleanup_protected_nodes.insert(target.id());
    for &(candidate, _) in &next_continuations {
        cleanup_protected_nodes.insert(candidate.id());
    }
    break;
}

// Case B: before search.
if unique_anchor != curr_dir && anchor_gap > dist.max_insert as i64 {
    stats.withheld_no_unique_anchor += 1;
    for &(candidate, _) in &valid_continuations {
        cleanup_protected_nodes.insert(candidate.id());
    }
    break;
}
```

Restore the pre-F23 search/probe signatures, `curr_dir` scoring, separations, and self-return condition. Keep anchor updates. Remove the existing empty-candidate counter increments and update its documentation so it counts horizon rejections only. Leave β, `classify_hop`, thresholds, and historical tests unchanged.