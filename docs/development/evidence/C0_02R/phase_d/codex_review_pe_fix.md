VERDICT:
**FAIL**

The orientation changes are correct for inward-facing paired-end reads. The new distance condition is only a coarse topology filter; it does not establish that the counted fragments support the selected path.

Line numbers below refer to commit `f241d9d3c9d21470f63bd85b5e1bd6eb0bfc8d30`. The three scoped files are unchanged between that commit and the current HEAD.

FINDINGS:

**P1 — Unitig-pair counts can support a path that is incompatible with the authentic fragments**

- **Location:** [exspander.rs:223](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/src/resolve/exspander.rs:223), lines 241–244; [mapper.rs:144](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/src/resolve/mapper.rs:144), lines 158–159.
- **Necessary conditions:** Authentic pairs connect the same directed flanks through another graph traversal. Their aggregate count wins, and the selected shortcut satisfies `gap <= max_insert`.
- **Minimal counterexample:** Use exact, reciprocal 30-base overlaps and this topology:

  ```text
  A  ──→ R ──→ B
         └──→ C ──→ B
  A2 ──→ R
  ```

  Set `|A|=2000`, `|R|=120`, `|C|=270`, `|B|=1000`; make C’s prefix and suffix match R’s suffix/B’s prefix. The two repeat occurrences genuinely traverse `A→R→C→B` and `A2→R→B`.

  Supply five distinct authentic 600-base fragments from the first traversal, with 150-base mates at A’s end and B’s start. They map to `(A,B)` and its reciprocal. On that authentic traversal:

  ```text
  fragment span = 150 + (120−30−30) + (270−30) + 150
                = 600
  ```

  The resolver selects the direct exit B, whose inferred span for those same read placements is:

  ```text
  fragment span = 150 + 60 + 150 = 360
  ```

  Nevertheless, support is five versus zero, and `gap=60` passes the cutoff. Exact sequence-overlap validation also passes. The fragments support the longer traversal, but their counts authorize the shortcut.

  This is production-reachable through `map_paired_reads`; it does not require manually forging matrix entries.
- **Existing tests:** **Not covered.** They inject counts directly and have no alternate traversal or read coordinates.
- **Fix direction:** Retain oriented read placements, offsets, lengths, and fragment identity or an equivalent distance-evidence representation. Score support against each candidate path’s implied fragment span before applying support and confidence thresholds.

**P2 — The empirical insert estimator uses k-mer positions as read positions**

- **Location:** [mapper.rs:116](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/src/resolve/mapper.rs:116), lines 123 and 147–153; the resulting bound is consumed at `exspander.rs:244`.
- **Necessary conditions:** At least 20 qualifying intra-unitig pairs are mapped, activating empirical estimation.
- **Minimal counterexample:** For a forward R1/reverse R2 pair with 150-base reads, `k=31`, true fragment span 450, and successful first anchors:

  ```text
  R1 anchor coordinate = S
  R2 anchor coordinate = S + 450 − 31
  reported insert      = (450 − 31) + 150 = 569
  ```

  Twenty such pairs produce `median=569`, `mad=20`, and `max_insert=639`. Correctly measured, the same estimator policy would give `median=450` and `max_insert=520`.

  For sampled read offsets `t1,t2`, this configuration instead gives:

  ```text
  reported insert = true span + read1_length − k − t1 − t2
  ```

  Thus anchor availability changes the apparent distribution. The intra-unitig branch also does not require inward-facing orientation or proper ordering.
- **Existing tests:** **Not covered.** None calls `map_paired_reads` or checks estimated insert sizes.
- **Fix direction:** Recover oriented read intervals from anchor coordinates and sampled offsets, validate pair orientation/order, and measure fragment endpoints. Estimate the distribution only from qualifying concordant pairs.

This estimator defect predates the commit, but the commit newly relies on its upper bound for rejecting repeat bridges.

A–B. Directed orientation and counting

At `mapper.rs:80–83` and `119–123`, the graph and read canonical-strand bits combine to recover the read’s direction relative to the stored unitig.

For an inward-facing pair, recording:

```text
(u1, twin(u2))
(u2, twin(u1))
```

at lines 158–159 correctly represents the path and its reverse complement. This also handles either physical strand carrying R1.

`get_support()` at lines 48–49 now performs an exact directed-key lookup:

- No forward/reverse orientation mixing.
- No normalization.
- No summation of reciprocal entries.
- No double counting caused by the two production insertions: inter-unitig endpoints have different node IDs, so those keys differ.

`PairLinkMatrix` counts input pairs; it does not deduplicate repeated input fragments or validate their distance. Those are separate limitations, not orientation mixing in `get_support()`.

C–D. What `gap` measures

For oriented `A→R→B`, let the overlaps be `o₁,o₂`. In the spelled path:

```text
start(B) = end(A) + |R| − o₁ − o₂
gap      = start(B) − end(A)
```

Therefore, `gap` is the **signed separation between A’s end and B’s start**. It is not the fragment span.

For a forward read starting at coordinate `s` in oriented A and its reverse mate ending at coordinate `e` in oriented B:

```text
fragment span = (|A| − s) + gap + e
```

The implementation discards the information needed to evaluate this expression.

Consequently:

1. **Repeat span versus fragment span:** `gap` measures boundary separation, not fragment span.
2. **False acceptance:** **Yes.** The P1 counterexample has authentic counted pairs, but none supports the shortcut at the actual library span.
3. **False rejection:** **Possible.** With fewer than 20 qualifying calibration pairs, the fallback upper bound is 700. A real 1100-base library can supply five valid pairs across `gap=800`, with 150 bases on each flank, yet the bridge is rejected. Conversely, with an accurate hard upper bound and valid placements, rejecting `gap > max_insert` is a sound coarse exclusion.
4. **Distribution and variance:** The cutoff uses only `max_insert`. MAD influences that value indirectly, but no candidate fragment span is tested against the distribution. Ignoring `min_insert` is problematic because actual fragment spans are never reconstructed; comparing **gap itself** against `min_insert` would also be incorrect.
5. **Arithmetic:** Negative gaps are legitimate. For example, `|R|=31` and two 30-base overlaps give `gap=-29`, meaning the flanking sequence intervals overlap. It does not imply a negative fragment span. Casting the `u32` length and `u8` overlaps before subtraction prevents unsigned underflow, and their ranges fit safely in `i64`. A `usize→i64` wrap requires an externally supplied enormous `max_insert`; production estimator/default values cannot reach it. No production-reachable overflow defect in this expression was established.

E. Test adequacy

**Inadequate for the stated risks.** The three tests cover exact lookup isolation, rejection of one oversized synthetic repeat, and sequence spelling for one accepted synthetic bridge.

They do not cover real read mapping, reciprocal orientation generation, fragment coordinates, alternate traversals, insert estimation, distribution boundaries, fallback false rejection, or negative gaps.

Verification: the existing cached test executable passed all three tests. Its build provenance is **UNVERIFIED**. Counterexample execution against a fresh build is **UNVERIFIED** because the environment is read-only; the examples above follow the inspected production control flow.

Final Question:
Does `gap <= max_insert` provide sufficient paired-end evidence for accepting the selected path?

**Answer: NO.**

It checks one topology-derived lower-bound component of fragment length. Acceptance requires authentic, correctly oriented fragments whose placements yield plausible spans **along the selected traversal**. Neither the scalar comparison nor aggregate directed-pair counts establishes that requirement.