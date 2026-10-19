# fastBact — External Model Audit Findings Ledger

This document maintains the persistent, machine-traceable record of all cross-model audits and external AI invocations according to the *fastBact External AI Model and Permission Routing Policy*.

---

## Invocation Schema Template

Every cross-model execution MUST record the following fields:

```markdown
### [AUDIT-YYYYMMDD-NNN] <Short Title>

- **Task ID**: <e.g., C0-02R Phase D>
- **Agent**: <Codex | Claude Code | Cursor Agent>
- **Model**: <Exact Model ID verified in CLI environment>
- **Reasoning Effort**: <low | medium | high | xhigh | none>
- **Permission Level**: <READ_ONLY | WORKSPACE_WRITE | HUMAN_APPROVAL>
- **Git SHA**: <40-char full commit SHA>
- **CLI Command**: `<exact command line array executed, redacted of tokens/keys>`
- **Exit Code**: <integer>
- **Output Path**: `<path to raw output or receipt>`
- **Findings**:
  - <Bulleted summary of architectural / correctness observations>
- **Next Action**:
  - <Specific implementation, rejection, or verification follow-up>
```

---

## Audit History

### [AUDIT-20261019-001] Codex Medium Read-Only Review for PE Fidelity & Insert-Size Bounds

- **Task ID**: `C0-02R-PE-FIX-AUDIT`
- **Agent**: Codex CLI (v0.163.0)
- **Model**: `gpt-6.1-sol`
- **Reasoning Effort**: `medium` (`-c model_reasoning_effort="medium"`)
- **Permission Level**: `READ_ONLY` (`-s read-only`)
- **Git SHA**: `f241d9d3c9d21470f63bd85b5e1bd6eb0bfc8d30`
- **CLI Command**:
  ```bash
  codex exec -s read-only -m gpt-6.1-sol -c model_reasoning_effort="medium" -o docs/development/evidence/C0_02R/phase_d/codex_review_pe_fix.md - < docs/development/evidence/C0_02R/phase_d/.prompt_codex_review.txt
  ```
- **Exit Code**: `0` (Completed; 63,169 tokens used)
- **Output Path**: [`docs/development/evidence/C0_02R/phase_d/codex_review_pe_fix.md`](./evidence/C0_02R/phase_d/codex_review_pe_fix.md)
- **Verdict**: **FAIL**
- **Findings**:
  - **Orientation Verification**: Confirmed that `PairLinkMatrix::get_support()` and `map_paired_reads()` in commit `f241d9d` correctly handle inward-facing read pairs without orientation mixing, normalization errors, or double counting.
  - **P1 (Unitig-pair counts can support a path incompatible with authentic fragments)**:
    - *Location*: `src/resolve/exspander.rs:223`, lines 241–244; `src/resolve/mapper.rs:144`, lines 158–159.
    - *Defect*: `PairLinkMatrix` only records counts of links between unitigs without preserving read offset/coordinates. When multiple alternate paths connect `A` and `B` (e.g. `A -> R -> C -> B` true 600 bp span vs shortcut `A -> R -> B`), read pairs originating from the longer traversal produce links between `(A, B)`. The direct exit `B` from `R` has `gap = 60 <= max_insert`, so the aggregate count selects the shortcut, even though the true fragments do not span 360 bp.
    - *Coverage*: Not covered by existing synthetic unit tests.
    - *Fix Direction*: Retain oriented read placements, offsets, and fragment lengths; evaluate candidate path support against implied fragment span.
  - **P2 (Empirical insert estimator uses k-mer anchor positions rather than read positions)**:
    - *Location*: `src/resolve/mapper.rs:116`, lines 123, 147–153.
    - *Defect*: Intra-unitig distance calculation `d = (pos2 - pos1).abs() + r1.len()` uses k-mer anchor positions instead of read endpoints, shifting apparent insert sizes by up to `read_len - k`.
    - *Coverage*: Not covered by existing unit tests.
    - *Fix Direction*: Recover oriented read intervals from anchor coordinates and sampled offsets; estimate distribution only from concordant inward-facing pairs.
  - **Final Question Answer**: `Does gap <= max_insert provide sufficient paired-end evidence for accepting the selected path?` **NO**. `gap` only measures net separation across the repeat node (`start(B) - end(A)`), not actual fragment span (`(|A| - s) + gap + e`). It is a necessary geometric lower-bound exclusion filter, but not sufficient evidence for path acceptance.
- **Phase E Empirical Verification (2026-10-19)**:
  - **P1 Status**: **`VERIFIED`**. Confirmed in production pipeline test `tests/test_c0_02r_p1_counterexample.rs` (commit `b4472a1`). exSPAnder erroneously selected shortcut $A \to R \to B$ (1038 bp) over true path $A \to R \to C \to B$ (1247 bp) because `PairLinkMatrix` only stored aggregate counts `(A, B) = 25`. Minimal design proposal drafted in [`docs/development/P1_MINIMAL_DESIGN_PROPOSAL.md`](./P1_MINIMAL_DESIGN_PROPOSAL.md).
  - **P2 Status**: **`VERIFIED` & `REPAIRED`** (in feature worktree). Confirmed in deterministic red test `tests/test_c0_02r_p2_estimator.rs` (+119 bp estimation error on unpatched code). Minimal repair implemented in `src/resolve/mapper.rs` (commit `480ff5f`), restoring true fragment span recovery across arbitrary seed offsets and inward strand orientations. All 106 tests green.
  - **R1-001 517 bp Contig Status**: **`VERIFIED` (False Hypothesis Eliminated)**. The extra 517 bp contig `contig_164` aligned 100.0% identity (0 mismatches, 0 indels) to `NC_000913.3`. Diagnostic QUAST (Slurm Job `2980935`) with `contig_164` removed yielded identical error rates (2.86 mismatches, 0.96 indels). The true root cause is consensus polishing base shifts across 10 other contigs (e.g. `contig_51`, `contig_105`) caused by anchor selection differences.
- **Next Action**:
  - Maintain C0-02R as **NOT ACCEPTED**.
  - Await user guidance on initiating Phase F (production implementation of P1 fragment span verification and consensus polishing stabilization).

---

## Prepared Audits (Awaiting User Approval)

*(None currently pending approval)*

