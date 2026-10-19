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
- **Next Action**:
  - Report findings to user. Await completion of ongoing Slurm benchmark job 2980439.
  - Do NOT modify Rust algorithms until authorized. Formulate minimal Antigravity patch retaining read offsets / concordant fragment span verification, or evaluate escalation to Claude for bioinfo algorithm review.

---

## Prepared Audits (Awaiting User Approval)

*(None currently pending approval)*

