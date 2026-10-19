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

*(No external model invocations have been launched yet. All preliminary development, deterministic refutations, and Slurm benchmark submissions have been executed directly by Antigravity CLI).*

---

## Prepared Audits (Awaiting User Approval)

### [PREPARED-20261019-001] Codex Medium Read-Only Review for PE Fidelity & Insert-Size Bounds

- **Task ID**: `C0-02R-PE-FIX-AUDIT`
- **Agent**: Codex CLI
- **Model**: `gpt-6.1-sol`
- **Reasoning Effort**: `medium` (`-c model_reasoning_effort="medium"`)
- **Permission Level**: `READ_ONLY` (`-s read-only`)
- **Git SHA**: `f241d9d3c9d21470f63bd85b5e1bd6eb0bfc8d30`
- **CLI Command**:
  ```bash
  codex exec -s read-only -m gpt-6.1-sol -c model_reasoning_effort="medium" -o docs/development/evidence/C0_02R/phase_d/codex_review_pe_fix.md "Review commit f241d9d in src/resolve/mapper.rs and src/resolve/exspander.rs regarding directed paired-end link support and insert size constraint gap <= dist.max_insert. Verify Rust correctness, integer overflow safety, and edge invariants."
  ```
- **Exit Code**: PENDING
- **Output Path**: `docs/development/evidence/C0_02R/phase_d/codex_review_pe_fix.md`
- **Findings**: PENDING_USER_APPROVAL
- **Next Action**: Execute upon user approval of first external CLI call.

