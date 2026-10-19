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
