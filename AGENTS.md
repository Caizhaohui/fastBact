# fastBact Agent Rules

## General policy (mandatory load)

Multi-AI Collaboration V3.2 is the only general autonomous-development policy:

`~/.config/ai-collaboration/POLICY.md`

Antigravity (`agy`), Codex (`codex`), OpenCode (`opencode`), and Cursor (`cursor-agent` / `agent`) do not auto-load that path. Before any work in this repository, read that file in full and follow it. If it cannot be read, stop and report `POLICY_UNREADABLE`. Do not reconstruct the policy from memory, and do not treat older V3.1 text as the current general policy.

When dispatching a worker, put this sentence in the worker prompt:

`Read and follow ~/.config/ai-collaboration/POLICY.md before acting. It is not auto-loaded.`

Do not `@include` the policy into this file. Antigravity truncates an expanded rule file at 24 KB, and the other entry points do not expand that syntax.

This repository does not copy the general policy. Project rules below do not replace it. If a project rule and the general policy conflict, keep the stricter safety or correctness requirement and record the conflict in [`docs/development/HANDOFF.md`](docs/development/HANDOFF.md).

Embedded Multi-AI Collaboration V3.1 role text in this file is superseded as of 2026-10-28.

## fastBact project contract

- Maintain: `F23_ALGORITHM_CONTRACT_BLOCKED`.
- Core unresolved issue: Contract divergence between F9(b) loop re-entrancy waiver and T34 legal self-loop bridging.
- Next action: Architect first clarifies the technical contract. Once defined, Developer implements sequentially, Sol conducts independent review, and Luna runs deterministic verification. C0-02R remains **NOT ACCEPTED**.
- Scientific status, evidence, and phase history live in [`docs/development/HANDOFF.md`](docs/development/HANDOFF.md) and [`docs/development/ACTIVE_TASK.md`](docs/development/ACTIVE_TASK.md). This migration does not adjudicate those records.

## fastBact execution constraints

- Before choosing GLM or Gemini for a task, check `python3 ~/.local/share/ai-routing/quota_manager.py status`. `UNKNOWN` means availability is unconfirmed.
- Before a Claude task, inspect `~/.local/state/ai-routing/claude-budget.json` and the task authorization. Opus standing authorization is defined only in the general policy.
- Read HANDOFF and `docs/development/ACTIVE_TASK.md` when entering active project work.
- HANDOFF is the fastBact handoff: current binding, verified CLI matrix, and project status. `docs/development/ACTIVE_TASK.md` is the technical task and phase status. General roles and autonomous execution are defined only in the general policy.
- Use isolated Git worktrees for write tasks. Exactly one writer per worktree; other agents on that worktree remain read-only.
- Run Rust compilation, tests, assembly pipelines, and QUAST only on Slurm compute nodes in the `qfree` partition.
- Preserve historical benchmark and audit artifacts. Never force-push or automatically merge `main`.
- The user remains `FINAL_AUTHORIZATION_AUTHORITY` for scientific-contract changes, merge, push, and release.
