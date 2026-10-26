# fastBact Agent Rules

This file is the F23 worktree entry point. Follow the shared correctness,
safety, and execution constraints in
[`docs/development/HANDOFF.md`](docs/development/HANDOFF.md).

- Default controller: Codex GPT-6 Luna Low. If Codex quota is exhausted,
  Antigravity Gemini 3.8 Flash takes over. If both are unavailable, preserve a
  HANDOFF and pause; Cursor is a development worker, never controller.
- OpenCode through the existing `zhipuai-coding-plan` subscription is the sole
  formal GLM route. Use `zhipuai-coding-plan/glm-5.3-flash` for routine work and
  `zhipuai-coding-plan/glm-5.3` for complex diagnosis and independent review.
  Reasoning variants, including `#max`, are unverified. Never silently fall
  back to another provider.
- Before selecting GLM/Gemini, read the shared status with
  `python3 ~/.local/share/ai-routing/quota_manager.py status`. UNKNOWN means
  availability is unconfirmed. This recorder does not parse real CLI output or
  perform automatic probes; see HANDOFF.md for event recording and retry rules.
- Antigravity handles tests, Slurm, benchmarks, bioinformatics, and auxiliary
  development. Cursor Pro Auto/Grok handles suitable worker tasks when needed,
  without preempting an active worker or taking controller duties.
- Claude Sonnet/Opus and Codex GPT-6.1 Sol retain their specialist roles in
  HANDOFF.md. DeepSeek API requires explicit user authorization.
- Read [`docs/development/ACTIVE_TASK.md`](docs/development/ACTIVE_TASK.md)
  when entering project work. C0-02R remains NOT ACCEPTED until separately
  supported by acceptance evidence.
- This F23 worktree is independent from the main worktree. Do not infer a
  merge or cherry-pick. Run Rust tests and assembly/QUAST workloads only on
  Slurm compute nodes in the `qfree` partition. Preserve untracked and
  historical artifacts.
