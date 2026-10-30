# fastBact Development Handoff

## F24 final technical acceptance and integration handoff — 2026-10-30

- CURRENT_PROJECT: fastBact
- CURRENT_MILESTONE: F24 — Provenance-Safe Contiguity Recovery
- CURRENT_STATUS: F24_ACCEPTED_WITH_P1-4_WAIVER; integrated into main and verified clean on Slurm qfree.
- USER_DECISION: Formal technical acceptance approved with conditions (2026-10-30).
- INTEGRATED_MAIN_SHA: `cfbd00503b3031853bb41441160717d1e5bc3e69`
- IMMUTABLE_FREEZE_TAG: `f24-accepted-2026-10-30`
- POST_INTEGRATION_BINARY_SHA256: `53998c592cfe2155fca21a57fa1f39b36f5fc5fb1aeb11db970e760a4742b8fd` (Slurm Job `3072360`)
- POST_INTEGRATION_TESTS: Slurm Job `3072360` on `vnode1` (`qfree`): 213/213 release tests passed, CLI smoke test OK.
- P1-4 STATIC FLOOR WAIVER RATIONALE:
  - Original static floor requirement: 46,800 bp for R1-003 NGA50 (retained in contract record).
  - Observed R1-003 NGA50: 44,702 bp (+155.3% over frozen F23 baseline 17,509 bp).
  - Qualification: 44,702 bp represents the *empirically observed safe contiguity boundary under the implemented conservative evidence isolation criteria*, and is explicitly not described as a mathematically proven maximum.
  - Waiver rationale: Significant contiguity improvement (+155.3%), substantial base-level accuracy gains (mismatches 0.13 vs 0.59, indels 0.00), zero new misassemblies, and conservative evidence isolation taking precedence over speculative repeat joining.
- RESIDUAL F23-ERA MISASSEMBLIES PRESERVED & DOCUMENTED:
  - R1-001 (*E. coli* K-12): Strictly 5 extensive misassemblies (identical to F23 baseline; AUD-01 2.1 Mb chimera remains 100% eliminated; 5th misassembly on `contig_48` is a QUAST min-alignment/ambiguity threshold filtering artifact over a pre-existing baseline precursor).
  - R1-003 (*S. aureus* MW2): Strictly 0 extensive misassemblies, 1 residual local misassembly (identical to baseline).
- EMPIRICAL DETERMINISM: $\sigma^2 = 0.000000$ across 3 paired computational replicates (Slurm Job `3066737`).
- PRESERVED BENCHMARK ARTIFACTS: Slurm Job `3052717` (F23 baseline), Job `3065856` (initial F24 paired), Job `3066737` (repaired F24 3-replicate benchmark).
- FULL ACCEPTANCE REPORT: [`FASTBACT_F24_FINAL_ACCEPTANCE_WITH_WAIVER_REPORT.md`](../../../../.gemini/antigravity-cli/brain/7e7479aa-7a0c-4a62-84bd-459360983831/FASTBACT_F24_FINAL_ACCEPTANCE_WITH_WAIVER_REPORT.md).

## F24 evidence closure handoff — 2026-10-29

- CURRENT_PROJECT: fastBact
- CURRENT_MILESTONE: F24 — Provenance-Safe Contiguity Recovery
- CURRENT_STATUS: BLOCKED_WITH_SPECIFIC_EVIDENCE; independent acceptance review by Codex GPT-6.1 Sol Medium returned REJECT / FINDINGS
- ACTIVE_WORKTREE: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f24-contiguity-worktree`
- CANDIDATE_SHA: `3d44889218de85f7674fe73609db283047762fc3`
- BASELINE_SHA: `68b19f06d22537cf199e8fe124ea6b62bb0bb50a` (peeled frozen tag `f23-accepted-2026-10-29`)
- CANDIDATE_BINARY_SHA256: `8348395b9de19f7a3d21bf1179737d3b10427ee782c6fa6461f9c0b3608f1820`
- BASELINE_BINARY_SHA256: `97adb3c2ae386d453ba3f50fa61920836fd6c37d83924c0e0bc04d6763caaeb7`
- Contract: [F24_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md](../../F24_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md)
- Evidence Closure Report: [F24_FINAL_ACCEPTANCE_EVIDENCE_CLOSURE.md](./F24_FINAL_ACCEPTANCE_EVIDENCE_CLOSURE.md)
- COMPLETED_GATES:
  - Invariants F24-I01 through F24-I12 sequentially implemented in `src/resolve/exspander.rs`
  - Slurm Job `3065181`: 106/106 tests passed across 9 frozen test suites
  - Slurm Job `3065341`: 212/212 full workspace release tests passed with CLI smoke test OK
  - Slurm Job `3065210`: Real-data A/B QUAST benchmark (R1-001 NGA50 +78.0%, R1-003 NGA50 +164.1%, 0 new misassemblies, AUD-01 chimera absent)
- OPEN_FINDINGS:
  - Base-error increase under contract §8.3 (R1-001 mismatches 2.24 $\to$ 2.54; R1-003 mismatches 0.59 $\to$ 0.68)
  - Attribution lattice defaults to Single without propagating inherited Shared state (`exspander.rs:407`)
  - Evidence filtering predicate checks $p > b$ instead of $p > b + H$ (`exspander.rs:865`)
  - R1-003 NGA50 (46,236 bp) misses static absolute floor (46,800 bp) by 564 bp
  - Unbounded `Vec::contains` in origin propagation
- REMAINING_BUDGET: Within allocated parameters; no further paid-provider calls required.
- NEXT_ACTION: Awaiting user / Principal Architect guidance on resolution options (Option A: algorithmic tightening vs Option B: contract ruling on variance). Preserve isolated worktree. Merge/push strictly blocked.

**Status conflict recorded:** root AGENTS.md still says `F23_ALGORITHM_CONTRACT_BLOCKED` and C0-02R NOT ACCEPTED. The dated acceptance record below and the current user's explicit accepted-baseline instruction establish F23 as the frozen comparator for this architectural task. Older text is preserved as history; this specification neither reverses F23 acceptance nor grants F24 acceptance. Scientific safety, isolated worktree, compute, merge/push/release restrictions remain in force.

**Architectural qualification:** the requested unconditional inference from `length >= max_insert` to a locus-specific exit conflicts with frozen shared-long-node protections. F24 retains the stricter safety requirement: physical isolation and independent attribution are separate certificates. Missing attribution fails closed, even if recovery targets are missed.

## General policy

Normative cross-project collaboration rules are Multi-AI Collaboration V3.2:

`~/.config/ai-collaboration/POLICY.md`

That file is the only general policy. This handoff does not restate it. Antigravity, Codex, OpenCode, and Cursor do not auto-load that path. The mandatory load instruction is in the repository [`AGENTS.md`](../../AGENTS.md). Read the policy before acting. If it cannot be read, stop and report `POLICY_UNREADABLE`.

On 2026-10-28 the embedded V3.1 role policy in this file was superseded. Project evidence, phase history, the verified CLI matrix, and the fastBact safety constraints below are retained. They are not a second collaboration policy.

If a fastBact constraint is stricter than the general policy, keep the stricter safety or correctness requirement and record the conflict in this file.

**Recorded relationship (2026-10-28):** Opus architecture calls that meet the standing authorization in the general policy do not need a new per-call user prompt. Other Claude calls, budget breaches, and unverifiable balances still follow the Claude budget section below. This does not authorize scientific-contract changes, merge, push, or release.

**F23 final acceptance and completion status (2026-10-29; authoritative over older status text):**
Status: **`F23_COMPLETED`** and **`C0-02R_ACCEPTED`** (Formally authorized by user on 2026-10-29 under `USER_FINAL_MERGE_AUTHORIZATION_GRANTED`).

1. **User Authorization & Formal Scientific Acceptance**:
   - On 2026-10-29, the user as `FINAL_AUTHORIZATION_AUTHORITY` formally granted final merge and acceptance authorization (`USER_FINAL_MERGE_AUTHORIZATION_GRANTED`).
   - Scientific contract `C0-02R` formally transitions to **`ACCEPTED`**.
   - Phase F23 formally transitions to **`F23_COMPLETED`**.
   - Integrated into `main` via integration branch `integration/f23-final-acceptance`.
2. **Definitive Provenance & Evidence Matrix**:
   - Candidate Final Source: `src/resolve/exspander.rs` SHA-256 `fd59e8bdc231dbb83468c711f3d87834c51638d75941f587eb0aa1c89298790c`.
   - Candidate Release Binary: `target/release/fastBact` SHA-256 `2f1b6229c241479e9e3480adc6d60ab009d96f7ca8f7ff78307412a74b1c9a11` (Slurm Job `3052676`).
   - Deterministic Gate Suite: 118/118 PASS, Exit `0:0`, `ALL_GATES_PASSED` (Slurm Job `3052557` on `vnode1`).
   - Full Workspace Release Tests: 194/194 PASS, Exit `0:0` (Slurm Job `3049336` on `vnode1`).
   - Real-Data A/B Benchmark (QUAST v5.2.0, Slurm Job `3052717` on `vnode1`):
     - R1-003 (*E. coli* BW25113): extensive misassemblies reduced from 2 to 0; indels reduced from 0.04 to 0.00/100kb; mismatches reduced by 61%.
     - R1-001 (*E. coli* K-12 MG1655): historical AUD-01 2.1 Mb `contig_16` chimera completely eliminated. The 5th misassembly (`contig_48`) confirmed as a QUAST ambiguity filter threshold artifact over a pre-existing baseline precursor.
3. **Independent Review Closure (Codex GPT-6.1 Sol Medium)**:
   - Sol R6 targeted review issued an explicit **`PASS`** on Anchor Absorption safety under ARG/PAG, confirming that evidence stored under $B$ or $X$ has zero deciding weight for walks originating from $A$.
4. **Execution Safety**:
   - Integrated locally into `main`.
   - Zero remote push, zero GitHub release, zero git tags. All execution strictly within local repository boundaries.

**Opus technical ruling (2026-10-27): `APPROVE_WITH_CONDITIONS`.** Selected design direction: track attributable/shared context from graph topology; do not reset the anchor inside shared context; exclude evidence keyed at shared-context nodes; withhold when no attributable in-horizon support exists; keep node length only as a secondary condition. Before production implementation, statically verify where fork evidence is keyed and whether unchanged T26/T27/T30/T34 contracts conflict. Then implement and independently review F1–F9, including shared-long intermediate/seed, short-predecessor bypass, adversarial weights, legal continuity, distance boundaries, reverse traversal, merge remap and re-entry; require F1–F4 to fail on current source. This is conditional design approval only: no code changed, no tests were run, and release/C0-02R acceptance remains blocked. Full ruling: /hpcfs/fhome/caizhh/.local/state/fastbact-f23-post-audit-tests/validation/F23_OPUS_FINAL_TECHNICAL_DECISION.md (SHA-256 `37255753b3a15112fc3c322a78c6ee256a9c62438bd0ccc6d36b5b6aa0985f13`). CLI reports `claude-opus-5-5` and estimates USD 0.36093875; upstream identity, effort, actual relay deduction and current balance remain unverified.

**Capability evidence (2026-10-27):** Codex local model cache lists Luna,
Terra and Sol; actual account availability may change. OpenCode lists both GLM
IDs and GLM-5.3 previously completed independent review with effort unverified.
Cursor CLI lists High Fast, XHigh Fast and Auto; Auto was used on F23 tests.
`agy models` lists `gemini-3.8-flash-low`, `-medium`, and `-high`. Medium
interactive sandbox displayed the requested Gemini and answered one no-tool
request; Low previously ran `pwd` after one-time approval. Neither proves
backend model/pool identity or headless Slurm capability. Headless command
permission was denied; do not bypass permissions or enable broad allow rules.
Claude Code 2.1.216 routes ordinary tasks through the OmniaKey launcher and
`api.omniakey.com`; `claude auth status` bypasses that launcher and reports
raw first-party status, so it does not establish whether relay tasks work. On
2026-10-27, no-tool requests succeeded through the launcher with explicit
`claude-haiku-5-5`, `claude-sonnet-5-5`, and `claude-opus-5-5`; CLI
`modelUsage` keys matched those requests. This is relay-reported identity, not
independent proof of the upstream backend. Sonnet and Haiku also completed
read-only temporary-file checks; Haiku stream output recorded a `Read` tool
call. Shell, editing, subagents, and write-worker permissions remain unverified.
Future paid tasks require task-specific authorization; verify actual provider,
model and effort and use least-privilege tools. Do not assume CLI cost estimates
match relay billing.

### Claude Code OmniaKey soft budget

The user reports USD 50 recharged to the OmniaKey relay. The shared local
reminder is `~/.local/state/ai-routing/claude-budget.json`. Its reported
recharge is not a verified remaining balance: until a relay statement is
recorded, `last_verified_balance_usd` and reliable cumulative spend are null.
Use relay deductions as authoritative; CLI `total_cost_usd` is an estimate.
This file is a manual reminder, not a billing limit, monitor or hard block.

- Before each Claude task, check the latest recorded balance and the current
  task's budget authorization. Unknown balance must be reported as unknown;
  do not claim the budget is sufficient. Paid Claude calls still need task authorization, except an Opus architecture call that already satisfies the standing authorization in `~/.config/ai-collaboration/POLICY.md`.
- After a task, update reliable cumulative spend and balance only from relay
  billing or another reliable source. Record the source and timestamp. Keep
  estimated CLI costs separate from verified deductions.
- Retain a suggested USD 10 reserve. Below USD 20, prefer GLM, Gemini, Cursor
  or Codex subscriptions. Below USD 10, request renewed explicit budget
  authorization for a Claude call. Below USD 5, stop routine Claude work
  unless the user explicitly approves it.
- Sonnet 5.5 is the main Claude expert for key algorithm design and independent
  review. Opus 5.5 handles major disputes and high-risk correctness under
  the standing authorization in the general policy. Sonnet and Haiku still
  require separate task authorization. Use Haiku 5.5 only for a proven speed or cost
  advantage. Routine development and tests prefer subscription resources.
- Do not recharge automatically, switch to another metered provider, or call
  Claude merely to spend the balance. The soft reminder does not enforce
  charges at the relay.

**F23 current evidence:** independent OpenCode audit `CONDITIONAL_PASS`; new
targeted tests, full Release suite, and R1-001/R1-003 A/B have not run.
C0-02R remains **NOT ACCEPTED**. A static prediction is not a test result.


**Project**: `fastBact`
**Primary Worktree Branch**: `fix/c0-02r-final-gate`
**Base Commit**: `b9727c3ab372f7b995e439d92b72c7a4bfd051cf` (C0-02 baseline)
**Main/Upstream**: `dac6a57` (`origin/main`, `origin/master`)
**Default Controller**: Antigravity CLI (Gemini 3.8 Flash, Medium Reasoning, Default Speed, Auth: Google AI Pro). Sole default execution controller.
**Development and Review Resources**: Role duties and provider order are defined only in `~/.config/ai-collaboration/POLICY.md`. Local verified IDs are in the CLI matrix below.
**Model Routing**: The general routing policy is `~/.config/ai-collaboration/POLICY.md`. This handoff records the fastBact binding, verified CLI matrix, and project status. Root [`AGENTS.md`](../../AGENTS.md) is the project entry and carries the mandatory load instruction. No separate `MODEL_ROUTING.md` exists.
**Active Task Board**: [`ACTIVE_TASK.md`](./ACTIVE_TASK.md) (mirrored at [`repo root`](../../ACTIVE_TASK.md))
**Specific Phase Ledger**: [`C0_02R_HANDOFF.md`](./C0_02R_HANDOFF.md)
**Audit Findings Registry**: [`AUDIT_FINDINGS.md`](./AUDIT_FINDINGS.md)

### Current F23 Handoff

- **Status**: independent F23 audit `CONDITIONAL_PASS`; new targeted tests and further acceptance gates remain pending. C0-02R remains **NOT ACCEPTED**.
- **F23 worktree**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f23-unique-anchor-worktree`, branch `feature/c0-02r-phase-f23-unique-anchor`, HEAD `7bdc6f962f93d54cda192d47fb410e046fd5a627`. The main worktree remains separate at `fix/c0-02r-final-gate`; do not infer merge or cherry-pick.
- **Implementation**: Horizon guard commit `08efc6b1aae6c7499bd79618e12f1269f556962a`; F23 evidence/diagnosis commit `7bdc6f962f93d54cda192d47fb410e046fd5a627`.
- **Regression history**: Job `3010914` reported 38/42 with T26, T27, T30, and T34 failing. Job `3011496` later reported 42/42 but remains **TESTED_SNAPSHOT_UNVERIFIED** because it recorded no source SHA.
- **Current-HEAD verification**: Job `3020949` ran on `vnode1` against full SHA `7bdc6f962f93d54cda192d47fb410e046fd5a627` (tree `33f8a148c94cd73b6cfbe15b4de6bf2628f88d4f`). Commands: `cargo test --release --test test_c0_02r_state_semantics_gate` (42 passed, 0 failed, 0 ignored; exit 0) and `cargo test --release --test test_c0_02r_unique_anchor_gate` (2 passed, 0 failed, 0 ignored; exit 0). Slurm exit `0:0`. Unpiped log: `/hpcfs/fhome/caizhh/.local/state/codex-with-chatgpt/f23-verify-7bdc6f962f93d54cda192d47fb410e046fd5a627-3020949.log`; SHA-256 `2c4f7ea7145be3d8f0f4842d79251f13be43779a2f028359be8edb71615a770c`.
- **Test infrastructure**: Jobs `3010866` and `3010889` are Slurm `COMPLETED` with exit 0 despite inner integration-test failures. Their scripts pipe Cargo output through `tail`/`grep` without `pipefail`, so pipeline exit status can mask test failure. Job `3020949` is the current-HEAD targeted verification and records Cargo and Slurm exit statuses without output pipelines.
- **Independent-review routing**: Prepare F23 review through OpenCode with explicit model `zhipuai-coding-plan/glm-5.3` in a fresh read-only session. The installed OpenCode model directory confirms this model ID. Reasoning variants, including `#max`, remain unverified; record effort as `EFFORT_UNVERIFIED`. Do not treat prior ZCode configuration or simulated adapter output as current CLI evidence. Do not start the review until the read-only execution profile is confirmed.
- **Next gate**: targeted test execution, independent evidence review, then Release and selective A/B only if their gates pass. Existing baseline tests are green; new targeted tests remain unrun and C0-02R acceptance remains pending.

---

## 1. Operating Rules & Engineering Principles

1. **Roles and autonomous execution**:
   - Defined only by `~/.config/ai-collaboration/POLICY.md` (Multi-AI Collaboration V3.2).
   - This file does not assign a second controller or developer priority.
   - **DeepSeek API**: explicit user authorization is required; never use it as an automatic fallback. This project constraint is additional to the general provider list.
2. **Shared Correctness and Safety Constraints**:
   - All agents follow the same deterministic correctness, evidence, permissions, and artifact-preservation rules in this handoff.
   - Write tasks require an isolated Git worktree. Assign one writer per worktree at a time; other agents sharing it remain read-only unless explicitly assigned as its sole writer.
   - Before an external CLI call, verify its actual parameters, model ID, authentication mode, and non-interactive behavior. Record the actual model and effort used. Never silently switch to a paid API or a higher-cost model.
   - If Antigravity Controller is unavailable, preserve HANDOFF and pause automatic controller work without changing default controller identity.
3. **Slurm Compute Exclusivity**:
   - All Rust compilation (`cargo build`, `cargo test`, `cargo check`), assembly pipelines, and QUAST evaluations MUST run on Slurm compute nodes (`qfree` partition, typically `vnode1`).
   - Login nodes must never execute compilation or assembly workloads.
4. **Git Protection & Historical Data Integrity**:
   - Never auto-merge `main` or force-push (`git push --force`).
   - Never execute destructive operations: `git reset --hard`, `git clean -fd`, `rm -rf`.
   - Never overwrite or delete historical benchmark results or frozen baseline artifacts.
5. **Deterministic Verification**:
   - All core correctness issues and bug fixes must have deterministic counter-examples proving red before green, followed by regression testing.

---

## 2. Execution Permissions and CLI Safety

### Permission Levels

| Permission Level | Description | Execution Constraints |
| :--- | :--- | :--- |
| **`READ_ONLY`** (Default) | Code analysis, algorithm inspection, and existing result review | Shares project source view. No workspace modifications permitted. |
| **`WORKSPACE_WRITE`** | Concrete authorized code edits | Allowed **only** inside an isolated Git worktree. Requires explicit user approval. Never modify active working branch. |
| **`HUMAN_APPROVAL`** | High-risk, destructive, or high-cost actions | Requires explicit user sign-off before invocation. |

### Mandatory Approval Triggers
The following require explicit user authorization before execution:
- Any Claude call outside the Opus 5.5 standing authorization in `~/.config/ai-collaboration/POLICY.md`, and any Opus call that would exceed that authorization.
- Any transition of an external CLI to **`WORKSPACE_WRITE`**.
- Any authentication mode that could incur additional API charges (e.g., fallback to pay-per-token API keys).
- Submitting large-scale or high-duration Slurm tasks.
- Installing external dependencies over the network or broadening filesystem access.

### Strictly Prohibited Flags
- `codex --dangerously-bypass-approvals-and-sandbox` / `--yolo`
- `claude --dangerously-skip-permissions` / `--allow-dangerously-skip-permissions`
- `agent --force` / `--yolo`
- `agy --dangerously-skip-permissions`

### CLI Capability, Verified Model & Flag Routing Matrix

Role duties and autonomous escalation live in `~/.config/ai-collaboration/POLICY.md`. The matrix below is the verified fastBact invocation record. It is not a second policy.

| Tool | Verified Binary Path | Current Auth Mode | Verified binding | Verified Model ID | Reasoning Effort | Default Permissions / CLI Invocation Pattern |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Antigravity** | `~/.local/bin/agy` | Google AI Pro | **Controller** (Task planning, sequencing, gate checks, handoffs) | `gemini-3.8-flash-medium` | Medium | `agy --model gemini-3.8-flash-medium --sandbox`; approve commands individually; never use `--dangerously-skip-permissions` |
| **OpenCode** | `~/.opencode/bin/opencode` | Existing Zhipu AI Coding Plan login | **Developer (Priority 1)**: Routine development | `zhipuai-coding-plan/glm-5.3-flash` | CLI default; variant unverified | `opencode run --model zhipuai-coding-plan/glm-5.3-flash --format json <prompt>`; never add `--auto` by default |
| **OpenCode** | `~/.opencode/bin/opencode` | Existing Zhipu AI Coding Plan login | **Developer (Priority 1)**: Complex development | `zhipuai-coding-plan/glm-5.3` | `EFFORT_UNVERIFIED`; `#max` not verified | `opencode run --model zhipuai-coding-plan/glm-5.3 --format json <prompt>`; prioritize consuming GLM quota |
| **Cursor Pro** | `~/.local/bin/cursor-agent` | Cursor Pro subscription | **Developer (Priority 2)**: Fallback implementation | `grok-4.7-high-fast` (default) / `grok-4.7-xhigh-fast` (complex) | High Fast / XHigh Fast | `cursor-agent -p --model <model> <prompt>`; strictly when GLM is exhausted; if CLI is constrained, use human handoff |
| **Codex** | `~/.local/bin/codex` | ChatGPT OAuth | **Reviewer**: Code, test, Oracle & regression review | `gpt-6.1-terra` | Medium | `codex exec -m gpt-6.1-terra -c model_reasoning_effort="medium" <prompt>`; strictly independent of implementer |
| **Codex** | `~/.local/bin/codex` | ChatGPT OAuth | **Reviewer**: Core algorithm, state machine & correctness | `gpt-6.1-sol` | Medium (High upon authorization) | `codex exec -m gpt-6.1-sol -c model_reasoning_effort="medium" <prompt>`; independent review |
| **Codex** | `~/.local/bin/codex` | ChatGPT OAuth | **Validator**: Test validation, Slurm/Cargo exit code & gate checks | `gpt-6-luna` | Low (higher if justified) | `codex exec -m gpt-6-luna -c model_reasoning_effort="low" <prompt>`; checks evidence, never replaces test evidence |
| **Codex** | `~/.local/bin/codex` | ChatGPT OAuth | **Architect (On Demand)**: High-risk design & scientific contracts | `gpt-6-astra` | Low / Medium / High | Verified in cache (`gpt-6-astra`); single invocation per decision; cannot override deterministic tests |
| **Claude Code** | `~/.local/share/omniakey/bin/claude` | User's third-party relay API | **Architect (On Demand)**: High-risk technical rulings | `claude-opus-5-5` | Verified per call | Requires explicit task budget authorization; check `claude-budget.json` reminder; reserve USD 10 |

The OpenCode model directory lists `zhipuai-coding-plan/glm-5.3`, `glm-5.3-flash`, and `glm-5.3-highspeed`. The CLI accepts a `#variant` syntax, but no reasoning variant has been verified for these models. Do not claim GLM-5.3 Max.

Antigravity 1.4.8 lists `gemini-3.8-flash-low`, `gemini-3.8-flash-medium`, and `gemini-3.8-flash-high`. Gemini 3.8 Flash Medium serves as the primary controller. Do not use `--dangerously-skip-permissions`.

Cursor CLI is available via `~/.local/bin/cursor-agent` (`grok-4.7-high-fast`, `grok-4.7-xhigh-fast`). If non-interactive CLI worker execution is restricted or requires interactive IDE handling, record as human handoff mode.

Codex models `gpt-6.1-terra`, `gpt-6.1-sol`, `gpt-6-luna`, and `gpt-6-astra` are verified in local model cache (`~/.codex/models_cache.json`).

### Execution Safety & Subprocess Rules
- When calling external CLIs via Python, always use argument arrays (`argv: list[str]`) with `shell=False`.
- Never perform unvetted string concatenation into shell commands.
- External model reviews are informational evidence; they do not automatically constitute verification PASS.
- Every invocation must be logged to [`AUDIT_FINDINGS.md`](./AUDIT_FINDINGS.md) with Task ID, Agent, Model, Effort, Permissions, Git SHA, sanitized CLI command, exit code, output path, findings, and next actions.

### Quota & Rate Limit Handling
- Before choosing a GLM/Gemini worker for a task, read the shared summary with `python3 ~/.local/share/ai-routing/quota_manager.py status`; use `route routine`, `route slurm`, or another supported task type as a recommendation. Do not read manager source/tests on every session start.
- State path: `~/.local/state/ai-routing/quota-state.json`; tool: `~/.local/share/ai-routing/quota_manager.py`. GLM and Gemini pools are independent. `UNKNOWN` means availability is unconfirmed, not exhausted or verified available.
- Record only an explicit quota exhaustion observation. Use `record-limit <pool_id> --observed-at <ISO-8601-UTC>` and add `--retry-after <ISO-8601-UTC>` only when the platform supplied a trustworthy time. Without one, the manager schedules the first retry about five hours after first limit and later retries about 45 minutes after another explicit limit. Expiry only makes a normal task retry eligible; it never changes status automatically.
- After a normal successful request, use `record-success <pool_id> --observed-at <ISO-8601-UTC> --identity-verified --pool-confirmed` only when provider, model, pool, and success are actually confirmed. Ordinary network, auth, permission, or model errors do not clear a limit. Preserve task and worktree state when handing work to Cursor; do not preempt a running worker or allow multiple writers in one worktree.
- This manager is a local state recorder, not a real CLI parser or monitor. Stage 3A adapter envelopes are simulation-only. Do not infer limits from arbitrary log text, run automatic probes, or assume a five-hour timer proves recovery.

---

## 3. Phase Execution Status Summary (FB-P0-02R)

| Phase | Scope | Baseline / Target | Status | Key Milestone / Commit |
| :--- | :--- | :--- | :---: | :--- |
| **Phase A** | Directed Edge Verification Contract | `b9727c3` | **PASS** | Commit `f500d0f`: Replaced weak path validation with 9-point contract; eliminated unverified `GraphEdge::new()`. Tests E1–E8. |
| **Phase B** | Graph Mutation Safety & Tri-Contract Boundary | `f500d0f` | **PASS** | Commit `0d13f99`: Formalized 3 contracts (snapshot rollback, protected rejected short nodes, isolated node cleanup). Tests G1–G6. |
| **Phase C** | Modular Alignment & Provenance Engine | `0d13f99` | **PASS** | Commit `b318061`: Replaced position diffing with `minimap2 --cs=short` event extraction; machine-verifiable manifest architecture. |
| **Phase D** | 3-Way A/B/C Benchmark & Final Acceptance | `b318061` | **FAIL (NOT ACCEPTED)** | Slurm benchmark job `2980439` completed. Invariant FASTA on S1-COV-150X/003, but gate failures on R1-001 (mismatches 2.14 $\to$ 2.86, indels 0.48 $\to$ 0.96) and R1-003 (mismatches 1.51 $\to$ 1.67). Independent Codex audit completed (`AUDIT-20261019-001`, Verdict: FAIL; identified P1 shortcut vulnerability and P2 estimator k-mer offset bias). |
| **Phase E** | Evidence, Root Cause & Minimal Repair | `f241d9d` | **COMPLETED (C0-02R NOT ACCEPTED)** | E1: Proved R1-001 517 bp contig had 0 errors; root cause is consensus polish shifts in 10 other contigs (`VERIFIED`). E2: Repaired P2 insert-size estimator in isolated worktree (`feature/c0-02r-p2-estimator`, Commit `480ff5f`, `VERIFIED`). E3: Constructed deterministic counterexample proving Codex P1 shortcut vulnerability in production pipeline (Commit `b4472a1`, `VERIFIED`). E4: Formulated non-invasive lightweight fragment coordinate design proposal (`docs/development/P1_MINIMAL_DESIGN_PROPOSAL.md`). |
| **Phase F11-B** | Checked Single-Use Resolver Contract | `678c930` | **PASS** | Formalized RAII SingleUseResolverGuard, atomic edge commitment, and error result API. |
| **Phase F17-D/I** | Unified Search Admissibility & $\beta(\text{curr})$ Bound | `ac445f6` | **IMPLEMENTED** | Claude Opus 5.5 design Option A implemented by Grok 4.7 xhigh fast. Unified classify_hop with evidence-aware $\beta$ bound. |
| **Phase F18** | Independent Codex Audit of F17-I | `c37ef6e` | **CONDITIONAL PASS** | Audit `AUDIT-20261022-013`: 0 P0, 0 P1, 5 P2 test findings. Commit `39f11ca`. |
| **Phase F19** | Codex-Led P2 Verification Repair | `a993047` | **AUDITED (4/5 CLOSED, 1 P2 OPEN)** | Repaired oracle and sensitivity tests. Audit `AUDIT-20261022-014`. Commit `0606be4`. 183 tests pass on `vnode1`. |
| **Phase F20** | Codex-Led Selective Real-Data A/B Benchmark (R1-001, R1-003) | `0606be4` vs `b9727c3` | **COMPLETED (REGRESSION DETECTED)** | Slurm job `3001561` on `vnode1`. R1-003 improved (-50% misassemblies, 0 indels, +10% NGA50); R1-001 regressed (+1 misassembly, indels 0.48 $\to$ 1.03). Audit `AUDIT-20261022-015`. C0-02R strictly NOT ACCEPTED. |
| **Phase F21-A** | Targeted Forensic Investigation of Real-Data Regression | `4b956a4` | **FORENSIC COMPLETED** | Proved R1-001 extra misassembly is real chimeric join (`contig_16`, 7 short nodes). Proved base errors are 84.6% boundary polishing artifacts. Audit `AUDIT-20261023-016`. C0-02R strictly NOT ACCEPTED. |
| **Phase F21-B** | Repeat-Join Decision Provenance & Mechanism Verification | `07e0e5f` | **AUDITED (FAIL / QUALIFIED)** | Static evidence reconstruction of the 7-node repeat chain decision. Audit `AUDIT-20261023-017`: derived evidence insufficient to certify mechanism (`u1476` locus specificity contradicted; pair-weight reconciliation failures). Commit `2ae4acc`. C0-02R strictly NOT ACCEPTED. |
| **Phase F22-D** | Runtime Repeat-Join Decision Trace (`resolver-trace` diagnostics) | `2ae4acc` | **DIAGNOSTIC COMPLETE** | Feature-gated runtime decision/pair tracing; bit-identical non-interference verified vs F20 Candidate B on R1-001/R1-003. Root cause: anchor loss onto shared repeat node `u2102` (`unitig_1476`) let ectopic 3.66 Mb PE evidence hijack exit choice at `DEC_000022` (16:6, ratio 2.819 > 1.2). Diagnostics commit `9c543a0` (not merged to production). Counterfactual matrix recommends F23 = Option A Unique Anchor Invariant (B/D/E rejected). Independent audit deferred to post-F23 benchmark. C0-02R strictly NOT ACCEPTED. |
| **Phase F23** | Provenance-Aware Repeat Resolution | `2ae4acc` vs `68b19f0` | **ACCEPTED (COMPLETED)** | Option A+B unique anchor retention; Sol Round 6 PASS; Slurm `3052717` eliminated AUD-01 chimera; user authorized merge to main on 2026-10-29. |
| **Phase F24** | Provenance-Safe Contiguity Recovery | `68b19f0` (baseline) | **READY_FOR_FINAL_ACCEPTANCE** | 12 architectural invariants implemented; Sol Round 11 ACCEPT; Slurm `3065181` 106/106 PASS; Slurm `3065210` R1-001 NGA50 +78.0%, R1-003 NGA50 +164.1%, 0 new misassemblies. |

---

## 4. HPC Cluster Environment & Compute Resources

- **Scheduler**: Slurm 22.05+
- **Default Work Allocation**: Partition `qfree`, Node `vnode1`, 8 CPUs (`-c 8`), 12 GB RAM (`--mem=12G`).
- **Python**: `/hpcfs/fpublic/app/miniforge3/conda/bin/python3` (Python 3.12+ with pytest, uv).
- **Rust Toolchain**: `rustc 1.83.0+` (`~/.cargo/bin`).
- **QUAST**: `~/.local/bin/quast.py`.
- **Minimap2**: System / local binary in `PATH`.
