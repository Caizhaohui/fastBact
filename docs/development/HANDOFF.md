# fastBact Development Handoff

**Project**: `fastBact`
**Main Worktree Branch**: `fix/c0-02r-final-gate`
**F23 Audit Worktree**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f23-unique-anchor-worktree`, branch `feature/c0-02r-phase-f23-unique-anchor`, source code baseline `7bdc6f962f93d54cda192d47fb410e046fd5a627`.
**Base Commit**: `b9727c3ab372f7b995e439d92b72c7a4bfd051cf` (C0-02 baseline)
**Main/Upstream**: `dac6a57` (`origin/main`, `origin/master`)
**Primary Controller**: Codex GPT-6 Luna Low. Antigravity Gemini 3.8 Flash takes over controller duties only when Codex quota is exhausted. If both are unavailable, save HANDOFF and pause; Cursor remains a worker.
**Formal GLM Route**: OpenCode CLI through the existing Zhipu AI Coding Plan only. Routine work uses `zhipuai-coding-plan/glm-5.3-flash`; complex diagnosis and independent review use `zhipuai-coding-plan/glm-5.3`. The installed model list confirms these IDs. Reasoning variants, including `#max`, are unverified and must not be claimed. Do not fall back to another provider without authorization.
**Development Workers**: Antigravity Gemini 3.8 Flash for tests, Slurm, benchmark, bioinformatics, and auxiliary development; Cursor Pro Auto/Grok for suitable work when GLM/Gemini is limited or unsuitable. Cursor never controls or preempts a running task. Claude Sonnet/Opus and Codex GPT-6.1 Sol keep their authorized specialist roles. DeepSeek API requires explicit authorization.
**Active Task Board**: [`ACTIVE_TASK.md`](./ACTIVE_TASK.md) (mirrored at [`repo root`](../../ACTIVE_TASK.md))
**Specific Phase Ledger**: [`C0_02R_HANDOFF.md`](./C0_02R_HANDOFF.md)
**Audit Findings Registry**: [`AUDIT_FINDINGS.md`](./AUDIT_FINDINGS.md)

---

## 1. Operating Rules & Engineering Principles

1. **Controller and Worker Assignment**:
   - Codex GPT-6 Luna Low is the default controller. Antigravity Gemini 3.8 Flash takes over controller duties when Codex quota is exhausted. If both are unavailable, preserve HANDOFF and pause; Cursor is never controller.
   - OpenCode with the existing Zhipu AI Coding Plan is the sole formal GLM route. Use `zhipuai-coding-plan/glm-5.3-flash` for routine development and `zhipuai-coding-plan/glm-5.3` for complex diagnosis and independent review. The installed model list confirms these IDs. Reasoning variants, including `#max`, are unverified; do not claim Max or silently select another provider.
   - Antigravity Gemini 3.8 Flash handles tests, Slurm, benchmarks, bioinformatics, and auxiliary development.
   - Cursor Pro Auto/Grok handles suitable development work when GLM/Gemini is limited or unsuitable. Preserve an active task and worktree; never preempt a running Cursor worker or assign two writers to one worktree.
   - Claude Sonnet 5.5 High handles key algorithm design and independent review. Claude Opus 5.5 High is reserved for major architecture and difficult correctness disputes. Third-party paid calls require explicit authorization.
   - Codex GPT-6.1 Sol Medium handles key Rust algorithms and high-risk defects. High requires explicit authorization.
   - DeepSeek API requires explicit user authorization and is never automatic fallback.
2. **Shared Correctness and Safety Constraints**:
   - All agents follow deterministic correctness, evidence, permissions, and artifact-preservation rules.
   - Write tasks use an isolated Git worktree with one writer. Do not infer merge or cherry-pick between independent worktrees.
   - Before any external CLI call, verify the actual model ID, provider, auth route, permission mode, and invocation behavior. Do not silently switch to a paid provider or an unverified model variant.
   - Preserve task, Git, and test state when changing workers. If Codex and Antigravity are both unavailable, write HANDOFF and pause controller work.
3. **Slurm Compute Exclusivity**:
   - All Rust compilation/tests, assembly pipelines, and QUAST evaluations run on Slurm compute nodes in the `qfree` partition.
4. **Git Protection & Historical Data Integrity**:
   - Never force-push or automatically merge. Never run destructive reset/clean commands. Preserve untracked files and historical benchmark/evidence artifacts.
5. **Deterministic Verification**:
   - Core correctness changes require deterministic counterexamples and evidence before acceptance.

---

## 2. External Model & Permission Routing Policy

### Permission Levels

| Permission Level | Description | Execution Constraints |
| :--- | :--- | :--- |
| **`READ_ONLY`** | Code and evidence inspection | No workspace modifications |
| **`WORKSPACE_WRITE`** | Authorized code edits | Isolated worktree; one writer |
| **`HUMAN_APPROVAL`** | High-cost or sensitive action | Explicit user authorization before invocation |

### Approval and Permission Rules
- Third-party paid Claude calls require explicit user authorization; Opus requires explicit authorization.
- Do not use `--auto`, `--yolo`, `--force`, or dangerous permission bypass flags for review.
- An OpenCode F23 review must use a verified read-only agent/profile. That execution profile is not yet confirmed; keep the audit prepared but do not start it.
- Never inspect or print authentication material. Do not invoke an unapproved provider or fallback.

### Current Model Routing

| Tool | Provider / Auth | Task | Model ID | Reasoning | Permission / Invocation |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **OpenCode** | Existing Zhipu AI Coding Plan | Routine implementation, tests, scripts | `zhipuai-coding-plan/glm-5.3-flash` | CLI default; variant unverified | Explicit model and provider; no automatic provider fallback |
| **OpenCode** | Existing Zhipu AI Coding Plan | Complex diagnosis and F23 independent review | `zhipuai-coding-plan/glm-5.3` | `EFFORT_UNVERIFIED`; `#max` unverified | F23 review requires verified read-only agent/profile |
| **Antigravity** | Existing Gemini subscription | Tests, Slurm, benchmarks, bioinformatics; fallback controller | Gemini 3.8 Flash | Use an installed selectable effort and record it | Preserve normal permission prompts; no bypass |
| **Cursor Pro** | Subscription | Suitable development worker | Auto / Grok 4.7 High or XHigh, if available | Actual mode must be checked | Manual handoff until execution entry is verified; never controller |
| **Claude Sonnet** | Authorized third-party/local session | Key algorithm design and independent review | Sonnet 5.5 | High | Explicit authorization for paid calls |
| **Claude Opus** | Authorized third-party/local session | Major architecture/correctness disputes | Opus 5.5 | High | Explicit authorization |
| **Codex** | ChatGPT OAuth | High-risk Rust algorithms and defects | GPT-6.1 Sol | Medium; High by authorization | Read-only for review; use isolated worktree for writes |

OpenCode 2.1.4 currently lists `zhipuai-coding-plan/glm-5.3`, `glm-5.3-flash`, and `glm-5.3-highspeed`. The CLI accepts `provider/model#variant` syntax, but no GLM reasoning variant has been verified. Do not claim GLM-5.3 Max. ZCode use in older logs is historical and is not a current execution route.

### Simplified Quota State
- Before selecting a GLM/Gemini worker, read `python3 ~/.local/share/ai-routing/quota_manager.py status`; use its `route <task_type>` output as a recommendation. Do not read the manager source/tests at every session start.
- Shared state: `~/.local/state/ai-routing/quota-state.json`. GLM and Gemini pools are independent. `UNKNOWN` means availability is unconfirmed.
- Record only an explicit limit with `record-limit <pool_id> --observed-at <ISO-8601-UTC>`. Add `--retry-after <ISO-8601-UTC>` only if the provider supplied a trustworthy recovery time. Otherwise the first retry is eligible about five hours after the first limit; a further explicit limit schedules another retry about 45 minutes later. Expiry does not automatically restore availability.
- Record a successful normal task with `record-success <pool_id> --observed-at <ISO-8601-UTC> --identity-verified --pool-confirmed` only when provider, model, pool, and execution success are confirmed. Network, auth, permission, and model errors do not clear a limit.
- This is a local state recorder, not a real CLI output parser or monitor. Stage 3A adapter envelopes are simulation-only. Do not infer limits from arbitrary text or run automatic probes.

### Execution Safety & Subprocess Rules
- Python subprocesses use argument arrays and `shell=False`; never build unvetted shell strings.
- External reviews are informational and do not automatically establish verification PASS.
- Log authorized external review invocations in `AUDIT_FINDINGS.md` with task ID, actual model/effort, permissions, source SHA, sanitized command, exit code, output location, findings, and next actions.

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
| **Phase F21-A** | Targeted Forensic Investigation of Real-Data Regression | `4b956a4` | **FORENSIC COMPLETED** | Proved R1-001 extra misassembly is real chimeric join (`contig_16`, 7 short nodes). Proved base errors are 84.6% boundary polishing artifacts. Audit `AUDIT-20261023-016`. C0-02R strictly NOT ACCEPTED. Commit `07e0e5f`. |
| **Phase F21-B** | Repeat-Join Decision Provenance & Mechanism Verification | `07e0e5f` | **AUDITED (FAIL / QUALIFIED)** | Provenance ledger & graph geometry verified; Baseline A pre-existence confirmed in R1-003; coverage loss reconciled to 25.5 kb. Audit `AUDIT-20261023-017`. Status: `PHASE_F21B_DIAGNOSTIC_EVIDENCE_REQUIRED`. C0-02R strictly NOT ACCEPTED. |

---

## 4. HPC Cluster Environment & Compute Resources

- **Scheduler**: Slurm 22.05+
- **Default Work Allocation**: Partition `qfree`, Node `vnode1`, 8 CPUs (`-c 8`), 12 GB RAM (`--mem=12G`).
- **Python**: `/hpcfs/fpublic/app/miniforge3/conda/bin/python3` (Python 3.12+ with pytest, uv).
- **Rust Toolchain**: `rustc 1.83.0+` (`~/.cargo/bin`).
- **QUAST**: `~/.local/bin/quast.py`.
- **Minimap2**: System / local binary in `PATH`.
