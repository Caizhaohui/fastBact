# fastBact Development Handoff

**Project**: `fastBact`  
**Current Active Branch**: `fix/c0-02r-final-gate`  
**Base Commit**: `b9727c3ab372f7b995e439d92b72c7a4bfd051cf` (C0-02 baseline)  
**Main/Upstream**: `dac6a57` (`origin/main`, `origin/master`)  
**Primary Coordinator & Implementer**: Antigravity CLI  
**Active Task Board**: [`ACTIVE_TASK.md`](./ACTIVE_TASK.md) (mirrored at [`repo root`](../../ACTIVE_TASK.md))  
**Specific Phase Ledger**: [`C0_02R_HANDOFF.md`](./C0_02R_HANDOFF.md)  
**Audit Findings Registry**: [`AUDIT_FINDINGS.md`](./AUDIT_FINDINGS.md)  

---

## 1. Operating Rules & Engineering Principles

1. **Self-Contained Primary Development**:
   - Everyday code implementation, bug fixing, test development, benchmark runs, documentation, and maintenance are handled directly by Antigravity CLI to conserve external quota.
   - External CLI dispatch is never invoked for ceremonial multi-agent interaction.
2. **Strict External CLI Dispatch & Escalation**:
   - **Codex CLI (`~/.local/bin/codex`)**: Reserved for complex Rust language correctness and memory/graph mutation invariant audits.
   - **Claude Code (`~/.local/share/omniakey/bin/claude`)**: Reserved for bioinformatics algorithms and biological/scientific correctness reviews.
   - **Cursor Agent (`~/.local/bin/agent`)**: Optional engineering assistance (available economic models only).
   - **Escalation Path**:
     $$\text{Antigravity Implementation} \longrightarrow \text{Codex Medium Review} \longrightarrow \text{Antigravity Fix} \longrightarrow \text{Claude Sonnet High Review (if needed)} \longrightarrow \text{Claude Opus High (critical unresolved disputes only)}$$
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

## 2. External Model & Permission Routing Policy

### Permission Levels

| Permission Level | Description | Execution Constraints |
| :--- | :--- | :--- |
| **`READ_ONLY`** (Default) | Code analysis, algorithm inspection, and existing result review | Shares project source view. No workspace modifications permitted. |
| **`WORKSPACE_WRITE`** | Concrete authorized code edits | Allowed **only** inside an isolated Git worktree. Requires explicit user approval. Never modify active working branch. |
| **`HUMAN_APPROVAL`** | High-risk, destructive, or high-cost actions | Requires explicit user sign-off before invocation. |

### Mandatory Approval Triggers
The following require explicit user authorization before execution:
- Every invocation of **Claude Opus**.
- **First-time invocation** of any external CLI tool.
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

| Tool | Verified Binary Path | Current Auth Mode | Task Type | Verified Model ID | Reasoning Effort | Default Permissions / CLI Invocation Pattern |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Codex** | `~/.local/bin/codex` | ChatGPT OAuth Tokens (`~/.codex/auth.json`) | Standard Rust Review | `gpt-6.1-sol` | Medium | `codex exec -s read-only -m gpt-6.1-sol -c model_reasoning_effort="medium" -o <out>` |
| **Codex** | `~/.local/bin/codex` | ChatGPT OAuth Tokens | Complex P0 Correctness | `gpt-6.1-sol` | High | `codex exec -s read-only -m gpt-6.1-sol -c model_reasoning_effort="high" -o <out>` |
| **Claude** | `~/.local/share/omniakey/bin/claude` | OmniaKey / Local Session | Standard Algorithm Review | `sonnet` | Medium | `claude -p --permission-mode plan --model sonnet "<prompt>"` |
| **Claude** | `~/.local/share/omniakey/bin/claude` | OmniaKey / Local Session | Complex Algorithm Review | `sonnet` | High | `claude -p --permission-mode plan --model sonnet "<prompt>"` (with high-effort prompt) |
| **Claude** | `~/.local/share/omniakey/bin/claude` | OmniaKey / Local Session | Critical Graph Disputes | `opus` | High | *Requires explicit user approval*. `claude -p --permission-mode plan --model opus "<prompt>"` |
| **Cursor** | `~/.local/bin/agent` | Local Cursor Agent CLI | Optional Engineering | `composer-2.5` / `cursor-grok-4.5-low` | Default | `agent -p --mode=ask --model composer-2.5 "<prompt>"` (*Note: `-p` is NOT read-only; `--mode=ask` enforces read-only*) |

### Execution Safety & Subprocess Rules
- When calling external CLIs via Python, always use argument arrays (`argv: list[str]`) with `shell=False`.
- Never perform unvetted string concatenation into shell commands.
- External model reviews are informational evidence; they do not automatically constitute verification PASS.
- Every invocation must be logged to [`AUDIT_FINDINGS.md`](./AUDIT_FINDINGS.md) with Task ID, Agent, Model, Effort, Permissions, Git SHA, sanitized CLI command, exit code, output path, findings, and next actions.

### Quota & Rate Limit Handling
- If Codex or Claude returns quota exceeded or rate limit:
  1. Immediately stop invocations of that CLI.
  2. Do not loop or retry infinitely.
  3. Do not switch automatically to pay-as-you-go API keys.
  4. Do not auto-upgrade to higher tier models.
  5. Record quota block in this document and continue unblocked Antigravity tasks.

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

---

## 4. HPC Cluster Environment & Compute Resources

- **Scheduler**: Slurm 22.05+
- **Default Work Allocation**: Partition `qfree`, Node `vnode1`, 8 CPUs (`-c 8`), 12 GB RAM (`--mem=12G`).
- **Python**: `/hpcfs/fpublic/app/miniforge3/conda/bin/python3` (Python 3.12+ with pytest, uv).
- **Rust Toolchain**: `rustc 1.83.0+` (`~/.cargo/bin`).
- **QUAST**: `~/.local/bin/quast.py`.
- **Minimap2**: System / local binary in `PATH`.
