# fastBact Development Handoff

**Project**: `fastBact`  
**Current Active Branch**: `fix/c0-02r-final-gate`  
**Base Commit**: `b9727c3ab372f7b995e439d92b72c7a4bfd051cf` (C0-02 baseline)  
**Main/Upstream**: `dac6a57` (`origin/main`, `origin/master`)  
**Primary Coordinator**: Antigravity CLI  
**Active Task**: [`ACTIVE_TASK.md`](../../ACTIVE_TASK.md)  
**Specific Phase Handoff**: [`C0_02R_HANDOFF.md`](./C0_02R_HANDOFF.md)  

---

## 1. Operating Rules & Engineering Principles

1. **Self-Contained Primary Development**: Routine development, debugging, and testing are carried out directly by Antigravity CLI to conserve external API quotas.
2. **Strict External CLI Dispatch**:
   - **Codex CLI (`~/.local/bin/codex`)**: Reserved for complex Rust language correctness and memory/graph mutation invariant audits.
   - **Claude Code (`~/.local/share/omniakey/bin/claude`)**: Reserved for bioinformatics algorithms and biological/scientific correctness reviews.
   - **Cursor Agent (`~/.local/bin/agent`)**: Optional engineering assistance.
   - **User Approval**: First-time invocation of any external CLI, any Opus model usage, or tasks with significant resource/quota cost must receive explicit user approval.
   - **Read-Only**: External CLIs are strictly read-only and never write directly to the working directory.
   - **Audit Trail**: Every cross-model invocation must record target Git SHA, exact command line, exit code, and structured audit outcome in `ACTIVE_TASK.md` and this document.
3. **Slurm Compute Exclusivity**:
   - All Rust compilation (`cargo build`, `cargo test`, `cargo check`), assembly runs, and QUAST evaluations MUST run on Slurm compute nodes (`qfree` partition, typically `vnode1`).
   - Login nodes must never run heavy compilation or assembly pipelines.
4. **Git Protection & Historical Data Integrity**:
   - Never auto-merge `main` or force-push (`push -f`).
   - Never overwrite or delete historical benchmark results or frozen baseline artifacts.
5. **Deterministic Verification**:
   - All core correctness issues and bug fixes must have deterministic counter-examples proving red before green, followed by regression testing.

---

## 2. Phase Execution Status Summary (FB-P0-02R)

| Phase | Scope | Baseline / Target | Status | Key Milestone / Commit |
| :--- | :--- | :--- | :---: | :--- |
| **Phase A** | Directed Edge Verification Contract | `b9727c3` | **PASS** | Commit `f500d0f`: Replaced weak path validation with 9-point contract; eliminated unverified `GraphEdge::new()`. Tests E1–E8. |
| **Phase B** | Graph Mutation Safety & Tri-Contract Boundary | `f500d0f` | **PASS** | Commit `0d13f99`: Formalized 3 contracts (snapshot rollback, protected rejected short nodes, isolated node cleanup). Tests G1–G6. |
| **Phase C** | Modular Alignment & Provenance Engine | `0d13f99` | **PASS** | Commit `b318061`: Replaced position diffing with `minimap2 --cs=short` event extraction; machine-verifiable manifest architecture. |
| **Phase D** | 3-Way A/B/C Benchmark & Final Acceptance | `b318061` | **IN PROGRESS** | Fixed residual polish anchor stealing (`fc926ec`); resolved PE link strand orientation leak and repeat length vs insert size constraint (`tests/test_c0_02r_pe_evidence.rs`, verified in Slurm job `2979536`). Full 3-way gate pending. |

---

## 3. HPC Cluster Environment & Compute Resources

- **Scheduler**: Slurm 22.05+
- **Default Work Allocation**: Partition `qfree`, Node `vnode1`, 8 CPUs (`-c 8`), 12 GB RAM (`--mem=12G`).
- **Python**: `/hpcfs/fpublic/app/miniforge3/conda/bin/python3` (Python 3.12+ with pytest, uv).
- **Rust Toolchain**: `rustc 1.83.0+` (`~/.cargo/bin`).
- **QUAST**: `~/.local/bin/quast.py`.
- **Minimap2**: System / local binary in `PATH`.

---

## 4. Current State & Pending Handoff Actions

See [`ACTIVE_TASK.md`](../../ACTIVE_TASK.md) for the active checklist and immediate execution roadmap.
