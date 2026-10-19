# ACTIVE TASK: fastBact C0-02R Phase D (Final Gate & Provenance Acceptance)

**Status**: BENCHMARK_RUNNING  
**Branch**: `fix/c0-02r-final-gate`  
**Current HEAD**: `97510e6`  
**Benchmark Target Commit**: `f241d9d` (sources/C worktree)  
**Target Node**: Slurm `vnode1` (`qfree`, `-c 8 --mem=12G`)  
**Slurm Benchmark Job**: `2980439` (Status: `RUNNING`, Started: 17:37:10, Elapsed: ~25 min)  
**Benchmark Output Dir**: `benchmark-results/C0_02R_pe_fix_20261019`  
**Phase D Acceptance Verdict**: **UNKNOWN** (Awaiting completion of Job 2980439)  
**Coordinator**: Antigravity CLI  
**Last Updated**: 2026-10-19

---

## 1. Task Objective

Complete C0-02R (Final Correctness and Evidence Recovery for FB-P0-02 repeat resolution):
1. Resolve all regressions observed on real isolates (R1-001, R1-003) in previous runs (`2978178`, `2978440`).
2. Run clean 3-way A/B/C benchmark (`submit_c0_02r_final_gate.sh`) across all 6 standard datasets on Slurm compute node `vnode1`.
3. Satisfy conservative acceptance criteria:
   - 0 increased extensive or local misassemblies vs A and B.
   - 0 increased mismatches or indels per 100kbp vs A and B.
   - 0 decreased genome fraction vs A and B.
   - Invariant bit-identical FASTA identity on `S1-COV-150X` and `S1-COV-003`.
4. Perform machine-verified provenance audit via `verify_c0_02r_manifest.py`.

---

## 2. Recent Progress & Root-Cause Analysis

### Previous Gate Failure (`2978440`, commit `fc926ec`)
- **R1-001**: mismatches increased vs B (2.14 $\to$ 2.86), indels increased vs B (0.48 $\to$ 0.96).
- **R1-003**: mismatches increased vs B (1.51 $\to$ 1.67).

### Root Causes Identified
1. **Strand & Orientation Leak in `src/resolve/mapper.rs`**:
   - `PairLinkMatrix::get_support()` took `max(direct, norm)`, normalizing both unitigs to forward strand regardless of actual orientation.
   - `map_paired_reads()` redundantly added non-twinned `(u1, u2)` alongside `(u1, u2_forward)`, corrupting physical inward-facing pair semantics.
2. **Missing Insert Size Constraint in `src/resolve/exspander.rs`**:
   - `resolve_repeats_exspander_with_stats()` accepted `_dist: &InsertSizeDist` but never inspected `dist.max_insert`.
   - `extend_one_direction()` traversed repeat bridges regardless of length, permitting read pairs to span repeats wider than physical library fragment distribution.
   - **Verification on Insert-Size Semantic Scope**: The insert-size check `gap <= dist.max_insert as i64` explicitly measures the span across the specific repeat node being bridged (`target_node.seq_len - overlaps`), and does NOT prematurely reject paths simply because cumulative path length exceeds insert size.

### Empirical Evidence & Verification (`2979536`)
- **Deterministic Test**: Added `tests/test_c0_02r_pe_evidence.rs`.
- **Red Run** (HEAD `fc926ec`): Failed 2/3 assertions (link support was direction-agnostic; repeat wider than insert size was falsely resolved).
- **Green Run** (Working tree): 3/3 passed; release full-regression 97+ tests passed (exit 0).
- **Committed**: Commit `f241d9d`.

---

## 3. Work in Progress / Next Steps

1. [x] Establish `ACTIVE_TASK.md` and `docs/development/HANDOFF.md`.
2. [x] Commit PE evidence fixes (`src/resolve/mapper.rs`, `src/resolve/exspander.rs`, `tests/test_c0_02r_pe_evidence.rs`, and manifest verification suite) $\to$ Commit `f241d9d`.
3. [x] Submit clean 3-way benchmark job via `submit_c0_02r_final_gate.sh` to Slurm compute node `vnode1` $\to$ Job `2980439` (target `benchmark-results/C0_02R_pe_fix_20261019`).
4. [ ] Monitor Slurm job 2980439 execution to completion (currently `PENDING (Resources)`).
5. [ ] Verify manifest integrity with `verify_c0_02r_manifest.py` on compute node upon job completion.
6. [ ] Review 6-dataset acceptance matrix; verify zero regression.
7. [ ] Archive Phase D receipts in `docs/development/evidence/C0_02R/phase_d/`.
8. [ ] Execute prepared Codex read-only audit upon user approval.

---

## 4. External Model Invocation Governance & Registry

### Routing & Security Status
- **Codex CLI**: Verified local binary (`~/.local/bin/codex`), Auth: ChatGPT token subscription. Regular execution configured with `-c model_reasoning_effort="medium"`.
- **Claude Code**: Verified wrapper (`~/.local/share/omniakey/bin/claude`). **Security Finding**: Routes to third-party API gateway `https://api.omniakey.com` via API key `~/.config/omniakey/claude-api-key` (not web subscription OAuth). **Status: Automatic invocation strictly PROHIBITED; requires explicit user approval.**
- **Cursor Agent**: Verified local binary (`~/.local/bin/agent`). Default `--mode=ask` enforced for read-only. Status: Capability-checked only.

### Executed Audits
- **[AUDIT-20261019-001] Codex Medium Read-Only Review for PE Fidelity & Insert-Size Bounds**:
  - **Task ID**: `C0-02R-PE-FIX-AUDIT`
  - **Agent**: Codex CLI (`gpt-6.1-sol`, medium reasoning effort, read-only)
  - **Exit Code**: 0 (63,169 tokens)
  - **Verdict**: **FAIL**
  - **Findings**:
    - Orientation handling verified correct; no mixing or double counting.
    - **P1**: Unitig-pair count aggregation allows shortcuts when alternate longer traversals share endpoints (`src/resolve/exspander.rs:223`, `src/resolve/mapper.rs:144`).
    - **P2**: Insert size estimator uses k-mer anchor coordinates instead of read endpoints (`src/resolve/mapper.rs:116`).
    - **Sufficiency of `gap <= max_insert`**: **NO** (scalar gap is only boundary separation, not fragment span; does not establish that authentic fragments span the selected traversal).
  - **Full Report**: [`docs/development/evidence/C0_02R/phase_d/codex_review_pe_fix.md`](./evidence/C0_02R/phase_d/codex_review_pe_fix.md)
  - **Ledger Record**: [`docs/development/AUDIT_FINDINGS.md`](./AUDIT_FINDINGS.md)

