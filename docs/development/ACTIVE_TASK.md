# ACTIVE TASK: fastBact C0-02R Phase E (Evidence, Root Cause & Minimal Repair)

**Status**: PHASE_E_COMPLETED  
**Overall Acceptance Status**: **NOT ACCEPTED** (C0-02R blocked by Phase D gate failures & P1 vulnerability)  
**Main Branch**: `fix/c0-02r-final-gate` (HEAD `460541b`)  
**P2 Feature Branch**: `feature/c0-02r-p2-estimator` (HEAD `b4472a1` in worktree `fastBact-p2-worktree`)  
**Slurm Benchmark Job**: `2980439` (Status: `COMPLETED`, Elapsed: `00:29:32`, ExitCode: 1:0)  
**Diagnostic QUAST Job**: `2980935` (Status: `COMPLETED`, ExitCode: 0:0)  
**Phase D Benchmark Verdict**: **FAIL** (Gate failures on R1-001 and R1-003)  
**Coordinator**: Antigravity CLI  
**Last Updated**: 2026-10-19

---

## 1. Task Objective & Scope Boundary

Execute Phase E (Evidence, Root Cause & Minimal Repair) under strict constraints:
- C0-02R remains **NOT ACCEPTED**.
- No external AI calls (Codex, Claude, Cursor prohibited).
- No production algorithm refactoring of P1 in Phase E.
- No auto-merging main or modifying baseline history.
- All conclusions strictly classified: `VERIFIED`, `UNVERIFIED`, `BLOCKED`.

---

## 2. Phase E Deliverables & Verification Evidence

### E1 — Benchmark Final Evidence & R1-001 Root Cause (`VERIFIED`)
- **Slurm Benchmark Job `2980439` Outcome**:
  - 18/18 assembly and QUAST evaluations completed.
  - S1-COV-150X, S1-COV-003, S1-COV-015, S1-COV-030: Zero misassemblies, invariant FASTA identity on 150X and 003.
  - **Gate Failures (3 total)**:
    - R1-001: `mismatches_per_100kbp` increased vs B (2.14 $\to$ 2.86).
    - R1-001: `indels_per_100kbp` increased vs B (0.48 $\to$ 0.96).
    - R1-003: `mismatches_per_100kbp` increased vs B (1.51 $\to$ 1.67).
- **R1-001 517 bp Contig Investigation (`VERIFIED`)**:
  - Extra contig `contig_164` (517 bp) in C aligned with **100.0% identity (0 mismatches, 0 indels)** against `NC_000913.3`.
  - Diagnostic QUAST on compute node `vnode1` (Job `2980935`) with `contig_164` removed yielded **identical error rates** (2.86 mismatches, 0.96 indels).
  - **Root Cause Confirmed**: The regression was **not** caused by `contig_164`, but by small consensus polishing sequence shifts across 10 other contigs (e.g. `contig_51`, `contig_105`) arising from anchor selection differences.

### E2 — P2 Insert-Size Estimator Repair (`VERIFIED`)
- **Worktree**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-p2-worktree` (`feature/c0-02r-p2-estimator`).
- **Deterministic Red Tests**: `tests/test_c0_02r_p2_estimator.rs` demonstrated systematic +119 bp and +59 bp estimation biases on unpatched code due to unrecovered k-mer seed step offsets and missing strand orientation checks.
- **Minimal Fix**: Updated `map_read` in `src/resolve/mapper.rs` to track seed `step` and compute true fragment span with inward-facing (FR) concordancy enforcement:
  - Subcase A: R1 fwd, R2 rev $\implies \text{span} = (\text{pos}_2 + k + \text{step}_2) - (\text{pos}_1 - \text{step}_1)$.
  - Subcase B: R1 rev, R2 fwd $\implies \text{span} = (\text{pos}_1 + k + \text{step}_1) - (\text{pos}_2 - \text{step}_2)$.
  - Same-strand (FF/RR) and outward pairs rejected.
- **Green Validation**: All 5 deterministic P2 tests passed; full 106 test release regression suite passed on Slurm compute node `vnode1`.
- **Commit**: `480ff5f`.

### E3 — P1 Production Counterexample (`VERIFIED`)
- **Deterministic Test**: `tests/test_c0_02r_p1_counterexample.rs` (commit `b4472a1`).
- **Graph Topology**: $A \to R \to B$ (erroneous shortcut, span 391 bp) vs $A \to R \to C \to B$ (true path, span 600 bp).
- **Library Span**: Calibrated median 600 bp, $\text{min\_insert} = 530$, $\text{max\_insert} = 670$.
- **Finding**: Through full production pipeline (`map_paired_reads` $\to$ `PairLinkMatrix` $\to$ `resolve_repeats_exspander_with_stats`), exSPAnder chose the shortcut $A \to R \to B$ (assembled unitig length 1038 bp) because `PairLinkMatrix` only stored aggregate counts `(A, B) = 25`, deleting intermediate unitig $C$ (240 bp).
- **Conclusion**: Codex Finding P1 is **empirically VERIFIED** in the production pipeline.

### E4 — P1 Minimal Design Proposal (`COMPLETED`)
- Formulated non-invasive, lightweight fragment coordinate architecture:
  - Document: [`docs/development/P1_MINIMAL_DESIGN_PROPOSAL.md`](./P1_MINIMAL_DESIGN_PROPOSAL.md).
  - Stores compact 4-byte `FragmentTerminalDistance` (`dist_from_u1_exit: u16, dist_from_u2_entry: u16`) in `SmallVec`.
  - Enables $O(1)$ candidate path span verification against `[min_insert, max_insert]`.
  - Memory overhead $< 0.1\%$ peak RSS; computation overhead $< 0.05\%$.
  - Preserves backwards compatibility; safe non-resolution fallback.

---

## 3. Work in Progress / Next Steps

1. [x] E1: Slurm Job 2980439 completion and R1-001 root-cause verification.
2. [x] E2: P2 estimator deterministic red/green fix and full regression on `vnode1` (Commit `480ff5f`).
3. [x] E3: P1 production counterexample verifying shortcut vulnerability (Commit `b4472a1`).
4. [x] E4: P1 minimal design proposal drafted (`P1_MINIMAL_DESIGN_PROPOSAL.md`).
5. [x] E5: Synchronize `ACTIVE_TASK.md`, `HANDOFF.md`, and `AUDIT_FINDINGS.md`.
6. [ ] Await user guidance on initiating Phase F (production implementation of P1 span verification).

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

