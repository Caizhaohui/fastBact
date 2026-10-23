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

### [AUDIT-20261019-001] Codex Medium Read-Only Review for PE Fidelity & Insert-Size Bounds

- **Task ID**: `C0-02R-PE-FIX-AUDIT`
- **Agent**: Codex CLI (v0.163.0)
- **Model**: `gpt-6.1-sol`
- **Reasoning Effort**: `medium` (`-c model_reasoning_effort="medium"`)
- **Permission Level**: `READ_ONLY` (`-s read-only`)
- **Git SHA**: `f241d9d3c9d21470f63bd85b5e1bd6eb0bfc8d30`
- **CLI Command**:
  ```bash
  codex exec -s read-only -m gpt-6.1-sol -c model_reasoning_effort="medium" -o docs/development/evidence/C0_02R/phase_d/codex_review_pe_fix.md - < docs/development/evidence/C0_02R/phase_d/.prompt_codex_review.txt
  ```
- **Exit Code**: `0` (Completed; 63,169 tokens used)
- **Output Path**: [`docs/development/evidence/C0_02R/phase_d/codex_review_pe_fix.md`](./evidence/C0_02R/phase_d/codex_review_pe_fix.md)
- **Verdict**: **FAIL**
- **Findings**:
  - **Orientation Verification**: Confirmed that `PairLinkMatrix::get_support()` and `map_paired_reads()` in commit `f241d9d` correctly handle inward-facing read pairs without orientation mixing, normalization errors, or double counting.
  - **P1 (Unitig-pair counts can support a path incompatible with authentic fragments)**:
    - *Location*: `src/resolve/exspander.rs:223`, lines 241–244; `src/resolve/mapper.rs:144`, lines 158–159.
    - *Defect*: `PairLinkMatrix` only records counts of links between unitigs without preserving read offset/coordinates. When multiple alternate paths connect `A` and `B` (e.g. `A -> R -> C -> B` true 600 bp span vs shortcut `A -> R -> B`), read pairs originating from the longer traversal produce links between `(A, B)`. The direct exit `B` from `R` has `gap = 60 <= max_insert`, so the aggregate count selects the shortcut, even though the true fragments do not span 360 bp.
    - *Coverage*: Not covered by existing synthetic unit tests.
    - *Fix Direction*: Retain oriented read placements, offsets, and fragment lengths; evaluate candidate path support against implied fragment span.
  - **P2 (Empirical insert estimator uses k-mer anchor positions rather than read positions)**:
    - *Location*: `src/resolve/mapper.rs:116`, lines 123, 147–153.
    - *Defect*: Intra-unitig distance calculation `d = (pos2 - pos1).abs() + r1.len()` uses k-mer anchor positions instead of read endpoints, shifting apparent insert sizes by up to `read_len - k`.
    - *Coverage*: Not covered by existing unit tests.
    - *Fix Direction*: Recover oriented read intervals from anchor coordinates and sampled offsets; estimate distribution only from concordant inward-facing pairs.
  - **Final Question Answer**: `Does gap <= max_insert provide sufficient paired-end evidence for accepting the selected path?` **NO**. `gap` only measures net separation across the repeat node (`start(B) - end(A)`), not actual fragment span (`(|A| - s) + gap + e`). It is a necessary geometric lower-bound exclusion filter, but not sufficient evidence for path acceptance.
- **Phase E Empirical Verification (2026-10-19)**:
  - **P1 Status**: **`VERIFIED`**. Confirmed in production pipeline test `tests/test_c0_02r_p1_counterexample.rs` (commit `b4472a1`). exSPAnder erroneously selected shortcut $A \to R \to B$ (1038 bp) over true path $A \to R \to C \to B$ (1247 bp) because `PairLinkMatrix` only stored aggregate counts `(A, B) = 25`. Minimal design proposal drafted in [`docs/development/P1_MINIMAL_DESIGN_PROPOSAL.md`](./P1_MINIMAL_DESIGN_PROPOSAL.md).
  - **P2 Status**: **`VERIFIED` & `REPAIRED`** (in feature worktree). Confirmed in deterministic red test `tests/test_c0_02r_p2_estimator.rs` (+119 bp estimation error on unpatched code). Minimal repair implemented in `src/resolve/mapper.rs` (commit `480ff5f`), restoring true fragment span recovery across arbitrary seed offsets and inward strand orientations. All 106 tests green.
  - **R1-001 517 bp Contig Status**: **`VERIFIED` (False Hypothesis Eliminated)**. The extra 517 bp contig `contig_164` aligned 100.0% identity (0 mismatches, 0 indels) to `NC_000913.3`. Diagnostic QUAST (Slurm Job `2980935`) with `contig_164` removed yielded identical error rates (2.86 mismatches, 0.96 indels). The true root cause is consensus polishing base shifts across 10 other contigs (e.g. `contig_51`, `contig_105`) caused by anchor selection differences.
- **Next Action**:
  - Maintain C0-02R as **NOT ACCEPTED**.
  - Await user guidance on initiating Phase F (production implementation of P1 fragment span verification and consensus polishing stabilization).

---

## Prepared Audits (Awaiting User Approval)

*(None currently pending approval)*


---

### [AUDIT-20261022-014] Codex Independent Incremental Audit of P2 Verification Repairs

- **Task ID**: `C0-02R Phase F19`
- **Agent**: Codex CLI (v0.163.0)
- **Model**: `gpt-6.1-sol`
- **Reasoning Effort**: `medium` (`-c model_reasoning_effort="medium"`)
- **Permission Level**: `READ_ONLY` (`-s read-only`)
- **Target Branch**: `feature/c0-02r-phase-f19-codex`
- **Worktree**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f19-codex-worktree`
- **Audit Target HEAD**: `a9930471a73612159a83694d3d661f060042efc8` (Commit `a993047`)
- **Baseline Commit**: `39f11ca90a5de13a5b3d5572d29c191c40d9c03d` (Commit `39f11ca`)
- **CLI Command**:
  ```bash
  codex exec -C /hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f19-codex-worktree \
    -s read-only -m gpt-6.1-sol -c model_reasoning_effort="medium" \
    -o /hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f19-codex-worktree/docs/development/evidence/C0_02R/phase_f19/codex_f19_incremental_audit.md \
    - < docs/development/evidence/C0_02R/phase_f19/.prompt_codex_incremental_audit.txt
  ```
- **Exit Code**: `0` (Completed; 89,794 tokens used)
- **Output Path**: [`docs/development/evidence/C0_02R/phase_f19/codex_f19_incremental_audit.md`](./evidence/C0_02R/phase_f19/codex_f19_incremental_audit.md)
- **Overall Verdict**: **FAIL** (on claim of 100% resolution of all 5 P2 findings; 4 resolved, 1 open P2)
- **Production Code Status**: **COMPLETELY IMMUTABLE & VERIFIED SOUND** (`git diff 39f11ca..HEAD -- src/` is 0 lines). Zero P0 and zero P1 confirmed blockers.
- **C0-02R Acceptance**: **NOT ACCEPTED** (Strictly maintained pending final gate resolution)
- **Finding Dispositions**:
  1. **F18-P2-001 (Oracle Reciprocal Selection)**: **CLOSED**. `independent_resolve_edge` now evaluates all reverse edges to `from.twin()`, removes premature filtering, enforces reciprocal support equality (`rev_max_support == chosen_support`), and enforces strict uniqueness on the reverse strand. All 5 controls verified sound.
  2. **F18-P2-002 (Oracle Initialization Long-Fork Reuse)**: **CLOSED**. `reference_oracle_search` rejects initial hop if target is a long fork node self-looping back to itself; preserves short repeat unrolling. Both cases verified in differential test.
  3. **F18-P2-003 (T30 Reachability & Sensitivity Control)**: **CLOSED**. Repaired junction sequence mismatch ($X$ starts with `j_b`), adjusted $C$ length to 129 bp so Decision 1 commits cleanly and reaches Decision 2 at $P_2$. Cloned 6-edge topology sensitivity control proves downstream evidence sensitivity.
  4. **F18-P2-004 (T37 $\beta$ Boundary Sensitivity)**: **CLOSED**. Self-loop repeat re-entry evaluates sharp mathematical boundary $\text{sep}_v + \beta == \text{max\_insert}$ ($L_P = 150 \implies \text{Probe} \implies \text{Withheld}$, $L_P = 151 \implies \text{SafePrune} \implies \text{Committed}$).
  5. **F18-P2-005 (T39 Queue Capacity Limit)**: **PARTIALLY RESOLVED (Open P2)**. Queue saturation at $N = 1025$ is reached, but state budget limit ($S_{\text{MAX}} = 256$) produces identical truncation if queue capacity guard is removed in black-box test. Requires a queue-specific sensitivity assertion (such as an internal search test) to close.
- **Next Action**:
  - Maintain C0-02R strictly as **NOT ACCEPTED**.
  - Await user guidance on whether to address the T39 queue-specific sensitivity assertion or proceed directly to selective real-data benchmarking (`R1-001`, `R1-003`).

---

### [AUDIT-20261022-015] Codex Independent Real-Data Benchmark Audit of Phase F20

- **Task ID**: `C0-02R Phase F20`
- **Agent**: Codex CLI (v0.163.0)
- **Model**: `gpt-6.1-sol`
- **Reasoning Effort**: `medium` (`-c model_reasoning_effort="medium"`)
- **Permission Level**: `READ_ONLY` (`-s read-only`)
- **Target Branch**: `feature/c0-02r-phase-f19-codex`
- **Worktree**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f19-codex-worktree`
- **Audited HEAD Commit**: `0606be423d187771f5925a83da2cc377c40c2a3d` (Code: `a9930471a73612159a83694d3d661f060042efc8`, Tree: `30967bf4940844aeb6b744732bcad6d2a72daab4`)
- **Baseline Commit**: `b9727c3ab372f7b995e439d92b72c7a4bfd051cf` (C0-02 baseline)
- **Slurm Benchmark Job ID**: `3001561` on `vnode1` (Exit code `0:0`, 00:13:13 elapsed)
- **Benchmark Scope**: Selective real-data A/B benchmark (4 cells: R1-001 A, R1-001 B, R1-003 B, R1-003 A)
- **Output Report Path**: [`docs/development/evidence/C0_02R/phase_f20/codex_independent_benchmark_audit.md`](./evidence/C0_02R/phase_f20/codex_independent_benchmark_audit.md)
- **Overall Verdict**: **FAIL (REGRESSION DETECTED)**
- **C0-02R Acceptance**: **Strictly NOT ACCEPTED**
- **Independent Acceptance Evaluation**:
  - Exactly 56 metrics extracted across all 4 assemblies.
  - Acceptance criteria: 7 comparisons PASS, 7 comparisons REGRESSION $\implies$ Status: `REGRESSION_DETECTED`.
  - **R1-003 (Marked Improvement)**:
    - Extensive misassemblies: 2 $\to$ 1 (-50%, PASS)
    - Mismatches/100kbp: 1.51 $\to$ 1.08 (-28.5%, PASS)
    - Indels/100kbp: 0.04 $\to$ 0.00 (-100%, PASS)
    - N50 bp: 67,452 $\to$ 70,868 (+5.06%, PASS)
    - NGA50 bp: 63,648 $\to$ 70,010 (+10.00%, PASS)
    - Local misassemblies: 1 $\to$ 1 (0%, PASS)
    - Genome fraction: 98.330% $\to$ 98.233% (-0.097 pp, REGRESSION)
  - **R1-001 (Quality Regression)**:
    - Extensive misassemblies: 4 $\to$ 5 (+25%, REGRESSION)
    - Mismatches/100kbp: 2.14 $\to$ 2.38 (+11.2%, REGRESSION)
    - Indels/100kbp: 0.48 $\to$ 1.03 (+114.6%, REGRESSION)
    - N50 bp: 58,318 $\to$ 57,462 (-1.47%, REGRESSION)
    - NGA50 bp: 58,318 $\to$ 57,127 (-2.04%, REGRESSION)
    - Genome fraction: 98.364% $\to$ 98.281% (-0.083 pp, REGRESSION)
    - Local misassemblies: 0 $\to$ 0 (PASS)
    - Active unitigs post-exspander: 315 $\to$ 1,728 (5.49×, conservative withholding leaves small unitigs unmerged below 500 bp cutoff)
- **Auditor Findings**:
  1. *Acceptance Blocker*: R1-001 fails 6 quality criteria (misassemblies, mismatches, indels, genome fraction, N50, NGA50).
  2. *Acceptance Blocker*: R1-003 fails genome fraction (-0.097 pp) despite significant improvements in all 5 other criteria.
  3. *Provenance & Execution Integrity*: All raw input FASTQs, reference genomes, and binaries match recorded SHA-256 digests bit-for-bit.
  4. *Reproducibility*: Baseline A assemblies reproduce historical benchmark outputs bit-for-bit.
- **Next Action**:
  - Stop automated benchmark expansion. Do not expand to full 18 datasets.
  - Maintain C0-02R strictly as **NOT ACCEPTED**.
  - Conduct focused forensic investigation into R1-001's additional misassembly breakpoint and indel loci before modifying any code.

---

### [AUDIT-20261023-016] Codex Independent Forensic Audit of Real-Data Regression

- **Task ID**: `C0-02R Phase F21-A`
- **Agent**: Codex CLI (v0.163.0)
- **Model**: `gpt-6.1-sol`
- **Reasoning Effort**: `medium` (`-c model_reasoning_effort="medium"`)
- **Permission Level**: `READ_ONLY` (`-s read-only`)
- **Target Branch**: `feature/c0-02r-phase-f21-forensics`
- **Worktree**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21-forensics-worktree`
- **Audited Commit Baseline**: `4b956a46935248779e17b8e3453f4bcb188f8e1f` (Phase F20 evidence HEAD)
- **Production Tree SHA**: `30967bf4940844aeb6b744732bcad6d2a72daab4` (**0 lines changed in `src/`**)
- **Output Report Path**: [`docs/development/evidence/C0_02R/phase_f21/codex_independent_forensic_audit.md`](./evidence/C0_02R/phase_f21/codex_independent_forensic_audit.md)
- **Overall Forensic Verdict**: **QUALIFIED PASS / AUDITED WITH CORRECTIONS** (All factual, structural, and numerical measurements verified; locus-specific decision-level causal attributions qualified as supported hypotheses)
- **C0-02R Acceptance**: **Strictly NOT ACCEPTED**
- **Key Forensic Findings Independently Replicated**:
  1. *Structural Misassembly Ground Truth*: Candidate B's additional extensive misassembly on R1-001 (`contig_16`, length 71,975 bp, relocation of -2,095,535 bp between 3.66 Mb and 1.57 Mb) is an **authentic chimeric join** formed during Stage 5 (`unitig_1631` bit-identical to `contig_16`). It connects `unitig_1927` to `unitig_169` through 7 short intermediate nodes (`u1926`, `u3275`, `u1475`, `u1476`, `u4349`, `u4465`, `u3919`, total 2,092 bp). Tag: `NEW_JOIN_ERROR_CONFIRMED`.
  2. *Base Error Pre/Post-Polish Attribution*: At Stage 5 (pre-polishing), Candidate B has **fewer mismatches (24 vs 29)** and **identical indels (4 vs 4, len 176 bp)** compared to Baseline A. The entire excess of mismatches and indels emerged during the Stage 5 to Stage 6 transition.
  3. *Boundary Concentration*: **84.6% of Candidate B's errors (132/156)** are located in the $101–350\text{ bp}$ terminal boundary zone. Deep inside contigs ($> 350\text{ bp}$), Candidate B has fewer errors than Baseline A (20 vs 22). Tag: `POLISHING_BOUNDARY_REGRESSION_SUPPORTED`.
  4. *Unitig Fate & Genome Fraction*: 1,561 short unitigs (totalling 235 kb, 90.3% of all unitigs in the graph) fell below the 500 bp contig output cutoff. At Stage 5, Candidate B's genome fraction was actually higher (+0.38 pp on R1-001, +0.41 pp on R1-003). Tag: `OUTPUT_FILTER_SEQUENCE_LOSS_SUPPORTED`.
  5. *R1-003 vs R1-001 Divergence*: Both isolates share the 1.56 Mb $\leftrightarrow$ 3.66 Mb repeat family. R1-003 improved because Candidate B eliminated a second false join (+590 kb relocation on `contig_4`), while R1-001 regressed because Baseline A had avoided the 2.1 Mb repeat whereas Candidate B bridged it.
- **Production Code Status**: **100% BIT-IDENTICAL & IMMUTABLE** (0 lines changed in `src/`).
- **Next Action**:
  - Maintain C0-02R strictly as **NOT ACCEPTED**.
  - Await user review of the forensic report before designing any repeat-family uniqueness filters or polishing boundary safeguards.

---

### [AUDIT-20261023-017] Codex Independent Trace Audit of Repeat-Join Decision Provenance & Mechanism Verification

- **Task ID**: `C0-02R Phase F21-B`
- **Agent**: Codex CLI (v0.163.0)
- **Model**: `gpt-6.1-sol`
- **Reasoning Effort**: `medium` (`-c model_reasoning_effort="medium"`)
- **Permission Level**: `READ_ONLY` (`-s read-only`)
- **Target Branch**: `feature/c0-02r-phase-f21b-trace`
- **Worktree**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21b-trace-worktree`
- **Audited Commit Baseline**: `07e0e5fd927cc3c8be96b22927e33b8ae6b9e138` (Phase F21-A forensic HEAD)
- **Production Tree SHA**: `30967bf4940844aeb6b744732bcad6d2a72daab4` (**0 lines changed in `src/`**)
- **Output Report Path**: [`docs/development/evidence/C0_02R/phase_f21b/codex_independent_trace_audit.md`](./evidence/C0_02R/phase_f21b/codex_independent_trace_audit.md)
- **Overall Forensic Verdict**: **FAIL / QUALIFIED FINDINGS** (Strict verification: verified production tree frozen at `30967bf...`, confirmed R1-003 Baseline A relocation, confirmed terminal polishing error concentration; qualified offline decision scores and required formal diagnostic trace instrumentation)
- **Phase Completion Status**: **`PHASE_F21B_DIAGNOSTIC_EVIDENCE_REQUIRED`**
- **C0-02R Acceptance**: **Strictly NOT ACCEPTED**
- **Key Findings & Reconciliations**:
  1. *Production Code Freeze*: Verified bit-identical production tree `30967bf4940844aeb6b744732bcad6d2a72daab4` with 0 lines changed in `src/`.
  2. *Chimeric Relocation Verified*: Structural misassembly on `contig_16` (-2,095,535 bp relocation between 3.66 Mb and 1.57 Mb across 7 short nodes) is 100% verified.
  3. *Cross-Isolate Baseline Pre-existence*: Corroborated that Baseline A in R1-003 already committed the exact same chimeric relocation on `contig_16` (-2,094,638 bp), confirming Candidate B reproduced existing graph behavior while resolving an independent 590 kb misassembly.
  4. *Terminal Polishing Artifacts*: 84.6% (132/156) of base errors occur in 101–350 bp terminal zone. Net increases from Stage 5 $\to$ 6 are +85 mismatches and +43 indels for B vs +69 mismatches and +18 indels for A.
  5. *Unique Reference Coverage Loss*: Reconciled component rows in `unique_coverage_loss.tsv` to 25,574 bp unique loss (mapped 32,601 bp - redundant 7,027 bp), accounting for the QUAST reference coverage drop of 25,715 bp (0.554%).
  6. *Resolver Decision Traversal Qualification*: While the graph topology and 0-spanning-pair reality are confirmed, the exact walker state transitions and decision weights are documented offline interpretations. Dedicated runtime diagnostic tracing instrumentation is required in Phase F22.
- **Next Action**:
  - Maintain C0-02R strictly as **NOT ACCEPTED**.
  - Update Phase F21-B governance records and report to `PHASE_F21B_DIAGNOSTIC_EVIDENCE_REQUIRED`.
  - Design Phase F22 targeted interventions: flanking unique-anchor requirement for repeat traversals, terminal polishing boundary guard, and runtime diagnostic logging.

