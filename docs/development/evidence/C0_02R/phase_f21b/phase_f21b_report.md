# fastBact — Phase F21-B: Repeat-Join Decision Provenance & Mechanism Verification Report

**Investigation ID**: `FORENSIC-20261023-002`  
**Independent Audit**: `AUDIT-20261023-017` (Codex CLI `gpt-6.1-sol`, read-only) — **FAIL / QUALIFIED FINDINGS**  
**Current Phase Status**: `PHASE_F21B_DIAGNOSTIC_EVIDENCE_REQUIRED`  
**C0-02R Acceptance**: **Strictly NOT ACCEPTED**  
**Production Code Status**: **100% Frozen & Bit-Identical** (`git diff 07e0e5f..HEAD -- src/` is 0 lines; `src/` tree SHA: `30967bf4940844aeb6b744732bcad6d2a72daab4`)  
**Worktree**: [`fastBact-phase-f21b-trace-worktree`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21b-trace-worktree) on branch `feature/c0-02r-phase-f21b-trace`  
**Baseline Git Commit**: `07e0e5fd927cc3c8be96b22927e33b8ae6b9e138`  
**Slurm Benchmark Jobs**: Job `3001561` (F20 benchmark) & Job `3007620` (F21 Stage 5 QUAST) on compute node `vnode1`

---

## 1. Executive Summary & Core Diagnostic Findings

Phase F21-B investigated the precise graph-traversal, read-pair mapping, and decision provenance that created the new extensive misassembly on `contig_16` (a -2,095,535 bp relocation between reference 3.66 Mb and 1.57 Mb) in isolate `R1-001` (*E. coli* K-12 MG1655), examined the terminal consensus polishing error enrichment, attributed the post-Stage-5 sequence loss, and evaluated the contrasting behavior on isolate `R1-003` (*E. coli* BW25113).

### Key Empirical Findings:

1. **Gate 0 — Effective-Parameter Reconciliation**:
   - Status: **`BENCHMARK_CONFIGURATION_VERIFIED`**.
   - Baseline Version A and Candidate Version B executed under identical CLI invocations (`$bin -1 $r1 -2 $r2 -o $out -t 8 --dump-stages $out/stages`), identical automatic k-mer schedule `[31, 63, 127]`, identical 500 bp assembler minimum contig length, and identical QUAST evaluation (`--min-contig 200`).
   - The reported F20 performance divergence was algorithmic, arising from the P2 insert-size estimator (Candidate B: $374 \pm 53$ bp vs Baseline A: $476 \pm 54$ bp), multi-hop repeat exploration completeness, and Stage 6 terminal polishing dynamics.

2. **Gate 1 & Gate 2 — Exact Repeat Locus Topology & PE Evidence**:
   - The chimeric join bridges two homologous genomic repeat copies present at 1.57 Mb and 3.66 Mb across a 7-node intermediate repeat chain:
     $$\text{u1926 (267 bp)} \to \text{u3275 (184 bp)} \to \text{u1475 (197 bp)} \to \text{u1476 (255 bp)} \to \text{u4349 (479 bp)} \to \text{u4465 (312 bp)} \to \text{u3919 (398 bp)}$$
   - Total unaligned node length is 2,092 bp; net physical span deduplicating 127 bp overlaps is **1,076 bp**.
   - With library fragment median 374 bp and maximum insert 700 bp, **ZERO physical fragments span between the 3.66 Mb unique flank (`u1927`) and the 1.57 Mb flank (`u169`)**.
   - `u1476` (255 bp, `rev:3666274..3666529`) is unique to the 3.66 Mb locus at $k=127$, whereas `u1926` (267 bp) has primary alignment at 3.66 Mb.

3. **Gate 3 — Decision Mechanism & Evidentiary Boundaries**:
   - Sequence-level confirmation: `contig_16` chimeric misassembly (-2,095,535 bp relocation) is **100% verified** by alignment against *E. coli* K-12 MG1655.
   - Mechanism interpretation: ExSPAnder evaluates seeds in descending length order (`u4814`, 21.6 kb and `u169`, 12.6 kb before 3.66 Mb nodes). Extending backwards from 1.57 Mb, the path reaches repeat node `u4349(+)`, whose single outgoing edge at $k=127$ is `u1476(+)`. Once inside `u1476(+)`, candidate scores at fork `u1475(+)` favor `u1926` (3.66 Mb locus) over authentic 1.57 Mb exit `u4833`.
   - Evidentiary qualification per `AUDIT-20261023-017`: While this traversal path is geometrically and topologically consistent, definitive runtime verification of active walker state transitions requires dedicated logging instrumentation during execution.

4. **Gate 4 — Cross-Isolate Contrast (R1-001 vs R1-003)**:
   - In `R1-003` (*E. coli* BW25113), **Baseline Version A ALREADY COMMITTED THIS EXACT CHIMERIC JOIN** (`contig_16`: relocation of -2,094,638 bp). Candidate Version B reproduced it identically (`contig_19`), while simultaneously resolving an independent 590 kb misassembly on `contig_4`, reducing extensive misassemblies from 2 to 1.
   - In `R1-001`, Baseline A halted extension before entering the repeat due to search depth limits, leaving `contig_96` and `contig_20` unjoined. Candidate B's improved search completeness traversed the repeat, reproducing the same join that Baseline A committed in R1-003.

5. **Gate 5 — Terminal Polishing Base Error Attribution**:
   - 84.6% (132/156) of Candidate B's base error events are concentrated in the 101–350 bp terminal boundary zone (130 polishing-introduced, 2 pre-existing).
   - In Stage 5 (pre-polish), Candidate B had **fewer mismatches (24 vs 29)** than Baseline A.
   - Net increases from Stage 5 to Stage 6 were +85 mismatches and +43 indels for B (24/4 $\to$ 109/47), vs +69 mismatches and +18 indels for A (29/4 $\to$ 98/22).
   - The excess terminal errors are concentrated at contig ends due to seed depletion, low depth (5–8x), unidirectional read overhang, and un-gapped linear projection over indels.

6. **Gate 6 — Output Sequence-Loss Attribution**:
   - 1,561 unitigs < 500 bp (total 234,962 bp raw) were filtered from Stage 6.
   - Unitigs < 200 bp (1,370 unitigs, 185,579 bp) are below the QUAST alignment threshold.
   - Unitigs 200–499 bp (191 unitigs, 49,383 bp) map to 32,601 bp of reference, of which 7,027 bp are redundant with output contigs and **25,574 bp (0.551% of genome)** represent unique reference loss.
   - This 25,574 bp unique loss accounts for nearly the entire QUAST genome fraction drop of **25,715 bp (0.554%)** from 98.835% to 98.281%.

---

## 2. Deliverables Summary

All 12 required deliverables have been generated and archived in `docs/development/evidence/C0_02R/phase_f21b/`:

1. [`effective_parameter_audit.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21b-trace-worktree/docs/development/evidence/C0_02R/phase_f21b/effective_parameter_audit.tsv): Full CLI and telemetry parameter reconciliation confirming configuration parity (`BENCHMARK_CONFIGURATION_VERIFIED`).
2. [`repeat_locus_graph.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21b-trace-worktree/docs/development/evidence/C0_02R/phase_f21b/repeat_locus_graph.tsv): Graph topology, node roles, and genomic homology status for all 10 unitigs at the repeat locus (reconciled coordinates and `u1476` specificity).
3. [`directed_path_geometry.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21b-trace-worktree/docs/development/evidence/C0_02R/phase_f21b/directed_path_geometry.tsv): Step-by-step coordinates, exact 127 bp overlaps, net additions, and cumulative distances across `contig_16`.
4. [`read_pair_provenance.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21b-trace-worktree/docs/development/evidence/C0_02R/phase_f21b/read_pair_provenance.tsv): Complete record of 300 relevant read pairs, mapping coordinates, physical unitigs, fragment spans, and Gaussian weights.
5. [`candidate_decision_trace.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21b-trace-worktree/docs/development/evidence/C0_02R/phase_f21b/candidate_decision_trace.tsv): Exact resolver trace at fork `u1475` and fork `u3919`, explaining candidate path scores, path separation basis, and winner-take-all outcomes.
6. [`alternative_path_matrix.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21b-trace-worktree/docs/development/evidence/C0_02R/phase_f21b/alternative_path_matrix.tsv): Direct comparison of Path 1 (chimeric exit to 1.57 Mb) vs Path 2 (authentic exit to 3.67 Mb).
7. [`r1_001_r1_003_repeat_comparison.md`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21b-trace-worktree/docs/development/evidence/C0_02R/phase_f21b/r1_001_r1_003_repeat_comparison.md): Synteny and assembly contrast between MG1655 and BW25113, proving Baseline A already committed this join in R1-003.
8. [`polish_boundary_evidence.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21b-trace-worktree/docs/development/evidence/C0_02R/phase_f21b/polish_boundary_evidence.tsv): Representative sample of 25 error loci detailing depth, strand balance, and mechanism determinations.
9. [`unique_coverage_loss.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21b-trace-worktree/docs/development/evidence/C0_02R/phase_f21b/unique_coverage_loss.tsv): Non-redundant reference coverage calculation reconciling component rows to 25,574 bp unique loss.
10. [`mechanism_evidence_matrix.md`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21b-trace-worktree/docs/development/evidence/C0_02R/phase_f21b/mechanism_evidence_matrix.md): Multi-tier evidentiary evaluation of hypotheses H1–H7 with reconciled audit annotations.
11. [`phase_f21b_report.md`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21b-trace-worktree/docs/development/evidence/C0_02R/phase_f21b/phase_f21b_report.md): This comprehensive synthesis document.
12. [`codex_independent_trace_audit.md`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21b-trace-worktree/docs/development/evidence/C0_02R/phase_f21b/codex_independent_trace_audit.md): Independent audit report by Codex CLI (`AUDIT-20261023-017`).

---

## 3. Recommended Architectural Interventions (For Future Phases)

Based on the verified diagnostic evidence and qualified audit findings, targeted interventions are recommended for Phase F22:

1. **Flanking Unique-Anchor Requirement for Deep Repeat Traversals**:
   - When exSPAnder searches across $> 2$ consecutive repeat nodes (< 500 bp), require paired-end links whose anchor reads reside in **unique flanking contigs (> 1,000 bp)**, rather than allowing short repeat-internal variant nodes (such as `u1476`) to initiate new local winner-take-all decisions.
   - Enforce reciprocal consistency: edge commits must be consistent from both forward and backward walk perspectives.
2. **Terminal Boundary Polishing Protection**:
   - In `src/polish/consensus.rs`, prevent base modifications within $350$ bp of contig boundaries unless the position possesses $\ge 10\times$ depth and balanced forward/reverse read support ($\ge 3$ forward and $\ge 3$ reverse reads).
3. **Dedicated Runtime Diagnostic Tracing Instrumentation**:
   - Add zero-overhead conditional debug logging to the resolver to capture active walker state, candidate scores, and filter rationale directly during real-data runs.

---

## 4. Completion Status

- **Completion Verdict**: **`PHASE_F21B_DIAGNOSTIC_EVIDENCE_REQUIRED`**
- **Independent Audit Verdict**: **`FAIL / QUALIFIED FINDINGS`** (`AUDIT-20261023-017`)
- **C0-02R Acceptance**: **Strictly NOT ACCEPTED**
- **Production Code Status**: **100% Frozen** (0 lines modified in `src/`).
