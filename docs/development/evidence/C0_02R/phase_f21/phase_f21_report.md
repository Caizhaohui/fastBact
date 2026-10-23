# fastBact — Phase F21-A Comprehensive Forensic Investigation Report

**Investigation ID**: `FORENSIC-20261023-001`  
**Phase Status**: `PHASE_F21_FORENSIC_COMPLETE_PENDING_DECISION`  
**C0-02R Acceptance**: **Strictly NOT ACCEPTED**  
**Evidence Worktree**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f21-forensics-worktree`  
**Branch**: `feature/c0-02r-phase-f21-forensics`  
**Base Commit**: `4b956a46935248779e17b8e3453f4bcb188f8e1f` (Phase F20 evidence HEAD)  
**Production Tree**: `30967bf4940844aeb6b744732bcad6d2a72daab4` (**0 lines changed in `src/`**)  
**Slurm Benchmark Jobs**: Job `3001561` (F20 benchmark) & Job `3007620` (F21 Stage 5 QUAST) on compute node `vnode1`

---

## 1. Executive Summary & Direct Answers to Core Questions

Phase F21-A conducted a multi-stage, reference-coordinate level forensic investigation into the assembly quality regressions observed on isolate `R1-001` during the Phase F20 selective real-data benchmark.

### Core Question 1: Does the extra extensive misassembly ($4 \to 5$) on R1-001 represent a real new error?
**YES (NEW_JOIN_ERROR_CONFIRMED)**.
- Baseline A and Candidate B share **4 identical common misassembly events** (two inversions at 1.20 Mb, one relocation at 1.87 Mb $\leftrightarrow$ 291 kb, and one relocation at 1.87 Mb $\leftrightarrow$ 1.97 Mb).
- The 5th misassembly is exclusively present on Candidate B's `contig_16` (length 71,975 bp).
- It is an authentic chimeric relocation of **-2,095,535 bp**, joining reference locus `3,666,175 bp` directly to `1,570,641 bp`.
- The join bridges **7 short intermediate nodes** (totalling 2,092 bp unaligned length: `u1926`, `u3275`, `u1475`, `u1476`, `u4349`, `u4465`, `u3919`) across an *E. coli* homologous repeat family. It exists bit-identically in Stage 5 `unitig_1631`, confirming that consensus polishing did not create the join.
- While the structural defect is confirmed, locus-specific causal attributions (e.g. cross-hybridizing repeat pairs vs search traversal reachability) are supported mechanistic hypotheses rather than proven facts.

### Core Question 2: Are the increased mismatches ($2.14 \to 2.38$) and indels ($0.48 \to 1.03$) primarily introduced by polishing?
**YES (POLISHING_BOUNDARY_REGRESSION_SUPPORTED)**.
- At Stage 5 (pre-polishing, `stage_05_exspander.fa`), Candidate B has **FEWER mismatches (24 vs 29)** and **IDENTICAL indels (4 vs 4, len 176 bp)** compared to Baseline A.
- **The net excess of mismatches and indels emerged entirely during the Stage 5 to Stage 6 transition.**
- **84.6% of Candidate B's errors (132/156) are concentrated within the $101–350\text{ bp}$ terminal boundary zone of contigs.**
- In deep internal sequence ($> 350\text{ bp}$ from contig ends), Candidate B has **FEWER total errors than Baseline A (20 vs 22: 18 mismatches, 2 indels vs 18 mismatches, 4 indels)**.
- The base error elevation is an artifact of consensus polishing read recruitment at newly created contig ends, not an algorithm degradation inside the assembled sequence.

### Core Question 3: What is the mechanistic link between residual active unitigs ($315 \to 1728$) and the genome fraction drop?
**OUTPUT_FILTER_SEQUENCE_LOSS_SUPPORTED**.
- When exSPAnder withholds ambiguous repeat joins, it leaves short unmerged unitigs in the graph.
- In Candidate B, **1,370 unitigs (185.6 kb)** are $< 200\text{ bp}$, and **191 unitigs (49.4 kb)** are $200–499\text{ bp}$. Because fastBact's default output cutoff is 500 bp, all 1,561 short unitigs (totalling 235 kb, 90.3% of all unitigs in the graph) were not emitted into final `contigs.fa`.
- At Stage 5 (where all contigs are evaluated), Candidate B actually achieved a **HIGHER genome fraction than Baseline A (98.835% vs 98.456%, +0.379 pp)**. The apparent loss in Stage 6 is an artifact of length threshold filtering applied to unmerged repetitive unitigs.

### Core Question 4: Why did R1-003 improve while R1-001 regressed?
- Both isolates possess the exact same homologous repeat family between ~1.56 Mb and ~3.66 Mb.
- In `R1-003`, Baseline A **already committed** this chimeric join (`contig_16`). Candidate B also committed it (`contig_19`), but Candidate B **eliminated a second chimeric join** (+590 kb relocation on `contig_4`), cutting misassemblies from 2 to 1 and dropping indels to 0.00.
- In `R1-001`, Baseline A avoided the join because its search did not traverse the repeat chain. Candidate B's multi-hop search generalization bridged the 7 short intermediate nodes, accepting the chimeric join for the first time on this isolate (+1 misassembly).

---

## 2. Evidence Artifact Inventory

| Artifact | File Path | Key Content / Finding |
| :--- | :--- | :--- |
| **Forensic Manifest** | [`forensic_manifest.tsv`](./forensic_manifest.tsv) | Bit-exact verification of all inputs, binaries, tools, Slurm jobs, and parameters. |
| **Misassembly Matrix** | [`misassembly_event_matrix.tsv`](./misassembly_event_matrix.tsv) | Single-base reference coordinate comparison of all 12 events across A and B. |
| **Locus Analysis** | [`misassembly_locus_analysis.md`](./misassembly_locus_analysis.md) | Dissection of `contig_16`, tiling of Stage 4 unitigs, and repeat homology. |
| **Stage 5 vs 6 Table** | [`stage5_stage6_comparison.tsv`](./stage5_stage6_comparison.tsv) | 8-condition QUAST comparison proving pre-polish vs post-polish error origin. |
| **Base Error Loci** | [`base_error_loci.tsv`](./base_error_loci.tsv) | 396 individual mismatch and indel coordinates with pre/post-polish attribution. |
| **Boundary Distribution** | [`boundary_error_distribution.tsv`](./boundary_error_distribution.tsv) | Bin-stratified error distribution proving 84.6% boundary concentration. |
| **Unitig Fate Summary** | [`unitig_fate_summary.tsv`](./unitig_fate_summary.tsv) | Fate accounting of all 1,728 unitigs proving length threshold sequence loss. |
| **R1-003 Contrast** | [`r1_003_contrast.md`](./r1_003_contrast.md) | Mechanistic contrast explaining the divergence between R1-001 and R1-003. |
| **Hypothesis Matrix** | [`causal_hypothesis_matrix.md`](./causal_hypothesis_matrix.md) | Evidentiary evaluation of H1–H5 with OBSERVED / SUPPORTED tags. |
| **Independent Audit** | [`codex_independent_forensic_audit.md`](./codex_independent_forensic_audit.md) | Codex CLI independent read-only verification report. |

---

## 3. Engineering Recommendations for Next Steps

1. **Maintain Zero-Change Policy on Production Code**:
   - Production code in `src/` remains completely frozen.
   - C0-02R acceptance remains strictly **NOT ACCEPTED**.
2. **Path to Eliminating the Chimeric Join (`contig_16`)**:
   - The root cause of EV-05 on R1-001 (and EV-02 on R1-003) is that exSPAnder permits multi-hop paths through many short nodes without verifying that the supporting read pairs do not cross-hybridize with homologous repeats elsewhere in the genome.
   - Future repair should introduce a **multi-hop repeat penalty** or a **repeat-family uniqueness filter** that requires unambiguous anchor support when traversing $> 4$ consecutive short nodes ($< 500\text{ bp}$).
3. **Path to Eliminating Polishing Boundary Errors**:
   - Future work on consensus polishing (`src/polish/consensus.rs`) should adjust terminal read trimming or require minimum terminal coverage ($\ge 10×$) before applying consensus base corrections within $\le 350\text{ bp}$ of contig ends.
