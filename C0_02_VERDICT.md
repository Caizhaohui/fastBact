# fastBact C0-02 — Junction Overlap Correctness Verdict

**Repository**: `https://github.com/Caizhaohui/fastBact`  
**Task ID**: C0-02 (FB-P0-02: exSPAnder Path Splicing Ignores Directed Edge Overlap Correctness Defect)  
**Baseline Git Commit**: `3950b2ab2584b24d39172b25bbbcb565232891b7` (`fix/c0-01-polish-orientation`)  
**Target Branch**: `fix/c0-02-junction-overlap`  
**Execution Node**: HPC Slurm cluster `vnode1` (partition `qfree`, `-c 8 --mem=12G`)  
**Date**: October 16, 2026  

---

## 1. Executive Summary

FB-P0-02 identified a critical correctness defect in the paired-end repeat resolution engine (`src/resolve/exspander.rs`):
1. `extend_one_direction()` returned only unitig identifiers (`DirUnitigId`), permanently discarding the traversed `GraphEdge` and its recorded overlap length.
2. Contig sequence assembly hardcoded `let overlap = (k - 1).min(u_seq.len())`, blindly assuming all junctions have an overlap of $k - 1$.
3. For any junction with an actual overlap of $k$, this hardcoded offset duplicated one base at the junction, injecting widespread 1-bp insertion indels and downstream coordinate shifts across assembled contigs.
4. When a candidate target had both $k$ and $k - 1$ overlap edges, exSPAnder treated the edges as competing candidate nodes, triggering self-competition (winner-take-all ratio $1.0 < 2.5$) and aborting path extension.
5. Path exploration and graph mutations were coupled in-place without pre-commit sequence validation, and post-exSPAnder cleanup unconditionally deleted candidate unitigs (<500 bp) from rejected paths.

The controlled code implementation in `src/resolve/exspander.rs`:
1. Introduced internal `PathStep` capturing directed unitig ID, chosen `GraphEdge`, validated overlap, and support evidence.
2. Implemented `resolve_edge_between()` with explicit, auditable edge selection: validates actual sequence suffix/prefix, selects higher support when unambiguous, and rejects ties as `Ambiguous`.
3. Separated path exploration from graph modification: validates all nodes, edge orientations, and base-by-base sequence overlaps in `validate_and_splice_path()` before atomic commit.
4. Added candidate node protection (`rejected_nodes`), ensuring no valid graph topology is destroyed when candidate paths are rejected.
5. Added comprehensive junction telemetry tracking `valid_k_joins`, `valid_k_minus_1_joins`, `rejected_ambiguous_joins`, and `rejected_invalid_joins`.

### Final Verdict: **PASS**

---

## 2. Evaluation Matrix & Categorization

| Dimension | Evaluation Category | Status | Empirical Evidence & Reference |
| :--- | :--- | :---: | :--- |
| **1. Edge Identity Preservation** | Traversed `GraphEdge` retained with exact overlap | **VERIFIED** | `PathStep` preserves edge identity across forward and backward extensions. Validated in T1–T10 unit tests. |
| **2. Overlap Correctness** | Both $k$ and $k - 1$ joins assembled base-accurately | **VERIFIED** | Full DNA sequence equality verified in T1 ($k-1$), T2 ($k$), and T5 (mixed). Indels dropped to 0.00 / 100 kbp on SYN-2.5M. |
| **3. Reverse-Orientation Correctness** | Stranded RC extraction & twin edge conversion | **VERIFIED** | Stranded sequence extraction and dual edge mapping verified in T4 (RC orientation) and T10 (backward twin reversal). |
| **4. Ambiguous-Edge Handling** | Auditable edge selection without arbitrary order | **VERIFIED** | T3 verifies higher support selection (overlap 31, support 20 > 10); T9 verifies safe rejection on tied support. |
| **5. Graph Mutation Safety** | Atomic commit & candidate node protection | **VERIFIED** | T6 rejects invalid overlap; T8 proves 100% graph invariance after failure; unresolved candidate nodes protected from cleanup. |
| **6. C0-01 Regression** | Reverse polishing & coordinate fix preserved | **VERIFIED** | All 29 C0-01 unit tests passed with 100% success; zero modifications to `src/polish/consensus.rs`; bit-identical output on invariant datasets. |
| **7. Assembly Quality** | No structural degradation, indels minimized | **VERIFIED** | 0 new misassemblies across all 6 datasets; indels dropped by 97.9%–99.7% on real isolates (R1-001, R1-003); indels reduced to 0.00 on synthetic. |
| **8. Benchmark Reproducibility** | Full replication from archived source & data | **VERIFIED** | Binaries built from clean commits; input FASTQ and reference SHA-256 hashes verified; execution logs and Slurm job IDs archived. |

---

## 3. Minimal Counterexample & Baseline Reproduction

Minimal unit test matrix `tests/test_c0_02_junction_overlap.rs` executed against unmodified baseline `3950b2a` on compute node `vnode1` (`-c 8 --mem=12G`):

- **Baseline Failures (8/10 tests failed)**:
  - `test_t2_two_unitigs_valid_overlap_k`: panicked at `Length must be 1169 bp, not 1170 bp` (1-bp duplicate base).
  - `test_t3_same_directed_pair_both_k_and_k_minus_1_edges`: panicked at `Should merge into exactly 1 active unitig` (self-competition abort).
  - `test_t4_reverse_complement_orientation`: panicked on RC sequence comparison ($k-1$ drift).
  - `test_t5_three_unitigs_mixing_k_and_k_minus_1`: panicked on mixed sequence comparison (1740 bp vs 1739 bp).
  - `test_t6_edge_overlap_mismatches_actual_sequence`: panicked at `Splicing must be rejected` (blindly merged corrupted sequence).
  - `test_t8_failure_preserves_graph_state_atomically`: panicked at `Node count unchanged` (small candidate node deleted by cleanup).
  - `test_t9_multiple_overlap_edges_explicit_selection_and_ambiguity_rejection`: panicked on duplicate edge resolution.
  - `test_t10_backward_path_reversal_references_correct_edge_and_overlap`: panicked on backward path twin overlap.
- **Fixed C0-02 Version**: **10/10 tests passed (100% success)**.

Complete failure traces archived in [`C0_02_BASELINE_REPRO.md`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/C0_02_BASELINE_REPRO.md) and [`C0_02_JUNCTION_UNIT_TESTS.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/C0_02_JUNCTION_UNIT_TESTS.tsv).

---

## 4. A/B Benchmark Results (6 Standard Datasets)

All assemblies and QUAST evaluations executed on Slurm compute node `vnode1` (`-c 8 --mem=12G`).

| Dataset ID | Organism / Description | Baseline Mismatches / 100kbp | Fixed Mismatches / 100kbp | Baseline Indels / 100kbp | Fixed Indels / 100kbp | Indel Reduction | k Joins | k-1 Joins | N50 (Base $\to$ Fix) | Status |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **S1-COV-150X** | Synthetic 200 kb (150x) | 0.00 | **0.00** | 0.00 | **0.00** | 0.0% | 0 | 0 | 199,983 $\to$ 199,983 | **VERIFIED** |
| **S1-COV-003** | Synthetic 200 kb (300x) | 0.00 | **0.00** | 0.00 | **0.00** | 0.0% | 0 | 0 | 199,999 $\to$ 199,999 | **VERIFIED** |
| **SYN-2.5M** | Synthetic 2.5 Mb (100x) | 0.00 | **0.00** | 0.40 | **0.00** | **-100.0%** | 10 | 0 | 2,499,977 $\to$ 2,499,967 | **VERIFIED** |
| **S2-REAL-001** | *E. coli* K-12 synthetic (100x) | 0.07 | **0.02** | 0.04 | **0.00** | **-100.0%** | 2 | 7 | 285,687 $\to$ 285,687 | **VERIFIED** |
| **R1-001** | Real *E. coli* MG1655 (78x) | 12.28 | **2.14** | 23.13 | **0.48** | **-97.9%** | 987 | 18 | 58,335 $\to$ 58,318 | **VERIFIED** |
| **R1-003** | Real *E. coli* BW25113 (132x) | 6.79 | **1.51** | 12.25 | **0.04** | **-99.7%** | 564 | 7 | 67,459 $\to$ 67,452 | **VERIFIED** |

### Structural Metrics Continuity
- **Genome Fraction**: Maintained at $\ge 98.33\%$ on real isolates, and $99.99\%$ on synthetic datasets.
- **Misassemblies**: 0 extensive misassemblies on all synthetic datasets; 4 on R1-001, 2 on R1-003 (strictly identical to baseline; 0 new misassemblies).
- **Contiguity**: N50/NGA50 changes reflect solely the elimination of duplicate 1-bp insertion indels at verified $k$-overlap junctions (e.g., exactly -10 bp on SYN-2.5M across 10 $k$-overlap junctions).
- **Invariant Datasets**: `S1-COV-150X` and `S1-COV-003` have 0 exSPAnder repeat stitches; their FASTA outputs are bit-identical (identical SHA-256).

---

## 5. Deliverable Inventory

- [`C0_02_BASELINE_REPRO.md`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/C0_02_BASELINE_REPRO.md): Baseline defect reproduction, failure traces, and mathematical analysis.
- [`C0_02_JUNCTION_UNIT_TESTS.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/C0_02_JUNCTION_UNIT_TESTS.tsv): Full T1–T10 test suite results and behavioral verification.
- [`C0_02_GRAPH_CONSERVATION.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/C0_02_GRAPH_CONSERVATION.tsv): Graph topology conservation, active unitigs, and atomic commit safety.
- [`C0_02_ASSEMBLY_AB.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/C0_02_ASSEMBLY_AB.tsv): 6-dataset A/B benchmark metrics (N50, NGA50, indels, mismatches, RSS, wall time).
- [`C0_02_JUNCTION_DELTAS.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/C0_02_JUNCTION_DELTAS.tsv): Curated junction breakpoint events and sequence delta forensics.
- [`C0_02_VERDICT.md`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/C0_02_VERDICT.md): Formal verdict and claims matrix.

All deliverables are present in the workspace root.  
Branch `fix/c0-02-junction-overlap` is ready for commit and push.
