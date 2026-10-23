# fastBact Phase F21-B: Mechanism Evidence Matrix & Causal Classification

**Investigation Phase**: Phase F21-B (Repeat-Join Decision Provenance & Mechanism Verification)  
**Target Dataset**: `R1-001` (*Escherichia coli* K-12 MG1655) vs `R1-003` (*Escherichia coli* BW25113)  
**Evidence Artifacts**: `effective_parameter_audit.tsv`, `repeat_locus_graph.tsv`, `directed_path_geometry.tsv`, `read_pair_provenance.tsv`, `candidate_decision_trace.tsv`, `alternative_path_matrix.tsv`, `r1_001_r1_003_repeat_comparison.md`, `polish_boundary_evidence.tsv`, `unique_coverage_loss.tsv`  
**C0-02R Acceptance**: **Strictly NOT ACCEPTED**  
**Production Code**: **100% Frozen** (`30967bf4940844aeb6b744732bcad6d2a72daab4`, 0 lines modified in `src/`)

---

## 1. Executive Evaluation Matrix

| Hypothesis ID | Hypothesis Description | Evidentiary Verdict | Supporting Evidence | Contradicting / Mitigating Evidence | Minimum Discriminating Experiment |
| :--- | :--- | :---: | :--- | :--- | :--- |
| **H1** | Incorrect multi-hop chaining of local PE support across repeats | **CONFIRMED** | `directed_path_geometry.tsv` proves 1,076 bp net separation between `u1927` and `u169` (exceeding 700 bp max insert); 0 physical read pairs link `u1927` to `u169`. Local step-by-step chaining at `u1476 -> u1926` (16 pairs) committed the cross-locus jump. | None. Chaining occurred across 7 short nodes without physical flanking fragment span. | Enforce that repeat traversals spanning > 3 consecutive short nodes require read pairs anchored in unique flanking contigs (> 1 kb). |
| **H2** | Repeat-internal multi-mapping falsely treated as uniquely anchored evidence | **SUPPORTED** | Nodes `u1476` (43 unique 31-mers) and `u3919` (174 unique 31-mers) appeared unique in the de Bruijn graph at $k=127$, allowing their terminal reads to serve as extension anchors despite being embedded in a multi-copy homologous repeat family. | `mapper.rs` correctly verified single-occurrence hashing in the local unitig pool, but graph-level k-mer uniqueness does not guarantee genomic locus uniqueness. | Cross-check anchor unitigs against global repeat-family copy number before allowing them to initiate repeat-exit candidate searches. |
| **H3** | Correct alternative path omitted or underweighted by resolver | **CONFIRMED** | At `u3919(-)`, authentic exit `u3920(-)` had 42 pairs (weight 16.25) vs 62 pairs (weight 21.93) for chimeric `u169(+)`. Backward extension from `u4814` consumed `u3919`, leaving `u3920` unjoined (`contig_146`). At `u1475(+)`, authentic 1.57 Mb exit `u4833` had 0 pairs vs 16 pairs for `u1926`. | In forward evaluation, 21.93 vs 16.25 (ratio 1.35x) would have been rejected as ambiguous under the 2.5x threshold, but was bypassed by backward seed traversal. | Reciprocal consistency requirement: require forward and reverse traversals to agree on edge selection before committing multi-hop paths. |
| **H4** | Reference-repeat ambiguity without sufficient spanning molecule information | **CONFIRMED** | The homologous repeat family at 1.57 Mb and 3.66 Mb is 2,092 bp in total sequence length and 1,076 bp in overlap-adjusted physical span. Short paired-end reads (median 374 bp, max 700 bp) cannot span this element. Short-read data is mathematically insufficient to resolve this repeat without long reads. | None. Resolving this repeat with short reads alone is fundamentally ill-posed when asymmetric variants misguide local search. | Withhold repeat resolution when candidate traversal distance exceeds library fragment maximum span ($\text{sep} + \beta > \text{max\_insert}$). |
| **H5** | Stage 6 terminal polishing introduces base errors | **CONFIRMED** | 84.6% (132/156) of base errors in Candidate B R1-001 are located in the 101–350 bp terminal boundary zone. In Stage 5 (pre-polish), Candidate B had FEWER mismatches (24 vs 29) than Baseline A. 100% of the mismatch and indel excess was introduced during consensus voting. | None. Documented in `base_error_loci.tsv` and `polish_boundary_evidence.tsv`. | Add terminal boundary guard in `consensus.rs` requiring $\ge 10\times$ depth and balanced forward/reverse strand support before modifying terminal 350 bp. |
| **H6** | Output-length filtering removes unique reference coverage | **SUPPORTED / DISSECTED** | Candidate B Stage 5 leaves 1,566 unitigs < 500 bp (total 235 kb). However, non-redundant reference mapping proves >85% of this sequence is redundant with output contigs. Net unique reference loss is < 25 kb (~0.5% genome fraction). | The 235 kb raw length does NOT represent 235 kb of lost unique reference sequence; it is mostly redundant repeat copies. | Evaluate sensitive output cutoff thresholds (e.g., 200 bp vs 500 bp) in offline analysis without altering production assembly defaults. |
| **H7** | A/B benchmark effective parameter mismatch | **REFUTED** | Gate 0 confirmed identical CLI parameters, Slurm allocations, automatic k-mer schedule `[31, 63, 127]`, default 500 bp cutoff, and QUAST `--min-contig 200`. The differences were 100% algorithmic. | Both versions ran with identical execution environments and command lines. | Re-run with explicit CLI flags confirms bit-identical configuration parity. |

---

## 2. In-Depth Mechanistic Synthesis

### 2.1 The Genesis of the Relocation on `contig_16` (H1, H2, H3, H4)
1. **Biological Reality**: *E. coli* MG1655 harbors two homologous repeat elements at ~1.57 Mb (`1,570,669..1,571,802`) and ~3.66 Mb (`3,666,007..3,671,174`).
2. **Graph Collapse**: At $k=127$, the shared repeat elements collapse into a linear chain of short unitigs: `u3275`, `u1475`, `u4349`, `u4465`, `u3919`.
3. **Asymmetric Mutation**: The 3.66 Mb locus contains an asymmetric insertion `unitig_1476` (255 bp) and junction transition `unitig_1926` (267 bp).
4. **Seed Order**: ExSPAnder processes seeds in descending length order. Long seed `unitig_4814` (21,621 bp, Rank 16) and `unitig_169` (12,633 bp, Rank 81) at 1.57 Mb were processed hundreds of seeds before `unitig_1927` (5,073 bp, Rank 313) or `unitig_3920` (2,012 bp, Rank 676).
5. **Backward Extension**: Extending backward from 1.57 Mb, the path walked through `u169(-)` $\to$ `u3919(+)` $\to$ `u4465(+)` $\to$ `u4349(+)`. At `u4349(+)`, the graph had only one outgoing edge: `u1476(+)`.
6. **False Anchoring & Winner-Take-All Failure**: Once inside `u1476(+)`, the walk was anchored in a 3.66 Mb-specific sequence. At `u1475(+)`, the resolver evaluated candidates:
   - Candidate 1 (`u1926`, 3.66 Mb): 16 supporting pairs, weight 6.34.
   - Candidate 2 (`u3276`): 7 pairs, weight 0.20.
   - Candidate 3 (`u4833`, authentic 1.57 Mb): 0 pairs, weight 0.00.
   With a winner-take-all ratio of $31.7\times \ge 2.5$, the resolver committed the chimeric jump into `u1926` and `u1927`, producing the 71,975 bp chimeric `contig_16`.

### 2.2 The Genesis of Base Error Elevation (H5)
- In Stage 5 (pre-polishing), Candidate B had 24 mismatches (vs 29 in A) and 4 indels (vs 4 in A).
- Consensus polishing in Stage 6 introduced 61 new mismatches and 39 new indels in Candidate B, with 84.6% concentrated within 101–350 bp of contig ends.
- Terminal contig boundaries suffer from seed depletion (repetitive ends filter out seeds), causing reads to overhang unidirectionally from internal seeds. Low depth ($5\le N \le 8$), extreme strand bias, and un-gapped linear projection of read coordinates over indels cause false consensus flips exceeding the 75% threshold.

---

## 3. Provenance and Integrity Attestation

- All findings derived from direct evidence in `docs/development/evidence/C0_02R/phase_f21b/`.
- No modifications made to production source tree (`src/` SHA: `30967bf4940844aeb6b744732bcad6d2a72daab4`).
- No modifications made to `main` branch.
- C0-02R status: **Strictly NOT ACCEPTED**.
