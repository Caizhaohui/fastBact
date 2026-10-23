# fastBact Phase F21-B: Cross-Isolate Repeat Comparison (R1-001 vs R1-003)

**Dataset R1-001**: *Escherichia coli* K-12 MG1655 (`SRR28349173`, Reference: `NC_000913.3`, 4,641,652 bp)  
**Dataset R1-003**: *Escherichia coli* BW25113 (`SRR28349164`, Reference: `CP009273.1`, 4,631,469 bp)  
**Evaluator**: QUAST v5.2.0 (`--min-contig 200`)  
**Investigation Scope**: Gate 4 Cross-Isolate Contrast of Homologous Repeat Loci (1.57 Mb $\leftrightarrow$ 3.66 Mb)

---

## 1. Executive Summary & Core Finding

In Phase F20, Candidate Version B exhibited contrasting assembly outcomes:
- **R1-003**: Extensive misassemblies **decreased** from 2 to 1 (-50%), mismatches improved from 1.51 to 1.08/100kbp, and indels dropped from 0.04 to 0.00.
- **R1-001**: Extensive misassemblies **increased** from 4 to 5 (+25%), accompanied by an increase in polished terminal base errors.

Our cross-isolate reference-coordinate dissection uncovers the biological and topological reason for this divergence:

1. **The ~2.09 Mb Chimeric Relocation is NOT Unique to Candidate B; it Was Already Committed by Baseline A in R1-003**:
   - In R1-003 Baseline A, `contig_16` contains the **exact same chimeric relocation** (QUAST: inconsistency = -2,094,638 bp between coordinates 1..13,066 and 13,067..75,502).
   - In R1-003 Candidate B, `contig_19` reproduces this exact chimeric join identically.
   - However, Candidate B simultaneously corrected an independent 590,017 bp relocation on `contig_4`, reducing R1-003 misassemblies from **2 to 1**.
2. **Why Baseline A Did Not Make the Join in R1-001**:
   - In R1-001, the repeat locus contains an extra node `unitig_4465` (312 bp) in the simplified de Bruijn graph at $k=127$, requiring a 7-node chain.
   - Baseline A's shallow candidate search budget failed to traverse this 7-node repeat chain, halting extension at `unitig_3275` (emitting unjoined `contig_96`, `contig_142`, and `contig_20`).
3. **Why Candidate B Committed the Join in R1-001**:
   - Candidate B's search completeness improvements (multi-hop traversal, short prefix reuse, and unified `classify_hop` with evidence-aware $\beta$ bound) enabled exSPAnder to traverse through the 7-node repeat chain in R1-001, triggering the exact same chimeric joining mechanism that had already occurred in R1-003 Baseline A.

---

## 2. Cross-Isolate Genomic and Graph Comparison

| Feature | R1-001 (*E. coli* MG1655) | R1-003 (*E. coli* BW25113) | Alignment Status |
| :--- | :--- | :--- | :--- |
| **Reference Genome** | `NC_000913.3` (4,641,652 bp) | `CP009273.1` (4,631,469 bp) | Isolate-specific reference |
| **Left Flank Locus** | `3,661,061..3,666,134` (3.66 Mb) | `3,654,421..3,659,494` (3.65 Mb) | Synteny conserved |
| **Right Flank Locus** | `1,571,675..1,584,308` (1.57 Mb) | `1,565,932..1,578,565` (1.57 Mb) | Synteny conserved |
| **Repeat Family Copy Number** | 2 genomic copies | 2 genomic copies | Bit-identical copies |
| **Asymmetric Node `u1476` (255 bp)** | Present ONLY at 3.66 Mb | Present ONLY at 3.65 Mb | Conserved asymmetry |
| **Shared Repeat Unitigs** | `u3275`, `u1475`, `u4349`, `u4465`, `u3919` | `u3275`, `u1475`, `u4349`, `u3919` | R1-001 has extra node `u4465` |
| **Baseline A Assembly Fate** | **Unjoined** (`contig_96` + `contig_20`) | **Chimeric Join Committed** (`contig_16`) | Divergent in Baseline A |
| **Candidate B Assembly Fate** | **Chimeric Join Committed** (`contig_16`) | **Chimeric Join Committed** (`contig_19`) | Identical in Candidate B |
| **QUAST Misassembly Classification** | Relocation (-2,095,535 bp) | Relocation (-2,094,638 bp) | Identical genomic event |

---

## 3. Detailed Mechanism Analysis

### 3.1 Graph Topology Divergence in Baseline A
In R1-003, the de Bruijn graph at $k=127$ formed a simpler 6-node repeat transition (`u1926 -> u3275 -> u1475 -> u1476 -> u4349 -> u3919`), omitting `u4465`. Baseline A's search depth was just sufficient to traverse this 6-node path, committing `contig_16` with the -2.09 Mb relocation.

In R1-001, `u4465` added an extra 185 bp net separation, pushing the total intermediate path to 7 nodes. Baseline A's search halted before reaching `u4349`, leaving `contig_96` (at 3.66 Mb) and `contig_20` (at 1.57 Mb) unjoined. Thus, Baseline A avoided the misassembly in R1-001 purely by virtue of search incompleteness (accidental avoidance through truncation).

### 3.2 Candidate B Traversal and Chimeric Commitment
Candidate B's search completeness enhancements successfully traversed the 7-node repeat path in R1-001. However, because `unitig_1476` (255 bp) physically exists only at 3.66 Mb, backward extension from the 1.57 Mb body entered `unitig_1476` and found 16 supporting paired-end links pointing exclusively to 3.66 Mb `unitig_1926` (weight 6.34) and 0 links pointing to 1.57 Mb `unitig_4833`. The resolver committed this path, producing the same misassembly that Baseline A had committed in R1-003.

### 3.3 Summary of Isolate Outcome Differences
The apparent "regression" on R1-001 and "improvement" on R1-003 arose because:
1. Candidate B's repeat-traversal logic behaves **consistently and deterministically** across both datasets, discovering the same path topology in both.
2. In R1-003, Baseline A had already committed this join, so Candidate B's behavior did not alter the locus, while Candidate B's superior repeat resolution elsewhere fixed `contig_4` (reducing misassemblies 2 $\to$ 1).
3. In R1-001, Baseline A had previously halted prematurely due to search depth limits, so Candidate B's traversal appeared as a new extensive misassembly (+1).

---

## 4. Evidence Reproducibility

- **R1-001 QUAST Report**: `benchmark-results/C0_02R_f20_selective_20261022/quast/R1-001/B/contigs_reports/contigs_report_contigs.mis_contigs.info`
- **R1-003 QUAST Report**: `benchmark-results/C0_02R_f20_selective_20261022/quast/R1-003/B/contigs_reports/contigs_report_contigs.mis_contigs.info`
- **Reference Fasta MG1655**: `NC_000913.3.fasta` (SHA-256: `6b195feda4c66140f6762742eb8b30c2652f02b45878b174f5b00ef85ecc95d7`)
- **Reference Fasta BW25113**: `CP009273.1.fasta` (SHA-256: `37eade01311f1f81659fad53bf907a93139075a75d115fd3aec2ef09beac2e90`)
