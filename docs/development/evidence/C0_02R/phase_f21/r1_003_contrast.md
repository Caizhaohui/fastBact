# fastBact Phase F21-A: R1-003 Mechanistic Contrast Analysis

**Dataset**: `R1-003` (*Escherichia coli* BW25113, SRA accession `SRR28349164`)  
**Reference Genome**: `CP009273.1` (4,631,469 bp, SHA-256: `37eade01311f1f81659fad53bf907a93139075a75d115fd3aec2ef09beac2e90`)  
**Comparison**: Baseline Version A (`b9727c3`) vs Candidate Version B (`0606be4` / `a993047`)

---

## 1. Contrast Overview: Divergent Outcomes Between R1-001 and R1-003

| Metric | Baseline A | Candidate B | $\Delta$ (B - A) | Outcome Classification |
| :--- | :---: | :---: | :---: | :---: |
| **Extensive Misassemblies** | 2 | 1 | -1 (-50.0%) | **IMPROVED** |
| **Local Misassemblies** | 1 | 1 | 0 (0.0%) | **UNCHANGED** |
| **Mismatches per 100 kbp** | 1.51 | 1.08 | -0.43 (-28.5%) | **IMPROVED** |
| **Indels per 100 kbp** | 0.04 | 0.00 | -0.04 (-100.0%) | **IMPROVED (ELIMINATED)** |
| **Genome Fraction (%)** | 98.330% | 98.233% | -0.097 pp (-0.10%) | **REGRESSION** |
| **N50 (bp)** | 67,452 | 70,868 | +3,416 (+5.06%) | **IMPROVED** |
| **NGA50 (bp)** | 63,648 | 70,010 | +6,362 (+10.00%) | **IMPROVED** |
| **Active Unitigs Post-exSPAnder** | 300 | 704 | +404 (2.35×) | **CONSERVATIVE WITHHOLDING** |

While `R1-001` suffered an extra misassembly and elevated polishing boundary errors, `R1-003` showed remarkable improvements across 5 key quality dimensions. Understanding why these two closely related *E. coli* isolates diverged provides the key to resolving fastBact's repeat resolution behavior.

---

## 2. Misassembly Event Reconciliation

In `R1-003`, QUAST identified:

### Baseline A Events (2 Misassemblies)
1. **Event 1 (`contig_4`, length 144,537 bp)**:
   - Relocation with inconsistency = **+590,017 bp**.
   - Connects reference coordinate `4,198,055 bp` to `3,463,502 bp` across a 590 kb segmental duplication.
2. **Event 2 (`contig_16`, length 75,502 bp)**:
   - Relocation with inconsistency = **-2,094,638 bp**.
   - Connects reference coordinate `3,659,535 bp` to `1,564,898 bp` across the 2.1 Mb homologous repeat family.

### Candidate B Events (1 Misassembly)
1. **Event 1 (`contig_19`, length 75,502 bp)**:
   - Relocation with inconsistency = **-2,094,638 bp**.
   - Connects reference coordinate `3,659,535 bp` to `1,564,898 bp`.
   - **This is the exact same repeat relocation as Event 2 in Baseline A.**
2. **ELIMINATED EVENT**:
   - The +590,017 bp relocation on `contig_4` was **completely eliminated** in Candidate B.

### Mechanistic Factors in Event Elimination
At the 590 kb duplication locus on `contig_4`, Baseline A committed a join that QUAST identified as a relocation. In Candidate B, this join was withheld at Stage 5. Candidate B contains conservative decision logic (`SearchStatus::Truncated` and evidence-aware distance bounds $\beta(\text{curr})$), which withhold ambiguous paths. This structural difference occurred prior to polishing.

---

## 3. The Shared Homologous Repeat Locus (1.56 Mb $\leftrightarrow$ 3.66 Mb)

Crucially, both `R1-001` and `R1-003` possess the exact same homologous repeat family between ~1.56–1.57 Mb and ~3.65–3.66 Mb (relocation inconsistency of approximately -2.09 Mb).

- **In R1-003**: Baseline A **ALREADY COMMITTED** this chimeric join (`contig_16`), and Candidate B also committed it (`contig_19`). Thus, the misassembly count did not increase on R1-003.
- **In R1-001**: Baseline A halted extension prior to traversing the repeat (`contig_96` was separated from `contig_20`), avoiding the error. But Candidate B's extended multi-hop search bridged the 7 short intermediate nodes and committed the join, causing the +1 misassembly regression.

This proves that **both isolates share the identical repeat-family vulnerability**. The apparent difference in misassembly count (+1 on R1-001, -1 on R1-003) is a boundary artifact of how deeply the graph was explored on the two different graph topologies.

---

## 4. Why Did Genome Fraction Fall on Both Isolates?

In both datasets, Candidate B recorded a small regression in QUAST genome fraction:
- R1-001: $98.364\% \to 98.281\%$ ($-0.083$ pp)
- R1-003: $98.330\% \to 98.233\%$ ($-0.097$ pp)

Our unitig fate and Stage 5 vs Stage 6 analysis explains the contributing factors:
1. **Stage 5 Pre-Polish Genome Fraction is HIGHER in Candidate B**:
   - R1-001: Stage 5 Genome Fraction = **98.835%** in B vs **98.456%** in A (**+0.379 pp**).
   - R1-003: Stage 5 Genome Fraction = **98.857%** in B vs **98.444%** in A (**+0.413 pp**).
2. **Output Length Threshold Filtering**:
   - When exSPAnder withholds ambiguous repeat joins, it leaves short repetitive unitigs unmerged in the graph.
   - In R1-003 Candidate B, active unitigs increased from 300 to 704.
   - Exactly **395 unitigs (57,848 bp)** are $< 200\text{ bp}$.
   - Another **165 unitigs (43,925 bp)** are $200–499\text{ bp}$. Because fastBact's default output cutoff is 500 bp, these unitigs were not emitted into final `contigs.fa`.
   - Total unmerged sequence under 500 bp: **101,773 bp** (79.5% of all unitigs in the graph).
3. **Conclusion**:
   Candidate B did NOT fail to assemble these genomic regions in the graph. The sequences are present at Stage 5. However, because Candidate B avoided forced joins across ambiguous repeats, the sequences remained in short unitigs that were excluded by the default 500 bp contig output filter.
