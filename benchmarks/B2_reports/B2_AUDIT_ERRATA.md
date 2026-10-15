# FB-BENCH-B2-C1b-3B: Formal Audit Errata & Ground Truth Clarification

- **Date:** 2026-10-15
- **Branch:** `fix/b2-c1b3b-var1`
- **Baseline Git SHA:** `dac6a57417b390cf7385ff07adf569e550cc8ca9`
- **Execution Target:** Slurm HPC Cluster Node `vnode1`

---

## 1. Executive Summary

During the `B2-C1b-3A-R` audit, three technical nuances were uncovered in the diagnostic scripts and preliminary reports. While these issues do not block the implementation of `Var 1` (CLI Subsampling Decoupling), they are formally documented here as **Audit Errata** to ensure absolute mathematical and biological rigor before proceeding to `Var 2` and `Var 3`.

---

## 2. Erratum 1: `to_u64()` Folding Value Duplication vs. Distinct K-mer Hash Collision

### Issue Description
In `src/bin/audit_depth_trace.rs`, the folding collision loop for `S2-REAL-001` was implemented as:
```rust
for i in 0..=(seq.len() - 127) {
    let can = ...;
    if !seen_127.insert(can.to_u64()) {
        col_127 += 1;
    }
}
```
This loop reported **38,937 folding collisions** on `S2-REAL-001` (4.65 Mb *E. coli*).

### Methodological Root Cause & Verification Requirement
The audit script inserted raw sliding-window k-mers directly into a `FxHashSet<u64>` without first deduplicating identical k-mers in the reference sequence.
- *E. coli* K-12 MG1655 (`S2-REAL-001`) contains seven 5kb ribosomal RNA (`rrn`) operons, insertion sequences, and repetitive mobile elements.
- When an identical 127-mer appears 7 times in the reference genome, its identical sequence naturally produces the identical `to_u64()` value 7 times.
- **Identical genomic repeat copies are NOT hash collisions.** A hash collision occurs if and only if **two distinct k-mers $k_1 \ne k_2$ produce the same folded hash $h(k_1) = h(k_2)$**.
- While repeats explain multiple observations of the same sequence, evaluating true 127-mer folded hash collisions requires an independent test specifically comparing distinct 256-bit canonical k-mers against their 64-bit folded projections.

### Corrective Action for Future Milestones
All future k-mer collision benchmarks must first insert into `FxHashSet<K>` (full 256-bit representation) to obtain the true unique distinct k-mer set, and only then evaluate whether distinct elements map to identical `to_u64()` keys.

---

## 3. Erratum 2: Scope of `final_contig_accuracy_pct` Metric & QUAST Interpretation

### Issue Description
In `B2_CMS_GRAPH_FATE.tsv`, the column `final_contig_accuracy_pct` reported an accuracy of **85.56%** for the un-subsampled 300X dataset (`S1-COV-003`), while QUAST reported **0 misassemblies**.

### Methodological Root Cause & Interpretive Boundary
- **Intermediate Unitig Concordance vs. Final FASTA**: The audit script inspected the intermediate simplified graph (`raw_graph`) after `clip_tips` and `pop_bubbles`, directly iterating over all active unitigs in memory.
  - In the full production CLI pipeline (`src/main.rs`), these intermediate unitigs pass through `curate_contigs` (stripping unitigs with length $< 500$ bp or coverage $< 2.0$) and `consensus.rs` (polishing against read pileups).
  - Residual low-coverage error unitigs present in the intermediate graph do not appear in the final curated FASTA.
- **QUAST Misassembly $\ne$ 100% Nucleotide Identity**:
  - A QUAST report of 0 extensive misassemblies indicates the absence of large-scale structural rearrangements (>1 kb inversions, relocations, or translocations).
  - 0 misassemblies does **NOT** equal 100% nucleotide identity; small indels, single-nucleotide mismatches, or local unaligned flanking regions may still exist. Accuracy claims must differentiate structural correctness from per-base sequence identity.

---

## 4. Erratum 3: `S2-REAL-001` Ground Truth Baseline Clarification

### Issue Description
In preliminary drafts, `S2-REAL-001` (4.65 Mb real *E. coli*) was inconsistently referenced against theoretical single-contig expectations or unverified experimental numbers (~4.08 Mb).

### Methodological Root Cause & Official Measured Baseline
- Bacterial short-read paired-end sequencing (150 bp PE, ~450 bp insert size) cannot bridge repeat structures longer than the library fragment size without long-read scaffolding.
- *E. coli* K-12 contains seven 5kb `rrn` ribosomal operons. In all short-read de novo assemblers (SPAdes, Shovill, Unicycler short-read mode), the assembly graph breaks at these repeats.
- **Official Archived Measured Baseline for `S2-REAL-001` (v0.1.0-baseline / `dac6a57` / C1b-1)**:
  - **Contigs (≥500 bp):** **28**
  - **N50:** **285,687 bp**
  - **NGA50:** **285,687 bp**
  - **Misassemblies:** 1 (5kb `rrn` repeat boundary)
  - **Genome Fraction:** 99.18%
- **Clarification**: Unverified figures such as 4.08 Mb must **NOT** be used as `baseline_measured`. The official archived baseline of **28 contigs, N50 = 285,687 bp** is the sole verified reference baseline for `S2-REAL-001`.
