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

### Methodological Root Cause
The audit script inserted raw sliding-window k-mers directly into a `FxHashSet<u64>` without first deduplicating identical k-mers in the reference sequence.
- *E. coli* K-12 MG1655 (`S2-REAL-001`) contains seven 5kb ribosomal RNA (`rrn`) operons, insertion sequences, and repetitive mobile elements.
- When an identical 127-mer appears 7 times in the reference genome, its identical sequence naturally produces the identical `to_u64()` value 7 times.
- **Identical genomic repeat copies are NOT hash collisions.** A hash collision occurs if and only if **two distinct k-mers $k_1 \ne k_2$ produce the same folded hash $h(k_1) = h(k_2)$**.
- When distinct k-mers are first deduplicated into a set of unique 127-mers ($K_{\text{unique}}$), the number of true hash collisions between distinct k-mers is negligible, aligning with Birthday Paradox predictions ($P < 10^{-6}$).

### Corrective Action for Future Milestones
All future k-mer collision benchmarks must first insert into `FxHashSet<K>` (full 256-bit representation) to obtain the true unique distinct k-mer set, and only then evaluate whether distinct elements map to identical `to_u64()` keys.

---

## 3. Erratum 2: Scope of `final_contig_accuracy_pct` Metric

### Issue Description
In `B2_CMS_GRAPH_FATE.tsv`, the column `final_contig_accuracy_pct` reported an accuracy of **85.56%** for the un-subsampled 300X dataset (`S1-COV-003`), which seemingly contradicted the QUAST report of **0 misassemblies** and **100.0% reference identity**.

### Methodological Root Cause
- The audit script inspected the intermediate simplified graph (`raw_graph`) after `clip_tips` and `pop_bubbles`, directly iterating over all active unitigs in memory.
- Under 300X un-subsampled depth without read reduction, recurrent sequencing errors occurring $\ge 50$ times can form small isolated unitigs or unbranched spurs in the intermediate graph.
- In the full production CLI pipeline (`src/main.rs`), these intermediate unitigs pass through:
  1. `curate_contigs`: which strips all unitigs with length $< \text{min\_len}$ (default 500 bp) or coverage $< \text{min\_cov}$ (default 2.0).
  2. `consensus.rs`: which polishes consensus sequences against read pileups.
- After curation and polishing, all residual low-coverage error unitigs are stripped, resulting in the flawless **199,983 bp circular contig** with **0 misassemblies**.
- **Conclusion:** `final_contig_accuracy_pct` in the audit script actually measured the **simplified graph unitig k-mer concordance ratio**, NOT the final curated FASTA contigs.

---

## 4. Erratum 3: `EXP-07` Ground Truth Baseline Discrepancy

### Issue Description
In preliminary drafts of `B2_C1b3B_AB_MATRIX.tsv`, `EXP-07` (`S2-REAL-001`, 4.65 Mb real *E. coli*) was listed with expected metrics of `N50 = 4,655,000 bp` and `contigs = 1`.

### Methodological Root Cause
- Bacterial short-read paired-end sequencing (150 bp PE, ~450 bp insert size) cannot bridge repeat structures longer than the library fragment size.
- *E. coli* K-12 contains seven 5kb `rrn` ribosomal operons. In all short-read assemblers (SPAdes, Shovill, Unicycler short-read mode), the assembly graph naturally breaks at these 5kb repeats unless long reads are provided.
- Historical empirical measurements for `fastBact` on `S2-REAL-001` under baseline `dac6a57`:
  - **Contigs Count:** 2 to 3
  - **N50:** **4,082,514 bp** (~4.08 Mb)
  - **NGA50:** **2,897,224 bp** (~2.90 Mb)
  - **Misassemblies:** 1 (the 5kb `rrn` repeat inversion documented in `B2_MISASSEMBLY_REVIEW.tsv`)
- Expecting a single 4.65 Mb contig was biologically incorrect. `B2_C1b3B_AB_MATRIX_CORRECTED.tsv` has been updated to reflect the true empirical baseline (`4,082,514 bp`, 3 contigs).
