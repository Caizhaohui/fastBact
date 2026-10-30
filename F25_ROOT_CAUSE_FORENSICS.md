# fastBact Phase F25 — Residual Misassembly Forensics Report

**Document ID**: `F25_ROOT_CAUSE_FORENSICS`  
**Phase**: F25 — Assembly Correctness & Residual Misassembly Forensics  
**Date**: 2026-10-30  
**Controller**: Antigravity CLI / Gemini 3.8 Flash Medium  
**Independent Reviewer**: Codex CLI / GPT-6.1 Sol Medium (`AUD-F25-03`, Sol Round 1 Review)  
**Policy File**: `~/.config/ai-collaboration/POLICY.md` (Multi-AI Collaboration V3.2)  
**Target Dataset**: R1-001 (*Escherichia coli* str. K-12 substr. MG1655, BioProject PRJNA1088182, Run SRR28349173)  
**Reference Genome**: `NC_000913.3.fasta` (4,641,652 bp, SHA-256 `6b195feda4c66140f6762742eb8b30c2652f02b45878b174f5b00ef85ecc95d7`)  
**Investigated Artifact**: F24 Accepted Baseline Assembly (`contigs.fa` from Slurm Job `3066737` / `3072360`)  

---

## 1. Forensic Summary & Event Classification

All 5 extensive misassemblies reported by QUAST v5.2.0 in the accepted F24 baseline for dataset R1-001 have been independently analyzed from raw FASTQ reads (`SRR28349173_1.fastq.gz` and `_2.fastq.gz`), through intermediate graph stages (Stage 1 cDBG, Stage 4 Multi-K, Stage 5 exSPAnder, Stage 6 Polish), to final contigs.

| Event # | QUAST Contig ID | Event Type | Reference Breakpoints | Inconsistency | Raw Read 70-bp Probe Support | Provenance Origin Stage | Biological & Algorithmic Root Cause | Formal Independent Classification |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---|:---:|
| **1** | `contig_53` | Inversion | 1,207,805 / 1,207,806 | N/A | **29 reads** (INV) vs 18 (WT) | Stage 4 (`unitig_4190`) | e14 prophage Pin recombinase (`pinE`) inversion boundary (`pixL` site); biological mixed population | **`EVALUATION_ARTIFACT`** |
| **2** | `contig_53` | Inversion | 1,209,602 / 1,209,603 | N/A | **37 reads** (INV) vs 14 (WT) | Stage 4 (`unitig_4190`) | e14 prophage Pin recombinase (`pinE`) inversion boundary (`pixR` site); biological mixed population | **`EVALUATION_ARTIFACT`** |
| **3** | `contig_63` | Relocation | 889,647 / 274,524 | 614,497 bp | **73 reads** | Stage 4 (`unitig_5348`) | Multi-copy IS5 transposon insertion; QUAST arbitrarily maps 626 bp repeat to 274 kb (out of 10 identical loci) | **`EVALUATION_ARTIFACT`** |
| **4** | `contig_200` | Relocation | 1,873,031 / 291,276 | 1,581,629 bp | **66 reads** | Stage 4 (`unitig_1974`) | Multi-copy IS1 transposon terminal repeat (126 bp = $k-1$); QUAST arbitrarily maps to 291 kb | **`EVALUATION_ARTIFACT`** |
| **5** | `contig_222` | Relocation | 1,872,615 / 1,978,503 | 105,463 bp | **51 reads** | Stage 4 (`unitig_3500`) | Multi-copy IS1 transposon terminal repeat (126 bp = $k-1$); QUAST arbitrarily maps to 1.98 Mb | **`EVALUATION_ARTIFACT`** |

### Verified Forensic Facts:
1. **Zero Chimeras Among Investigated Events**: None of the 5 reported extensive misassemblies represents an erroneous joining of distant genomic loci by repeat resolution. The historical AUD-01 2.1 Mb chimera remains 100% eliminated.
2. **Biological Heterogeneity of Inversion (Events 1 & 2)**: The inversion in `contig_53` is an authentic genomic rearrangement mediated by the Pin site-specific recombinase in the e14 prophage element. Independent raw read scans by Codex Sol confirm both inverted (29 and 37 reads) and wild-type reference (18 and 14 reads) junctions in the sequencing library, establishing a biological mixed population in the sequenced culture. The assembled sequence correctly represents the major inverted isoform.
3. **Multi-Mapping Repeat Ambiguity (Events 3, 4, 5)**: The reported "relocations" are evaluation artifacts caused by QUAST's alignment heuristic when multi-copy mobile genetic elements (IS5 and IS1) are present at contig boundaries. Because these elements exist in multiple identical copies across the chromosome (10 copies for IS5, 6-7 copies for IS1), QUAST's aligner assigns the terminal repeat segment to an arbitrary distant copy, misinterpreting an authentic repeat entry as an extensive relocation.

---

## 2. Event-by-Event Forensic Evidence

### 2.1 Events 1 & 2: Pin Recombinase Inversion (`contig_53`)

- **Contig Geometry**: Length 30,191 bp.
  - Prefix: 1..15,240 (15,240 bp) $\to$ Ref `NC_000913.3` 1,192,566..1,207,805 (Forward orientation).
  - Core: 15,241..17,037 (1,797 bp) $\to$ Ref `NC_000913.3` 1,207,806..1,209,602 (Inverted orientation).
  - Suffix: 17,038..30,191 (13,154 bp) $\to$ Ref `NC_000913.3` 1,209,603..1,222,756 (Forward orientation).
- **Genomic Locus & Biology**:
  - The locus corresponds to the *E. coli* K-12 e14 prophage tail fiber assembly operon (`tfaP`, `tfaE`, `stfE`), immediately adjacent to the `pinE` site-specific recombinase gene (`NC_000913.3:1209684..1210239`).
  - Pin invertase catalyzes reversible DNA inversion of the 1,797 bp P-segment between inverted repeat recognition sites (`pixL` and `pixR`).
- **Assembly Pipeline Tracking**:
  - Stage 1: `unitig_21534` (163 bp) contains the junction matching coordinate 16,891 of `contig_53`.
  - Stage 4 (`stage_04_multik_127_simplified.fa`): Assembled into `unitig_4190` (1,829 bp), representing the complete inverted core unitig.
  - Stage 5 (`stage_05_exspander.fa`): Extended into `unitig_1784` (30,191 bp) by connecting through verified unique corridors.
  - Raw Read Support: Codex Sol independent verification established 29 reads for the left inverted junction and 37 reads for the right inverted junction, confirming physical biological existence in the sequencing library.
- **Verdict**: **`EVALUATION_ARTIFACT`**. The assembly is structurally true to the sequenced biological genome; the divergence from `NC_000913.3` is due to static reference phase selection over an actively inverting bacterial locus.

---

### 2.2 Event 3: Transposon IS5 Multi-Copy Mapping (`contig_63`)

- **Contig Geometry**: Length 26,211 bp.
  - Main Body: 1..25,585 (25,585 bp) $\to$ Ref `NC_000913.3` 915,231..889,647 (in gene `ybjL`).
  - Terminal Segment: 25,586..26,211 (626 bp) $\to$ Ref `NC_000913.3` 274,524..275,149.
- **Repeat Context**:
  - The 626 bp terminal sequence is 100.0% identical to the insertion sequence IS5 (`insH1` transposase).
  - IS5 is present in **10 identical copies** across the *E. coli* K-12 genome:
    - 274,524 bp (reverse)
    - 575,160 bp (reverse)
    - 688,420 bp (reverse)
    - 1,299,499 bp (forward)
    - 1,396,044 bp (forward)
    - 2,102,318 bp (reverse)
    - 2,289,488 bp (reverse)
    - 3,130,146 bp (forward)
    - 3,366,125 bp (reverse)
    - 3,652,605 bp (reverse)
- **Assembly Pipeline Tracking**:
  - The junction between the unique `ybjL` gene and the IS5 element is present in Stage 4 as `unitig_5348` (2,930 bp, covering `contig_63` coordinates 22,778–25,707).
  - Raw Read Support: Codex Sol verification confirmed **73 reads** in `SRR28349173` spanning the 70-bp junction probe.
  - exSPAnder Behavior: exSPAnder extended `unitig_5348` upstream across unique corridors (from coordinate 22,777 to 1) and also extended 504 bp further into the multi-copy IS5 sequence (from coordinate 25,707 to 26,211). Because this extension remains entirely within the homogeneous IS5 element, it represents repeat continuation rather than a chimeric joining to another unique locus.
- **Verdict**: **`EVALUATION_ARTIFACT`**. The 626 bp terminal segment is an authentic IS5 entry stub. QUAST's aligner arbitrarily selected the first reference copy at 274 kb, flagging an extensive relocation of 614,497 bp.

---

### 2.3 Events 4 & 5: Transposon IS1 Terminal Inverted Repeats (`contig_200` & `contig_222`)

- **Contig Geometries**:
  - `contig_200`: Length 1,285 bp. Flank: 1..1,159 (Ref `NC_000913.3` 1,874,189..1,873,031, `yeaK`). Terminal Stub: 1,160..1,285 (126 bp, Ref 291,276..291,401, `insA3`).
  - `contig_222`: Length 551 bp. Flank: 1..425 (Ref `NC_000913.3` 1,872,615..1,873,039, `dgcJ`). Terminal Stub: 426..551 (126 bp, Ref 1,978,503..1,978,628, `insA5`).
- **Biological Context & Exact $k-1$ Artifact**:
  - In the reference genome, `dgcJ` (1,872 kb) and `yeaK` (1,873 kb) flank coordinate 1,873,030 bp.
  - In isolate SRR28349173, an IS1 transposable element is inserted at coordinate ~1.873 Mb.
  - The assembler enters the IS1 element from both sides. Because $k=127$, the unitig entering the repeat carries exactly $k-1 = 126$ bp of the repeat sequence before terminating at the repeat junction.
  - The 126 bp sequence corresponds to the multi-copy IS1 terminal repeat (found at 6 distinct reference loci: 20 kb, 258 kb, 279 kb, 291 kb, 1.98 Mb, 3.58 Mb).
- **Assembly Pipeline Tracking**:
  - `contig_200` corresponds to `unitig_1974` (1,285 bp) in Stage 4.
  - `contig_222` corresponds to `unitig_3500` (551 bp) in Stage 4.
  - Both unitigs are **100% bit-identical between Stage 4 and Stage 5** (`stage_05_exspander.fa`). exSPAnder did not modify them.
  - Raw Read Support: Codex Sol verification confirmed **66 reads** spanning the `contig_200` junction and **51 reads** spanning the `contig_222` junction.
- **Verdict**: **`EVALUATION_ARTIFACT`**. Both events represent natural entry into a multi-copy IS1 transposition site. The 126 bp ($k-1$) repeat overlap aligns to arbitrary distant IS1 reference copies in QUAST, producing spurious relocation calls.

---

## 3. Conclusions and Correctness Contract

1. **Integrity of Repeat Resolution**: None of the 5 residual extensive misassembly events in R1-001 was caused by chimeric or erroneous join decisions in exSPAnder. The historical AUD-01 2.1 Mb chimera remains 100% eliminated.
2. **True Assembly Correctness**: fastBact's assembly correctly separates unique corridors from multi-copy repeats and preserves authentic genomic structures supported by high-depth raw read evidence.
3. **Strict Non-Regression**: Under the frozen F24 non-regression contract, fastBact preserves all verified contiguity gains while maintaining $\le 5$ extensive misassemblies on R1-001 (all 5 proven to be evaluation artifacts) and strictly 0 extensive misassemblies on R1-003.
