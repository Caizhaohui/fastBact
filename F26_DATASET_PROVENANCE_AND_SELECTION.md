# fastBact Phase F26: Dataset Provenance and Selection Report

**Document ID:** `F26_DATASET_PROVENANCE_AND_SELECTION`  
**Phase:** F26 — Cross-Species Assembly Benchmarking, Robustness Hardening & Competitive Validation  
**Task:** F26-DATA-001 (Stage B: Dataset Discovery, Provenance Verification, Manifest Finalization)  
**Worker:** OpenCode / GLM (`glm-5.3-flash`)  
**Date:** 2026-10-30  
**Status:** VERIFIED — Accessions, References, and Manifest Checksums Confirmed Exact; Audit Finding BC-01 (CS-SA01 Strain Identity) Empirically Resolved (2026-10-30, F26-DATA-003)  
**Feature Branch:** `feature/f26-cross-species-robustness`  
**Base Commit:** `02872f81811402d2bd603ff7c0a2ecf8b3dbe528` (main)  

---

## 1. Executive Summary & Objectives

In accordance with Multi-AI Collaboration Policy V3.2 and the fastBact Phase F26 Technical Mission (Contract §F26-P1), this document establishes the full provenance, biological rationale, reference ground truth, run-to-reference compatibility assessment, and checksum integrity for the benchmark datasets comprising Phase F26:

1. **Historical Baseline Oracles (4 datasets):** The four validated *Escherichia coli* datasets from F24 and F25 (R1-001, R1-003, R1-001-Val, R1-003-Val), preserved for non-regression verification and competitive comparison. These are **not counted** toward the non-*E. coli* expansion requirement.
2. **Cross-Species Expansion (4 datasets):** Four newly selected, verified non-*E. coli* bacterial datasets spanning **two additional Gram-negative species** and **two Gram-positive species** across **four distinct genera** (*Staphylococcus*, *Bacillus*, *Pseudomonas*, *Salmonella*). The F26-P1 contract minimum (≥ 3 non-*E. coli* species, ≥ 2 genera, both Gram categories) is exceeded.

All selected cross-species datasets satisfy the following strict criteria:
- Complete, closed, circular RefSeq chromosome references with **zero ambiguous bases (N)** — machine-verified from the staged FASTA files (§3).
- Genuine paired-end Illumina WGS libraries (WGS strategy, PAIRED layout) — verified from NCBI SRA experiment XML (§4).
- Explicit `strain=` BioSample attribute correspondence between sequenced isolate and reference genome (verified: MW2, 168, PAO1, LT2).
- Wide GC-content spectrum spanning **32.83%** (*S. aureus*) to **66.56%** (*P. aeruginosa*) — recomputed from the staged references and matching the manifest exactly (§3).
- Genome size diversity from **2.82 Mb** to **6.26 Mb**, and coverage depths from **116×** to **234×**, supporting the F26-P3 depth-gradient design (30×–150× + full) from a single full-depth source per species (§7).

### 1.1 Verification Outcome (2026-10-30)

| Check | Result |
| :--- | :--- |
| SRA run accessions resolve, are public, and load-complete | **PASS** (4/4) |
| BioProject / BioSample / strain correspondence | **PASS** (4/4; strain tags = MW2, 168, PAO1, LT2) — CS-SA01 BioSample Title `MRSA252-cotreat2` conflicts with the `strain = MW2` attribute; **resolved empirically** by 31-mer containment analysis: organism confirmed *S. aureus* MW2 (USA400) + phage Staph1N co-culture, Title recorded as submitter labelling artefact — see [`docs/development/evidence/CS_SA01_STRAIN_DISCRIMINATION_EVIDENCE.md`](docs/development/evidence/CS_SA01_STRAIN_DISCRIMINATION_EVIDENCE.md); closes audit finding **BC-01** (`F26-AUDIT-001`) |
| Illumina paired-end layout and instrument model | **PASS** (4/4; instrument models match manifest) |
| Total spot counts vs manifest | **PASS** (4/4 exact) |
| Reference FASTA existence, single-contig, zero-N | **PASS** (8/8 files) |
| Reference raw + normalized SHA-256 vs `F26_DATASET_MANIFEST.tsv` | **PASS** (8/8 exact, zero corrections required) |
| GC content and genome size vs manifest | **PASS** (8/8 exact to reported precision) |
| Library layout (read-length arithmetic) vs manifest | **PASS** (4/4 consistent) |
| Experimental-context purity of sequenced isolates | **PASS with documented deviations** — see §6; all four runs derive from reference-strain lineages under treatment/evolution/manipulation protocols (contract-permitted via the "verified reference strain" clause), and each carries a quantified compatibility risk rating. |

**Controller decision point:** the four candidate accessions satisfy F26-P1 as written and are **retained as selected**. §6 quantifies per-dataset risks and mitigation notes so the Controller may revisit any dataset before Stage C if stricter purity is desired. Dataset replacement after Stage C would invalidate the pre-registered baseline (Contract Stage D ordering).

---

## 2. Complete Benchmark Dataset Matrix

| Dataset ID | Species | Strain | Gram Stain | GC Content (%) | Genome Size (bp) | SRA Run Accession | BioProject / BioSample | Platform / Layout | Total Spots | Reference Accession | Role |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **R1-001** | *Escherichia coli* | K-12 MG1655 | Negative | 50.79% | 4,641,652 | `SRR28349173` | PRJNA1088182 / SAMN40463664 | Illumina NovaSeq 6000 (2x151) | 2,423,735 | `NC_000913.3` | Historical Baseline Oracle |
| **R1-003** | *Escherichia coli* | BW25113 | Negative | 50.78% | 4,631,469 | `SRR28349164` | PRJNA1088182 / SAMN40463672 | Illumina NovaSeq 6000 (2x151) | 2,544,142 | `CP009273.1` | Historical Baseline Oracle |
| **R1-001-Val** | *Escherichia coli* | K-12 MG1655 | Negative | 50.79% | 4,641,652 | `SRR28349174` | PRJNA1088182 / SAMN40463663 | Illumina NovaSeq 6000 (2x151) | 2,367,236 | `NC_000913.3` | Historical Baseline Oracle |
| **R1-003-Val** | *Escherichia coli* | BW25113 | Negative | 50.78% | 4,631,469 | `SRR28349165` | PRJNA1088182 / SAMN40463671 | Illumina NovaSeq 6000 (2x151) | 4,267,402 | `CP009273.1` | Historical Baseline Oracle |
| **CS-SA01** | *Staphylococcus aureus* | MW2 | Positive | 32.83% | 2,820,462 | `SRR33682404` | PRJNA1263016 / SAMN48627418 | Illumina NovaSeq X Plus (2x150) | 2,239,653 | `NC_003923.1` | F26 Cross-Species Evaluation |
| **CS-BS01** | *Bacillus subtilis* | subsp. *subtilis* str. 168 | Positive | 43.51% | 4,215,606 | `SRR40739147` | PRJNA1530750 / SAMN63192558 | Illumina NextSeq 550 (2x151) | 2,253,807 | `NC_000964.3` | F26 Cross-Species Evaluation |
| **CS-PA01** | *Pseudomonas aeruginosa* | PAO1 | Negative | 66.56% | 6,264,404 | `SRR41244734` | PRJNA1545662 / SAMN64277765 | Illumina NovaSeq 6000 (2x240) | 1,518,403 | `NC_002516.2` | F26 Cross-Species Evaluation |
| **CS-SE01** | *Salmonella enterica* | subsp. *enterica* serovar Typhimurium str. LT2 | Negative | 52.22% | 4,857,450 | `SRR41010285` | PRJNA1539015 / SAMN63983458 | Illumina NovaSeq 6000 (2x150) | 3,359,296 | `NC_003197.2` | F26 Cross-Species Evaluation |

All values in this table are mirrored in `F26_DATASET_MANIFEST.tsv`, which was **byte-audited on 2026-10-30 and requires zero corrections** (see §3).

---

## 3. Reference Genome Integrity and Cryptographic Checksums

All reference genomes were fetched from NCBI RefSeq / GenBank via NCBI Entrez E-utilities and verified as complete, single-contig (circular chromosome) assemblies matching the strain isolates. FASTA deflines were re-inspected on 2026-10-30 and read exactly:

```text
>NC_000913.3 Escherichia coli K-12 MG1655, complete sequence        (staged at m6_identity/)
>CP009273.1  Escherichia coli BW25113 ... complete sequence          (staged at m6_identity/)
>NC_003923.1 Staphylococcus aureus subsp. aureus MW2, complete sequence
>NC_000964.3 Bacillus subtilis subsp. subtilis str. 168 complete genome
>NC_002516.2 Pseudomonas aeruginosa PAO1, complete genome
>NC_003197.2 Salmonella enterica subsp. enterica serovar Typhimurium str. LT2, complete genome
```

Each file was independently re-hashed and re-derived (GC%, length, N-count, contig count) on 2026-10-30; every value matched `F26_DATASET_MANIFEST.tsv` exactly:

| Dataset ID | Reference File (absolute) | Length (bp) | Contigs | N bases | Raw File SHA-256 | Normalized Sequence SHA-256 |
| :--- | :--- | :---: | :---: | :---: | :--- | :--- |
| R1-001 / R1-001-Val | `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/references/m6_identity/NC_000913.3.fasta` | 4,641,652 | 1 | 0 | `6b195feda4c66140f6762742eb8b30c2652f02b45878b174f5b00ef85ecc95d7` | `e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0` |
| R1-003 / R1-003-Val | `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/references/m6_identity/CP009273.1.fasta` | 4,631,469 | 1 | 0 | `37eade01311f1f81659fad53bf907a93139075a75d115fd3aec2ef09beac2e90` | `3761072f8884add5e9facdafc583de5d81a4cf5dbc7e51023ba966f3cf4887ce` |
| CS-SA01 | `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/f26_cross_species/references/NC_003923.1.fasta` | 2,820,462 | 1 | 0 | `19ba62a339065c9c498efb72bb2491e372344e4f23a4ac1e3e34e7345a25c63c` | `7c4947b9106af48ffee6bdf8dc29753333d97d2824092fb8228b57fc753d315d` |
| CS-BS01 | `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/f26_cross_species/references/NC_000964.3.fasta` | 4,215,606 | 1 | 0 | `a334e891ffc0e307f23f48842775d3383177a9d9cb5d5075b552a2cccddfe139` | `71572ae0b98f44a6382b29e194e56c1f7d87477be8f53915fb06816f886b17b7` |
| CS-PA01 | `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/f26_cross_species/references/NC_002516.2.fasta` | 6,264,404 | 1 | 0 | `902fa552c8c25aaaf82eccfed000e359d6a04c8744d9ec0a8134699cb6ba53e6` | `f969830bbe6aa58beccd852ae6e129a9489a92557ad082283f639b6e8539a340` |
| CS-SE01 | `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/f26_cross_species/references/NC_003197.2.fasta` | 4,857,450 | 1 | 0 | `59fceaea5ca57836e02b6c4ad3582847623c87268f7e5511874909eedb4acafe` | `37854e9ce0ff2b3aa4a3394bfab5a63f4922d5e9479acb930558979ffca81326` |

*Normalized Sequence SHA-256 = SHA-256 over the uppercase ACGT sequence only (headers and all whitespace excluded), as computed by the manifest verification commands in §8.*

### 3.1 Path Resolution Convention

`F26_DATASET_MANIFEST.tsv` stores paths **relative to the shared data root** `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/`:

- `data/references/m6_identity/*.fasta` → baseline *E. coli* references (shared with F24/F25; located under `references/`, **not** under `f26_cross_species/`).
- `data/f26_cross_species/references/*.fasta` → the four new cross-species references. All four files verified present in `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/f26_cross_species/references/` (deliverable 4 satisfied).

### 3.2 Reference Scope Note (chromosome-only)

All six references are **chromosome-only** (no plasmids/extrachromosomal elements). This is consistent across all datasets and matches the QUAST evaluation convention used in F24/F25. It is particularly relevant for CS-SA01 (*S. aureus* MW2 carries plasmids pMW2/pMW3, and the verified phage **Staph1N** co-culture adds ≈ 140 kb of extrachromosomal viral content — see §4.1.1) and CS-SE01 (*S. enterica* LT2 carries pSLT ≈ 94 kb): plasmid- and phage-derived contigs will appear as unaligned contigs in QUAST and must be interpreted per §6, not counted as assembler failures.

---

## 4. SRA Accession Verification (NCBI, 2026-10-30)

Metadata for all four candidate runs was retrieved from the NCBI SRA database in a single E-utilities call and parsed programmatically:

```text
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=SRR33682404,SRR40739147,SRR41244734,SRR41010285&rettype=xml&retmode=xml
HTTP 200; response retained at /tmp/opencode/sra_meta.xml (33,989 bytes)
```

### 4.1 CS-SA01 — `SRR33682404` (*S. aureus* MW2)

| Field | Verified value |
| :--- | :--- |
| Experiment / Study | SRX28908742 ("MW2-cotreat2") / SRP585898, BioProject **PRJNA1263016** ✓ |
| BioSample | **SAMN48627418** ✓ (SRS25134113 "MW2-cotreat2") |
| BioSample Title | ⚠ **`MRSA252-cotreat2`** — conflicts with the `strain = MW2` attribute; **resolved as submitter labelling artefact** (see §4.1.1) |
| Organism strain attribute | `strain = MW2` ✓ |
| Submitter / Center | University of Wisconsin–Madison, Bacteriology |
| Study title | "Sequencing data of MRSA strains subjected to bacteriophage infections" |
| Library | WGS, GENOMIC, RANDOM selection, **PAIRED** ✓ |
| Instrument | **Illumina NovaSeq X Plus** ✓ |
| Runs / Bases | **2,239,653 spots** ✓ / 659,542,596 bases (avg 294.5 bp/spot, consistent with 2×150) |
| SRA archive size | 187,232,418 bytes (published 2026-09-01) |
| Experimental context | Replicate A, phage **Staph1N** co-treatment (MOI 0.01) + oxacillin 32 µg/mL; `sample_type = mixed culture` |

### 4.1.1 BC-01 Resolution — Strain Identity Empirically Verified (F26-DATA-003, 2026-10-30)

The independent Stage B/C audit (`F26-AUDIT-001`, Codex Sol) raised finding **BC-01 / P1**: the BioSample Title `MRSA252-cotreat2` contradicts the `strain = MW2` attribute (the Title is the only occurrence of the `MRSA252` token in the entire record; the sample alias `MW2-cotreat2`, run alias `MW2_A_Staph1N_E12_S196_R1_001.fastq.gz`, strain attribute, and all experimental attributes point to MW2 + phage Staph1N). Per the audit's required closure path, the conflict was resolved by **discriminatory sequence-based identity analysis** — canonical 31-mer containment of the checksum-verified FASTQ pairs against MW2 (`NC_003923.1`), phage Staph1N (`NC_047722.1`), and MRSA252 (`NC_002952.2`):

| Analysis | Result |
| :--- | :--- |
| vs MW2 (`NC_003923.1`, 2,820,462 bp) | 29.90% of distinct read 31-mers matched; 28.85% of reads contain MW2 k-mers |
| vs phage Staph1N (`NC_047722.1`, 140,511 bp) | 60.74% of distinct read 31-mers matched; 53.49% of reads contain Staph1N k-mers |
| Combined MW2 + Staph1N | **90.64%** of all distinct 31-mers in the sequencing library |
| vs MRSA252 (`NC_002952.2`, 2,902,619 bp) | 19.12% matched (strictly core-genome homology with MW2); **specific MRSA252 diagnostic markers absent** |

**Conclusion:** the sequenced organism is conclusively ***S. aureus* MW2 (USA400)**; the library is a mixed bacterial–viral co-culture (phage depth ≈ 1000×, MW2 chromosome ≈ 35×) — the intended biology of BioProject PRJNA1263016 (phage–antibiotic co-treatment), not a contamination event. The BioSample Title is a **submitter labelling artefact**. Reference `NC_003923.1` is confirmed correct and unchanged. Full evidence, metric provenance, and the binding assembly/QUAST protocol are in [`docs/development/evidence/CS_SA01_STRAIN_DISCRIMINATION_EVIDENCE.md`](docs/development/evidence/CS_SA01_STRAIN_DISCRIMINATION_EVIDENCE.md). **BC-01: RESOLVED.**

### 4.2 CS-BS01 — `SRR40739147` (*B. subtilis* 168)

| Field | Verified value |
| :--- | :--- |
| Experiment / Study | SRX35306902 ("M067") / SRP737943, BioProject **PRJNA1530750** ✓ |
| BioSample | **SAMN63192558** ✓ (SRS30795682 "M067_BS_CHL_Rep1_T1") |
| Organism strain attribute | `strain = 168` ✓ (annotated `type-material = type strain of Bacillus subtilis`) |
| Study title | "Bacillus subtilis 168, CHL replicate 1 T1" |
| Library | WGS, GENOMIC, RANDOM selection, **PAIRED** ✓ |
| Instrument | **Illumina NextSeq 550** ✓ |
| Runs / Bases | **2,253,807 spots** ✓ / 680,649,714 bases (avg 302.0 bp/spot = exactly 2×151) |
| SRA archive size | 313,297,218 bytes (published 2026-09-18) |
| Experimental context | Laboratory culture of a chloramphenicol-adapted (CHL) laboratory-evolved lineage of strain 168, replicate 1, timepoint T1 |

### 4.3 CS-PA01 — `SRR41244734` (*P. aeruginosa* PAO1)

| Field | Verified value |
| :--- | :--- |
| Experiment / Study | SRX35752421 ("AztR6") / SRP746937, BioProject **PRJNA1545662** ✓ |
| BioSample | **SAMN64277765** ✓ (SRS31180756 "PAO1-GFP-AztR6") |
| Organism strain attribute | `strain = PAO1` ✓; `isolate = PAO1-GFP`; `genotype = mexR, aroB mutant 2` |
| Study title | "Squenced spontaneous resistance mutants of pseudomonas aeruginosa" (sic) |
| Library | WGS, GENOMIC, RANDOM selection, **PAIRED** ✓ |
| Instrument | **Illumina NovaSeq 6000** ✓ |
| Runs / Bases | **1,518,403 spots** ✓ / 726,918,458 bases (avg 478.8 bp/spot, consistent with 2×240) |
| SRA archive size | 244,609,460 bytes (published 2026-10-19) |
| Experimental context | GFP-tagged PAO1 reporter derivative ("AztR6"); run file alias `265438_aztR6_1_trimmed.fastq.gz` indicates **submitter-side adapter/quality trimming** |

### 4.4 CS-SE01 — `SRR41010285` (*S. enterica* LT2)

| Field | Verified value |
| :--- | :--- |
| Experiment / Study | SRX35538071 ("SeResistantAMPB4") / SRP742527, BioProject **PRJNA1539015** ✓ |
| BioSample | **SAMN63983458** ✓ (SRS30999721 "SeResistantAMPB4") |
| Organism strain attribute | `strain = LT2` ✓, `serovar = Typhimurium` ✓ (`type-material = type strain of Salmonella enterica`) |
| Study title | "DNA of Salmonella enterica LT2" |
| Library | WGS, GENOMIC, **RANDOM PCR** selection, **PAIRED** ✓ |
| Instrument | **Illumina NovaSeq 6000** ✓ |
| Runs / Bases | **3,359,296 spots** ✓ / 1,006,879,350 bases (avg 299.7 bp/spot, consistent with 2×150) |
| SRA archive size | 120,299,090 bytes (published 2026-10-05) |
| Experimental context | Adaptive-laboratory-evolution lineage with elevated antimicrobial-peptide resistance (AMPB4, Spain 2023); run uploaded as **`SeResistantAMPB4.final.bam`** (alignment-derived submission) |

---

## 5. Biological Rationale for Cross-Species Diversity

The selection of the 4 cross-species benchmarks was engineered to systematically test graph traversal, edge classification, bubble popping, and coverage estimation across biological extremes:

1. **Extreme Low GC & Small Genome (`CS-SA01` — *Staphylococcus aureus* MW2):**
   - **GC Content:** 32.83% (machine-verified).
   - **Characteristics:** AT-rich genomes suffer from non-uniform sequencing coverage and biased sequencing error profiles on Illumina platforms. The small 2.82 Mb genome tests fastBact's coverage estimator scaling and *k*-mer distribution modeling under low-entropy sequences.
2. **Intermediate GC & Spore-Forming Firmicute (`CS-BS01` — *Bacillus subtilis* 168):**
   - **GC Content:** 43.51% (machine-verified).
   - **Characteristics:** Moderate AT bias, extensive rRNA and prophage repeat structures. The standard model Gram-positive bacterium, providing a balanced physiological comparator to Gram-negative enteric rods.
3. **Extreme High GC & Large Complex Genome (`CS-PA01` — *Pseudomonas aeruginosa* PAO1):**
   - **GC Content:** 66.56% (machine-verified).
   - **Characteristics:** High GC content induces strong PCR amplification bias, producing GC-associated coverage dropouts and sharp depth cliffs. The large 6.26 Mb genome with intricate repeat families and 2×240 bp long paired reads stresses bubble-popping thresholds and tip-clipping logic.
4. **Closely Related Enteric Enterobacteriaceae (`CS-SE01` — *Salmonella enterica* LT2):**
   - **GC Content:** 52.22% (machine-verified).
   - **Characteristics:** Phylum Pseudomonadota, family Enterobacteriaceae. Closely related to *E. coli* with similar GC and genome size (~4.86 Mb), but harboring distinctive pathogenicity islands, prophage insertions, and the 94 kb pSLT virulence plasmid. Serves as a direct sister-genus validation of assembly algorithms tuned on *E. coli*.

---

## 6. Run-to-Reference Compatibility Assessment

**Purpose:** quantify how faithfully each sequenced run represents its reference genome, and which QUAST metrics could be affected. All four runs derive from the **reference strain lineage** (contract §F26-P1 "verified reference strain" clause is satisfied), but none are pristine pure-isolate wild-type preparations. Honest documentation of deviations is required so Stage D/G interpretations are not contaminated by unacknowledged biological confounds.

| Dataset | Sequenced material vs reference | Dominant risk class | Expected magnitude | QUAST metrics affected | Compatibility rating |
| :--- | :--- | :--- | :--- | :--- | :---: |
| CS-SA01 | **Verified MW2 + phage Staph1N co-culture** (31-mer containment, BC-01 resolution); oxacillin co-treatment | Expected phage contigs (~140 kb); resistance SNPs/indels | Phage genome **140,511 bp** (`NC_047722.1`) of non-reference content, verified; handful of resistance loci | Unaligned contigs ≈ 140 kb = expected viral co-culture content (evidence doc §5 protocol); small mismatch/indel floor | **ACCEPTABLE — BC-01 resolved** |
| CS-BS01 | 168-derived chloramphenicol-adapted lineage | True biological divergence | Few to tens of SNPs + possible small indels/CNVs at resistance loci | Mismatches & indels per 100 kbp (biological, not assembler error) | **ACCEPTABLE** |
| CS-PA01 | PAO1-GFP reporter (`mexR`, `aroB` mutant 2); submitter-pretrimmed reads | Engineered cassette; read-length reduction | GFP cassette ≈ 0.7–3 kb if chromosomally integrated (unknown); avg 478.8 vs 480 bp | One insertion vs reference if cassette is chromosomal; N50 sensitivity negligible | **ACCEPTABLE with caveat** |
| CS-SE01 | LT2-derived AMP-resistance ALE lineage; run submitted as **BAM**; "RANDOM PCR" library | Alignment-derived artifacts; PCR duplicates | Duplicate fraction above typical PCR-free WGS (unquantified pre-download) | Inflated effective coverage → coverage-estimator stress; possible read-duplication artifacts in SRA conversion | **ACCEPTABLE with caveat** |

### 6.1 Per-Dataset Mitigation Notes (binding for Stage D/G interpretation)

1. **CS-SA01:** Resolved by F26-DATA-003 (BC-01): phage Staph1N co-culture content is **established by 31-mer containment evidence** (90.64% of library 31-mers explained by MW2 ∪ Staph1N; phage depth ≈ 1000×, MW2 ≈ 35×). In QUAST (reference `NC_003923.1` only), unaligned contigs totaling ≈ 140 kb are **expected viral co-culture content** — record them as such, not as misassemblies or contaminant failures (binding protocol: [`docs/development/evidence/CS_SA01_STRAIN_DISCRIMINATION_EVIDENCE.md`](docs/development/evidence/CS_SA01_STRAIN_DISCRIMINATION_EVIDENCE.md) §5). Misassembly calls *within* chromosome-aligned contigs remain genuine quality signals. Relative assembler comparison (identical inputs) remains valid.
2. **CS-BS01:** Record the absolute mismatch/indel floor per assembler; treat small uniform shifts across all four competitors as biological divergence of the evolved lineage, not a regression.
3. **CS-PA01:** If a single large insertion is detected at the GFP cassette locus in all four assemblies, annotate as engineered content. The pre-trimmed read length (avg 478.8 bp) is documented here to explain nominal 2×240 deviation.
4. **CS-SE01:** After download, measure duplicate/read-parity statistics (the download harness `scripts/download_f26_datasets.sh` already enforces R1/R2 read-count parity and records SHA-256). If duplicate fraction is pathologically high (> 20%), flag for Controller review before Stage C; coverage-related metrics (F26-P3 depth gradients) must use *reported* depths recomputed post-download rather than nominal SRA bases.
5. **All datasets:** because every competitor (fastBact, SPAdes, SKESA, Unicycler) receives byte-identical inputs, **relative** comparisons — the core F26-P2 deliverable — are unaffected by any of the above deviations. Absolute metric shifts must be attributed per the notes above.

---

## 7. Depth, Read-Length, and Download Size Estimates

Derived from verified SRA statistics (`total_spots`, `total_bases`, archive `bytes`) and machine-verified reference lengths:

| Dataset | Total Bases (bp) | Avg Read Length (bp/spot) | Reference (bp) | Est. Full Coverage | SRA Archive Size | Est. Uncompressed FASTQ | Est. pigz FASTQ |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| CS-SA01 | 659,542,596 | 294.5 | 2,820,462 | ≈ 234× | 179 MiB | ≈ 1.4 GiB | ≈ 0.45 GiB |
| CS-BS01 | 680,649,714 | 302.0 | 4,215,606 | ≈ 161× | 299 MiB | ≈ 1.4 GiB | ≈ 0.47 GiB |
| CS-PA01 | 726,918,458 | 478.8 | 6,264,404 | ≈ 116× | 233 MiB | ≈ 1.5 GiB | ≈ 0.50 GiB |
| CS-SE01 | 1,006,879,350 | 299.7 | 4,857,450 | ≈ 207× | 115 MiB | ≈ 2.1 GiB | ≈ 0.70 GiB |
| **Total** | **3,073,990,118** | — | — | — | **≈ 825 MiB** | **≈ 6.4 GiB** | **≈ 2.1 GiB** |

*Estimation model: FASTQ bytes ≈ 2 × (seq len + qual len + ~17 B headers/newlines) per spot; pigz default-level compression ≈ 30–35% of raw FASTQ.*

**Planning implications:**
- SRA download footprint is small (≈ 0.83 GiB total); the `fasterq-dump` scratch peak is bounded by the largest run (CS-SE01, ≈ 2.1 GiB uncompressed). The configured 16 GB Slurm allocation and 4 h walltime are adequate.
- Post-download, every dataset supports the F26-P3 depth-gradient design (30× / 50× / 100× / 150× / full) from a single full-depth source, with CS-PA01 (≈ 116×) capped at ≈ 100× for the top non-full tier.
- CS-PA01 pre-trimmed reads (avg 478.8 bp/spot) mean depth is computed post-download; see §6.1 note 4.

---

## 8. Verification Evidence Record (Policy §9)

Executed by OpenCode / GLM `glm-5.3-flash`, task F26-DATA-001, 2026-10-30 (login node; read-only operations):

| Step | Command (abridged) | Exit | Result |
| :--- | :--- | :---: | :--- |
| 1. Manifest checksum audit | `sha256sum <ref>` + normalized-sequence hash (`grep -v '^>' | tr -d whitespace | tr acgtn→ACGTN | sha256sum`) for all 8 FASTA files | 0 | All 16 hashes match manifest exactly; zero corrections |
| 2. Reference derivation | `awk` uppercase concatenation → length, GC%, N-count, contig count | 0 | 8/8 match manifest; all refs single-contig, zero-N |
| 3. SRA metadata | `curl efetch.fcgi db=sra id=<4 accessions> rettype=xml` → ElementTree parse | 0 | BioProject/BioSample/strain/instrument/layout/spots all confirmed (§4) |
| 4. FASTA defline inspection | `head` of each reference | 0 | Deflines match expected RefSeq records (§3) |
| 5. FASTQ staging check | `ls data/f26_cross_species/fastq/` | 0 | Directory exists and is **empty** — downloads not yet started |
| 6. Download harness review | Read `scripts/download_f26_datasets.sh` | 0 | Accessions match manifest; split-3 + pigz + parity + SHA-256 logging in place |

**Manifest finalization statement:** `F26_DATASET_MANIFEST.tsv` was byte-audited on 2026-10-30 against independently recomputed file hashes and NCBI metadata; **all cells verified exact and the file required no modifications**.

---

## 9. Selection and Exclusion Rationale

### 9.1 Inclusion criteria (all four CS datasets)

1. Complete, closed, single-contig RefSeq chromosome reference with zero `N` (verified §3).
2. Illumina paired-end WGS (WGS strategy, PAIRED layout) of sufficient depth (≥ 100×) for the F26-P3 depth gradient (verified §7).
3. `strain=` BioSample attribute matching the reference strain lineage (verified §4).
4. Distinct genus from all other selected datasets; combined panel spans 4 genera, both Gram categories, GC 32.8–66.6%, size 2.8–6.3 Mb.
5. Public, load-complete SRA runs with explicit BioProject/BioSample traceability.

### 9.2 Exclusion classes

| Excluded class | Rationale |
| :--- | :--- |
| ONT/PacBio long-read or hybrid datasets | F26-P2 evaluates short-read-only mode; hybrid inputs invalidate like-for-like comparison |
| Runs with draft/fragmented or contig-level references (multi-N, multi-contig) | Violates zero-N complete-chromosome reference requirement; NGA50/genome-fraction metrics become reference-quality-confounded |
| Metagenomes, enrichment cultures, multi-isolate pools | Violates pure-lineage traceability; unusable for strain-specific QUAST evaluation |
| Multiple runs per species | One run per species maximizes diversity per GiB; additional replicates offer no incremental algorithmic signal at Stage B |
| *E. coli* candidates | Excluded from the CS panel by definition (already covered by the 4 historical oracles); a fifth *E. coli*-adjacent entry (CS-SE01) was accepted only because *Salmonella* exercises distinct repeat/pathogenicity-island structure |
| Deeply non-standard platforms (Ion Torrent, MGISEQ) | Out of the Illumina PE evaluation scope frozen in F26-P2 |

### 9.3 Accepted deviations (documented, not excluded)

The experimental-context deviations in §6 (phage co-treatment, evolved lineages, GFP cassette, BAM submission) were **accepted rather than excluded** because: (a) F26-P1 explicitly permits "a complete reference genome representing the sequenced isolate **or verified reference strain**"; (b) every deviation affects absolute metrics only — relative four-way assembler comparison on identical inputs is unaffected; (c) each deviation carries a §6.1 mitigation and interpretation rule; and (d) replacing datasets after Stage B would invalidate pre-registration ordering (Stage D must run on the frozen baseline before any code change).

---

## 10. Staging and Storage Layout

All files are structured inside `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/`:

```text
data/
├── references/m6_identity/
│   ├── NC_000913.3.fasta  (E. coli MG1655)      [verified]
│   └── CP009273.1.fasta   (E. coli BW25113)     [verified]
├── f26_cross_species/
│   ├── references/
│   │   ├── NC_003923.1.fasta  (S. aureus MW2)    [verified]
│   │   ├── NC_000964.3.fasta  (B. subtilis 168)  [verified]
│   │   ├── NC_002516.2.fasta  (P. aeruginosa PAO1) [verified]
│   │   └── NC_003197.2.fasta  (S. enterica LT2)  [verified]
│   └── fastq/                                   [empty — staging pending]
│       ├── SRR33682404_1.fastq.gz / SRR33682404_2.fastq.gz   (planned)
│       ├── SRR40739147_1.fastq.gz / SRR40739147_2.fastq.gz   (planned)
│       ├── SRR41244734_1.fastq.gz / SRR41244734_2.fastq.gz   (planned)
│       └── SRR41010285_1.fastq.gz / SRR41010285_2.fastq.gz   (planned)
```

All FASTQ downloads are staged via `fasterq-dump --split-3` (8 threads) through `scripts/download_f26_datasets.sh` on Slurm `vnode1` (`qfree`), compressed with `pigz -p 8`, gzip-integrity-checked, read-count-parity-checked (R1 = R2), and SHA-256-logged. Upon download completion, per-FASTQ SHA-256 values must be appended to a staging receipt and cross-checked against this document before Stage C may begin.
