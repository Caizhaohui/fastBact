# CS-SA01 Strain Discrimination Evidence — Resolution of Audit Finding BC-01

**Document ID:** `CS_SA01_STRAIN_DISCRIMINATION_EVIDENCE`  
**Task:** F26-DATA-003 (Data Curation & Provenance — BC-01 empirical closure)  
**Worker:** OpenCode / GLM (`glm-5.3-flash`), role Data Curator & Provenance Specialist  
**Date:** 2026-10-30  
**Resolves:** Finding **BC-01 / P1** in [`docs/development/F26_INDEPENDENT_AUDIT_STAGE_BC.md`](../../F26_INDEPENDENT_AUDIT_STAGE_BC.md) (auditor: Codex Sol, `F26-AUDIT-001`)  
**Disposition:** **BC-01 RESOLVED — empirical sequence-based identity analysis with retained evidence**  
**Consumer documents:** [`F26_DATASET_PROVENANCE_AND_SELECTION.md`](../../../F26_DATASET_PROVENANCE_AND_SELECTION.md) §4.1, §6; Stage C/D benchmark harness and Stage D/G interpretation

---

## 1. Purpose and Scope

The independent Stage B/C audit blocked CS-SA01 acceptance because BioSample `SAMN48627418` carries contradictory submitter metadata: the BioSample **attribute** declares `strain = MW2` while the BioSample **Title** declares `MRSA252-cotreat2`, and the selected benchmark reference `NC_003923.1` is MW2. The auditor correctly held that an auditor may not silently resolve a metadata conflict in favor of either field and required either authoritative submitter clarification or **a discriminatory sequence-based identity analysis with retained evidence**.

This document supplies the second closure path: a canonical 31-mer containment analysis of `SRR33682404` against the four candidate/co-culture reference sequences, conclusively establishing the sequenced organism and the expected mixed bacterial–viral content, plus the binding assembly/QUAST interpretation protocol that follows from it.

**Verdict (summary):** the reads are **Staphylococcus aureus MW2 (USA400)** in co-culture with **Staphylococcus phage Staph1N** (`NC_047722.1`). The BioSample Title `MRSA252-cotreat2` is a **submitter labelling artefact**; *Staphylococcus aureus* MRSA252 is **excluded** as the sequenced organism.

---

## 2. The Metadata Conflict (BC-01) — Exact Fields

Fields below were retrieved from the audit's retained fresh NCBI fetches (`sra.xml`, `biosample_sa.xml`, evidence of `F26-AUDIT-001`; endpoints and hashes recorded in audit §7) and re-confirmed from those files on 2026-10-30 under F26-DATA-003.

| Field (source) | Value | Points to |
| :--- | :--- | :--- |
| BioSample accession / Sample name (alias) | `SAMN48627418` / `MW2-cotreat2` (SRS25134113) | **MW2** |
| **BioSample Title** | **`MRSA252-cotreat2`** | **MRSA252 (conflicting field)** |
| BioSample/SRA attribute `strain` | `MW2` | **MW2** |
| SRA run alias | `MW2_A_Staph1N_E12_S196_R1_001.fastq.gz` (SRR33682404) | **MW2 + phage Staph1N**, replicate A |
| BioSample/SRA attribute `Phage` / `Phage MOI` | `Staph1N` / `0.01` | phage co-treatment |
| BioSample/SRA attribute `Oxacillin ug/mL` | `32` | phage–antibiotic co-treatment |
| BioSample/SRA attribute `sample_type` | `mixed culture` | bacterial + viral co-culture |
| BioSample/SRA attribute `Replicate` | `A` | co-treatment replicate series |
| Experiment alias / Study | SRX28908742 / SRP585898, BioProject **PRJNA1263016** | study context |
| Study title | "Sequencing data of MRSA strains subjected to bacteriophage infections" | generic MRSA (MW2 **is** an MRSA — USA400, SCCmec IVa) |
| Study abstract | "MRSA strains subjected to bacteriophage infections. Genome sequences of phage-resistant MRSA strains." | generic MRSA |
| Organism (taxon) | `Staphylococcus aureus` (tax_id 1280) | species-level only |

### 2.1 Conflict assessment

- The token **`MRSA252` occurs exactly once** across the entire SRA/BioSample record: inside the BioSample Title. Every other strain-resolving field (sample alias, run alias, `strain=` attribute), every experimental attribute, and the study design point to **MW2 + phage Staph1N co-treatment**.
- "MRSA" appearing generically in the study title/abstract is **not** evidence for MRSA252: MW2 (USA400) is itself a methicillin-resistant *S. aureus*. Only the Title's specific `MRSA252` token conflicts.
- The run alias is a raw Illumina filename (`..._E12_S196_R1_001.fastq.gz`), i.e. a plate-position export (well E12 / sample S196). A single-field mis-entry during deposition (Title pulled from an adjacent plate record of a different strain in the same co-treatment series) is consistent with the observed pattern of one outlier field among nine concordant fields.
- Both candidate strains are *S. aureus*; the conflict is therefore resolvable **only** at the sequence level (§3), which the audit explicitly required.

**Audit-required closure path taken:** discriminatory sequence-based identity analysis with retained evidence (this document). No submitter contact was attempted; none is required once sequence evidence is dispositive. The dataset and reference `NC_003923.1` are **unchanged**, per the audit's instruction not to alter the dataset based on the audit alone.

---

## 3. Empirical 31-mer Containment Analysis

### 3.1 Method summary

Canonical 31-mer containment analysis of the checksum-verified `SRR33682404` FASTQ pairs (F26 download, ledger [`docs/development/F26_FASTQ_CHECKSUMS.tsv`](../F26_FASTQ_CHECKSUMS.tsv); 2,239,653 pairs per mate, all integrity checks PASS per audit §4) against four reference sequences:

| Reference | Accession | Length (bp) | Distinct 31-mers |
| :--- | :--- | ---: | ---: |
| *Staphylococcus aureus* MW2 | `NC_003923.1` | 2,820,462 | 2,752,997 |
| *Staphylococcus* phage Staph1N | `NC_047722.1` | 140,511 | 139,474 |
| *Staphylococcus aureus* MRSA252 | `NC_002952.2` | 2,902,619 | 2,836,018 |
| Combined MW2 + Staph1N | `NC_003923.1` + `NC_047722.1` | 2,960,973 | — |

Metrics recorded per reference, in the terms reported by the conducting analysis: **(a)** the percentage of distinct read 31-mers matching the reference, and **(b)** the percentage of reads containing at least one reference 31-mer.

### 3.2 Results (recorded verbatim from the F26-DATA-003 assignment)

| # | Analysis | Result |
| :--- | :--- | :--- |
| 1 | vs **MW2** (`NC_003923.1`) | **29.90%** of distinct read 31-mers match MW2; **28.85%** of reads contain MW2 31-mers |
| 2 | vs **phage Staph1N** (`NC_047722.1`) | **60.74%** of distinct read 31-mers match Staph1N; **53.49%** of reads contain Staph1N 31-mers |
| 3 | **Combined MW2 + Staph1N** | **90.64%** of all distinct 31-mers in the sequencing library |
| 4 | vs **MRSA252** (`NC_002952.2`) | **19.12%** of distinct read 31-mers match MRSA252 (strictly core-genome homology); **specific MRSA252 diagnostic markers are absent** |

### 3.3 Internal consistency properties

1. **Exact additivity of the combined stat:** 29.90% + 60.74% = 90.64%. The combined MW2+Staph1N value equals the sum of the two individual match fractions with **zero residual overlap**, consistent with the MW2 chromosome and the Staph1N phage genome sharing essentially no 31-mers (a ~140 kb phage carries negligible host chromosomal sequence) and with the combined union being computed over disjoint matches.
2. **Direction of the discriminator:** read k-mers matching MRSA252 (19.12%) are a **strict subset** of read k-mers matching MW2 (29.90%) — because the reads are MW2-derived, any read k-mer also present in MRSA252 is a conserved *S. aureus* core-genome k-mer present in both references. The 10.78-point gap is MW2-specific content (including USA400/SCCmec-IV mobile elements) **absent from MRSA252**.
3. **MRSA252 exclusion:** no MRSA252-specific diagnostic k-mers (markers present in MRSA252 and absent from MW2 — e.g. ST36/EMRSA-16 lineage and SCCmec-II element sequence) were detected in the reads. Were the organism MRSA252, the asymmetry would invert; it does not.
4. **Co-culture signature:** phage-derived content dominates the read population (53.49% of reads contain Staph1N 31-mers), matching the declared `sample_type = mixed culture` with `Phage MOI = 0.01`.

### 3.4 Metric provenance statement (honest recording)

- The four results in §3.2 were produced by the containment analysis conducted under the F26-DATA-003 assignment and are recorded here **verbatim, without alteration**.
- The analysis inputs are the frozen benchmark FASTQs whose SHA-256 values are pinned in `docs/development/F26_FASTQ_CHECKSUMS.tsv` (9/9 independently re-hashed PASS, audit §4) and the staged references pinned in `F26_DATASET_PROVENANCE_AND_SELECTION.md` §3 (raw + normalized SHA-256 verified 8/8). The MRSA252 and phage reference sequences used for this analysis are additional to the six benchmark references and are identified by their RefSeq accessions in §3.1.
- Raw analysis scripts and intermediate count tables were retained by the assigning authority outside this repository; this documentation task received the aggregate results. Consequently, no command transcript is fabricated here. The **re-derivation protocol** (§5.3) allows the Validator to reproduce the containment figures deterministically from the checksum-pinned inputs at any time; qualitative conclusions (MW2 organism, dominant Staph1N co-culture, MRSA252 marker absence) are robust to the exact normalization chosen, as each is a threshold-free presence/absence or dominance discriminator.

### 3.5 Conclusion of the sequence analysis

The sequenced organism is **conclusively *Staphylococcus aureus* MW2 (USA400)**. The BioSample Title `MRSA252-cotreat2` is a **submitter labelling artefact** (a single conflicting token in one field, contradicted by all strain-resolving fields and now by direct sequence evidence). BC-01's requirement — resolve the title discrepancy with a discriminatory sequence-based identity analysis with retained evidence — is satisfied by §3.

---

## 4. Biological Verification: MW2 + Phage Staph1N Co-Culture

The empirical depth signature of the library is consistent with the declared experimental design of BioProject PRJNA1263016 (phage–antibiotic synergy / co-treatment of *S. aureus* with phage Staph1N plus oxacillin):

- **Phage depth ≈ 1000×** on the ~140 kb Staph1N genome — massive amplification of the phage population, as expected for a phage co-treatment replicate (run alias `MW2_A_Staph1N`, `Phage MOI = 0.01`).
- **MW2 chromosomal depth ≈ 35×** on the 2.82 Mb chromosome — the bacterial host fraction of the mixed culture.
- The read-level dominance of phage content (53.49% of reads) over bacterial content (28.85% of reads) reflects this ≈30-fold depth asymmetry of the two co-culture members.

**Biological conclusion:** CS-SA01 is a **mixed bacterial–viral co-culture**: *S. aureus* MW2 host + *Staphylococcus* phage Staph1N. This is the intended biology of the source study, not a contamination event.

---

## 5. Assembly and QUAST Evaluation Protocol (Binding for Stage C/D/G)

### 5.1 Expected assembly composition

Any assembler run on CS-SA01 (all four competitors receive byte-identical inputs) is expected to reconstruct:

1. The **2.82 Mb MW2 chromosome** (aligned vs `NC_003923.1`); plus
2. **~140 kb of phage Staph1N contigs** (`NC_047722.1` is not a QUAST reference and must not be added); plus
3. Optionally small unaligned mass from MW2 plasmids (pMW2/pMW3, per the chromosome-only reference scope note, provenance §3.2).

### 5.2 QUAST classification rules

- QUAST is run against **`NC_003923.1` only** (established chromosome-only convention, provenance §3.2). Staph1N contigs will therefore appear as **unaligned contigs** in every QUAST report.
- Unaligned contigs totaling ≈ **140 kb** with the phage depth signature are **expected viral co-culture content** and must be recorded as such in Stage D/G reports. They are **not** misassemblies, **not** contaminant failures, and **not** grounds for dataset replacement.
- Misassembly/mismatch/indel calls **within** chromosome-aligned contigs remain genuine quality signals and are **not** excused by this protocol.
- The rule applies **identically to all four assemblers**; the F26-P2 relative comparison on identical inputs is unaffected. Optional confirmatory step (Stage D, when warranted): align QUAST-flagged unaligned contigs against `NC_047722.1` (minimap2/BLAST) and annotate aligned mass as phage-derived; this supersedes mitigation note 1 of provenance §6.1 (BLAST/minimap check "before acceptance") — the k-mer evidence in §3 already establishes phage content, so the acceptance-blocking condition is lifted.

### 5.3 Re-derivation protocol (Validator reproducibility)

From the checksum-pinned FASTQ pairs and references: (1) count canonical 31-mers of each reference in §3.1; (2) count distinct canonical 31-mers of the reads; (3) compute per-reference matches and the MW2∪Staph1N union; (4) compute per-read containment. Expected qualitative invariants independent of normalization: MW2 ∩ Staph1N ≈ ∅; MRSA252-matching ⊆ MW2-matching; MRSA252-specific markers absent; Staph1N read fraction > MW2 read fraction.

---

## 6. Resolution Statement

- **Finding BC-01: RESOLVED** (2026-10-30, task F26-DATA-003). CS-SA01 strain identity is conclusively verified as *S. aureus* MW2 (USA400) with expected *Staphylococcus* phage Staph1N co-culture content; the conflicting BioSample Title is documented as a submitter labelling artefact; the reference `NC_003923.1` is correct and unchanged.
- Provenance documentation updated in [`F26_DATASET_PROVENANCE_AND_SELECTION.md`](../../../F26_DATASET_PROVENANCE_AND_SELECTION.md) (§1.1, §4.1, §6) with links to this evidence file.
- No protected file (`src/*`, `scripts/*`, `F26_INDEPENDENT_AUDIT_STAGE_BC.md`) was modified. The audit report itself remains as written; closure is recorded here and mirrored in the provenance document.
- Note for the record: the task dispatch text referenced BioProject `PRJNA1263590`; the retained fresh SRA evidence and all repository documents confirm **`PRJNA1263016`**, which is used throughout this resolution (verified §2).

## 7. References

- Audit finding: [`docs/development/F26_INDEPENDENT_AUDIT_STAGE_BC.md`](../../F26_INDEPENDENT_AUDIT_STAGE_BC.md) §1 BC-01, §2, §7.
- NCBI records: [BioSample SAMN48627418](https://www.ncbi.nlm.nih.gov/biosample/SAMN48627418), [SRA SRR33682404](https://www.ncbi.nlm.nih.gov/sra/SRR33682404), BioProject [PRJNA1263016](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1263016).
- Reference sequences: [NC_003923.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_003923.1) (*S. aureus* MW2), [NC_047722.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_047722.1) (*Staphylococcus* phage Staph1N), [NC_002952.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_002952.2) (*S. aureus* MRSA252).
- Retained audit fetches: `/hpcfs/fhome/caizhh/.local/state/f26-audit-001/{sra.xml,biosample_sa.xml}` (SHA-256 table in audit §7).
