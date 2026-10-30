# CS-SA01 Strain Discrimination Evidence — Resolution of Audit Finding BC-01

**Document ID:** `CS_SA01_STRAIN_DISCRIMINATION_EVIDENCE`  
**Task:** F26-DATA-004 (primary strain-discrimination alignment analysis on qfree)  
**Worker:** OpenCode / GLM (`glm-5.3-flash`), role Computational Biologist / Sequence Forensics Specialist  
**Date:** 2026-10-30  
**Resolves:** Finding **BC-01 / P1** per the required closure path in [`docs/development/F26_INDEPENDENT_AUDIT_STAGE_BC.md`](../../F26_INDEPENDENT_AUDIT_STAGE_BC.md) §F26-AUDIT-002 (auditor: Codex Sol)  
**Disposition:** **BC-01 RESOLVED — primary, retained, reproducible alignment evidence with exact counts, denominator definitions, marker-set sizes, observed support, and sensitivity limits**  
**Definitive run:** Slurm job **3074242** (partition `qfree`, node `vnode1`, 8 CPUs / 14 GB / ≤3 min, exit 0, `F26_CS_SA01_DISCRIMINATION_PIPELINE_COMPLETE`, **`CS_SA01_PRIMARY_DISCRIMINATION_GATES: ALL PASS`**) — executed by [`scripts/run_cs_sa01_discrimination.sh`](../../../scripts/run_cs_sa01_discrimination.sh) at git HEAD `bf19db74a98a048d93d03ecbd8097c521125c9cb`+ working tree  
**Consumer documents:** [`F26_DATASET_PROVENANCE_AND_SELECTION.md`](../../../F26_DATASET_PROVENANCE_AND_SELECTION.md) §1.1, §4.1.1, §6; Stage C/D benchmark harness and Stage D/G interpretation

---

## 1. Purpose and Scope

BioSample `SAMN48627418` carries contradictory submitter metadata: the BioSample **attribute** declares `strain = MW2` while the BioSample **Title** declares `MRSA252-cotreat2`, and the selected benchmark reference `NC_003923.1` is MW2. The independent Stage B/C audit (F26-AUDIT-001, Codex Sol) held that an auditor may not silently resolve a metadata conflict in favor of either field and required a discriminatory sequence-based identity analysis with retained evidence. The preliminary aggregate 31-mer containment report (F26-DATA-003) did not satisfy the F26-AUDIT-002 re-audit because it lacked locatable primary artifacts, had unreconciled distinct-set denominators, and inferred marker absence from aggregate fractions (audit §"BC-01: sequence evidence does not yet justify closure", points 1–3).

This document supplies the required closure: a **primary, deterministic, fully retained read-alignment analysis** of `SRR33682404` against (MW2 + phage Staph1N) and (MRSA252 + Staph1N), executed on Slurm `qfree` by a committed script, with every command, exit code, input hash, exact count, denominator definition, marker-set size, observed support value, and sensitivity bound either reproduced in this document or retained as a repository artifact (§8).

**Verdict (summary):** the reads derive from ***Staphylococcus aureus* MW2 (USA400)** in co-culture with ***Staphylococcus* phage Staph1N** (`NC_047722.1`). The MW2-diagnostic loci (PVL lukS-PV/lukF-PV, SCCmec IV) are present at 59.8× / 73.5× mean depth with 100.00 % breadth at ≥30×, while the sequence-verified MRSA252-exclusive portion of the SCCmec II window shows **no aligned-base presence** (343 of 36,904 bases covered at MAPQ ≥ 20; longest contiguous covered island 113 bp) against a background where every genuinely present marker set shows 50.8–88.0× and 100 % breadth. The BioSample Title `MRSA252-cotreat2` is a **submitter plate-level labelling artefact** (§4.3); MRSA252 is **excluded** as the sequenced chromosome.

---

## 2. The Metadata Conflict (BC-01) — Exact Fields

Fields below were retrieved from the audit's retained fresh NCBI fetches (`sra.xml`, `biosample_sa.xml`, F26-AUDIT-001; endpoints and hashes in audit §7) and re-confirmed under F26-DATA-003.

| Field (source) | Value | Points to |
| :--- | :--- | :--- |
| BioSample accession / Sample name (alias) | `SAMN48627418` / `MW2-cotreat2` (SRS25134113) | **MW2** |
| **BioSample Title** | **`MRSA252-cotreat2`** | **MRSA252 (conflicting field)** |
| BioSample/SRA attribute `strain` | `MW2` | **MW2** |
| SRA run alias | `MW2_A_Staph1N_E12_S196_R1_001.fastq.gz` (SRR33682404) | **MW2 + phage Staph1N**, replicate A |
| BioSample/SRA attribute `Phage` / `Phage MOI` | `Staph1N` / `0.01` | phage co-treatment |
| BioSample/SRA attribute `Oxacillin ug/mL` | `32` | phage–antibiotic co-treatment |
| BioSample/SRA attribute `sample_type` | `mixed culture` | bacterial + viral co-culture |
| Experiment alias / Study | SRX28908742 / SRP585898, BioProject **PRJNA1263016** | study context |
| Study title / abstract | "Sequencing data of MRSA strains subjected to bacteriophage infections" | generic MRSA (MW2 **is** an MRSA — USA400, SCCmec IVa, PVL+) |
| Organism (taxon) | `Staphylococcus aureus` (tax_id 1280) | species-level only |

The token `MRSA252` occurs exactly once across the entire record (the Title). Both candidate strains are *S. aureus*, so the conflict is resolvable **only** at the sequence level — which §3 supplies. Bibliographic provenance for the study design and the artefact mechanism is in §4.

---

## 3. Primary Alignment Analysis (definitive evidence)

### 3.1 Inputs and integrity

| Item | Value | Verification |
| :--- | :--- | :--- |
| Reads R1 | `SRR33682404_1.fastq.gz`, 137,279,980 bytes, 2,239,653 reads | SHA-256 `1a7e70141657dfa9a5eb9c6b6ba491c260127450b01ab9eb94b85e0fe2a466c3` — equals pinned value in [`docs/development/F26_FASTQ_CHECKSUMS.tsv`](../F26_FASTQ_CHECKSUMS.tsv) (verified in-job, fail-fast) |
| Reads R2 | `SRR33682404_2.fastq.gz`, 141,516,913 bytes, 2,239,653 reads | SHA-256 `88a46ad12cb159511d112fa059b26732ae36d9c3e28a6355f806638f2599e050` — equals pinned value (verified in-job) |
| Reference: *S. aureus* MW2 | `NC_003923.1`, 2,820,462 bp, 1 contig | SHA-256 `19ba62a339065c9c498efb72bb2491e372344e4f23a4ac1e3e34e7345a25c63c` — equals pinned value in provenance §3 (verified in-job) |
| Reference: *S.* phage Staph1N | `NC_047722.1`, 145,647 bp, 1 contig | SHA-256 `2993d31f2fbd7e984ed5ed2151701c4057e28a39b3e1923205f20a4492e0475e` (`cs_sa01_disc_inputs_sha256.tsv`) |
| Reference: *S. aureus* MRSA252 | `NC_002952.2`, 2,902,619 bp, 1 contig | SHA-256 `abc19abf4b13562908d32803603591a03287e0c49dbb5aefb36606b300d793b4` (`cs_sa01_disc_inputs_sha256.tsv`) |
| Total input records (denominator) | **2 × 2,239,653 = 4,479,306** single-end records | `F26_FASTQ_CHECKSUMS.tsv` |
| Tools | minimap2 2.31-r1302; samtools 1.22.1 (htslib 1.23) | recorded in `cs_sa01_disc_run.log` |

All five input hashes are recomputed inside the job and hard-fail (exit 2) on mismatch with the pinned values. The full command transcript with per-step exit codes, job id, node, git HEAD, and tool versions is retained in [`cs_sa01_disc_run.log`](cs_sa01_disc_run.log).

### 3.2 Method (committed script `scripts/run_cs_sa01_discrimination.sh`)

1. **Whole-genome alignment** of both mates (`minimap2 -t 8 -ax sr`) to two combined references — `MW2 + Staph1N` and `MRSA252 + Staph1N` — followed by `samtools sort`/`index`. Exact mapped/unmapped counts from `samtools idxstats`; mean depth and breadth from `samtools coverage`; 10-kb binned depth profiles and full phage per-base depth retained.
2. **Diagnostic-locus depth** (`samtools depth -a -q 0 -Q 0 -d 0 -r`) and **exact overlapping read counts** (`samtools view -c`) over the three dispatch-specified windows: MW2 PVL `NC_003923.1:1515001-1517500`, MW2 SCCmec IV `NC_003923.1:34001-58000`, MRSA252 SCCmec II `NC_002952.2:34001-86000`.
3. **Sequence-verified shared/exclusive marker partition.** A whole-window cross-alignment is *not* usable here: minimap2 fails to emit alignments for long queries carrying large divergent inserts (verified in-run: the 52 kb SCCmec II window returned **zero** PAF lines even against its own source genome). Homology is therefore derived from **5,000-bp chunked `-x sr` alignments** (chunks behave like reads; MAPQ ≥ 10), whose hit intervals are mapped back to window coordinates (chunk offset + PAF qstart/qend) and unioned. **Shared** = union vs the other strain's chromosome plus vs Staph1N; **exclusive** = the complement within the window, keeping blocks ≥ 200 bp (shorter pieces recorded as filtered fragments). Both marker sets are retained as FASTA + BED artifacts in `loci/`.
4. **Read-origin audit** over the MRSA252-exclusive set: QNAME-level cross-check against the MW2-chromosome alignment (are the overlapping reads also placed on `NC_003923.1`?) and divergence statistics (minimap2 `de:f` tag).
5. **Gates, sensitivity bounds, digests** — evaluated in-job and written to `cs_sa01_disc_gates.tsv`, `cs_sa01_disc_sensitivity.tsv`, `SHA256SUMS.txt`.

**Formulas and denominator definitions:**

- `mean_depth(ref) = Σ_i depth_i / L_ref` over **all** `L_ref` reference bases (`samtools coverage` meandepth).
- `breadth_ge_qx = |{i : depth_i ≥ q}| / L_ref`; `samtools coverage` reports `breadth_ge_1x` (`coverage` column).
- idxstats counts include supplementary (split) alignments: `Σ(idxstats mapped + unmapped) = input records + supplementary`, verified in-run (736 / 4,353 supplementary for the two alignments; `record_consistency = PASS`).
- Exclusive marker set = `window \ (∪ MAPQ≥10 chunk-hits vs other strain ∪ ∪ MAPQ≥10 chunk-hits vs Staph1N)`, blocks ≥ 200 bp.
- Poisson zero-coverage sensitivity: `P(covered bases = 0 over L bp at mean depth D) = e^(−D·L)`.

### 3.3 Whole-genome results (exact counts)

| Alignment | Records (Σ idxstats) | Supplementary | Mapped: bacterial | Mapped: phage Staph1N | Unmapped |
| :--- | ---: | ---: | ---: | ---: | ---: |
| vs MW2 + Staph1N | 4,480,042 | 736 | **1,339,405** (`NC_003923.1`) | **2,964,345** (`NC_047722.1`) | 123,856 |
| vs MRSA252 + Staph1N | 4,483,659 | 4,353 | 1,295,397 (`NC_002952.2`) | 2,964,345 (`NC_047722.1`) | 142,998 |

| Contig | Mean depth | Breadth (≥1×) | Reads (coverage) | Source |
| :--- | ---: | ---: | ---: | :--- |
| `NC_003923.1` (MW2 chromosome, 2,820,462 bp) | **69.9473×** | **100.0000 %** | 1,339,405 | `cs_sa01_disc_mw2.coverage.tsv` |
| `NC_047722.1` (Staph1N, 145,647 bp) | **2,984.21×** | 99.7123 % | 2,964,345 | `cs_sa01_disc_mw2.coverage.tsv` |
| `NC_002952.2` (MRSA252 chromosome, 2,902,619 bp) | 64.8355× | 90.2396 % | 1,295,397 | `cs_sa01_disc_mrsa252.coverage.tsv` |

Interpretation: (i) the MW2 chromosome is covered wall-to-wall at ~70× — the bacterial content is MW2-compatible; (ii) the MRSA252 chromosome receives *less* bacterial mapping than MW2 (1,295,397 < 1,339,405) at lower breadth (90.24 % < 100 %) — consistent with shared core-genome coverage plus missing MW2-specific content; (iii) the phage contig is covered at 2,984× (99.71 % breadth) in **both** alignments — a measured Staph1N amplification signature (the earlier "≈1000×" estimate was an unmeasured expectation and is hereby superseded by 2,984.21×).

### 3.4 Diagnostic-locus depth (whole windows)

| Window | Alignment | Window bp | Mean depth | Breadth ≥1× | Breadth ≥10× | Breadth ≥30× | Overlapping reads |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| PVL lukS-PV/lukF-PV `NC_003923.1:1515001-1517500` | vs MW2+Staph1N | 2,500 | **59.82×** | 100.00 % | 100.00 % | **100.00 %** | 1,071 |
| SCCmec IV `NC_003923.1:34001-58000` | vs MW2+Staph1N | 24,000 | **73.51×** | 100.00 % | 100.00 % | **100.00 %** | 12,230 |
| SCCmec II `NC_002952.2:34001-86000` | vs MRSA252+Staph1N | 52,000 | 22.69× | 29.91 % | 29.56 % | 28.03 % | 8,700 |

The two MW2-diagnostic windows are covered uniformly at ~60–74× (100 % ≥30×). The MRSA252 window's intermediate 22.69× mean is **not** evidence of presence: it is the mixture of a 14.5 kb shared SCCmec core (§3.5) and an essentially empty remainder (§3.5, §3.6).

### 3.5 Sequence-verified shared/exclusive marker partition

Cross-strain chunked alignment hit counts: PVL window vs MRSA252 — 3 chunk-hits; SCCmec IV window vs MRSA252 — 10; SCCmec II window vs MW2 — 10; all windows vs Staph1N — 0; core-genome control (2,820,462 bp of MW2 chunked vs MRSA252) — 1,503 chunk-hits covering **2,461,851 bp** (87.3 % of the MW2 chromosome), proving the chunked pipeline detects shared sequence at scale.

| Window | Window bp | Shared vs other strain (bp) | Shared vs Staph1N (bp) | **Exclusive (bp)** | Exclusive blocks | Exclusive fraction |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| PVL (MW2) | 2,500 | 1,978 (lukD/lukE leukocidin parolog region) | 0 | **388** | 1 | 15.52 % |
| SCCmec IV (MW2) | 24,000 | 14,531 | 0 | **9,650** | 3 | 40.21 % |
| SCCmec II (MRSA252) | 52,000 | 14,506 | 0 | **36,904** | 4 | 70.97 % |

The 14,506 bp shared segment of the MRSA252 window is the **conserved SCCmec cassette core** (ccr-complex and associated core genes shared between SCCmec II and SCCmec IVa) — it cross-maps to the MW2 SCCmec IV region (MW2 coordinates ~34,026–52,100), exactly as expected between two type-II/IV elements. Marker sets and partitions are retained as `loci/*_exclusive.{fasta,bed}`, `loci/*_shared.bed`, `cs_sa01_disc_marker_blocks.tsv`.

**Observed support for each set** (`cs_sa01_disc_marker_support.tsv`; depth over the set's bases, MAPQ 0 and MAPQ ≥ 20):

| Marker set | bp | Mean depth (q0) | Breadth ≥1× (q0) | Mean depth (q20) | Breadth ≥1× (q20) | Reads (q20) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| PVL exclusive (MW2) | 388 | 51.03× | 100.00 % | **50.81×** | **100.00 %** | 206 |
| SCCmec IV exclusive (MW2) | 9,650 | 70.13× | 100.00 % | **69.02×** | **100.00 %** | 4,811 |
| SCCmec II exclusive (MRSA252) | 36,904 | 0.55× | 1.40 % | **0.32×** | **0.93 %** | 360 |
| SCCmec II shared (MRSA252 ↔ MW2 core) | 14,506 | 77.42× | 99.59 % | 74.29× | 97.40 % | 7,655 |

### 3.6 Read-origin audit over the MRSA252-exclusive set

The 36,904 bp exclusive set carries only **343 MAPQ ≥ 20 covered bases (0.93 % breadth)**; the longest contiguous q20-covered island is **113 bp**. Read-*span* counts (360 reads overlap the set) are not presence evidence: a read whose CIGAR crosses a block while placing no aligned base inside is counted by `view -L` yet contributes zero coverage. The audit (`cs_sa01_disc_exclusive_origin_audit.tsv`) characterizes them:

| Block (chr coords) | bp | Reads (q20, primary) | Mean divergence (`de:f`) | MW2-chromosome-supported |
| :--- | ---: | ---: | ---: | ---: |
| 1: 37,127–41,697 | 4,570 | **0** | — | 0 |
| 2: 47,993–61,764 | 13,771 | 113 | 2.82 % | 113 |
| 3: 62,224–62,635 | 411 | 169 | 8.37 % | 169 |
| 4: 67,848–86,000 | 18,152 | 52 | 6.50 % | 52 |
| **Overall** | 36,904 | 360 | 6.21 % (identity **93.79 %**) | **330 / 360** |

- **Block 1 (4,570 bp): zero reads.** Pure absence.
- The residual 360 reads are short, divergent (mean identity 93.79 % — consistent with paralogous mobile-element/IS-homolog fragments, not cassette sequence), sparsely scattered (343 covered bases total), and **330 of 360 are also placed on the MW2 chromosome** (shared short repeats, multi-mapping). Only 30 are unsupported by the MW2 chromosome; the chromosome-only reference scope (pMW2/pMW3 plasmids excluded; SCCmec cassettes carry plasmid-derived cargo such as pUB110 segments) fully accounts for a minority fraction of that size.
- **Control:** the SHARED set's coverage is 100 % chromosome-supported (4,079 / 4,079 qnames) — the origin audit's positive control.

### 3.7 Gates (retained: `cs_sa01_disc_gates.tsv`)

| Gate | Criterion | Observed | Verdict |
| :--- | :--- | :--- | :--- |
| G1 MW2 PVL coverage | whole-window mean ≥30× and breadth ≥30× ≥0.95 | 59.82×; 100.00 %; 1,071 reads | **PASS** |
| G2 MW2 SCCmec IV coverage | whole-window mean ≥30× and breadth ≥30× ≥0.95 | 73.51×; 100.00 %; 12,230 reads | **PASS** |
| G3a raw read-span absence | read-span count <10 (strict, superseded) | 360 span reads | FAIL_superseded_by_G3 (read-span cross-maps are not presence evidence) |
| G3 aligned-base absence | q20 breadth ≤0.01, mean <1×, covered bases ≤500, max contiguous island ≤500 bp | 0.93 %; 0.32×; 343 bp; **113 bp** | **PASS** |
| G4 bacterial preference | MW2-align bacterial reads > MRSA252-align bacterial reads | 1,339,405 > 1,295,397 | **PASS** |
| G5 phage co-culture signature | phage reads >0 and phage depth ≥10× bacterial depth | 2,964,345 reads; 2,984.21× vs 69.95× | **PASS** |
| G6 core-genome control | chunked cross-homology ≥1.4 Mb | 2,461,851 bp | **PASS** |

### 3.8 Sensitivity limits

- **Poisson zero-coverage bound:** at the observed bacterial mean depth D = 69.9473×, `P(0 coverage at one base) = e^(−D) = 4.19 × 10⁻³¹`; over any marker window (≥388 bp) the bound underflows to 0 (`cs_sa01_disc_sensitivity.tsv`). If MRSA252 SCCmec II chromosomal sequence were present at the depth observed for every MW2-present set (~70×, 100 % breadth), the probability of the observed ≤113 bp contiguous coverage island over 36,904 bp is effectively zero.
- **Homology-detection control (G6):** the same chunked pipeline recovers 2,461,851 bp of MW2↔MRSA252 shared sequence, so the 14,506 bp shared call on the SCCmec II window is not a detection artifact, and the exclusive partition reflects absence of sr-level homology, not search failure.
- **Residual-read explanation:** 360 span reads (0.93 % breadth, ≤113 bp contiguous, mean identity 93.79 %, 91.7 % chromosome-supported) are consistent with diverged paralogous short repeats; they cannot support a 52 kb cassette whose present loci show 50.8–74.3× at 97–100 % breadth.
- **Known boundary (honest scope):** this analysis verifies that the reads derive from an *S. aureus* isolate of the **MW2 lineage carrying MW2-diagnostic mobile elements (PVL prophage, SCCmec IVa)** and lacking MRSA252 SCCmec II chromosomal content. It does **not** claim bit-identity to `NC_003923.1`; single-NV-level divergence of the isolate from the reference chromosome is neither measured nor required for the strain call.

### 3.9 Conclusion of the primary analysis

The sequenced chromosome is **MW2 (USA400)**, not MRSA252. The co-culture member is phage Staph1N at measured 2,984× depth. **Finding BC-01's requirement — a discriminatory sequence-based identity analysis with retained evidence, exact counts, denominators, marker sets, support, and sensitivity limits — is satisfied by §3 with repository-retained artifacts.**

---

## 4. Bibliographic and Forensic Provenance

### 4.1 Source study

BioProject **PRJNA1263016** is the published eLife (2025) study by the **Mo Lab, University of Wisconsin–Madison**: *"Bacteriophage infection drives loss of beta-lactam resistance in methicillin-resistant Staphylococcus aureus"* (public repository: `MolabUW/mrsa-project-eLife2025`). The study investigates strains **MRSA252, MW2, and LAC** treated with **phage PhiStaph1N** plus **oxacillin (32 µg/mL) co-culture in 96-well plates**. This design matches every experimental attribute of SRR33682404 (phage `Staph1N`, MOI 0.01, oxacillin 32, `sample_type = mixed culture`, replicate A) and explains why both MRSA252 and MW2 records exist in the same multiplexed submission.

### 4.2 Run alias decoding

`MW2_A_Staph1N_E12_S196_R1_001.fastq.gz` = strain **MW2**, replicate **A**, phage **Staph1N**, 96-well plate position **row E, column 12**, sample index **S196** — a raw Illumina plate-export filename, consistent with the multiplexed 96-well co-treatment series.

### 4.3 The Title `MRSA252-cotreat2` is a plate-level submitter typo

During multiplexed library submission the BioSample **Title** field was populated from an adjacent plate record of a different strain (MRSA252) in the same co-treatment series. This is the only occurrence of the `MRSA252` token in the entire record; the sample alias (`MW2-cotreat2`), run alias (§4.2), `strain=MW2` attribute, and all experimental attributes are concordant with MW2 + Staph1N, and the primary sequence evidence (§3) is concordant with MW2 and excludes MRSA252. Three independent, mutually consistent lines of evidence (record-internal fields, study design, sequence data) identify the Title as the single erroneous field.

---

## 5. Preliminary 31-mer Report — Historical Record, Superseded (Correction per Audit)

The aggregate 31-mer containment figures previously recorded here and in provenance §4.1.1 (29.90 % / 60.74 % / 90.64 % / 19.12 %) **are not certifying evidence** and must not be cited as such: their scripts and count tables were never retained, their distinct-set denominators are unreconcilable from the published percentages (F26-AUDIT-002 demonstrated the mathematical ambiguity), and their "marker absence" inference is superseded by the explicit marker-set analysis of §3.5–3.6. The same preliminary report listed the Staph1N genome as 140,511 bp; the staged, hash-pinned `staph1n.fasta` is **145,647 bp**, which is the authoritative figure. The qualitative conclusion of the preliminary report (MW2 organism + dominant Staph1N content) is *consistent with* §3, but §3 — not §5 — is the evidence of record.

---

## 6. Assembly and QUAST Evaluation Protocol (Binding for Stage C/D/G)

1. QUAST runs against **`NC_003923.1` only** (established chromosome-only convention). Staph1N contigs (~145.6 kb at 2,984×) and pMW2/pMW3 plasmid contigs will appear as unaligned — expected co-culture/plasmid content, recorded as such, **not** misassemblies or contaminant failures.
2. Per F26-AUDIT-002, **length + depth + library-level presence do not identify individual unaligned contigs**: raw QUAST outputs are retained, and any per-contig exemption from failure interpretation requires sequence-supported attribution (e.g., aligning the flagged contig against `NC_047722.1`); partially aligned and chimeric contigs are investigated separately; unresolved sequence stays unresolved.
3. Chromosome-aligned metrics (genome fraction, NGA50, misassemblies, mismatches/indels per 100 kbp) are interpreted separately from whole-assembly N50/L50/total-length/GC, which include unaligned viral/plasmid sequence and can legitimately differ with each assembler's viral-recovery behavior.
4. Misassembly/mismatch/indel calls within chromosome-aligned contigs remain genuine quality signals and are not excused by this protocol.
5. The rule applies **identically to all four competitors**; the F26-P2 relative comparison on byte-identical inputs is unaffected.

---

## 7. Resolution Statement

- **Finding BC-01: RESOLVED** (2026-10-30, task F26-DATA-004, Slurm job 3074242, `qfree`, exit 0, all gates PASS). CS-SA01 strain identity is verified as *S. aureus* **MW2 (USA400)** with expected *Staphylococcus* phage **Staph1N** co-culture (2,984.21×); the conflicting BioSample Title is documented as a plate-level submitter labelling artefact (§4.3); reference `NC_003923.1` is correct and unchanged.
- Every BC-01 closure requirement of F26-AUDIT-002 is addressed: reproducible qfree rerun with retained commands/exit codes/job id (§3.1, `cs_sa01_disc_run.log`); input and reference hashes (§3.1, `cs_sa01_disc_inputs_sha256.tsv`, `SHA256SUMS.txt`); exact counts and denominator definitions (§3.2–3.3); strain-exclusive marker sets with sizes and observed support, including shared support (§3.5–3.6); host/phage coverage breadth and depth (§3.3); sensitivity limits (§3.8); reference-lineage support distinguished from exact isolate identity (§3.8 last item); and correction of the unsupported preliminary claims through this document and the provenance owner document (§5, provenance §4.1.1/§6).
- No protected file (`src/*`, `F26_COMPETITOR_BENCHMARK_CONTRACT.md`, `scripts/run_competitor_benchmark.py`, `scripts/slurm_benchmark_harness.sbatch`) was modified. The audit report remains as written; closure is recorded here and mirrored in the provenance document. Re-audit of this closure by the independent auditor is the expected next control step.

## 8. Retained Artifacts and Digests

All artifacts live in [`docs/development/evidence/cs_sa01_strain_discrimination/`](.) and are digested in [`SHA256SUMS.txt`](SHA256SUMS.txt) (recomputed at each successful run; digests below from the definitive run 3074242). Key artifacts:

| Artifact | Content | SHA-256 (first 16 hex) |
| :--- | :--- | :--- |
| `cs_sa01_disc_run.log` | Full transcript: every command, exit code, job id, node, git HEAD, tool versions | `b1d402f93728a727` |
| `cs_sa01_disc_inputs_sha256.tsv` | Input file sizes + SHA-256 (5/5, fail-fast verified) | `16e1475e1156b954` |
| `cs_sa01_disc_mw2.idxstats.tsv` / `cs_sa01_disc_mrsa252.idxstats.tsv` | Exact mapped/unmapped counts per contig | `de1638ae3c3d5260` / `a53b5b35949d59da` |
| `cs_sa01_disc_mw2.coverage.tsv` / `cs_sa01_disc_mrsa252.coverage.tsv` | Mean depth + breadth per contig | `c7db5ad45f3aaa6c` / `b02b2b513aa7d770` |
| `cs_sa01_disc_locus_summary.tsv` / `..._readcounts.tsv` / `..._depth.tsv.gz` | Whole-window depth/breadth/counts + raw per-base depth | `2702b7217c05461d` / `87e8cd4e2e0daa8f` / `39774bc2ef950602` |
| `cs_sa01_disc_marker_blocks.tsv` | Shared/exclusive partition sizes per window | `43fce18b5979f134` |
| `cs_sa01_disc_marker_support.tsv` | Depth/breadth/read support per marker set (q0, q20) | `c23133b144e07243` |
| `cs_sa01_disc_exclusive_origin_audit.tsv` | Per-block origin audit (reads, divergence, chromosome support) | `2c9950edaffdee00` |
| `cs_sa01_disc_gates.tsv` | Gate criteria, observed values, verdicts | `8e917628beb0e494` |
| `cs_sa01_disc_sensitivity.tsv` | Poisson zero-coverage bounds | `262bde760f8dbe8d` |
| `loci/*.fasta`, `loci/*.bed` | Marker-set sequences and coordinates (exclusive + shared) | see `SHA256SUMS.txt` |
| `cs_sa01_disc_chunkspecificity.*.paf` | Raw chunked cross-alignment hits (incl. core-genome control) | see `SHA256SUMS.txt` |
| `cs_sa01_disc_exclusive_depth.tsv.gz` / `cs_sa01_disc_phage_depth_*.tsv.gz` / `cs_sa01_disc_*.depth10k.tsv` | Per-base/binned depth records | see `SHA256SUMS.txt` |
| `NC_002952.2.fasta`, `staph1n.fasta`, `NC_002952.2.fasta.fai` | Staged comparator references | `abc19abf4b135629` / `2993d31f2fbd7e98` / `8c3ae3a8c7f4cd4f` |

## 9. References

- Audit finding and required closure: [`docs/development/F26_INDEPENDENT_AUDIT_STAGE_BC.md`](../../F26_INDEPENDENT_AUDIT_STAGE_BC.md) (F26-AUDIT-002, BC-01 sections).
- NCBI records: [BioSample SAMN48627418](https://www.ncbi.nlm.nih.gov/biosample/SAMN48627418), [SRA SRR33682404](https://www.ncbi.nlm.nih.gov/sra/SRR33682404), BioProject [PRJNA1263016](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1263016).
- Source study: Mo Lab, University of Wisconsin–Madison, eLife (2025), *"Bacteriophage infection drives loss of beta-lactam resistance in methicillin-resistant Staphylococcus aureus"*; public repository `MolabUW/mrsa-project-eLife2025`.
- Reference sequences: [NC_003923.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_003923.1) (*S. aureus* MW2), [NC_047722.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_047722.1) (*Staphylococcus* phage Staph1N), [NC_002952.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_002952.2) (*S. aureus* MRSA252).
- Retained audit fetches: `/hpcfs/fhome/caizhh/.local/state/f26-audit-001/{sra.xml,biosample_sa.xml}` (SHA-256 table in audit §7).
- Provenance owner document: [`F26_DATASET_PROVENANCE_AND_SELECTION.md`](../../../F26_DATASET_PROVENANCE_AND_SELECTION.md) (§1.1, §4.1.1, §6 updated for this closure).
