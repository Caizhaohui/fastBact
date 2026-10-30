# fastBact Phase F25 — Architecture and Acceptance Contract

**Document ID**: `F25_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT`  
**Phase**: F25 — Assembly Correctness, Residual Misassembly Forensics, and Cross-Dataset Robustness  
**Date**: 2026-10-30  
**Controller**: Antigravity CLI / Gemini 3.8 Flash Medium  
**Policy File**: `~/.config/ai-collaboration/POLICY.md` (Multi-AI Collaboration V3.2)  
**Primary Worktree**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f25-correctness-worktree`  
**Feature Branch**: `feature/f25-assembly-correctness`  
**Frozen F24 Baseline Commit**: `92be30ae42e856185231c5b96f1d6c2d0ae54d2f` (tag: `f24-accepted-2026-10-30` at `cfbd00503b3031853bb41441160717d1e5bc3e69`)  
**Frozen F24 Baseline Binary SHA-256**: `53998c592cfe2155fca21a57fa1f39b36f5fc5fb1aeb11db970e760a4742b8fd`  

---

## 1. Executive Summary & Purpose

Phase F25 builds upon the formally accepted Phase F24 (*Provenance-Safe Contiguity Recovery*) foundation. F24 achieved substantial contiguity gains (+155.3% NGA50 on R1-003, +66.8% NGA50 on R1-001) while introducing zero new extensive or local misassemblies and significantly reducing base errors.

Phase F25 focuses on three primary scientific objectives:
1. **F25-P1: Residual Misassembly Root-Cause Investigation**: Rigorous forensic characterization and evidence collection for the 5 residual extensive misassembly events observed in dataset R1-001.
2. **F25-P2: Repeat Resolution and Evidence Diagnostics**: Verification and hardening of repeat boundary transitions, preservation of conservative evidence isolation ($p > b + H$), and diagnostic logging.
3. **F25-P3: Cross-Dataset Robustness**: Verification of algorithmic generalization across independent biological sequencing runs with verified reference provenance.

---

## 2. Frozen F24 Non-Regression Contract

All evaluations on the frozen acceptance datasets (R1-001 and R1-003) must meet or exceed the F24 accepted baseline under identical parameters:

| Metric | R1-001 F24 Baseline | R1-003 F24 Baseline | Non-Regression Invariant |
| :--- | :---: | :---: | :--- |
| **NGA50 (bp)** | 45,092 | 44,702 | No contiguity regression ($NGA50 \ge Baseline$) |
| **Genome Fraction (%)** | 98.115% | 98.060% | $\ge 98.0\%$ |
| **Extensive Misassemblies** | 5 | 0 | Strictly **$\le 5$** on R1-001; strictly **0** on R1-003 |
| **Local Misassemblies** | 0 | 1 | Strictly **0** on R1-001; strictly **$\le 1$** on R1-003 |
| **Mismatches / 100 kbp** | 0.85 | 0.13 | $\le 1.00$ / 100 kbp |
| **Indels / 100 kbp** | 0.53 | 0.00 | $\le 0.60$ / 100 kbp |
| **Computational Variance** | $\sigma^2 = 0.000000$ | $\sigma^2 = 0.000000$ | Strict bit-for-bit determinism across replicates |

### Invariant Rules:
1. **Never sacrifice structural accuracy for contiguity**: speculative graph joining across unverified repeat boundaries is strictly forbidden.
2. **Preserve $p > b + H$ evidence threshold**: conservative evidence isolation in repeat resolution remains normative.
3. **Preserve deterministic assembly**: execution under identical parameters must yield identical graph walks and sequences.
4. **No suppressed or fabricated results**: all QUAST outputs must reflect actual execution against ground truth references.

---

## 3. Scope of Forensic Investigation (F25-P1)

### Residual Misassemblies in R1-001 (F24 Final State):
1. **Events 1 & 2 (Inversion)**: `contig_53` (Ref 1,207,805 / 1,209,602, 1,797 bp inverted core).
2. **Event 3 (Relocation)**: `contig_63` (Ref 889,647 / 274,524, 626 bp IS5 repeat stub).
3. **Event 4 (Relocation)**: `contig_200` (Ref 1,874,189 / 291,276, 126 bp IS1 repeat stub).
4. **Event 5 (Relocation)**: `contig_222` (Ref 1,872,615 / 1,978,503, 126 bp IS1 repeat stub).

Each event will be analyzed across raw reads, stage-level FASTA files (Stage 1 through Stage 6), graph topology, and reference genomic context, and classified formally as:
- `CONFIRMED_ASSEMBLY_ERROR`: An erroneous join introduced by the assembler that is refuted by the read data.
- `EVALUATION_ARTIFACT`: A valid assembly of the sequenced sample that differs from the static linear reference due to natural phase variation/inversion or multi-mapping repeat alignment ambiguity in QUAST.
- `AMBIGUOUS_UNRESOLVED`: Insufficient evidence to establish definitive provenance.

---

## 4. Multi-AI Collaboration & Quality Gates

The phase progresses through standard V3.2 Quality Gates:
- **Gate A**: Static analysis, Rust compilation, zero compiler warnings.
- **Gate B**: Unit tests and regression suites (213/213 release tests minimum).
- **Gate C**: End-to-end synthetic assembly integration tests.
- **Gate D**: Independent algorithm correctness audit (Codex CLI).
- **Gate E**: High-performance Slurm compute execution on `qfree` partition.
- **Gate F**: Real-data paired benchmark evaluations on R1-001, R1-003, and cross-dataset validation.
- **Gate G**: Final scientific acceptance report and handoff documentation.

Final integration into `main` remains reserved for the user (`FINAL_AUTHORIZATION_AUTHORITY`).
