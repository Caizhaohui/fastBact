# fastBact Phase F26 — Architecture and Scientific Acceptance Contract

**Document ID:** `FASTBACT-F26-CONTRACT-01`  
**Date:** 2026-10-30  
**Phase:** F26 — Cross-Species Assembly Benchmarking, Robustness Hardening & Competitive Validation  
**Controller:** Antigravity CLI / Gemini 3.8 Flash Medium  
**Governing Policy:** Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)  
**Base Commit:** `02872f81811402d2bd603ff7c0a2ecf8b3dbe528` (main)  
**Feature Branch:** `feature/f26-cross-species-robustness`  
**Worktree:** `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f26-robustness-worktree`  

---

## 1. Mission and Scientific Objectives

The overarching mission of Phase F26 is to transition fastBact from an *Escherichia coli*-validated assembler into a robust, cross-species bacterial genome assembler systematically evaluated against industry-standard assemblers (SPAdes, SKESA, and Unicycler).

### F26-P1: Cross-Species Benchmark Dataset Expansion (Priority P0)
- Expand the benchmark evaluation beyond *E. coli* to at least three non-*E. coli* bacterial species spanning at least two distinct genera, encompassing both Gram-positive and Gram-negative representatives.
- Candidate genera / species:
  - Gram-positive: *Staphylococcus aureus*, *Bacillus subtilis*
  - Gram-negative: *Pseudomonas aeruginosa*, *Salmonella enterica*
- Rigorous Provenance Requirement:
  - Each dataset must possess verified SRA/ENA accession, BioSample metadata, Illumina paired-end layout (2x150 bp or equivalent), and a complete reference genome representing the sequenced isolate or verified reference strain.
  - Full traceability: record raw FASTA/FASTQ hashes and normalized sequence hashes.

### F26-P2: Reproducible Competitive Benchmarking (Priority P0)
- Evaluate fastBact head-to-head against established bacterial assemblers:
  - **SPAdes v4.2.0** (`/hpcfs/fhome/caizhh/.conda/envs/BactGenome/bin/spades.py`)
  - **SKESA v2.5.1** (`/hpcfs/fhome/caizhh/.conda/envs/LineaBact/bin/skesa`)
  - **Unicycler v0.5.1** (`/hpcfs/fhome/caizhh/.conda/envs/BactGenome/bin/unicycler`)
- Short-read only mode (Illumina PE); identical CPU allocations (8 threads) and memory caps.
- Consistent evaluation harness: QUAST v5.2.0 (`--min-contig 200`) against strain-specific reference genomes.
- Metrics collected: N50, NGA50, Genome Fraction (%), Raw Extensive Misassemblies, Local Misassemblies, Mismatches/100kb, Indels/100kb, Wall Time, Peak RSS.

### F26-P3: Coverage Sensitivity and Robustness (Priority P1)
- Evaluate assembly behavior across depth gradients (e.g., 30x, 50x, 100x, 150x, and full coverage).
- Pair-aware deterministic subsampling with recorded random seeds.
- Characterize graph fragmentation, repeat resolution behavior, memory scaling, and runtime efficiency.

### F26-P4: Evidence-Driven Algorithm Hardening (Priority P1)
- Identify reproducible assembly bottlenecks or defects discovered during cross-species benchmarking.
- **Mandatory Invariant**: Zero benchmark-driven speculative patching. Any algorithmic repair must:
  1. Identify a reproducible failure case.
  2. Isolate the responsible assembly stage.
  3. Formulate a minimal failing unit/regression test.
  4. Preserve conservative evidence isolation ($p > b + H$) and loop re-entrancy invariants.
  5. Undergo independent audit before merge.

---

## 2. Frozen Baseline Non-Regression Contract

The accepted F25 baseline results across the four historical datasets are strictly frozen as the non-regression oracle:

| Dataset | Strain | SRA Accession | Reference FASTA | F25 NGA50 | F25 Genome Fraction | Raw Misassemblies | Classified Artifacts | Local Misasm | Mismatches/100kb | Indels/100kb |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| **R1-001** | *E. coli* MG1655 | `SRR28349173` | `NC_000913.3` | 45,092 bp | 98.115% | 5 | 2 bio + 3 repeat | 0 | 0.85 | 0.53 |
| **R1-003** | *E. coli* BW25113 | `SRR28349164` | `CP009273.1` | 44,702 bp | 98.060% | 0 | 0 | 1 | 0.13 | 0.00 |
| **R1-001-Val** | *E. coli* MG1655 | `SRR28349174` | `NC_000913.3` | 54,554 bp | 98.233% | 0 | 0 | 0 | 0.81 | 0.24 |
| **R1-003-Val** | *E. coli* BW25113 | `SRR28349165` | `CP009273.1` | 39,314 bp | 98.097% | 0 | 0 | 1 | 0.26 | 0.00 |

### Non-Regression Invariants:
1. Under identical inputs and software environment, F25 baseline datasets must retain $\Delta_{\text{metric}} = 0$ and $\sigma^2 = 0.000000$ unless an intentional, scientifically justified algorithm improvement is audited and accepted.
2. Zero new confirmed extensive or local misassemblies.
3. No silent regression in base-level accuracy (mismatches / indels).
4. No unreviewed loosening of repeat-joining guardrails ($p > b + H$).

---

## 3. Stage Gates and Execution Protocol

- **Stage A**: Baseline and Environment Verification (Contract, Worktree, Tool Discovery).
- **Stage B**: Dataset Discovery, Provenance Verification, and Pre-Registration.
- **Stage C**: Unified Benchmark Harness Construction (fastBact, SPAdes, SKESA, Unicycler).
- **Stage D**: Preliminary Cross-Species Baseline Evaluation (Evaluating frozen baseline on new species).
- **Stage E**: Coverage Robustness & Sensitivity Analysis.
- **Stage F**: Evidence-Driven Algorithmic Hardening (if defects or justified improvements are discovered).
- **Stage G**: Final Multi-Tool Cross-Species A/B Benchmark & Determinism Verification.
- **Stage H**: Independent Scientific Acceptance Audit (Codex Sol).
- **Stage I**: Handoff Preparation and User Authorization Request (`F26_READY_FOR_USER_ACCEPTANCE`).
