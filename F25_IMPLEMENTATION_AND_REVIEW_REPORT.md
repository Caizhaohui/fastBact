# fastBact Phase F25 — Implementation and Review Report

**Document ID:** `FASTBACT-F25-IMPL-REV-01`  
**Date:** 2026-10-30  
**Status:** COMPLETE / AUDITED  
**Feature Branch:** `feature/f25-assembly-correctness`  
**Base Commit:** `92be30a38b5779c1da73f27f05c31718ef5dc059`  
**Target Release Binary:** SHA-256 `7c70e3572996ec777b62b8be0080cd8a0a1528ab01f03d0d4b9f3f56640a8112`  
**Author / Controller:** Antigravity CLI / Gemini 3.8 Flash Medium  
**Independent Acceptance Reviewer:** Codex CLI / `gpt-6.1-sol` (Medium Reasoning)  
**Governing Policy:** Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)  

---

## 1. Executive Summary

Phase F25 focuses on **Assembly Correctness, Residual Misassembly Forensics, and Cross-Dataset Robustness**. Following the successful resolution of the 2.1 Mb macro-chimera (AUD-01) in Phase F24 via conservative read-pair evidence isolation ($p > b + H$) and loop re-entrancy invariants, five residual extensive misassemblies remained in dataset R1-001.

In accordance with Phase F25 objectives:
1. **Forensic Investigation (F25-P1)**: The five residual misassemblies were traced across all assembly stages (Stage 1 De Bruijn graph unitigs through Stage 6 scaffold output) and evaluated against raw FASTQ reads and reference annotations. All five events were determined to be `EVALUATION_ARTIFACT`:
   - Two events (`contig_53` Left/Right) are caused by biological phase variation (active Pin-mediated DNA inversion in the e14 prophage) physically coexisting in the sequencing library.
   - Three events (`contig_63`, `contig_200`, `contig_222`) are caused by contigs entering multi-copy insertion sequences (IS5 and IS1) and safely terminating at the repeat boundary without false bridging. QUAST arbitrarily aligns the repeat-terminal stubs to distant genomic copies, triggering artificial extensive relocation misassembly calls.
2. **Codebase Stability & Preservation (F25-P2)**: The repeat resolution invariants implemented in F24 were verified to be strictly intact. Zero algorithmic regressions or speculative heuristic loosenings were introduced.
3. **Independent Acceptance Review**: An independent scientific code and forensic audit was conducted by Codex CLI (`gpt-6.1-sol`). Codex confirmed all forensic findings through direct raw read-pair probe scanning on compute node `vnode2`, verified the biological reality of the e14 inversion, validated the repeat stub mechanics, and concurred with the `EVALUATION_ARTIFACT` classifications.
4. **Deterministic Verification**: Full test suite (213/213 tests) passed cleanly under release mode on Slurm compute node `vnode1` (`qfree`).

---

## 2. Forensic Analysis & Classification Summary

| Event ID | Contig | Event Type | Coordinates | Reference Target | Direct Cause | Classification | Codex Finding |
|---|---|---|---|---|---|---|---|
| **E1** | `contig_53` | Inversion | 22,233–22,242 | 1,195,655 | Active *pin* inversion in e14 prophage | `EVALUATION_ARTIFACT` | Confirmed by raw reads |
| **E2** | `contig_53` | Relocation | 24,015–24,028 | 1,197,450 | Opposite junction of e14 inversion | `EVALUATION_ARTIFACT` | Confirmed by raw reads |
| **E3** | `contig_63` | Relocation | 12,042–12,043 | 888,970 | Terminal 626 bp enters IS5 (10 copies) | `EVALUATION_ARTIFACT` | Repeat stub artifact |
| **E4** | `contig_200`| Relocation | 2,752–2,753 | 4,008,790 | Terminal 126 bp ($k-1$) enters IS1 (6 copies) | `EVALUATION_ARTIFACT` | Repeat stub artifact |
| **E5** | `contig_222`| Relocation | 2,752–2,753 | 1,980,578 | Terminal 126 bp ($k-1$) enters IS1 (6 copies) | `EVALUATION_ARTIFACT` | Repeat stub artifact |

### 2.1 Biological Phase Variation (`contig_53`)
The reference *E. coli* K-12 MG1655 (`NC_000913.3`) represents one static isomeric state of the invertible P segment in the e14 prophage. fastBact faithfully assembled the alternative isomeric state present in the bacterial culture. Raw read pair scanning confirmed that both isomers are physically abundant in SRR28349173:
- Inverted junction: 29 reads (left), 37 reads (right).
- Reference junction: 18 reads (left), 14 reads (right).
Assembly of a true physical biological variant does not constitute an algorithmic misassembly.

### 2.2 Repeat Boundary Stubs (`contig_63`, `contig_200`, `contig_222`)
In all three cases:
- The assembler correctly traversed unique sequence into an insertion sequence repeat.
- Upon encountering copy-number ambiguity ($N_{\text{copies}} \ge 5$) with no unambiguous spanning mate pairs, the assembler terminated traversal safely, preserving contig integrity.
- QUAST evaluated the terminal stub ($k-1$ or repeat fragment) by arbitrarily selecting one reference copy. When this copy differed from the upstream genomic context, QUAST reported an extensive relocation.
- Neither false cross-repeat joins nor chimeric contigs occurred.

---

## 3. Independent Audit Review (Codex CLI / `gpt-6.1-sol`)

Codex CLI conducted an independent audit on 2026-10-30 under session `01a18a77-3e33-7221-a477-dff453eb5f23`.

### 3.1 Verification Methodology
Codex executed standalone Python forensic probes on Slurm compute node `vnode2`, scanning all 2.4 million read pairs of `SRR28349173`:
1. Extracted 70-bp junction probes flanking both the reference and inverted orientations of the e14 prophage.
2. Counted exact read matches for both isoforms.
3. Searched for terminal repeat sequences of `contig_63`, `contig_200`, and `contig_222` across `NC_000913.3`.
4. Traced stage-by-stage assembly logs from Stage 1 through Stage 6.

### 3.2 Audit Findings & Controller Resolution
Codex endorsed the scientific classification of all five events as `EVALUATION_ARTIFACT`, while raising three precision corrections on the initial report:
1. **Finding A (`contig_63` Extension Provenance)**:
   - *Audit observation*: `contig_63` gained 504 bp at its terminal IS5 repeat end during Stage 4 exSPAnder repeat resolution, in addition to 11,460 bp at its upstream unique end.
   - *Resolution*: Forensics report Section 3.3 updated to reflect the exact 504 bp extension and verify that the extension remained within the IS5 sequence without false bridging.
2. **Finding B (`contig_53` Precursor Citation)**:
   - *Audit observation*: Stage 1 `unitig_85574` was initially cited for the left inversion junction; however, Stage 4 `unitig_4190` was the true precursor spanning the junction.
   - *Resolution*: Forensics report Section 3.1 corrected with exact Stage 4 precursor `unitig_4190`.
3. **Finding C (Scope of Chimera Claim)**:
   - *Audit observation*: The claim "zero chimeric joins exist anywhere" exceeded the scope of the five investigated misassemblies.
   - *Resolution*: Scoped to: "None of the five residual misassembly events in R1-001 represents a chimeric join or an algorithm failure."

All required corrections were integrated into [`F25_ROOT_CAUSE_FORENSICS.md`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f25-correctness-worktree/F25_ROOT_CAUSE_FORENSICS.md).

---

## 4. Test Suite and Determinism Verification

On Slurm compute node `vnode1` (`qfree`), execution of the full test suite yielded:
- **Command**: `cargo test --release -- --nocapture`
- **Result**: `test result: ok. 213 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out`
- **Elapsed time**: 35.81s
- **Synthetic integration test**: 10,000 bp synthetic genome assembled with 100% genome fraction, 0 misassemblies, elapsed 1.29s.
- **CLI smoke test**: Validated `--help`, parameter parsing, and single-contig output.
- **Binary Checksum**: `7c70e3572996ec777b62b8be0080cd8a0a1528ab01f03d0d4b9f3f56640a8112`

---

## 5. Non-Regression and Safety Invariant Compliance

1. **Conservative Evidence Isolation**: Invariant $p > b + H$ is strictly preserved in `src/scaffolding/exspander.rs`.
2. **Deterministic Traversal**: All tie-breaking, graph exploration, and path selections use deterministic sorting and canonical k-mer keys.
3. **Zero Weakening of Guardrails**: No filters or thresholds were relaxed to artificially suppress misassembly counts.
4. **Contract Status**: All requirements of `F25_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md` are fulfilled.
