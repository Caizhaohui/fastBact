# fastBact Phase F26 — Task Ledger

**Phase:** F26 — Cross-Species Assembly Benchmarking, Robustness Hardening & Competitive Validation  
**Controller:** Antigravity CLI / Gemini 3.8 Flash Medium  
**Worktree:** `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f26-robustness-worktree`  
**Feature Branch:** `feature/f26-cross-species-robustness`  
**Base Commit:** `02872f81811402d2bd603ff7c0a2ecf8b3dbe528`  
**Policy:** Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)  

---

## Task Registry

| Task ID | Stage & Goal | Responsible Agent | Status | Inputs & Context | Primary Artifacts | Review / Audit Status | Blockers / Notes |
|:---|:---|:---|:---:|:---|:---|:---|:---|
| **F26-T01** | Stage A: Baseline & Environment Verification | Gemini Flash (Controller) | `DONE` | Main `02872f8`, tools check | Worktree created, `F26_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md` | Self-verified | None |
| **F26-T02** | Stage B: Dataset Discovery & Pre-Registration | Gemini Flash (Controller) | `RUNNING` | SRA / ENA public datasets | `F26_DATASET_PROVENANCE_AND_SELECTION.md`, `F26_DATASET_MANIFEST.tsv` | Pending | Target $\ge$ 3 non-*E. coli* species |
| **F26-T03** | Stage C: Benchmark Harness Development | Gemini Flash (Controller) | `TODO` | SPAdes, SKESA, Unicycler, QUAST | `F26_COMPETITOR_BENCHMARK_CONTRACT.md`, run scripts | Pending | Standardize 8 threads, 16GB, QUAST --min-contig 200 |
| **F26-T04** | Stage D: Preliminary Cross-Species Baseline Evaluation | Gemini Flash / Slurm | `TODO` | Frozen F25 binary on new datasets | Raw evaluation logs, baseline cross-species metrics | Pending | Must run before code changes |
| **F26-T05** | Stage E: Coverage Sensitivity & Robustness | Gemini Flash / Slurm | `TODO` | Subsampled depth gradients | `F26_COVERAGE_ROBUSTNESS_REPORT.md` | Pending | 30x, 50x, 100x, 150x, full |
| **F26-T06** | Stage F: Evidence-Driven Hardening (if needed) | Cursor / Grok High Fast | `TODO` | Trace logs from Stage D/E | Code changes & regression tests | Claude / Codex | Strictly evidence-backed |
| **F26-T07** | Stage G: Final Competitive & A/B Benchmark | Gemini Flash / Slurm | `TODO` | fastBact vs SPAdes vs SKESA vs Unicycler | `F26_COMPETITIVE_BENCHMARK_REPORT.md` | Pending | Machine-readable metrics TSV |
| **F26-T08** | Stage H: Independent Acceptance Audit | Codex CLI (`gpt-6.1-sol`) | `TODO` | Candidate diff, logs, datasets | Codex Audit Report & Verdict | Codex Sol High | Mandatory gate |
| **F26-T09** | Stage I: Final Scientific Acceptance & Handoff | Gemini Flash (Controller) | `TODO` | All deliverables | `F26_FINAL_SCIENTIFIC_ACCEPTANCE_REPORT.md` | Pending | Stop at `F26_READY_FOR_USER_ACCEPTANCE` |
