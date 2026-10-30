# fastBact Phase F25 — Task Ledger

**Phase**: F25 — Assembly Correctness, Residual Misassembly Forensics, and Cross-Dataset Robustness  
**Controller**: Antigravity CLI / Gemini 3.8 Flash Medium  
**Policy**: Multi-AI Collaboration V3.2  
**Worktree**: `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f25-correctness-worktree`  
**Feature Branch**: `feature/f25-assembly-correctness`  
**Base Commit**: `92be30a`  

---

## Task Registry

| Task ID | Goal and Scope | Assigned Agent & Model | Status | Commit / Diff Ref | Test & Review Results | Evidence Path | Remaining Blockers | Next Action |
| :--- | :--- | :--- | :---: | :--- | :--- | :--- | :--- | :--- |
| **F25-T01** | Stage A: Baseline audit & acceptance contract definition | Gemini 3.8 Flash Medium (Controller) | `DONE` | `F25_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md` | Non-regression matrix locked; 213/213 release tests verified on baseline | [F25_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f25-correctness-worktree/F25_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md) | None | Complete |
| **F25-T02** | Stage B: Forensic investigation of 5 residual misassemblies on R1-001 | Gemini 3.8 Flash Medium (Controller) | `DONE` | `F25_ROOT_CAUSE_FORENSICS.md` | All 5 events traced to Stage 1 raw reads; 0 chimeric joins from exSPAnder | [F25_ROOT_CAUSE_FORENSICS.md](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f25-correctness-worktree/F25_ROOT_CAUSE_FORENSICS.md) | None | Complete |
| **F25-T03** | Stage B Audit: Independent review of forensics & classification | Codex CLI (GPT-6.1 Sol Medium) | `DONE` | Audit Session `01a18a77-...` | Evaluated raw FASTQ probes (29/37/73/66/51 reads); confirmed all 5 `EVALUATION_ARTIFACT`; required & verified 2 wording refinements | [F25_AGENT_AUDIT_LOG.md](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f25-correctness-worktree/F25_AGENT_AUDIT_LOG.md) | None | Complete |
| **F25-T04** | Stage C/E: Release test suite & CLI verification on Slurm `qfree` | Gemini 3.8 Flash / Slurm | `DONE` | Slurm Job `3072530` | 213/213 release tests passed; CLI smoke tests OK; binary `7c70e357...` | `slurm-f25-release-3072530.out` | None | Complete |
| **F25-T05** | Stage E: Cross-dataset benchmark validation on Slurm `qfree` | Gemini 3.8 Flash / Slurm | `DONE` | Slurm Job `3072539` | 4 datasets complete: R1-001 & R1-003 identical to F24; R1-001-Val achieves NGA50 54,554 bp with 0 misassemblies; R1-003-Val achieves 0 extensive misassemblies | `assembly_metrics.tsv` | None | Complete |
| **F25-T06** | Stage F: Final acceptance reporting & documentation handoff | Gemini 3.8 Flash Medium (Controller) | `DONE` | Deliverable documents | `F25_FINAL_SCIENTIFIC_ACCEPTANCE_REPORT.md` produced; disposition `F25_READY_FOR_USER_ACCEPTANCE` | Acceptance reports | None | Handoff to User |
