# fastBact Phase F26: Worker Dispatch & Delegation Ledger

**Document ID:** `F26_WORKER_DISPATCH_LEDGER`  
**Phase:** F26 — Cross-Species Assembly Benchmarking, Robustness Hardening & Competitive Validation  
**Operating Mode:** `DELEGATION_FIRST_STRICT` (Hotfix Applied)  
**Controller:** Antigravity CLI / Gemini 3.8 Flash Medium (Orchestration & Coordination Only)  
**Policy Location:** `~/.config/ai-collaboration/POLICY.md`  
**Worktree:** `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f26-robustness-worktree`  
**Branch:** `feature/f26-cross-species-robustness`  
**Base Commit:** `02872f81811402d2bd603ff7c0a2ecf8b3dbe528` (main)

---

## 1. Operating Policy & Controller Firewall

Under Multi-AI Collaboration Policy V3.2 and the F26 Hotfix directive, the Gemini Controller is strictly constrained to orchestration, dependency tracking, worker dispatch, and verification.
All substantial activities—including dataset discovery, curation, validation scripts, Rust code, benchmark harnesses, and audits—are performed exclusively by external CLI workers:
- **OpenCode CLI** (`zhipuai-coding-plan/glm-5.3-flash`): Primary for dataset curation, metadata provenance, download/validation scripts.
- **Cursor CLI (`cursor-agent`)** (`grok-4.7-high` / `grok-4.7-high-fast`): Primary for benchmark harness implementation and Rust production/algorithm engineering.
- **Codex CLI** (`gpt-6.1-sol`): Primary for independent code review and scientific acceptance audits.
- **Claude CLI**: Fallback for high-complexity architecture design and difficult repairs.

---

## 2. Retrospective Direct Controller Action Log (Transition Record)

Prior to the F26 Hotfix enforcement, the Controller performed the following initial actions directly:
1. Environment verification, worktree initialization, and baseline commit freeze verification (`02872f8`).
2. Competitor executable discovery on HPC (`spades.py`, `unicycler`, `skesa`, `quast.py`).
3. Initial Entrez query draft and initial draft reference downloads.
4. Drafting of architecture contract and task ledgers.

*Transition Status:* No transfers are running. All subsequent research, data pipeline development, harness code, and audits are delegated to external workers under the pre-dispatch protocol below.

---

## 3. Worker Dispatch Registry

| Dispatch ID | Task ID | Worker CLI | Model | Role / Objective | Owned Files | Status | Session / Log Identifier |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **DISP-F26-001** | `F26-DATA-001` | `opencode` | `glm-5.3-flash` | Cross-Species Dataset Discovery & Provenance Audit | `F26_DATASET_PROVENANCE_AND_SELECTION.md`, `F26_DATASET_MANIFEST.tsv` | **COMPLETED** (Commit `87e54c1`) | `logs/disp_f26_001_opencode.log` |
| **DISP-F26-002** | `F26-BM-001` | `cursor-agent` | `grok-4.7-high` | Unified Competitor Benchmark Harness Implementation | `scripts/run_competitor_benchmark.py`, `scripts/slurm_benchmark_harness.sbatch`, `F26_COMPETITOR_BENCHMARK_CONTRACT.md` | **COMPLETED** (Commit `e11b88c`) | `logs/disp_f26_002_cursor.log` |
| **DISP-F26-003** | `F26-DATA-002` | `opencode` | `glm-5.3-flash` | Robust SRA Download & Checksum Validation Pipeline | `scripts/download_f26_datasets.sh`, `scripts/validate_fastq_pairs.py` | **COMPLETED** (Commit `a2a50c1`) | `logs/disp_f26_003_opencode.log` |
| **DISP-F26-004** | `F26-EXEC-001`| `slurm` (Job 3073557) | HPC `vnode1` | Execute Cross-Species FASTQ Downloads (`scripts/download_f26_datasets.sh`) | `data/f26_cross_species/fastq/*` | **COMPLETED** (Slurm 0:0, 100% spots) | `logs/download_f26_3073557.log` |
| **DISP-F26-005** | `F26-VAL-001` | `scripts/validate_fastq_pairs.py` | Local execution | Post-Download Pair Synchronization & Checksum Audit | `data/f26_cross_species/fastq/validation_report.json` | **COMPLETED** (4/4 PASS, 0 mismatches) | `data/f26_cross_species/fastq/validation_report.json` |
| **DISP-F26-006** | `F26-AUDIT-001`| `codex` | `gpt-6.1-sol` | Independent Audit of F26 Contracts, Datasets & Pipelines | `docs/development/F26_INDEPENDENT_AUDIT_STAGE_BC.md` | **COMPLETED** (Commit `cd332b5`, Verdict: BLOCKED) | `docs/development/F26_INDEPENDENT_AUDIT_STAGE_BC.md` |
| **DISP-F26-007** | `F26-BM-002` | `cursor-agent` | `grok-4.7-high` | Repair Harness Findings BC-02 (Wall-time) & BC-03 (Determinism) | `scripts/run_competitor_benchmark.py`, `scripts/slurm_benchmark_harness.sbatch`, `F26_COMPETITOR_BENCHMARK_CONTRACT.md` | **COMPLETED** (Commit `ae6ca39`, Self-Test PASS) | `logs/disp_f26_007_cursor.log` |
| **DISP-F26-008** | `F26-DATA-003` | `opencode` | `glm-5.3-flash` | Resolve Finding BC-01 via Empirical 31-mer Containment Evidence | `F26_DATASET_PROVENANCE_AND_SELECTION.md`, `docs/development/evidence/CS_SA01_STRAIN_DISCRIMINATION_EVIDENCE.md` | **COMPLETED** (Commit `d6f21b9`, BC-01 Resolved) | `logs/disp_f26_008_opencode.log` |
| **DISP-F26-009** | `F26-AUDIT-002`| `codex` | `gpt-6.1-sol` | Stage B-C Re-Audit and Verification of BC-01, BC-02, BC-03 | `docs/development/F26_INDEPENDENT_AUDIT_STAGE_BC.md` | **COMPLETED** (Commit `bf19db7`, BC-02/03 CLOSED, BC-01 retained) | `docs/development/F26_INDEPENDENT_AUDIT_STAGE_BC.md` |
| **DISP-F26-010** | `F26-DATA-004` | `opencode` | `glm-5.3-flash` | Primary Alignment Forensics & Diagnostic Marker Resolution for BC-01 | `scripts/run_cs_sa01_discrimination.sh`, `docs/development/evidence/CS_SA01_STRAIN_DISCRIMINATION_EVIDENCE.md`, `docs/development/evidence/cs_sa01_strain_discrimination/*` | **COMPLETED** (Commit `f7560a6`, Exit 0, Gates G1-G6 PASS) | `docs/development/evidence/cs_sa01_strain_discrimination/` |
| **DISP-F26-011** | `F26-AUDIT-003`| `codex` | `gpt-6.1-sol` | Stage B-C Final Audit Verification & Certification of BC-01 Closure | `docs/development/F26_INDEPENDENT_AUDIT_STAGE_BC.md` | **COMPLETED** (Commit `1d050b8`, Verdict: STAGE_BC_AUDIT_PASS) | `docs/development/F26_INDEPENDENT_AUDIT_STAGE_BC.md` |
| **DISP-F26-012** | `F26-BM-003` | `cursor-agent` | `grok-4.7-high` | Stage D Batch Benchmark Submitter & Results Collation Utility | `scripts/collate_benchmark_results.py`, `scripts/submit_stage_d_benchmarks.sh` | **COMPLETED** (Commit `968b06f`, Self-Test PASS) | `logs/disp_f26_013_cursor.log` |
| **DISP-F26-013** | `F26-EXEC-002`| `slurm` (Jobs 3074506-3074509) | HPC `qfree` | fastBact Stage D Replicate 1 Benchmarks (CS-SA01, CS-BS01, CS-PA01, CS-SE01) | `benchmarks/f26_results/raw/*_fastbact_r1/*` | **COMPLETED** (4/4 PASS, Genome Frac >98.5%, Misassemblies <=2) | `benchmarks/f26_results/summary.tsv` |
| **DISP-F26-014** | `F26-EXEC-003`| `slurm` (Jobs 3074516-3074527) | HPC `qfree` | Competitor Stage D Replicate 1 Benchmarks (SKESA, SPAdes, Unicycler) | `benchmarks/f26_results/raw/*` | **RUNNING** | Slurm Queue (Jobs 3074516-3074527) |
| **DISP-F26-015** | `F26-BM-005` | `cursor-agent` | `grok-4.7-high` | Stage E Deterministic FASTQ Pair Subsampling Utility | `scripts/subsample_fastq_pairs.py` | **COMPLETED** (Commit `8ecfac7`, Self-Test PASS) | `logs/disp_f26_015_cursor.log` |
| **DISP-F26-016** | `F26-BM-006` | `cursor-agent` | `grok-4.7-high` | Stage E Subsample Generation Pipeline & Manifest Builder | `scripts/generate_stage_e_subsamples.py` | **COMPLETED** (Commit `06c1732`, Self-Test PASS) | `logs/disp_f26_016_cursor.log` |
| **DISP-F26-017** | `F26-EXEC-004`| `slurm` (Job 3075001) | HPC `qfree` | Execute Stage E Subsample Generation & Manifest Generation | `data/f26_cross_species/subsamples/*`, `F26_STAGE_E_MANIFEST.tsv` | **RUNNING** | `logs/f26_stage_e_subsamples_3075001.out` |
| **DISP-F26-018** | `F26-REP-001`  | `cursor-agent` | `grok-4.7-high` | Author Stage D Cross-Species Baseline Benchmark Report | `docs/development/F26_CROSS_SPECIES_BASELINE_REPORT.md` | **DISPATCHED** | `logs/disp_f26_018_cursor.log` |

---

## 4. Dispatch Details & Execution Evidence

### DISP-F26-001: F26-DATA-001 (Dataset Discovery & Provenance Audit)
- **Task ID:** `F26-DATA-001`
- **Worker:** OpenCode (`~/.opencode/bin/opencode`)
- **Model:** `zhipuai-coding-plan/glm-5.3-flash`
- **Objective:** Perform rigorous scientific and bibliographic audit of the cross-species candidate datasets (CS-SA01, CS-BS01, CS-PA01, CS-SE01), verify exact strain identity, biological relevance (GC content, Gram status, genome size), reference accession compatibility, and SRA library layout.
- **Inputs:** `F26_DATASET_MANIFEST.tsv`, `F26_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md`, NCBI Entrez.
- **Owned Files:** `F26_DATASET_PROVENANCE_AND_SELECTION.md`, `F26_DATASET_MANIFEST.tsv`
- **Protected Files:** `src/*`, historical baseline references in `data/references/m6_identity/`.
- **Acceptance Criteria:** Document contains verified SRA accessions, exact organism strains, platform/read length, spot counts, reference genome accession numbers, download size estimates, and exclusion rationale for disqualified candidates.

### DISP-F26-002: F26-BM-001 (Unified Benchmark Harness Construction)
- **Task ID:** `F26-BM-001`
- **Worker:** Cursor Agent (`~/.local/bin/cursor-agent`)
- **Model:** `grok-4.7-high`
- **Objective:** Implement an automated, manifest-driven benchmark harness executing fastBact (`02872f8`), SPAdes v4.2.0, SKESA v2.5.1, and Unicycler v0.5.1, followed by QUAST v5.2.0 evaluation (`--min-contig 200`), with resource tracking (cgroup/time/RSS) and JSON/TSV metric collation.
- **Inputs:** `F26_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md`, competitor paths (`/hpcfs/fhome/caizhh/.conda/envs/BactGenome/bin/spades.py`, etc.).
- **Owned Files:** `scripts/run_competitor_benchmark.py`, `scripts/slurm_benchmark_harness.sbatch`, `F26_COMPETITOR_BENCHMARK_CONTRACT.md`
- **Protected Files:** `src/*`, baseline reference files, `data/*`.
- **Acceptance Criteria:** Runnable python/sbatch harness accepting `--manifest`, `--dataset`, `--tools`, outputting structured metrics to `benchmarks/f26_results/`, handling Slurm resource constraints gracefully.
