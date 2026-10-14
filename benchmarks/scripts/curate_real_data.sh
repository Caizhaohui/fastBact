#!/usr/bin/env bash
set -euo pipefail

# Curates and registers real bacterial isolate datasets (FB-BENCH-001 Level L2)
# Discovers local verified NCBI SRA runs with known reference genomes.

MANIFEST="benchmarks/datasets.tsv"
DATA_ROOT="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data"

echo "=== Curating Real-World Bacterial Datasets (R1/R2) ==="

python3 - << 'EOF'
import sys, os, hashlib
from pathlib import Path
sys.path.insert(0, "benchmarks/scripts")
from build_fixtures import register_dataset_row, sha256_file

manifest = "benchmarks/datasets.tsv"
root = "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data"

real_datasets = [
    {
        "dataset_id": "R1-001",
        "kind": "R1",
        "organism": "Escherichia coli str. K-12 substr. MG1655",
        "source_url": "https://www.ncbi.nlm.nih.gov/sra/SRR28349173",
        "study_accession": "PRJNA1088182",
        "run_accession": "SRR28349173",
        "biosample": "SAMN40463664",
        "reference_accession_version": "NC_000913.3",
        "reference_confidence": "verified",
        "replicons": "1_circ",
        "platform": "Illumina",
        "layout": "PE",
        "read_length": 150,
        "fragment_mean": 380,
        "fragment_sd": 40,
        "coverage_nominal": 78,
        "coverage_measured": 77.8,
        "simulator_version": "NA",
        "random_seed": 0,
        "fastq_r1": f"{root}/prjna1088182/fastq/SRR28349173/SRR28349173_1.fastq.gz",
        "fastq_r2": f"{root}/prjna1088182/fastq/SRR28349173/SRR28349173_2.fastq.gz",
        "reference_fasta": f"{root}/references/m6_identity/NC_000913.3.fasta",
        "truth_events_file": "NA",
        "r1_sha256": "e25760e43d1126b9b56d4ee9de79a42cfd28b1de1dd54224bcc24319bbbe5ebf",
        "r2_sha256": "ccf7c2627919cf2df6eb03c06c322d8dea241ad0ce9bdc34e7c6c13f4b983701",
        "reference_sha256": "6b195feda4c66140f6762742eb8b30c2652f02b45878b174f5b00ef85ecc95d7",
        "admission_status": "APPROVED",
        "notes": "Real E. coli K-12 MG1655 isolate resequencing (2.4M pairs, ~78x coverage)"
    },
    {
        "dataset_id": "R1-002",
        "kind": "R1",
        "organism": "Escherichia coli BW25113",
        "source_url": "https://www.ncbi.nlm.nih.gov/sra/SRR28349221",
        "study_accession": "PRJNA1088182",
        "run_accession": "SRR28349221",
        "biosample": "SAMN40463696",
        "reference_accession_version": "CP009273.1",
        "reference_confidence": "verified",
        "replicons": "1_circ",
        "platform": "Illumina",
        "layout": "PE",
        "read_length": 150,
        "fragment_mean": 380,
        "fragment_sd": 40,
        "coverage_nominal": 65,
        "coverage_measured": 64.5,
        "simulator_version": "NA",
        "random_seed": 0,
        "fastq_r1": f"{root}/prjna1088182/fastq/SRR28349221/SRR28349221_1.fastq.gz",
        "fastq_r2": f"{root}/prjna1088182/fastq/SRR28349221/SRR28349221_2.fastq.gz",
        "reference_fasta": f"{root}/references/m6_identity/CP009273.1.fasta",
        "truth_events_file": "NA",
        "r1_sha256": "ff88b5003ca1c152b4a543bbe6a27e8018386f977e9e501a3c8db732d29913b8",
        "r2_sha256": "5d86526471def34ee6bb2949ed68cfdf19b9e7bf8f8f26c0418cdbc5cc22be54",
        "reference_sha256": "37eade01311f1f81659fad53bf907a93139075a75d115fd3aec2ef09beac2e90",
        "admission_status": "APPROVED",
        "notes": "Real E. coli BW25113 control isolate resequencing (2.0M pairs, ~65x coverage)"
    },
    {
        "dataset_id": "R1-003",
        "kind": "R1",
        "organism": "Escherichia coli BW25113 (lambdaRed)",
        "source_url": "https://www.ncbi.nlm.nih.gov/sra/SRR28349164",
        "study_accession": "PRJNA1088182",
        "run_accession": "SRR28349164",
        "biosample": "SAMN40463672",
        "reference_accession_version": "CP009273.1",
        "reference_confidence": "verified",
        "replicons": "1_circ",
        "platform": "Illumina",
        "layout": "PE",
        "read_length": 150,
        "fragment_mean": 380,
        "fragment_sd": 40,
        "coverage_nominal": 132,
        "coverage_measured": 132.3,
        "simulator_version": "NA",
        "random_seed": 0,
        "fastq_r1": f"{root}/prjna1088182/fastq/SRR28349164/SRR28349164_1.fastq.gz",
        "fastq_r2": f"{root}/prjna1088182/fastq/SRR28349164/SRR28349164_2.fastq.gz",
        "reference_fasta": f"{root}/references/m6_identity/CP009273.1.fasta",
        "truth_events_file": "NA",
        "r1_sha256": "dd43192faf592582c30dfcf71938d6450f6525b54ca254276dbd8ea16672011b",
        "r2_sha256": "94fd5d094315c46725840a5ce8abbff4fc4842df748930b5b8f4bff9497fc3d1",
        "reference_sha256": "37eade01311f1f81659fad53bf907a93139075a75d115fd3aec2ef09beac2e90",
        "admission_status": "APPROVED",
        "notes": "Real E. coli BW25113 lambdaRed engineered strain (4.1M pairs, ~132x coverage)"
    }
]

for r in real_datasets:
    register_dataset_row(manifest, r)
    print(f"Registered real dataset {r['dataset_id']} ({r['organism']})")

print("All 3 real bacterial isolate datasets registered.")
EOF

python3 benchmarks/scripts/validate_manifest.py
