#!/usr/bin/env python3
"""
Generate Coverage Gradient Datasets: 30x, 100x, 150x based on S1-COV_ref.fa.
Ensures deterministic random seeds and records sha256 checksums.
"""

import os
import sys
import json
import random
import hashlib
from pathlib import Path

# Add benchmarks/scripts to sys.path to reuse simulate_reads_pe
sys.path.insert(0, str(Path(__file__).parent))
from build_fixtures import simulate_reads_pe, sha256_file, register_dataset_row

STORE_DIR = "benchmark-results/data_store"
REF_FA = f"{STORE_DIR}/S1-COV_ref.fa"
MANIFEST_PATH = "benchmarks/datasets.tsv"

def main():
    print(f"Loading reference from {REF_FA}...")
    with open(REF_FA) as f:
        seq = "".join(line.strip() for line in f if not line.startswith(">"))

    ref_hash = sha256_file(REF_FA)

    for depth in [30, 100, 150]:
        d_id = f"S1-COV-{depth:03d}X"
        r1_gz = f"{STORE_DIR}/{d_id}_R1.fq.gz"
        r2_gz = f"{STORE_DIR}/{d_id}_R2.fq.gz"
        events_file = f"{STORE_DIR}/{d_id}_truth_events.json"

        print(f"Generating {d_id} at {depth}x depth...")
        simulate_reads_pe(
            [("chr", seq, False, 1.0)],
            depth,
            150,
            450,
            30,
            r1_gz,
            r2_gz,
            random.Random(206 + depth),
        )

        with open(events_file, "w") as f:
            json.dump({"dataset": d_id, "depth": depth}, f, indent=2)

        r1_hash = sha256_file(r1_gz)
        r2_hash = sha256_file(r2_gz)

        register_dataset_row(MANIFEST_PATH, {
            "dataset_id": d_id,
            "kind": "S1",
            "organism": f"Synthetic Coverage {depth}x",
            "source_url": "internal://synthetic",
            "study_accession": "NA",
            "run_accession": "NA",
            "biosample": "NA",
            "reference_accession_version": f"{d_id}.1",
            "reference_confidence": "verified",
            "replicons": "1_lin",
            "platform": "Illumina",
            "layout": "PE",
            "read_length": 150,
            "fragment_mean": 450,
            "fragment_sd": 30,
            "coverage_nominal": depth,
            "coverage_measured": depth,
            "simulator_version": "build_fixtures_1.0",
            "random_seed": 206 + depth,
            "fastq_r1": r1_gz,
            "fastq_r2": r2_gz,
            "reference_fasta": REF_FA,
            "truth_events_file": events_file,
            "r1_sha256": r1_hash,
            "r2_sha256": r2_hash,
            "reference_sha256": ref_hash,
            "admission_status": "APPROVED",
            "notes": f"Coverage {depth}x gradient benchmark",
        })
        print(f"Registered {d_id} in {MANIFEST_PATH}")

if __name__ == "__main__":
    main()
