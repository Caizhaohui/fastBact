#!/usr/bin/env python3
"""
Manifest Validator for fastBact Automated Benchmark Suite (FB-BENCH-001)
Validates schema, column contracts, checksums, and admission status of benchmarks/datasets.tsv.
"""

import sys
import os
import hashlib
from pathlib import Path

MANDATORY_COLUMNS = [
    "dataset_id", "kind", "organism", "source_url", "study_accession",
    "run_accession", "biosample", "reference_accession_version", "reference_confidence",
    "replicons", "platform", "layout", "read_length", "fragment_mean",
    "fragment_sd", "coverage_nominal", "coverage_measured", "simulator_version",
    "random_seed", "fastq_r1", "fastq_r2", "reference_fasta", "truth_events_file",
    "r1_sha256", "r2_sha256", "reference_sha256", "admission_status", "notes"
]

VALID_KINDS = {"S0", "S1", "S2", "R1", "R2", "N"}
VALID_CONFIDENCES = {"verified", "qualified", "unavailable"}
VALID_ADMISSIONS = {"APPROVED", "PENDING", "REJECTED"}

def compute_sha256(filepath: str) -> str:
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()

def validate_manifest(tsv_path: str, check_files: bool = True) -> bool:
    if not os.path.exists(tsv_path):
        print(f"[ERROR] datasets.tsv not found at: {tsv_path}")
        return False

    with open(tsv_path, "r", encoding="utf-8") as f:
        lines = [l.strip() for l in f if l.strip()]

    if not lines:
        print(f"[ERROR] {tsv_path} is empty.")
        return False

    header = lines[0].split("\t")
    if header != MANDATORY_COLUMNS:
        print("[ERROR] Header mismatch!")
        print(f"  Expected ({len(MANDATORY_COLUMNS)}): {MANDATORY_COLUMNS}")
        print(f"  Found ({len(header)}): {header}")
        diff = set(MANDATORY_COLUMNS) ^ set(header)
        print(f"  Differences: {diff}")
        return False

    has_errors = False
    row_count = 0

    for idx, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != len(MANDATORY_COLUMNS):
            print(f"[ERROR Line {idx}] Expected {len(MANDATORY_COLUMNS)} columns, got {len(fields)}")
            has_errors = True
            continue

        row = dict(zip(MANDATORY_COLUMNS, fields))
        row_count += 1
        d_id = row["dataset_id"]

        # Validate Kind
        if row["kind"] not in VALID_KINDS:
            print(f"[ERROR Line {idx} ({d_id})] Invalid kind '{row['kind']}'. Must be one of {VALID_KINDS}")
            has_errors = True

        # Validate reference_confidence
        if row["reference_confidence"] not in VALID_CONFIDENCES:
            print(f"[ERROR Line {idx} ({d_id})] Invalid reference_confidence '{row['reference_confidence']}'")
            has_errors = True

        # Validate admission_status
        if row["admission_status"] not in VALID_ADMISSIONS:
            print(f"[ERROR Line {idx} ({d_id})] Invalid admission_status '{row['admission_status']}'")
            has_errors = True

        # Check numeric NA rules
        for num_col in ["read_length", "fragment_mean", "fragment_sd", "coverage_nominal"]:
            val = row[num_col]
            if val != "NA":
                try:
                    float(val)
                except ValueError:
                    print(f"[ERROR Line {idx} ({d_id})] Column '{num_col}' must be numeric or 'NA', got '{val}'")
                    has_errors = True

        # Check file existence and checksums if check_files enabled
        if check_files and row["admission_status"] == "APPROVED":
            for file_col, hash_col in [("fastq_r1", "r1_sha256"), ("fastq_r2", "r2_sha256"), ("reference_fasta", "reference_sha256")]:
                f_path = row[file_col]
                exp_hash = row[hash_col]
                if f_path != "NA":
                    if not os.path.exists(f_path):
                        print(f"[WARN Line {idx} ({d_id})] File {file_col} not found: {f_path}")
                    elif exp_hash != "NA":
                        actual_hash = compute_sha256(f_path)
                        if actual_hash != exp_hash:
                            print(f"[ERROR Line {idx} ({d_id})] SHA256 mismatch for {file_col} ({f_path})!")
                            print(f"  Expected: {exp_hash}")
                            print(f"  Actual:   {actual_hash}")
                            has_errors = True

    if has_errors:
        print(f"[FAIL] Manifest validation failed for {tsv_path}")
        return False

    print(f"[PASS] Successfully validated {row_count} datasets in {tsv_path}")
    return True

if __name__ == "__main__":
    tsv = sys.argv[1] if len(sys.argv) > 1 else "benchmarks/datasets.tsv"
    check_f = "--no-files" not in sys.argv
    success = validate_manifest(tsv, check_files=check_f)
    sys.exit(0 if success else 1)
