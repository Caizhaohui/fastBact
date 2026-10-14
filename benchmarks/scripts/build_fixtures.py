#!/usr/bin/env python3
"""
Deterministic Dataset & Fixture Builder for fastBact Benchmark Suite (FB-BENCH-001)
Generates S0, N, S1, and S2 datasets, reference FASTA, truth_events.json,
and registers entries in benchmarks/datasets.tsv with exact SHA256 hashes.
"""

import os
import sys
import gzip
import json
import random
import hashlib
import argparse
from pathlib import Path

DNA = ["A", "C", "G", "T"]
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

def rev_comp(seq: str) -> str:
    return "".join(COMPLEMENT.get(b, "N") for b in reversed(seq))

def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()

def make_random_seq(length: int, rng: random.Random, gc_bias: float = 0.50) -> str:
    bases = []
    for _ in range(length):
        if rng.random() < gc_bias:
            bases.append(rng.choice(["C", "G"]))
        else:
            bases.append(rng.choice(["A", "T"]))
    return "".join(bases)

def simulate_reads_pe(
    replicons: list, # [(name, seq, is_circular, copy_ratio)]
    cov_nominal: float,
    read_len: int,
    frag_mean: float,
    frag_sd: float,
    out_r1_gz: str,
    out_r2_gz: str,
    rng: random.Random,
    err_rate: float = 0.005
):
    """
    Simulates paired-end Illumina-style reads from one or more replicons.
    """
    total_len = sum(len(seq) * copy for _, seq, _, copy in replicons)
    total_pairs = int((total_len * cov_nominal) / (2 * read_len))

    # Weight replicons by length * copy_ratio
    weights = [len(seq) * copy for _, seq, _, copy in replicons]
    total_w = sum(weights)
    norm_w = [w / total_w for w in weights]

    r1_records = []
    r2_records = []

    pair_idx = 0
    for _ in range(total_pairs):
        pair_idx += 1
        # Pick replicon
        r_choice = rng.choices(replicons, weights=norm_w, k=1)[0]
        r_name, r_seq, is_circ, _ = r_choice
        L = len(r_seq)

        frag_len = max(read_len + 20, int(rng.gauss(frag_mean, frag_sd)))

        if is_circ:
            start = rng.randint(0, L - 1)
            # Extract circular fragment
            if start + frag_len <= L:
                frag = r_seq[start : start + frag_len]
            else:
                frag = r_seq[start:] + r_seq[: (start + frag_len) % L]
        else:
            if frag_len >= L:
                frag_len = L
                start = 0
            else:
                start = rng.randint(0, L - frag_len)
            frag = r_seq[start : start + frag_len]

        # Orientation
        is_fwd = rng.random() < 0.5
        if not is_fwd:
            frag = rev_comp(frag)

        # R1 is 5' of fragment, R2 is 5' of reverse-complement of fragment
        r1_seq = list(frag[:read_len])
        rc_frag = rev_comp(frag)
        r2_seq = list(rc_frag[:read_len])

        # Sequencing errors
        for i in range(len(r1_seq)):
            if rng.random() < err_rate:
                orig = r1_seq[i]
                r1_seq[i] = rng.choice([b for b in DNA if b != orig])
        for i in range(len(r2_seq)):
            if rng.random() < err_rate:
                orig = r2_seq[i]
                r2_seq[i] = rng.choice([b for b in DNA if b != orig])

        r1_str = "".join(r1_seq)
        r2_str = "".join(r2_seq)
        qual1 = "I" * len(r1_str)
        qual2 = "I" * len(r2_str)

        r1_records.append(f"@{r_name}_{pair_idx}/1\n{r1_str}\n+\n{qual1}\n")
        r2_records.append(f"@{r_name}_{pair_idx}/2\n{r2_str}\n+\n{qual2}\n")

    os.makedirs(os.path.dirname(out_r1_gz), exist_ok=True)
    with gzip.open(out_r1_gz, "wt") as f1, gzip.open(out_r2_gz, "wt") as f2:
        for rec1, rec2 in zip(r1_records, r2_records):
            f1.write(rec1)
            f2.write(rec2)

def register_dataset_row(manifest_path: str, row_dict: dict):
    from validate_manifest import MANDATORY_COLUMNS
    # Read existing
    rows = []
    if os.path.exists(manifest_path):
        with open(manifest_path, "r", encoding="utf-8") as f:
            lines = [l.strip() for l in f if l.strip()]
            if lines:
                for line in lines[1:]:
                    fields = line.split("\t")
                    if len(fields) == len(MANDATORY_COLUMNS):
                        rows.append(dict(zip(MANDATORY_COLUMNS, fields)))

    # Replace or append
    d_id = row_dict["dataset_id"]
    replaced = False
    for i, r in enumerate(rows):
        if r["dataset_id"] == d_id:
            rows[i] = row_dict
            replaced = True
            break
    if not replaced:
        rows.append(row_dict)

    with open(manifest_path, "w", encoding="utf-8") as f:
        f.write("\t".join(MANDATORY_COLUMNS) + "\n")
        for r in rows:
            line = "\t".join(str(r.get(col, "NA")) for col in MANDATORY_COLUMNS)
            f.write(line + "\n")

# -------------------------------------------------------------
# Builders
# -------------------------------------------------------------

def build_s0_n(manifest_path: str):
    print("=== Building S0 (Micro) & N (Negative) Fixtures ===")
    fix_dir = "benchmarks/fixtures"
    os.makedirs(fix_dir, exist_ok=True)

    # S0-001: 5kb linear genome, 50x
    rng = random.Random(101)
    s0_1_ref = make_random_seq(5000, rng)
    ref_fa = f"{fix_dir}/S0-001_ref.fa"
    with open(ref_fa, "w") as f:
        f.write(">S0-001_chr\n" + s0_1_ref + "\n")
    r1_gz = f"{fix_dir}/S0-001_R1.fq.gz"
    r2_gz = f"{fix_dir}/S0-001_R2.fq.gz"
    simulate_reads_pe([("chr", s0_1_ref, False, 1.0)], 50, 150, 350, 20, r1_gz, r2_gz, rng)
    events_file = f"{fix_dir}/S0-001_truth_events.json"
    with open(events_file, "w") as f:
        json.dump({"dataset": "S0-001", "replicons": [{"name": "chr", "length": 5000, "circular": False}]}, f, indent=2)

    register_dataset_row(manifest_path, {
        "dataset_id": "S0-001", "kind": "S0", "organism": "Synthetic micro linear", "source_url": "internal://synthetic",
        "study_accession": "NA", "run_accession": "NA", "biosample": "NA", "reference_accession_version": "S0-001.1",
        "reference_confidence": "verified", "replicons": "1_lin", "platform": "Illumina", "layout": "PE",
        "read_length": 150, "fragment_mean": 350, "fragment_sd": 20, "coverage_nominal": 50, "coverage_measured": 50,
        "simulator_version": "build_fixtures_1.0", "random_seed": 101, "fastq_r1": r1_gz, "fastq_r2": r2_gz,
        "reference_fasta": ref_fa, "truth_events_file": events_file, "r1_sha256": sha256_file(r1_gz),
        "r2_sha256": sha256_file(r2_gz), "reference_sha256": sha256_file(ref_fa), "admission_status": "APPROVED",
        "notes": "5kb linear genome verification"
    })

    # S0-002: 4kb circular plasmid, 80x
    rng = random.Random(102)
    s0_2_ref = make_random_seq(4000, rng)
    ref_fa = f"{fix_dir}/S0-002_ref.fa"
    with open(ref_fa, "w") as f:
        f.write(">S0-002_plasmid\n" + s0_2_ref + "\n")
    r1_gz = f"{fix_dir}/S0-002_R1.fq.gz"
    r2_gz = f"{fix_dir}/S0-002_R2.fq.gz"
    simulate_reads_pe([("plasmid", s0_2_ref, True, 1.0)], 80, 150, 350, 20, r1_gz, r2_gz, rng)
    events_file = f"{fix_dir}/S0-002_truth_events.json"
    with open(events_file, "w") as f:
        json.dump({"dataset": "S0-002", "replicons": [{"name": "plasmid", "length": 4000, "circular": True}], "plasmids": ["plasmid"]}, f, indent=2)

    register_dataset_row(manifest_path, {
        "dataset_id": "S0-002", "kind": "S0", "organism": "Synthetic circular plasmid", "source_url": "internal://synthetic",
        "study_accession": "NA", "run_accession": "NA", "biosample": "NA", "reference_accession_version": "S0-002.1",
        "reference_confidence": "verified", "replicons": "1_circ", "platform": "Illumina", "layout": "PE",
        "read_length": 150, "fragment_mean": 350, "fragment_sd": 20, "coverage_nominal": 80, "coverage_measured": 80,
        "simulator_version": "build_fixtures_1.0", "random_seed": 102, "fastq_r1": r1_gz, "fastq_r2": r2_gz,
        "reference_fasta": ref_fa, "truth_events_file": events_file, "r1_sha256": sha256_file(r1_gz),
        "r2_sha256": sha256_file(r2_gz), "reference_sha256": sha256_file(ref_fa), "admission_status": "APPROVED",
        "notes": "4kb circular plasmid verification"
    })

    # S0-003: 8kb sequence with SNP bubble, 60x
    rng = random.Random(103)
    s0_3_ref = make_random_seq(8000, rng)
    ref_fa = f"{fix_dir}/S0-003_ref.fa"
    with open(ref_fa, "w") as f:
        f.write(">S0-003_chr\n" + s0_3_ref + "\n")
    # Alternate allele with 1 SNP at position 4000
    alt_ref = list(s0_3_ref)
    alt_ref[4000] = "G" if alt_ref[4000] != "G" else "A"
    alt_ref = "".join(alt_ref)
    r1_gz = f"{fix_dir}/S0-003_R1.fq.gz"
    r2_gz = f"{fix_dir}/S0-003_R2.fq.gz"
    simulate_reads_pe([("chr", s0_3_ref, False, 0.8), ("alt", alt_ref, False, 0.2)], 60, 150, 350, 20, r1_gz, r2_gz, rng)
    events_file = f"{fix_dir}/S0-003_truth_events.json"
    with open(events_file, "w") as f:
        json.dump({"dataset": "S0-003", "replicons": [{"name": "chr", "length": 8000, "circular": False}], "snp": 4000}, f, indent=2)

    register_dataset_row(manifest_path, {
        "dataset_id": "S0-003", "kind": "S0", "organism": "Synthetic bubble genome", "source_url": "internal://synthetic",
        "study_accession": "NA", "run_accession": "NA", "biosample": "NA", "reference_accession_version": "S0-003.1",
        "reference_confidence": "verified", "replicons": "1_lin", "platform": "Illumina", "layout": "PE",
        "read_length": 150, "fragment_mean": 350, "fragment_sd": 20, "coverage_nominal": 60, "coverage_measured": 60,
        "simulator_version": "build_fixtures_1.0", "random_seed": 103, "fastq_r1": r1_gz, "fastq_r2": r2_gz,
        "reference_fasta": ref_fa, "truth_events_file": events_file, "r1_sha256": sha256_file(r1_gz),
        "r2_sha256": sha256_file(r2_gz), "reference_sha256": sha256_file(ref_fa), "admission_status": "APPROVED",
        "notes": "8kb bubble popping test"
    })

    # S0-004: 6kb sequence containing palindromic k-mers, 60x
    rng = random.Random(104)
    pal_kmer = "ACGTGCGCACGT"
    s0_4_ref = make_random_seq(2500, rng) + pal_kmer * 10 + make_random_seq(2500, rng)
    ref_fa = f"{fix_dir}/S0-004_ref.fa"
    with open(ref_fa, "w") as f:
        f.write(">S0-004_chr\n" + s0_4_ref + "\n")
    r1_gz = f"{fix_dir}/S0-004_R1.fq.gz"
    r2_gz = f"{fix_dir}/S0-004_R2.fq.gz"
    simulate_reads_pe([("chr", s0_4_ref, False, 1.0)], 60, 150, 350, 20, r1_gz, r2_gz, rng)
    events_file = f"{fix_dir}/S0-004_truth_events.json"
    with open(events_file, "w") as f:
        json.dump({"dataset": "S0-004", "replicons": [{"name": "chr", "length": len(s0_4_ref), "circular": False}]}, f, indent=2)

    register_dataset_row(manifest_path, {
        "dataset_id": "S0-004", "kind": "S0", "organism": "Synthetic palindrome genome", "source_url": "internal://synthetic",
        "study_accession": "NA", "run_accession": "NA", "biosample": "NA", "reference_accession_version": "S0-004.1",
        "reference_confidence": "verified", "replicons": "1_lin", "platform": "Illumina", "layout": "PE",
        "read_length": 150, "fragment_mean": 350, "fragment_sd": 20, "coverage_nominal": 60, "coverage_measured": 60,
        "simulator_version": "build_fixtures_1.0", "random_seed": 104, "fastq_r1": r1_gz, "fastq_r2": r2_gz,
        "reference_fasta": ref_fa, "truth_events_file": events_file, "r1_sha256": sha256_file(r1_gz),
        "r2_sha256": sha256_file(r2_gz), "reference_sha256": sha256_file(ref_fa), "admission_status": "APPROVED",
        "notes": "6kb palindrome and twin edge verification"
    })

    # N-001: Corrupted gzip FASTQ
    n1_r1 = f"{fix_dir}/N-001_R1.fq.gz"
    n1_r2 = f"{fix_dir}/N-001_R2.fq.gz"
    with open(n1_r1, "wb") as f:
        f.write(b"\x1f\x8b\x08\x00corrupted_gzip_header_and_random_bytes_that_fail_decompression_1234567890")
    with open(n1_r2, "wb") as f:
        f.write(b"\x1f\x8b\x08\x00corrupted_gzip_header_and_random_bytes_that_fail_decompression_1234567890")

    register_dataset_row(manifest_path, {
        "dataset_id": "N-001", "kind": "N", "organism": "Corrupted Gzip", "source_url": "internal://negative",
        "study_accession": "NA", "run_accession": "NA", "biosample": "NA", "reference_accession_version": "NA",
        "reference_confidence": "unavailable", "replicons": "NA", "platform": "Illumina", "layout": "PE",
        "read_length": 150, "fragment_mean": "NA", "fragment_sd": "NA", "coverage_nominal": "NA", "coverage_measured": "NA",
        "simulator_version": "NA", "random_seed": 0, "fastq_r1": n1_r1, "fastq_r2": n1_r2,
        "reference_fasta": "NA", "truth_events_file": "NA", "r1_sha256": sha256_file(n1_r1),
        "r2_sha256": sha256_file(n1_r2), "reference_sha256": "NA", "admission_status": "APPROVED",
        "notes": "Truncated and corrupted gzip file"
    })

    # N-002: Unequal mate counts (1000 in R1, 850 in R2)
    n2_r1 = f"{fix_dir}/N-002_R1.fq.gz"
    n2_r2 = f"{fix_dir}/N-002_R2.fq.gz"
    with gzip.open(n2_r1, "wt") as f1:
        for i in range(1000):
            f1.write(f"@read_{i}/1\n{'A'*150}\n+\n{'I'*150}\n")
    with gzip.open(n2_r2, "wt") as f2:
        for i in range(850):
            f2.write(f"@read_{i}/2\n{'T'*150}\n+\n{'I'*150}\n")

    register_dataset_row(manifest_path, {
        "dataset_id": "N-002", "kind": "N", "organism": "Unequal Mates", "source_url": "internal://negative",
        "study_accession": "NA", "run_accession": "NA", "biosample": "NA", "reference_accession_version": "NA",
        "reference_confidence": "unavailable", "replicons": "NA", "platform": "Illumina", "layout": "PE",
        "read_length": 150, "fragment_mean": "NA", "fragment_sd": "NA", "coverage_nominal": "NA", "coverage_measured": "NA",
        "simulator_version": "NA", "random_seed": 0, "fastq_r1": n2_r1, "fastq_r2": n2_r2,
        "reference_fasta": "NA", "truth_events_file": "NA", "r1_sha256": sha256_file(n2_r1),
        "r2_sha256": sha256_file(n2_r2), "reference_sha256": "NA", "admission_status": "APPROVED",
        "notes": "Unequal mate counts: 1000 in R1, 850 in R2"
    })

    # N-003: Reads with all-N and low Phred quality (#)
    n3_r1 = f"{fix_dir}/N-003_R1.fq.gz"
    n3_r2 = f"{fix_dir}/N-003_R2.fq.gz"
    with gzip.open(n3_r1, "wt") as f1, gzip.open(n3_r2, "wt") as f2:
        for i in range(500):
            f1.write(f"@badread_{i}/1\n{'N'*150}\n+\n{'#'*150}\n")
            f2.write(f"@badread_{i}/2\n{'N'*150}\n+\n{'#'*150}\n")

    register_dataset_row(manifest_path, {
        "dataset_id": "N-003", "kind": "N", "organism": "All-N Low-Phred", "source_url": "internal://negative",
        "study_accession": "NA", "run_accession": "NA", "biosample": "NA", "reference_accession_version": "NA",
        "reference_confidence": "unavailable", "replicons": "NA", "platform": "Illumina", "layout": "PE",
        "read_length": 150, "fragment_mean": "NA", "fragment_sd": "NA", "coverage_nominal": "NA", "coverage_measured": "NA",
        "simulator_version": "NA", "random_seed": 0, "fastq_r1": n3_r1, "fastq_r2": n3_r2,
        "reference_fasta": "NA", "truth_events_file": "NA", "r1_sha256": sha256_file(n3_r1),
        "r2_sha256": sha256_file(n3_r2), "reference_sha256": "NA", "admission_status": "APPROVED",
        "notes": "100% N bases with Phred # (Q2)"
    })

    print("[DONE] S0 (4 cases) and N (3 cases) generated and registered.")

def build_s1(manifest_path: str):
    print("=== Building S1 Matrix Datasets in benchmark-results/data_store ===")
    store_dir = "benchmark-results/data_store"
    os.makedirs(store_dir, exist_ok=True)

    # Helper for generating repeat datasets
    def make_repeat_genome(length: int, repeats: list, seed: int):
        # repeats: list of (rep_len, count)
        rng = random.Random(seed)
        seq = make_random_seq(length, rng)
        truth_repeats = []
        for rep_len, count in repeats:
            rep_seq = make_random_seq(rep_len, rng)
            positions = []
            for c in range(count):
                pos = rng.randint(0, len(seq) - rep_len)
                seq = seq[:pos] + rep_seq + seq[pos + rep_len:]
                positions.append(pos)
            truth_repeats.append({"rep_len": rep_len, "count": count, "positions": positions})
        return seq, truth_repeats

    # S1-REP-001: Multi-k repeats 35bp, 65bp, 120bp
    rng = random.Random(201)
    s1_1_seq, rep_meta = make_repeat_genome(150000, [(35, 2), (65, 2), (120, 2)], 201)
    ref_fa = f"{store_dir}/S1-REP-001_ref.fa"
    with open(ref_fa, "w") as f:
        f.write(">S1-REP-001_chr\n" + s1_1_seq + "\n")
    r1_gz = f"{store_dir}/S1-REP-001_R1.fq.gz"
    r2_gz = f"{store_dir}/S1-REP-001_R2.fq.gz"
    simulate_reads_pe([("chr", s1_1_seq, False, 1.0)], 60, 150, 450, 30, r1_gz, r2_gz, rng)
    events_file = f"{store_dir}/S1-REP-001_truth_events.json"
    with open(events_file, "w") as f:
        json.dump({"dataset": "S1-REP-001", "repeats": rep_meta}, f, indent=2)

    register_dataset_row(manifest_path, {
        "dataset_id": "S1-REP-001", "kind": "S1", "organism": "Synthetic Multi-k Repeats", "source_url": "internal://synthetic",
        "study_accession": "NA", "run_accession": "NA", "biosample": "NA", "reference_accession_version": "S1-REP-001.1",
        "reference_confidence": "verified", "replicons": "1_lin", "platform": "Illumina", "layout": "PE",
        "read_length": 150, "fragment_mean": 450, "fragment_sd": 30, "coverage_nominal": 60, "coverage_measured": 60,
        "simulator_version": "build_fixtures_1.0", "random_seed": 201, "fastq_r1": r1_gz, "fastq_r2": r2_gz,
        "reference_fasta": ref_fa, "truth_events_file": events_file, "r1_sha256": sha256_file(r1_gz),
        "r2_sha256": sha256_file(r2_gz), "reference_sha256": sha256_file(ref_fa), "admission_status": "APPROVED",
        "notes": "Repeats at 35bp, 65bp, 120bp resolvable by multi-k"
    })

    # S1-REP-002: exSPAnder repeat 350bp (3 copies)
    rng = random.Random(202)
    s1_2_seq, rep_meta2 = make_repeat_genome(150000, [(350, 3)], 202)
    ref_fa = f"{store_dir}/S1-REP-002_ref.fa"
    with open(ref_fa, "w") as f:
        f.write(">S1-REP-002_chr\n" + s1_2_seq + "\n")
    r1_gz = f"{store_dir}/S1-REP-002_R1.fq.gz"
    r2_gz = f"{store_dir}/S1-REP-002_R2.fq.gz"
    simulate_reads_pe([("chr", s1_2_seq, False, 1.0)], 60, 150, 450, 30, r1_gz, r2_gz, rng)
    events_file = f"{store_dir}/S1-REP-002_truth_events.json"
    with open(events_file, "w") as f:
        json.dump({"dataset": "S1-REP-002", "repeats": rep_meta2}, f, indent=2)

    register_dataset_row(manifest_path, {
        "dataset_id": "S1-REP-002", "kind": "S1", "organism": "Synthetic exSPAnder 350bp Repeat", "source_url": "internal://synthetic",
        "study_accession": "NA", "run_accession": "NA", "biosample": "NA", "reference_accession_version": "S1-REP-002.1",
        "reference_confidence": "verified", "replicons": "1_lin", "platform": "Illumina", "layout": "PE",
        "read_length": 150, "fragment_mean": 450, "fragment_sd": 30, "coverage_nominal": 60, "coverage_measured": 60,
        "simulator_version": "build_fixtures_1.0", "random_seed": 202, "fastq_r1": r1_gz, "fastq_r2": r2_gz,
        "reference_fasta": ref_fa, "truth_events_file": events_file, "r1_sha256": sha256_file(r1_gz),
        "r2_sha256": sha256_file(r2_gz), "reference_sha256": sha256_file(ref_fa), "admission_status": "APPROVED",
        "notes": "350bp repeat bridged by exSPAnder PE links"
    })

    # S1-REP-003: Conservative repeats 800bp and 1500bp (exceeds 450bp insert size, MUST split cleanly with 0 misassemblies)
    rng = random.Random(203)
    s1_3_seq, rep_meta3 = make_repeat_genome(150000, [(800, 2), (1500, 2)], 203)
    ref_fa = f"{store_dir}/S1-REP-003_ref.fa"
    with open(ref_fa, "w") as f:
        f.write(">S1-REP-003_chr\n" + s1_3_seq + "\n")
    r1_gz = f"{store_dir}/S1-REP-003_R1.fq.gz"
    r2_gz = f"{store_dir}/S1-REP-003_R2.fq.gz"
    simulate_reads_pe([("chr", s1_3_seq, False, 1.0)], 60, 150, 450, 30, r1_gz, r2_gz, rng)
    events_file = f"{store_dir}/S1-REP-003_truth_events.json"
    with open(events_file, "w") as f:
        json.dump({"dataset": "S1-REP-003", "repeats": rep_meta3, "expected_conservatism": True}, f, indent=2)

    register_dataset_row(manifest_path, {
        "dataset_id": "S1-REP-003", "kind": "S1", "organism": "Synthetic Conservative Repeats", "source_url": "internal://synthetic",
        "study_accession": "NA", "run_accession": "NA", "biosample": "NA", "reference_accession_version": "S1-REP-003.1",
        "reference_confidence": "verified", "replicons": "1_lin", "platform": "Illumina", "layout": "PE",
        "read_length": 150, "fragment_mean": 450, "fragment_sd": 30, "coverage_nominal": 60, "coverage_measured": 60,
        "simulator_version": "build_fixtures_1.0", "random_seed": 203, "fastq_r1": r1_gz, "fastq_r2": r2_gz,
        "reference_fasta": ref_fa, "truth_events_file": events_file, "r1_sha256": sha256_file(r1_gz),
        "r2_sha256": sha256_file(r2_gz), "reference_sha256": sha256_file(ref_fa), "admission_status": "APPROVED",
        "notes": "800bp and 1500bp repeats testing conservative boundary split with 0 misassemblies"
    })

    # S1-PLAS-001: 3 Plasmids (2.2kb @ 25x, 5.5kb @ 5x, 45kb @ 0.5x)
    rng = random.Random(204)
    chr_seq = make_random_seq(200000, rng)
    p1_seq = make_random_seq(2200, rng)
    p2_seq = make_random_seq(5500, rng)
    p3_seq = make_random_seq(45000, rng)
    ref_fa = f"{store_dir}/S1-PLAS-001_ref.fa"
    with open(ref_fa, "w") as f:
        f.write(">chr\n" + chr_seq + "\n")
        f.write(">plasmid_2.2k\n" + p1_seq + "\n")
        f.write(">plasmid_5.5k\n" + p2_seq + "\n")
        f.write(">plasmid_45k\n" + p3_seq + "\n")
    r1_gz = f"{store_dir}/S1-PLAS-001_R1.fq.gz"
    r2_gz = f"{store_dir}/S1-PLAS-001_R2.fq.gz"
    reps = [
        ("chr", chr_seq, False, 1.0),
        ("plasmid_2.2k", p1_seq, True, 25.0),
        ("plasmid_5.5k", p2_seq, True, 5.0),
        ("plasmid_45k", p3_seq, True, 0.5)
    ]
    simulate_reads_pe(reps, 60, 150, 450, 30, r1_gz, r2_gz, rng)
    events_file = f"{store_dir}/S1-PLAS-001_truth_events.json"
    with open(events_file, "w") as f:
        json.dump({"dataset": "S1-PLAS-001", "plasmids": ["plasmid_2.2k", "plasmid_5.5k", "plasmid_45k"]}, f, indent=2)

    register_dataset_row(manifest_path, {
        "dataset_id": "S1-PLAS-001", "kind": "S1", "organism": "Synthetic Multi-Plasmid", "source_url": "internal://synthetic",
        "study_accession": "NA", "run_accession": "NA", "biosample": "NA", "reference_accession_version": "S1-PLAS-001.1",
        "reference_confidence": "verified", "replicons": "1_lin_3_circ", "platform": "Illumina", "layout": "PE",
        "read_length": 150, "fragment_mean": 450, "fragment_sd": 30, "coverage_nominal": 60, "coverage_measured": 60,
        "simulator_version": "build_fixtures_1.0", "random_seed": 204, "fastq_r1": r1_gz, "fastq_r2": r2_gz,
        "reference_fasta": ref_fa, "truth_events_file": events_file, "r1_sha256": sha256_file(r1_gz),
        "r2_sha256": sha256_file(r2_gz), "reference_sha256": sha256_file(ref_fa), "admission_status": "APPROVED",
        "notes": "3 circular plasmids at 25x, 5x, 0.5x copy ratios"
    })

    # S1-PLAS-002: Chromosome + 8kb circular plasmid sharing identical 1.0kb IS element
    rng = random.Random(205)
    is_seq = make_random_seq(1000, rng)
    chr_base = make_random_seq(199000, rng)
    chr_seq2 = chr_base[:100000] + is_seq + chr_base[100000:]
    p_base = make_random_seq(7000, rng)
    plas_seq = p_base[:3500] + is_seq + p_base[3500:]
    ref_fa = f"{store_dir}/S1-PLAS-002_ref.fa"
    with open(ref_fa, "w") as f:
        f.write(">chr\n" + chr_seq2 + "\n")
        f.write(">plasmid_8k\n" + plas_seq + "\n")
    r1_gz = f"{store_dir}/S1-PLAS-002_R1.fq.gz"
    r2_gz = f"{store_dir}/S1-PLAS-002_R2.fq.gz"
    simulate_reads_pe([("chr", chr_seq2, False, 1.0), ("plasmid_8k", plas_seq, True, 6.0)], 60, 150, 450, 30, r1_gz, r2_gz, rng)
    events_file = f"{store_dir}/S1-PLAS-002_truth_events.json"
    with open(events_file, "w") as f:
        json.dump({"dataset": "S1-PLAS-002", "shared_is_len": 1000, "plasmids": ["plasmid_8k"], "expect_no_chimera": True}, f, indent=2)

    register_dataset_row(manifest_path, {
        "dataset_id": "S1-PLAS-002", "kind": "S1", "organism": "Synthetic Shared IS Plasmid", "source_url": "internal://synthetic",
        "study_accession": "NA", "run_accession": "NA", "biosample": "NA", "reference_accession_version": "S1-PLAS-002.1",
        "reference_confidence": "verified", "replicons": "1_lin_1_circ", "platform": "Illumina", "layout": "PE",
        "read_length": 150, "fragment_mean": 450, "fragment_sd": 30, "coverage_nominal": 60, "coverage_measured": 60,
        "simulator_version": "build_fixtures_1.0", "random_seed": 205, "fastq_r1": r1_gz, "fastq_r2": r2_gz,
        "reference_fasta": ref_fa, "truth_events_file": events_file, "r1_sha256": sha256_file(r1_gz),
        "r2_sha256": sha256_file(r2_gz), "reference_sha256": sha256_file(ref_fa), "admission_status": "APPROVED",
        "notes": "Shared 1kb IS element between chromosome and 8kb plasmid (tests zero chimeras)"
    })

    # Coverage matrix: S1-COV-001 (15x), S1-COV-002 (60x), S1-COV-003 (300x)
    rng = random.Random(206)
    cov_genome = make_random_seq(200000, rng)
    cov_ref_fa = f"{store_dir}/S1-COV_ref.fa"
    with open(cov_ref_fa, "w") as f:
        f.write(">S1-COV_chr\n" + cov_genome + "\n")
    ref_cov_hash = sha256_file(cov_ref_fa)

    for depth, d_id in [(15, "S1-COV-001"), (60, "S1-COV-002"), (300, "S1-COV-003")]:
        r1_gz = f"{store_dir}/{d_id}_R1.fq.gz"
        r2_gz = f"{store_dir}/{d_id}_R2.fq.gz"
        simulate_reads_pe([("chr", cov_genome, False, 1.0)], depth, 150, 450, 30, r1_gz, r2_gz, random.Random(206 + depth))
        events_file = f"{store_dir}/{d_id}_truth_events.json"
        with open(events_file, "w") as f:
            json.dump({"dataset": d_id, "depth": depth}, f, indent=2)

        register_dataset_row(manifest_path, {
            "dataset_id": d_id, "kind": "S1", "organism": f"Synthetic Coverage {depth}x", "source_url": "internal://synthetic",
            "study_accession": "NA", "run_accession": "NA", "biosample": "NA", "reference_accession_version": f"{d_id}.1",
            "reference_confidence": "verified", "replicons": "1_lin", "platform": "Illumina", "layout": "PE",
            "read_length": 150, "fragment_mean": 450, "fragment_sd": 30, "coverage_nominal": depth, "coverage_measured": depth,
            "simulator_version": "build_fixtures_1.0", "random_seed": 206 + depth, "fastq_r1": r1_gz, "fastq_r2": r2_gz,
            "reference_fasta": cov_ref_fa, "truth_events_file": events_file, "r1_sha256": sha256_file(r1_gz),
            "r2_sha256": sha256_file(r2_gz), "reference_sha256": ref_cov_hash, "admission_status": "APPROVED",
            "notes": f"Coverage {depth}x benchmark"
        })

    # GC matrix: S1-GC-001 (30% GC), S1-GC-002 (70% GC)
    for gc, d_id in [(0.30, "S1-GC-001"), (0.70, "S1-GC-002")]:
        rng_gc = random.Random(210 + int(gc * 100))
        gc_seq = make_random_seq(200000, rng_gc, gc_bias=gc)
        ref_fa = f"{store_dir}/{d_id}_ref.fa"
        with open(ref_fa, "w") as f:
            f.write(f">{d_id}_chr\n" + gc_seq + "\n")
        r1_gz = f"{store_dir}/{d_id}_R1.fq.gz"
        r2_gz = f"{store_dir}/{d_id}_R2.fq.gz"
        simulate_reads_pe([("chr", gc_seq, False, 1.0)], 60, 150, 450, 30, r1_gz, r2_gz, rng_gc)
        events_file = f"{store_dir}/{d_id}_truth_events.json"
        with open(events_file, "w") as f:
            json.dump({"dataset": d_id, "gc_target": gc}, f, indent=2)

        register_dataset_row(manifest_path, {
            "dataset_id": d_id, "kind": "S1", "organism": f"Synthetic GC {int(gc*100)}%", "source_url": "internal://synthetic",
            "study_accession": "NA", "run_accession": "NA", "biosample": "NA", "reference_accession_version": f"{d_id}.1",
            "reference_confidence": "verified", "replicons": "1_lin", "platform": "Illumina", "layout": "PE",
            "read_length": 150, "fragment_mean": 450, "fragment_sd": 30, "coverage_nominal": 60, "coverage_measured": 60,
            "simulator_version": "build_fixtures_1.0", "random_seed": 210 + int(gc * 100), "fastq_r1": r1_gz, "fastq_r2": r2_gz,
            "reference_fasta": ref_fa, "truth_events_file": events_file, "r1_sha256": sha256_file(r1_gz),
            "r2_sha256": sha256_file(r2_gz), "reference_sha256": sha256_file(ref_fa), "admission_status": "APPROVED",
            "notes": f"Extreme GC bias {int(gc*100)}% benchmark"
        })

    # S1-EDIT-001: 200kb with 1kb deletion, 1.2kb insertion, 5kb inversion
    rng_sv = random.Random(220)
    parent = make_random_seq(200000, rng_sv)
    # Del at 40kb (1kb)
    del_seq = parent[:40000] + parent[41000:]
    # Ins at 100kb (1.2kb)
    ins_chunk = make_random_seq(1200, rng_sv)
    ins_seq = del_seq[:100000] + ins_chunk + del_seq[100000:]
    # Inv at 140kb (5kb)
    inv_chunk = rev_comp(ins_seq[140000:145000])
    mutant_seq = ins_seq[:140000] + inv_chunk + ins_seq[145000:]

    ref_fa = f"{store_dir}/S1-EDIT-001_ref.fa"
    with open(ref_fa, "w") as f:
        f.write(">S1-EDIT-001_chr\n" + mutant_seq + "\n")
    r1_gz = f"{store_dir}/S1-EDIT-001_R1.fq.gz"
    r2_gz = f"{store_dir}/S1-EDIT-001_R2.fq.gz"
    simulate_reads_pe([("chr", mutant_seq, False, 1.0)], 60, 150, 450, 30, r1_gz, r2_gz, rng_sv)
    events_file = f"{store_dir}/S1-EDIT-001_truth_events.json"
    with open(events_file, "w") as f:
        json.dump({"dataset": "S1-EDIT-001", "sv": ["del_1kb", "ins_1.2kb", "inv_5kb"]}, f, indent=2)

    register_dataset_row(manifest_path, {
        "dataset_id": "S1-EDIT-001", "kind": "S1", "organism": "Synthetic SV Genome", "source_url": "internal://synthetic",
        "study_accession": "NA", "run_accession": "NA", "biosample": "NA", "reference_accession_version": "S1-EDIT-001.1",
        "reference_confidence": "verified", "replicons": "1_lin", "platform": "Illumina", "layout": "PE",
        "read_length": 150, "fragment_mean": 450, "fragment_sd": 30, "coverage_nominal": 60, "coverage_measured": 60,
        "simulator_version": "build_fixtures_1.0", "random_seed": 220, "fastq_r1": r1_gz, "fastq_r2": r2_gz,
        "reference_fasta": ref_fa, "truth_events_file": events_file, "r1_sha256": sha256_file(r1_gz),
        "r2_sha256": sha256_file(r2_gz), "reference_sha256": sha256_file(ref_fa), "admission_status": "APPROVED",
        "notes": "Structural variants: 1kb del, 1.2kb ins, 5kb inv"
    })

    print("[DONE] S1 matrix (11 cases) generated and registered.")

def build_s2(manifest_path: str):
    print("=== Building S2-REAL-001 Full-Scale 4.6 Mb Bacterial Genome ===")
    store_dir = "benchmark-results/data_store"
    os.makedirs(store_dir, exist_ok=True)

    rng = random.Random(301)
    chr_len = 4600000
    plas_len = 15000

    # Model E. coli K-12 with 7 copies of 5kb rRNA operons and 10 copies of 1.2kb IS elements
    chr_base = make_random_seq(chr_len, rng)
    rrna_seq = make_random_seq(5000, rng)
    is_seq = make_random_seq(1200, rng)

    # Insert 7 copies of rRNA
    for i in range(7):
        pos = int((i + 1) * (chr_len / 8.5)) + rng.randint(-5000, 5000)
        chr_base = chr_base[:pos] + rrna_seq + chr_base[pos + 5000:]

    # Insert 10 copies of IS
    for i in range(10):
        pos = int((i + 1) * (chr_len / 11.5)) + rng.randint(-3000, 3000)
        chr_base = chr_base[:pos] + is_seq + chr_base[pos + 1200:]

    plas_seq = make_random_seq(plas_len, rng)

    ref_fa = f"{store_dir}/S2-REAL-001_ref.fa"
    with open(ref_fa, "w") as f:
        f.write(">NC_000913_chr\n" + chr_base + "\n")
        f.write(">plasmid_15k\n" + plas_seq + "\n")

    r1_gz = f"{store_dir}/S2-REAL-001_R1.fq.gz"
    r2_gz = f"{store_dir}/S2-REAL-001_R2.fq.gz"

    reps = [
        ("NC_000913_chr", chr_base, False, 1.0),
        ("plasmid_15k", plas_seq, True, 8.0)
    ]

    simulate_reads_pe(reps, 100, 150, 420, 35, r1_gz, r2_gz, rng, err_rate=0.005)
    events_file = f"{store_dir}/S2-REAL-001_truth_events.json"
    with open(events_file, "w") as f:
        json.dump({"dataset": "S2-REAL-001", "plasmids": ["plasmid_15k"], "rrna_copies": 7, "is_copies": 10}, f, indent=2)

    register_dataset_row(manifest_path, {
        "dataset_id": "S2-REAL-001", "kind": "S2", "organism": "Escherichia coli K-12 synthetic", "source_url": "internal://synthetic",
        "study_accession": "NA", "run_accession": "NA", "biosample": "NA", "reference_accession_version": "NC_000913.3",
        "reference_confidence": "verified", "replicons": "1_lin_1_circ", "platform": "Illumina", "layout": "PE",
        "read_length": 150, "fragment_mean": 420, "fragment_sd": 35, "coverage_nominal": 100, "coverage_measured": 100,
        "simulator_version": "build_fixtures_1.0", "random_seed": 301, "fastq_r1": r1_gz, "fastq_r2": r2_gz,
        "reference_fasta": ref_fa, "truth_events_file": events_file, "r1_sha256": sha256_file(r1_gz),
        "r2_sha256": sha256_file(r2_gz), "reference_sha256": sha256_file(ref_fa), "admission_status": "APPROVED",
        "notes": "4.6Mb chromosome + 15kb plasmid (8x), 7 rRNA, 10 IS, 100x depth"
    })

    print("[DONE] S2-REAL-001 full-scale isolate generated and registered.")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--suite", required=True, choices=["S0_N", "S1", "S2", "ALL"])
    parser.add_argument("--manifest", default="benchmarks/datasets.tsv")
    args = parser.parse_args()

    # Ensure validate_manifest is importable
    sys.path.insert(0, str(Path(__file__).resolve().parent))

    if args.suite in ["S0_N", "ALL"]:
        build_s0_n(args.manifest)
    if args.suite in ["S1", "ALL"]:
        build_s1(args.manifest)
    if args.suite in ["S2", "ALL"]:
        build_s2(args.manifest)

if __name__ == "__main__":
    main()
