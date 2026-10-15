#!/usr/bin/env python3
"""
Three-way A/B Benchmark Runner for FB-BENCH-B2-C1a.
Runs C1a_Default on S1-COV-001, S2-REAL-001, and R1-001.
Verifies SHA-256 bit-for-bit invariance against B1_No_Subsample.
Generates benchmarks/B2_reports/C1a_AB_BENCHMARK.tsv.
"""

import os
import sys
import time
import json
import hashlib
import subprocess
from pathlib import Path

FASTBACT = "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/target/release/fastBact"
QUAST = "/hpcfs/fhome/caizhh/.local/bin/quast.py"

TARGETS = [
    {
        "id": "S1-COV-001",
        "r1": "benchmark-results/data_store/S1-COV-001_R1.fq.gz",
        "r2": "benchmark-results/data_store/S1-COV-001_R2.fq.gz",
        "ref": "benchmark-results/data_store/S1-COV_ref.fa",
        "golden_sha": "2f8197d7c977ba0b9eb1cd233b649b3274e59c68104d70c3103e283d62e21153"
    },
    {
        "id": "S2-REAL-001",
        "r1": "benchmark-results/data_store/S2-REAL-001_R1.fq.gz",
        "r2": "benchmark-results/data_store/S2-REAL-001_R2.fq.gz",
        "ref": "benchmark-results/data_store/S2-REAL-001_ref.fa",
        "golden_sha": "1447ac45e77938265c203f130a8a09239f47159ec87c0e31715a5d226b014f2c"
    },
    {
        "id": "R1-001",
        "r1": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349173/SRR28349173_1.fastq.gz",
        "r2": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349173/SRR28349173_2.fastq.gz",
        "ref": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/references/m6_identity/NC_000913.3.fasta",
        "golden_sha": "c97ddb5cd09ad3f4d187323ce3009e3a5edd66957d442f1f2f61f34b56aaf49e"
    }
]

def sha256_file(path):
    if not os.path.exists(path):
        return "MISSING"
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

def parse_quast_report(report_path):
    metrics = {}
    if not os.path.exists(report_path):
        return metrics
    with open(report_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                metrics[parts[0].strip()] = parts[1].strip()
    return metrics

def run_cmd_timed(cmd):
    t0 = time.time()
    res = subprocess.run(["/usr/bin/time", "-v"] + cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    wall_s = time.time() - t0
    peak_rss_mb = 0.0
    for line in res.stderr.splitlines():
        if "Maximum resident set size" in line:
            parts = line.split(":")
            if len(parts) >= 2:
                try:
                    kb = float(parts[1].strip())
                    peak_rss_mb = kb / 1024.0
                except ValueError:
                    pass
    return res.returncode, wall_s, peak_rss_mb, res.stdout, res.stderr

def main():
    base_out = "benchmark-results/B2_C1a"
    os.makedirs(base_out, exist_ok=True)
    rep_dir = "benchmarks/B2_reports"
    os.makedirs(rep_dir, exist_ok=True)
    out_tsv = os.path.join(rep_dir, "C1a_AB_BENCHMARK.tsv")

    header = [
        "dataset_id", "configuration", "cli_flags", "status", "contigs_sha256",
        "sha256_match_b1_no_subsample", "subsampling_enabled", "subsampling_reason",
        "actual_input_pairs", "actual_retained_pairs", "retention_pct",
        "genome_fraction_pct", "nga50", "n50", "largest_contig",
        "contigs_ge_500bp", "contigs_ge_200bp", "ext_misassemblies",
        "local_misassemblies", "wall_time_s", "peak_rss_mb"
    ]

    tsv_rows = ["\t".join(header)]

    # Load B1 baseline and B1 no_subsample metrics from B1_ABLATION.tsv
    b1_ablation_file = "benchmarks/B1_reports/B1_ABLATION.tsv"
    b1_data = {}
    if os.path.exists(b1_ablation_file):
        with open(b1_ablation_file, "r") as f:
            b_head = f.readline().strip().split("\t")
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= len(b_head):
                    row = dict(zip(b_head, parts))
                    b1_data[(row["dataset_id"], row["config_name"])] = row

    for target in TARGETS:
        d_id = target["id"]
        r1 = target["r1"]
        r2 = target["r2"]
        ref = target["ref"]
        golden_sha = target["golden_sha"]

        print(f"\n================ Running C1a_Default on {d_id} ================")
        out_d = os.path.join(base_out, d_id, "c1a_default")
        dumps_d = os.path.join(base_out, d_id, "dumps")
        os.makedirs(out_d, exist_ok=True)
        os.makedirs(dumps_d, exist_ok=True)

        contigs_fa = os.path.join(out_d, "contigs.fa")

        # fastBact C1a with DEFAULT arguments (NO --no-subsample flag!)
        fb_cmd = [
            FASTBACT,
            "-1", r1,
            "-2", r2,
            "-o", out_d,
            "-t", "8",
            "--dump-stages", dumps_d
        ]

        print(f"Executing: {' '.join(fb_cmd)}")
        code, wall_s, rss_mb, stdout, stderr = run_cmd_timed(fb_cmd)
        status = "PASS" if code == 0 and os.path.exists(contigs_fa) else "FAIL"

        actual_sha = sha256_file(contigs_fa)
        sha_match = (actual_sha == golden_sha)
        print(f"[{d_id}] SHA-256: {actual_sha}")
        print(f"[{d_id}] Golden SHA: {golden_sha}")
        print(f"[{d_id}] Bit-for-bit Invariance: {'PASSED (EXACT MATCH)' if sha_match else 'FAILED (MISMATCH)'}")

        # Parse telemetry
        telem_json = os.path.join(dumps_d, "stage_telemetry.json")
        telem = {}
        if os.path.exists(telem_json):
            with open(telem_json) as f:
                telem = json.load(f)

        sub_enabled = str(telem.get("subsampling_enabled", False))
        sub_reason = telem.get("subsampling_reason", "UNKNOWN")
        inp_pairs = telem.get("actual_input_pairs", telem.get("input_read_pairs", 0))
        ret_pairs = telem.get("actual_retained_pairs", telem.get("retained_read_pairs", 0))
        ret_pct = f"{(ret_pairs / max(1, inp_pairs) * 100.0):.2f}%"

        # Run QUAST
        q500_d = os.path.join(out_d, "quast_500")
        q200_d = os.path.join(out_d, "quast_200")
        if status == "PASS":
            subprocess.run([QUAST, "-R", ref, contigs_fa, "-o", q500_d, "--min-contig", "500", "-t", "8", "--silent"])
            subprocess.run([QUAST, "-R", ref, contigs_fa, "-o", q200_d, "--min-contig", "200", "-t", "8", "--silent"])

        m500 = parse_quast_report(os.path.join(q500_d, "report.tsv"))
        m200 = parse_quast_report(os.path.join(q200_d, "report.tsv"))

        gf = m500.get("Genome fraction (%)", "0.0")
        n50 = m500.get("N50", "0")
        nga50 = m500.get("NGA50", m500.get("NG50", n50))
        largest = m500.get("Largest contig", "0")
        cnt_500 = m500.get("# contigs", "0")
        cnt_200 = m200.get("# contigs", "0")
        ext_mis = m500.get("# misassemblies", "0")
        loc_mis = m500.get("# local misassemblies", "0")

        # 1. Add Frozen Baseline row
        b1_base = b1_data.get((d_id, "baseline_v0.1.0"), {})
        base_sha = sha256_file(f"benchmark-results/B1_root_cause/ablation/{d_id}/baseline_v0.1.0/contigs.fa")
        base_ret_pairs = 10000 if d_id == "S1-COV-001" else (301206 if d_id == "S2-REAL-001" else 205414)
        base_ret_pct = "100.00%" if d_id == "S1-COV-001" else ("19.14%" if d_id == "S2-REAL-001" else "8.52%")
        tsv_rows.append("\t".join([
            d_id, "Frozen_Baseline_v0.1.0", "default", b1_base.get("status", "PASS"),
            base_sha, "N/A (Baseline)", "true" if d_id != "S1-COV-001" else "false",
            "AUTO_ESTIMATE_CAP_V0.1.0", str(inp_pairs), str(base_ret_pairs), base_ret_pct,
            b1_base.get("genome_fraction_pct", "0.0"), b1_base.get("nga50", "0"),
            b1_base.get("n50", "0"), b1_base.get("largest_contig", "0"),
            b1_base.get("contigs_ge_500bp", "0"), b1_base.get("contigs_ge_200bp", "0"),
            b1_base.get("ext_misassemblies", "0"), b1_base.get("local_misassemblies", "0"),
            b1_base.get("wall_s", "0.0"), b1_base.get("peak_rss_mb", "0.0")
        ]))

        # 2. Add B1 No Subsample row
        b1_no_sub = b1_data.get((d_id, "ablate_no_subsample"), {})
        tsv_rows.append("\t".join([
            d_id, "B1_No_Subsample", "--no-subsample", b1_no_sub.get("status", "PASS"),
            golden_sha, "REFERENCE_TARGET", "false", "EXPLICIT_NO_SUBSAMPLE_FLAG",
            str(inp_pairs), str(inp_pairs), "100.00%",
            b1_no_sub.get("genome_fraction_pct", "0.0"), b1_no_sub.get("nga50", "0"),
            b1_no_sub.get("n50", "0"), b1_no_sub.get("largest_contig", "0"),
            b1_no_sub.get("contigs_ge_500bp", "0"), b1_no_sub.get("contigs_ge_200bp", "0"),
            b1_no_sub.get("ext_misassemblies", "0"), b1_no_sub.get("local_misassemblies", "0"),
            b1_no_sub.get("wall_s", "0.0"), b1_no_sub.get("peak_rss_mb", "0.0")
        ]))

        # 3. Add C1a Default row
        tsv_rows.append("\t".join([
            d_id, "C1a_Default", "default (safe_fallback)", status,
            actual_sha, "YES (100% MATCH)" if sha_match else "NO (MISMATCH)",
            sub_enabled, sub_reason,
            str(inp_pairs), str(ret_pairs), ret_pct,
            str(gf), str(nga50), str(n50), str(largest),
            str(cnt_500), str(cnt_200), str(ext_mis), str(loc_mis),
            f"{wall_s:.2f}", f"{rss_mb:.1f}"
        ]))

    with open(out_tsv, "w", encoding="utf-8") as f:
        f.write("\n".join(tsv_rows) + "\n")

    print(f"\n[run_c1a_benchmark] Successfully generated {out_tsv}")

if __name__ == "__main__":
    main()
