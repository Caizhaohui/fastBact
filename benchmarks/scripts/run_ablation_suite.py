#!/usr/bin/env python3
"""
Single-Factor Ablation Experiment Runner (FB-BENCH-B1 Phase B1-3)
Executes 9 configurations on S1-COV-001, S2-REAL-001, and R1-001, runs QUAST,
and generates benchmarks/B1_reports/B1_ABLATION.tsv.
"""

import os
import sys
import time
import subprocess
from pathlib import Path

FASTBACT = "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/target/release/fastBact"
QUAST = "/hpcfs/fhome/caizhh/.local/bin/quast.py"

TARGETS = [
    {
        "id": "S1-COV-001",
        "r1": "benchmark-results/data_store/S1-COV-001_R1.fq.gz",
        "r2": "benchmark-results/data_store/S1-COV-001_R2.fq.gz",
        "ref": "benchmark-results/data_store/S1-COV_ref.fa"
    },
    {
        "id": "S2-REAL-001",
        "r1": "benchmark-results/data_store/S2-REAL-001_R1.fq.gz",
        "r2": "benchmark-results/data_store/S2-REAL-001_R2.fq.gz",
        "ref": "benchmark-results/data_store/S2-REAL-001_ref.fa"
    },
    {
        "id": "R1-001",
        "r1": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349173/SRR28349173_1.fastq.gz",
        "r2": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349173/SRR28349173_2.fastq.gz",
        "ref": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/references/m6_identity/NC_000913.3.fasta"
    }
]

CONFIGS = [
    ("baseline_v0.1.0", []),
    ("ablate_no_trim", ["--no-trim"]),
    ("ablate_no_merge", ["--no-merge"]),
    ("ablate_no_subsample", ["--no-subsample"]),
    ("ablate_no_tip_clip", ["--no-tip-clip"]),
    ("ablate_no_bubble_pop", ["--no-bubble-pop"]),
    ("ablate_no_multik", ["--no-multik"]),
    ("ablate_no_exspander", ["--no-exspander"]),
    ("ablate_no_polish", ["--no-polish"]),
]

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
    return res.returncode, wall_s, peak_rss_mb

def main():
    base_dir = "benchmark-results/B1_root_cause/ablation"
    os.makedirs(base_dir, exist_ok=True)
    out_tsv = "benchmarks/B1_reports/B1_ABLATION.tsv"
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)

    header = [
        "dataset_id", "config_name", "cli_flags", "status", "ext_misassemblies",
        "local_misassemblies", "genome_fraction_pct", "nga50", "n50",
        "largest_contig", "contigs_ge_500bp", "contigs_ge_200bp",
        "mismatches_per_100kb", "indels_per_100kb", "wall_s", "peak_rss_mb",
        "delta_nga50_vs_base", "delta_gf_pct_vs_base"
    ]

    tsv_rows = ["\t".join(header)]

    for target in TARGETS:
        d_id = target["id"]
        r1 = target["r1"]
        r2 = target["r2"]
        ref = target["ref"]

        baseline_nga50 = None
        baseline_gf = None

        for cfg_name, flags in CONFIGS:
            out_dir = os.path.join(base_dir, d_id, cfg_name)
            os.makedirs(out_dir, exist_ok=True)
            contigs_fa = os.path.join(out_dir, "contigs.fa")

            cmd = [
                FASTBACT,
                "-1", r1,
                "-2", r2,
                "-o", out_dir,
                "-t", "8"
            ] + flags

            flags_str = " ".join(flags) if flags else "default"
            print(f"[ablation] Running {cfg_name} on {d_id} ({flags_str})...")
            code, wall_s, rss_mb = run_cmd_timed(cmd)

            status = "PASS" if code == 0 and os.path.exists(contigs_fa) else "FAIL"

            # Run QUAST
            q500_dir = os.path.join(out_dir, "quast_500")
            q200_dir = os.path.join(out_dir, "quast_200")

            if status == "PASS":
                subprocess.run([QUAST, "-R", ref, contigs_fa, "-o", q500_dir, "--min-contig", "500", "-t", "8", "--silent"])
                subprocess.run([QUAST, "-R", ref, contigs_fa, "-o", q200_dir, "--min-contig", "200", "-t", "8", "--silent"])

            m500 = parse_quast_report(os.path.join(q500_dir, "report.tsv"))
            m200 = parse_quast_report(os.path.join(q200_dir, "report.tsv"))

            ext_mis = m500.get("# misassemblies", "0")
            loc_mis = m500.get("# local misassemblies", "0")
            gf = float(m500.get("Genome fraction (%)", "0.0") or 0.0)
            n50 = int(float(m500.get("N50", "0") or 0))
            nga50_val = m500.get("NGA50", m500.get("NG50", str(n50)))
            try:
                nga50 = int(float(nga50_val or 0))
            except ValueError:
                nga50 = n50
            largest = m500.get("Largest contig", "0")
            cnt_500 = m500.get("# contigs", "0")
            cnt_200 = m200.get("# contigs", "0")
            mm_100k = m500.get("# mismatches per 100 kbp", "0.0")
            ind_100k = m500.get("# indels per 100 kbp", "0.0")

            if cfg_name == "baseline_v0.1.0":
                baseline_nga50 = nga50
                baseline_gf = gf
                delta_nga50 = "0"
                delta_gf = "0.00"
            else:
                delta_nga50 = f"{nga50 - (baseline_nga50 or 0):+d}"
                delta_gf = f"{gf - (baseline_gf or 0.0):+.2f}"

            row = [
                d_id, cfg_name, flags_str, status, str(ext_mis),
                str(loc_mis), f"{gf:.2f}", str(nga50), str(n50),
                str(largest), str(cnt_500), str(cnt_200),
                str(mm_100k), str(ind_100k), f"{wall_s:.2f}", f"{rss_mb:.1f}",
                delta_nga50, delta_gf
            ]
            tsv_rows.append("\t".join(row))

    with open(out_tsv, "w", encoding="utf-8") as f:
        f.write("\n".join(tsv_rows) + "\n")
    print(f"[run_ablation_suite] Successfully emitted {out_tsv}")

if __name__ == "__main__":
    main()
