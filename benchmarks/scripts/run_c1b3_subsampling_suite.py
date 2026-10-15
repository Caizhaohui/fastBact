#!/usr/bin/env python3
"""
C1b-3 Subsampling & Dynamic CMS Benchmark Runner.
Evaluates Coverage Gradients (15x to 300x) and Real Isolates under:
1. Default (Subsampling Disabled, 100% data retention)
2. Auto-Subsample (--enable-auto-subsample, default target depth 120x)
Verifies:
- Safe interception on low coverage (15x S1-COV-001)
- Safe bypass when coverage <= 120x
- Safe and beneficial subsampling when coverage > 120x (150x, 300x)
- Retention of assembly continuity, completeness, and 0 misassemblies
- Significant speedup and memory reduction at 300x
Outputs benchmarks/B2_reports/C1b3_SUBSAMPLING_BENCHMARK.tsv
"""

import os
import sys
import time
import hashlib
import subprocess
from pathlib import Path

FASTBACT = "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/target/release/fastBact"
QUAST = "/hpcfs/fhome/caizhh/.local/bin/quast.py"

DATASETS = [
    {
        "id": "S1-COV-001",
        "r1": "benchmark-results/data_store/S1-COV-001_R1.fq.gz",
        "r2": "benchmark-results/data_store/S1-COV-001_R2.fq.gz",
        "ref": "benchmark-results/data_store/S1-COV_ref.fa",
        "true_cov": 15.0,
    },
    {
        "id": "S1-COV-030X",
        "r1": "benchmark-results/data_store/S1-COV-030X_R1.fq.gz",
        "r2": "benchmark-results/data_store/S1-COV-030X_R2.fq.gz",
        "ref": "benchmark-results/data_store/S1-COV_ref.fa",
        "true_cov": 30.0,
    },
    {
        "id": "S1-COV-002",
        "r1": "benchmark-results/data_store/S1-COV-002_R1.fq.gz",
        "r2": "benchmark-results/data_store/S1-COV-002_R2.fq.gz",
        "ref": "benchmark-results/data_store/S1-COV_ref.fa",
        "true_cov": 60.0,
    },
    {
        "id": "S1-COV-100X",
        "r1": "benchmark-results/data_store/S1-COV-100X_R1.fq.gz",
        "r2": "benchmark-results/data_store/S1-COV-100X_R2.fq.gz",
        "ref": "benchmark-results/data_store/S1-COV_ref.fa",
        "true_cov": 100.0,
    },
    {
        "id": "S1-COV-150X",
        "r1": "benchmark-results/data_store/S1-COV-150X_R1.fq.gz",
        "r2": "benchmark-results/data_store/S1-COV-150X_R2.fq.gz",
        "ref": "benchmark-results/data_store/S1-COV_ref.fa",
        "true_cov": 150.0,
    },
    {
        "id": "S1-COV-003",
        "r1": "benchmark-results/data_store/S1-COV-003_R1.fq.gz",
        "r2": "benchmark-results/data_store/S1-COV-003_R2.fq.gz",
        "ref": "benchmark-results/data_store/S1-COV_ref.fa",
        "true_cov": 300.0,
    },
    {
        "id": "S2-REAL-001",
        "r1": "benchmark-results/data_store/S2-REAL-001_R1.fq.gz",
        "r2": "benchmark-results/data_store/S2-REAL-001_R2.fq.gz",
        "ref": "benchmark-results/data_store/S2-REAL-001_ref.fa",
        "true_cov": 100.0,
    },
    {
        "id": "R1-001",
        "r1": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349173/SRR28349173_1.fastq.gz",
        "r2": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349173/SRR28349173_2.fastq.gz",
        "ref": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/references/m6_identity/NC_000913.3.fasta",
        "true_cov": 77.8,
    },
    {
        "id": "R1-003",
        "r1": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349164/SRR28349164_1.fastq.gz",
        "r2": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349164/SRR28349164_2.fastq.gz",
        "ref": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/references/m6_identity/CP009273.1.fasta",
        "true_cov": 132.3,
    },
]

def sha256_file(path):
    if not os.path.exists(path):
        return "MISSING"
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

def run_cmd_timed(cmd):
    t0 = time.time()
    wrapped_cmd = ["/usr/bin/time", "-v"] + cmd
    res = subprocess.run(wrapped_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
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

def parse_fastbact_log(log_text):
    stats = {}
    for line in log_text.splitlines():
        if "exact_read_pairs:" in line:
            stats["exact_read_pairs"] = int(line.split("exact_read_pairs:")[1].strip())
        elif "total_bases:" in line:
            parts = line.split("total_bases:")[1].strip().split()[0]
            stats["total_bases"] = int(parts)
        elif "estimated_genome_size:" in line:
            parts = line.split("estimated_genome_size:")[1].strip().split()[0]
            stats["estimated_genome_size"] = int(parts)
        elif "estimated_coverage:" in line:
            parts = line.split("estimated_coverage:")[1].strip().rstrip("x")
            stats["estimated_coverage"] = float(parts)
        elif "estimator_confidence:" in line:
            stats["estimator_confidence"] = line.split("estimator_confidence:")[1].strip()
        elif "subsampling_enabled:" in line:
            stats["subsampling_enabled"] = line.split("subsampling_enabled:")[1].strip()
        elif "subsampling_reason:" in line:
            stats["subsampling_reason"] = line.split("subsampling_reason:")[1].strip()
        elif "subsample_factor:" in line:
            stats["subsample_factor"] = float(line.split("subsample_factor:")[1].strip())
        elif "actual_input_pairs:" in line:
            stats["actual_input_pairs"] = int(line.split("actual_input_pairs:")[1].strip())
        elif "actual_retained_pairs:" in line:
            stats["actual_retained_pairs"] = int(line.split("actual_retained_pairs:")[1].strip())
    return stats

def main():
    base_out = "benchmark-results/B2_C1b/c1b3_subsampling_suite"
    os.makedirs(base_out, exist_ok=True)
    rep_dir = "benchmarks/B2_reports"
    os.makedirs(rep_dir, exist_ok=True)
    out_tsv = os.path.join(rep_dir, "C1b3_SUBSAMPLING_BENCHMARK.tsv")

    header = [
        "dataset_id",
        "mode",
        "true_cov",
        "est_genome_size",
        "est_cov",
        "estimator_confidence",
        "subsample_factor",
        "subsampling_enabled",
        "subsampling_reason",
        "actual_input_pairs",
        "actual_retained_pairs",
        "retention_pct",
        "contigs_sha256",
        "genome_fraction_pct",
        "nga50",
        "n50",
        "largest_contig",
        "contigs_ge_500bp",
        "ext_misassemblies",
        "local_misassemblies",
        "wall_time_s",
        "peak_rss_mb",
    ]

    tsv_rows = ["\t".join(header)]

    configs = [
        ("Default_100pct", []),
        ("Auto_Subsample", ["--enable-auto-subsample"]),
    ]

    for d in DATASETS:
        d_id = d["id"]
        for cfg_name, extra_flags in configs:
            run_name = f"{d_id}_{cfg_name}"
            out_sub = os.path.join(base_out, run_name)
            os.makedirs(out_sub, exist_ok=True)
            contigs_fa = os.path.join(out_sub, "contigs.fa")

            print(f"\n=======================================================")
            print(f"Running fastBact {cfg_name} on {d_id}...")
            print(f"=======================================================")

            log_file = os.path.join(out_sub, "fastbact.log")
            if os.path.exists(contigs_fa) and os.path.exists(os.path.join(out_sub, "quast_out", "report.tsv")) and os.path.exists(log_file):
                print(f"Skipping execution for {run_name} (already completed), loading cached metrics...")
                with open(log_file, "r", encoding="utf-8") as f:
                    combined_log = f.read()
                wall = 0.0
                rss = 0.0
                for line in combined_log.splitlines():
                    if "Total Wall Time:" in line:
                        try:
                            # format: Total Wall Time: 1.14s
                            parts = line.split("Total Wall Time:")[1].strip().rstrip("s?")
                            wall = float(parts)
                        except:
                            pass
                    elif "Peak Memory (RSS):" in line:
                        try:
                            parts = line.split("Peak Memory (RSS):")[1].strip().split()[0]
                            rss = float(parts)
                        except:
                            pass
            else:
                cmd = [FASTBACT, "-1", d["r1"], "-2", d["r2"], "-o", out_sub, "-t", "8"] + extra_flags
                code, wall, rss, stdout, stderr = run_cmd_timed(cmd)
                combined_log = stdout + "\n" + stderr
                with open(log_file, "w", encoding="utf-8") as f:
                    f.write(combined_log)
                assert code == 0, f"fastBact failed on {run_name}: {stderr}"

            measured_sha = sha256_file(contigs_fa)
            fb_stats = parse_fastbact_log(combined_log)

            # Run QUAST
            quast_out = os.path.join(out_sub, "quast_out")
            quast_cmd = [
                QUAST,
                "-R", d["ref"],
                contigs_fa,
                "-o", quast_out,
                "--threads", "8",
                "--min-contig", "200",
                "--silent"
            ]
            subprocess.run(quast_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            report_tsv = os.path.join(quast_out, "report.tsv")
            qm = parse_quast_report(report_tsv)

            gf = qm.get("Genome fraction (%)", "N/A")
            nga50 = qm.get("NGA50", qm.get("NG50", "N/A"))
            n50 = qm.get("N50", "N/A")
            largest = qm.get("Largest contig", "N/A")
            cnt_500 = qm.get("# contigs (>= 500 bp)", "N/A")
            ext_mis = qm.get("# misassemblies", "N/A")
            loc_mis = qm.get("# local misassemblies", "N/A")

            input_pairs = fb_stats.get("actual_input_pairs", fb_stats.get("exact_read_pairs", 0))
            retained_pairs = fb_stats.get("actual_retained_pairs", input_pairs)
            ret_pct = f"{(retained_pairs / input_pairs * 100.0):.2f}%" if input_pairs > 0 else "100.00%"

            row = [
                d_id,
                cfg_name,
                f"{d['true_cov']:.1f}",
                str(fb_stats.get("estimated_genome_size", "N/A")),
                f"{fb_stats.get('estimated_coverage', 0.0):.1f}",
                str(fb_stats.get("estimator_confidence", "N/A")),
                f"{fb_stats.get('subsample_factor', 1.0):.4f}",
                str(fb_stats.get("subsampling_enabled", "false")),
                str(fb_stats.get("subsampling_reason", "N/A")),
                str(input_pairs),
                str(retained_pairs),
                ret_pct,
                measured_sha,
                gf,
                nga50,
                n50,
                largest,
                cnt_500,
                ext_mis,
                loc_mis,
                f"{wall:.2f}",
                f"{rss:.1f}",
            ]
            tsv_rows.append("\t".join(row))
            print(f"Finished {run_name}: ret={ret_pct}, wall={wall:.1f}s, rss={rss:.1f}MB, N50={n50}, mis={ext_mis}")

    with open(out_tsv, "w", encoding="utf-8") as f:
        f.write("\n".join(tsv_rows) + "\n")
    print(f"\nWritten benchmark report to {out_tsv}")

if __name__ == "__main__":
    main()
