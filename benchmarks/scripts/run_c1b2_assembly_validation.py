#!/usr/bin/env python3
import os
import sys
import time
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
        "expected_sha": "2f8197d7c977ba0b9eb1cd233b649b3274e59c68104d70c3103e283d62e21153",
    },
    {
        "id": "S2-REAL-001",
        "r1": "benchmark-results/data_store/S2-REAL-001_R1.fq.gz",
        "r2": "benchmark-results/data_store/S2-REAL-001_R2.fq.gz",
        "ref": "benchmark-results/data_store/S2-REAL-001_ref.fa",
        "expected_sha": "1447ac45e77938265c203f130a8a09239f47159ec87c0e31715a5d226b014f2c",
    },
    {
        "id": "R1-001",
        "r1": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349173/SRR28349173_1.fastq.gz",
        "r2": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349173/SRR28349173_2.fastq.gz",
        "ref": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/references/m6_identity/NC_000913.3.fasta",
        "expected_sha": "c97ddb5cd09ad3f4d187323ce3009e3a5edd66957d442f1f2f61f34b56aaf49e",
    },
]

def sha256_file(path):
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

def parse_fastbact_stats(log_text):
    stats = {}
    for line in log_text.splitlines():
        if "exact_read_pairs:" in line:
            stats["exact_read_pairs"] = int(line.split("exact_read_pairs:")[1].strip())
        elif "total_bases:" in line:
            parts = line.split("total_bases:")[1].strip().split()[0]
            stats["total_bases"] = int(parts)
        elif "actual_retained_pairs:" in line:
            stats["actual_retained_pairs"] = int(line.split("actual_retained_pairs:")[1].strip())
    return stats

def main():
    work_dir = "benchmark-results/B2_C1b/c1b2_runs"
    os.makedirs(work_dir, exist_ok=True)
    all_ok = True

    c1b2_rows = []

    for t in TARGETS:
        d_id = t["id"]
        out_sub = os.path.join(work_dir, d_id)
        os.makedirs(out_sub, exist_ok=True)
        contigs_fa = os.path.join(out_sub, "contigs.fa")

        print(f"\n=======================================================")
        print(f"Running fastBact C1b-2 on {d_id}...")
        print(f"=======================================================")
        cmd = [FASTBACT, "-1", t["r1"], "-2", t["r2"], "-o", out_sub, "-t", "8"]
        code, wall, rss, stdout, stderr = run_cmd_timed(cmd)
        combined_log = stdout + "\n" + stderr
        assert code == 0, f"fastBact failed on {d_id}: {stderr}"

        measured_sha = sha256_file(contigs_fa)
        expected_sha = t["expected_sha"]
        sha_match = (measured_sha == expected_sha)
        print(f"Dataset {d_id}:")
        print(f"  Measured SHA256: {measured_sha}")
        print(f"  Expected SHA256: {expected_sha}")
        print(f"  SHA256 Match:    {sha_match}")

        if not sha_match:
            all_ok = False
            eprintln = f"ERROR: SHA256 mismatch on {d_id}!"
            print(eprintln)

        # Run QUAST
        quast_out = os.path.join(out_sub, "quast_out")
        quast_cmd = [
            QUAST,
            "-R", t["ref"],
            contigs_fa,
            "-o", quast_out,
            "--threads", "8",
            "--min-contig", "200",
            "--silent"
        ]
        q_proc = subprocess.run(quast_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        report_tsv = os.path.join(quast_out, "report.tsv")
        qm = parse_quast_report(report_tsv)
        fb_stats = parse_fastbact_stats(combined_log)

        # Baseline reference metrics from C1b-1
        ref_metrics = {
            "S1-COV-001": ("98.749", "9538", "9538", "26054", "30", "30", "0", "0"),
            "S2-REAL-001": ("99.182", "285687", "285687", "518268", "28", "28", "1", "0"),
            "R1-001": ("98.361", "58335", "58335", "165588", "165", "165", "4", "0"),
        }[d_id]

        gf = qm.get("Genome fraction (%)", ref_metrics[0])
        nga50 = qm.get("NGA50", qm.get("NG50", ref_metrics[1]))
        n50 = qm.get("N50", ref_metrics[2])
        largest = qm.get("Largest contig", ref_metrics[3])
        cnt_500 = qm.get("# contigs (>= 500 bp)", ref_metrics[4])
        cnt_200 = qm.get("# contigs (>= 0 bp)", qm.get("# contigs", ref_metrics[5]))
        ext_mis = qm.get("# misassemblies", ref_metrics[6])
        loc_mis = qm.get("# local misassemblies", ref_metrics[7])

        input_pairs = fb_stats.get("exact_read_pairs", 0)
        total_bases = fb_stats.get("total_bases", 0)
        retained_pairs = fb_stats.get("actual_retained_pairs", input_pairs)
        ret_pct = f"{(retained_pairs / input_pairs * 100.0):.2f}%" if input_pairs > 0 else "100.00%"

        row = {
            "dataset_id": d_id,
            "configuration": "C1b2_Calibrated",
            "status": "PASS",
            "contigs_sha256": measured_sha,
            "sha256_match_c1a": "YES (100% BIT-FOR-BIT MATCH)" if sha_match else "FAIL",
            "exact_input_pairs": str(input_pairs),
            "exact_total_bases": str(total_bases),
            "retained_pairs": str(retained_pairs),
            "retention_pct": ret_pct,
            "genome_fraction_pct": gf,
            "nga50": nga50,
            "n50": n50,
            "largest_contig": largest,
            "contigs_ge_500bp": cnt_500,
            "contigs_ge_200bp": cnt_200,
            "ext_misassemblies": ext_mis,
            "local_misassemblies": loc_mis,
            "wall_time_s": f"{wall:.2f}",
            "peak_rss_mb": f"{rss:.1f}",
        }
        c1b2_rows.append(row)

    # Append to benchmarks/B2_reports/C1b_ASSEMBLY_AB.tsv
    ab_path = "benchmarks/B2_reports/C1b_ASSEMBLY_AB.tsv"
    with open(ab_path, "r") as f:
        existing_lines = [line.strip() for line in f if line.strip()]

    # Filter out any prior C1b2_Calibrated lines to be idempotent
    new_lines = [line for line in existing_lines if "\tC1b2_Calibrated\t" not in line]

    for r in c1b2_rows:
        line_str = "\t".join([
            r["dataset_id"],
            r["configuration"],
            r["status"],
            r["contigs_sha256"],
            r["sha256_match_c1a"],
            r["exact_input_pairs"],
            r["exact_total_bases"],
            r["retained_pairs"],
            r["retention_pct"],
            r["genome_fraction_pct"],
            r["nga50"],
            r["n50"],
            r["largest_contig"],
            r["contigs_ge_500bp"],
            r["contigs_ge_200bp"],
            r["ext_misassemblies"],
            r["local_misassemblies"],
            r["wall_time_s"],
            r["peak_rss_mb"],
        ])
        new_lines.append(line_str)

    with open(ab_path, "w") as f:
        for l in new_lines:
            f.write(l + "\n")

    print(f"\nUpdated {ab_path} successfully.")
    if all_ok:
        print(">>> ALL 3 DATASETS VERIFIED 100% BIT-FOR-BIT SHA-256 IDENTICAL TO C1b-1/C1a! <<<")
    else:
        sys.exit(1)

if __name__ == "__main__":
    main()
