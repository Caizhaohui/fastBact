#!/usr/bin/env python3
"""
FB-BENCH-B2-C1b-3B Var 2: High-k CMS Genome Size Propagation A/B Benchmark.
Compares Var 1 (Legacy prior contigs sum sizing) vs Var 2 (Explicit genome size propagation)
across the 6 required datasets:
- S1-COV-150X
- S1-COV-003 (300X)
- SYN-2.5M
- S2-REAL-001
- R1-001
- R1-003
Generates B2_C1b3B_VAR2_CMS_CONFIG.tsv and B2_C1b3B_VAR2_ASSEMBLY_AB.tsv.
"""

import os
import sys
import json
import time
import hashlib
import subprocess
from pathlib import Path

FASTBACT = "./target/release/fastBact"
QUAST = os.path.expanduser("~/.local/bin/quast.py")
WORK_DIR = "benchmark-results/B2_C1b/var2_ab_runs"
REPORTS_DIR = "benchmarks/B2_reports"

DATASETS = [
    {
        "id": "S1-COV-150X",
        "r1": "benchmark-results/data_store/S1-COV-150X_R1.fq.gz",
        "r2": "benchmark-results/data_store/S1-COV-150X_R2.fq.gz",
        "ref": "benchmark-results/data_store/S1-COV_ref.fa",
    },
    {
        "id": "S1-COV-003",
        "r1": "benchmark-results/data_store/S1-COV-003_R1.fq.gz",
        "r2": "benchmark-results/data_store/S1-COV-003_R2.fq.gz",
        "ref": "benchmark-results/data_store/S1-COV_ref.fa",
    },
    {
        "id": "SYN-2.5M",
        "r1": "benchmark-results/data_store/SYN-2.5M_R1.fq.gz",
        "r2": "benchmark-results/data_store/SYN-2.5M_R2.fq.gz",
        "ref": "benchmark-results/data_store/SYN-2.5M_ref.fa",
    },
    {
        "id": "S2-REAL-001",
        "r1": "benchmark-results/data_store/S2-REAL-001_R1.fq.gz",
        "r2": "benchmark-results/data_store/S2-REAL-001_R2.fq.gz",
        "ref": "benchmark-results/data_store/S2-REAL-001_ref.fa",
    },
    {
        "id": "R1-001",
        "r1": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349173/SRR28349173_1.fastq.gz",
        "r2": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349173/SRR28349173_2.fastq.gz",
        "ref": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/references/m6_identity/NC_000913.3.fasta",
    },
    {
        "id": "R1-003",
        "r1": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349164/SRR28349164_1.fastq.gz",
        "r2": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349164/SRR28349164_2.fastq.gz",
        "ref": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/references/m6_identity/CP009273.1.fasta",
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

def parse_multik_telemetry_from_log(log_text):
    records = []
    # [MULTI-K TELEMETRY] k=63: executed=true, early_stop=None, CMS source=..., size_used=... bp, width=..., bytes=... MB, unitigs: 719 -> 2
    for line in log_text.splitlines():
        if "[MULTI-K TELEMETRY]" in line:
            parts = line.split("[MULTI-K TELEMETRY]")[1].strip()
            # parse key-values
            # k=63: executed=true, early_stop=None, CMS source=..., size_used=... bp, width=..., bytes=... MB, unitigs: ...
            tokens = parts.split(":")
            k_val = int(tokens[0].split("=")[1])
            rest = ":".join(tokens[1:]).strip()
            records.append((k_val, rest))
    return records

def main():
    os.makedirs(WORK_DIR, exist_ok=True)
    os.makedirs(REPORTS_DIR, exist_ok=True)

    cms_config_rows = []
    cms_config_header = [
        "dataset",
        "variant",
        "stage_name",
        "k",
        "genome_size_source",
        "genome_size_used",
        "cms_width",
        "cms_allocated_bytes",
        "cms_allocated_mb",
        "active_unitigs",
        "active_edges",
        "upgrade_executed",
        "early_stop_reason",
    ]
    cms_config_rows.append("\t".join(cms_config_header))

    assembly_ab_rows = []
    assembly_ab_header = [
        "dataset",
        "variant",
        "n50_bp",
        "nga50_bp",
        "contigs_ge_500bp",
        "genome_fraction_pct",
        "extensive_misassemblies",
        "local_misassemblies",
        "peak_rss_mb",
        "wall_time_s",
        "k31_cms_mb",
        "k63_cms_mb",
        "k127_cms_mb",
        "k31_unitigs",
        "k63_unitigs",
        "k127_unitigs",
        "contigs_sha256",
        "sha256_match_vs_var1",
        "behavior_delta_note",
    ]
    assembly_ab_rows.append("\t".join(assembly_ab_header))

    results_by_dataset = {}

    for d in DATASETS:
        d_id = d["id"]
        results_by_dataset[d_id] = {}
        print(f"\n=======================================================")
        print(f"BENCHMARKING DATASET: {d_id}")
        print(f"=======================================================")

        for variant in ["Var1_Legacy", "Var2_Propagated"]:
            out_sub = os.path.join(WORK_DIR, f"{d_id}_{variant}")
            os.makedirs(out_sub, exist_ok=True)
            dump_sub = os.path.join(out_sub, "dumps")
            os.makedirs(dump_sub, exist_ok=True)

            cmd = [
                FASTBACT,
                "-1", d["r1"],
                "-2", d["r2"],
                "-o", out_sub,
                "--dump-stages", dump_sub,
                "-t", "8",
            ]
            if variant == "Var1_Legacy":
                cmd.append("--legacy-multik-cms")

            print(f"\nRunning {d_id} under {variant}...")
            code, wall_s, peak_rss, stdout, stderr = run_cmd_timed(cmd)
            combined_log = stdout + "\n" + stderr
            if code != 0:
                print(f"ERROR: {d_id} {variant} failed:\n{stderr}")
                sys.exit(1)

            contigs_fa = os.path.join(out_sub, "contigs.fa")
            sha256 = sha256_file(contigs_fa)
            print(f"  {variant} completed in {wall_s:.2f}s (Peak RSS: {peak_rss:.1f} MB), SHA256: {sha256}")

            # Run QUAST
            quast_out = os.path.join(out_sub, "quast_out")
            quast_cmd = [
                QUAST,
                "-R", d["ref"],
                contigs_fa,
                "-o", quast_out,
                "--threads", "8",
                "--min-contig", "200",
                "--silent",
            ]
            q_code, q_wall, q_rss, q_out, q_err = run_cmd_timed(quast_cmd)
            qm = parse_quast_report(os.path.join(quast_out, "report.tsv"))

            n50 = qm.get("N50", "N/A")
            nga50 = qm.get("NGA50", qm.get("NG50", "N/A"))
            contigs_ge_500 = qm.get("# contigs (>= 500 bp)", qm.get("# contigs", "N/A"))
            gf = qm.get("Genome fraction (%)", "N/A")
            mis = qm.get("# misassemblies", "0")
            local_mis = qm.get("# local misassemblies", "0")

            # Parse stage_telemetry.json
            telem_file = os.path.join(dump_sub, "stage_telemetry.json")
            stages_data = []
            if os.path.exists(telem_file):
                with open(telem_file, "r") as f:
                    try:
                        t_json = json.load(f)
                        stages_data = t_json.get("stages", [])
                    except Exception as e:
                        print(f"Warning: could not parse {telem_file}: {e}")

            k_cms = {}
            k_u = {}
            for st in stages_data:
                s_name = st.get("stage_name", "")
                k_val = st.get("k", 0)
                g_source = st.get("genome_size_source", "UNKNOWN")
                g_used = st.get("genome_size_used", 0)
                cms_w = st.get("cms_width", 0)
                cms_b = st.get("cms_allocated_bytes", 0)
                cms_mb = cms_b / (1024 * 1024) if cms_b > 0 else 0
                u = st.get("active_unitigs", 0)
                e = st.get("active_edges", 0)
                exec_flag = st.get("upgrade_executed", False)
                stop_reason = st.get("early_stop_reason", "NONE")

                # Track per-k stats for assembly table
                if k_val not in k_cms and cms_mb > 0:
                    k_cms[k_val] = cms_mb
                if k_val not in k_u:
                    k_u[k_val] = u

                # Record CMS config row
                if "raw" in s_name or s_name == "stage_01_raw_k1_cdbg":
                    cms_config_rows.append("\t".join([
                        d_id,
                        variant,
                        s_name,
                        str(k_val),
                        g_source,
                        str(g_used),
                        str(cms_w),
                        str(cms_b),
                        f"{cms_mb:.1f}",
                        str(u),
                        str(e),
                        str(exec_flag),
                        stop_reason,
                    ]))

            results_by_dataset[d_id][variant] = {
                "n50": n50,
                "nga50": nga50,
                "contigs_ge_500": contigs_ge_500,
                "gf": gf,
                "mis": mis,
                "local_mis": local_mis,
                "peak_rss": peak_rss,
                "wall_s": wall_s,
                "sha256": sha256,
                "k_cms": k_cms,
                "k_u": k_u,
            }

        # Compare Var 1 vs Var 2
        r1 = results_by_dataset[d_id]["Var1_Legacy"]
        r2 = results_by_dataset[d_id]["Var2_Propagated"]

        sha_match = (r1["sha256"] == r2["sha256"])
        delta_note = "BIT_IDENTICAL" if sha_match else "ALGORITHM_BEHAVIOR_DELTA"

        for variant in ["Var1_Legacy", "Var2_Propagated"]:
            res = results_by_dataset[d_id][variant]
            match_str = "BASE_REFERENCE" if variant == "Var1_Legacy" else ("100%_BIT_IDENTICAL" if sha_match else "DIFFERENT_SHA256")
            note_str = "LEGACY_BASELINE" if variant == "Var1_Legacy" else delta_note

            assembly_ab_rows.append("\t".join([
                d_id,
                variant,
                str(res["n50"]),
                str(res["nga50"]),
                str(res["contigs_ge_500"]),
                str(res["gf"]),
                str(res["mis"]),
                str(res["local_mis"]),
                f"{res['peak_rss']:.1f}",
                f"{res['wall_s']:.2f}",
                f"{res['k_cms'].get(31, 0):.1f} MB",
                f"{res['k_cms'].get(63, 0):.1f} MB",
                f"{res['k_cms'].get(127, 0):.1f} MB",
                str(res["k_u"].get(31, 'N/A')),
                str(res["k_u"].get(63, 'N/A')),
                str(res["k_u"].get(127, 'N/A')),
                res["sha256"],
                match_str,
                note_str,
            ]))

    cms_tsv_path = os.path.join(REPORTS_DIR, "B2_C1b3B_VAR2_CMS_CONFIG.tsv")
    with open(cms_tsv_path, "w") as f:
        f.write("\n".join(cms_config_rows) + "\n")
    print(f"\nSaved {cms_tsv_path}")

    assembly_tsv_path = os.path.join(REPORTS_DIR, "B2_C1b3B_VAR2_ASSEMBLY_AB.tsv")
    with open(assembly_tsv_path, "w") as f:
        f.write("\n".join(assembly_ab_rows) + "\n")
    print(f"Saved {assembly_tsv_path}")

    print("\n=======================================================")
    print("B2-C1b-3B Var 2 A/B Benchmark Execution Completed Successfully!")
    print("=======================================================")

if __name__ == "__main__":
    main()
