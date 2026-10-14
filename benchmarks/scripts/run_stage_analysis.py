#!/usr/bin/env python3
"""
Stage-by-Stage Forensics & Telemetry Extractor (FB-BENCH-B1 Phase B1-2)
Runs fastBact with --dump-stages on target datasets, runs QUAST on every intermediate
stage FASTA, and generates benchmarks/B1_reports/B1_STAGE_METRICS.tsv.
"""

import os
import sys
import glob
import json
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

def run_cmd(cmd, desc=""):
    print(f"[run_stage_analysis] {desc}: {' '.join(cmd)}")
    res = subprocess.run(cmd)
    if res.returncode != 0:
        print(f"[WARN] Command failed with exit code {res.returncode}")
    return res.returncode

def main():
    base_dir = "benchmark-results/B1_root_cause/stages"
    os.makedirs(base_dir, exist_ok=True)
    out_tsv = "benchmarks/B1_reports/B1_STAGE_METRICS.tsv"
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)

    header = [
        "dataset_id", "stage_id", "stage_name", "k_value", "input_pairs",
        "retained_pairs", "merged_pairs", "effective_cov", "unique_kmers",
        "solid_kmers", "solid_retention_pct", "cms_overestimate_pct",
        "active_unitigs", "active_edges", "tips_clipped_count", "tips_clipped_bp",
        "bubbles_popped_count", "bubbles_popped_bp", "contigs_ge_200bp",
        "contigs_ge_500bp", "total_len_500bp", "genome_fraction_pct",
        "n50", "nga50", "largest_contig", "ext_misassemblies", "local_misassemblies"
    ]

    tsv_rows = ["\t".join(header)]

    for target in TARGETS:
        d_id = target["id"]
        r1 = target["r1"]
        r2 = target["r2"]
        ref = target["ref"]

        d_dir = os.path.join(base_dir, d_id)
        dumps_dir = os.path.join(d_dir, "dumps")
        tool_out = os.path.join(d_dir, "tool_out")
        os.makedirs(dumps_dir, exist_ok=True)
        os.makedirs(tool_out, exist_ok=True)

        # 1. Run fastBact with --dump-stages if not already dumped
        telem_file = os.path.join(dumps_dir, "stage_telemetry.json")
        if not os.path.exists(telem_file):
            fb_cmd = [
                FASTBACT,
                "-1", r1,
                "-2", r2,
                "-o", tool_out,
                "-t", "8",
                "--dump-stages", dumps_dir
            ]
            run_cmd(fb_cmd, f"Running fastBact --dump-stages on {d_id}")

        # 2. Read stage_telemetry.json
        telemetry = {}
        if os.path.exists(telem_file):
            with open(telem_file, "r") as f:
                try:
                    telemetry = json.load(f)
                except Exception as e:
                    print(f"[WARN] Error reading telemetry {telem_file}: {e}")

        inp_pairs = telemetry.get("input_read_pairs", 0)
        ret_pairs = telemetry.get("retained_read_pairs", 0)
        mrg_pairs = telemetry.get("merged_pairs_count", 0)
        eff_cov = f"{telemetry.get('effective_coverage', 0.0):.2f}"
        uniq_kmers = telemetry.get("total_unique_kmers", 0)
        solid_kmers = telemetry.get("solid_kmers_retained", 0)
        solid_ret_pct = f"{telemetry.get('solid_retention_rate', 0.0) * 100.0:.2f}"
        cms_ov_pct = f"{telemetry.get('cms_vs_exact_agreement', {}).get('overestimation_rate', 0.0) * 100.0:.2f}"

        stage_telem_map = {}
        for st in telemetry.get("stages", []):
            stage_telem_map[st["stage_name"]] = st

        # 3. Find and sort dumped stage fastas
        stage_files = sorted(glob.glob(os.path.join(dumps_dir, "stage_*.fa")))
        print(f"[{d_id}] Found {len(stage_files)} stage files to evaluate")

        for idx, s_fa in enumerate(stage_files, start=1):
            s_base = os.path.splitext(os.path.basename(s_fa))[0]
            s_id = f"S{idx:02d}"

            # QUAST 500bp
            q500_dir = os.path.join(d_dir, f"quast_500_{s_base}")
            r500_tsv = os.path.join(q500_dir, "report.tsv")
            if not os.path.exists(r500_tsv):
                run_cmd([QUAST, "-R", ref, s_fa, "-o", q500_dir, "--min-contig", "500", "-t", "8", "--silent"],
                        f"QUAST 500bp on {s_base}")
            m500 = parse_quast_report(r500_tsv)

            # QUAST 200bp
            q200_dir = os.path.join(d_dir, f"quast_200_{s_base}")
            r200_tsv = os.path.join(q200_dir, "report.tsv")
            if not os.path.exists(r200_tsv):
                run_cmd([QUAST, "-R", ref, s_fa, "-o", q200_dir, "--min-contig", "200", "-t", "8", "--silent"],
                        f"QUAST 200bp on {s_base}")
            m200 = parse_quast_report(r200_tsv)

            st = stage_telem_map.get(s_base, {})
            k_val = st.get("k", 0)
            u_cnt = st.get("active_unitigs", m500.get("# contigs", "0"))
            e_cnt = st.get("active_edges", 0)
            t_cnt = st.get("deleted_tips_count", 0)
            t_bp = st.get("deleted_tips_bp", 0)
            b_cnt = st.get("popped_bubbles_count", 0)
            b_bp = st.get("popped_bubbles_bp", 0)

            cnt_200 = m200.get("# contigs", "0")
            cnt_500 = m500.get("# contigs", "0")
            tot_len = m500.get("Total length", "0")
            gf = m500.get("Genome fraction (%)", "0.0")
            n50 = m500.get("N50", "0")
            nga50 = m500.get("NGA50", m500.get("NG50", n50))
            largest = m500.get("Largest contig", "0")
            ext_mis = m500.get("# misassemblies", "0")
            loc_mis = m500.get("# local misassemblies", "0")

            row = [
                d_id, s_id, s_base, str(k_val), str(inp_pairs),
                str(ret_pairs), str(mrg_pairs), str(eff_cov), str(uniq_kmers),
                str(solid_kmers), str(solid_ret_pct), str(cms_ov_pct),
                str(u_cnt), str(e_cnt), str(t_cnt), str(t_bp),
                str(b_cnt), str(b_bp), str(cnt_200), str(cnt_500),
                str(tot_len), str(gf), str(n50), str(nga50),
                str(largest), str(ext_mis), str(loc_mis)
            ]
            tsv_rows.append("\t".join(row))

    with open(out_tsv, "w", encoding="utf-8") as f:
        f.write("\n".join(tsv_rows) + "\n")
    print(f"[run_stage_analysis] Successfully emitted {out_tsv}")

if __name__ == "__main__":
    main()
