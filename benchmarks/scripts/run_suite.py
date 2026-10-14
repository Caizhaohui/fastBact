#!/usr/bin/env python3
"""
Benchmark Suite Runner for fastBact (FB-BENCH-001)
Runs a set of datasets (filtered by suite/kind or dataset_id) across specified tools,
executing bench_compute.sh for each run. Designed to run on compute node.
"""

import os
import sys
import json
import subprocess
import argparse
from pathlib import Path

def parse_manifest(tsv_path):
    rows = []
    with open(tsv_path, "r", encoding="utf-8") as f:
        header = [h.strip() for h in f.readline().split("\t")]
        for line in f:
            if not line.strip(): continue
            parts = [p.strip() for p in line.split("\t")]
            if len(parts) == len(header):
                rows.append(dict(zip(header, parts)))
    return rows

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign", default="FB-BENCH-v0.1.0")
    parser.add_argument("--suite", choices=["S0_N", "S1", "S2", "R", "ALL"], default="S0_N")
    parser.add_argument("--datasets", nargs="*", default=[])
    parser.add_argument("--tools", nargs="+", default=["fastBact"])
    parser.add_argument("--rep", type=int, default=1)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--manifest", default="benchmarks/datasets.tsv")
    parser.add_argument("--results-dir", default="benchmark-results")
    args = parser.parse_args()

    manifest_rows = parse_manifest(args.manifest)

    # Filter datasets
    target_rows = []
    for r in manifest_rows:
        if r["admission_status"] != "APPROVED":
            continue
        if args.datasets and r["dataset_id"] not in args.datasets:
            continue
        if not args.datasets:
            if args.suite == "S0_N" and r["kind"] not in ["S0", "N"]:
                continue
            elif args.suite == "S1" and r["kind"] != "S1":
                continue
            elif args.suite == "S2" and r["kind"] != "S2":
                continue
            elif args.suite == "R" and r["kind"] not in ["R1", "R2"]:
                continue
        target_rows.append(r)

    print(f"[run_suite] Selected {len(target_rows)} datasets for tools={args.tools}, rep={args.rep}")

    # Determine git commit
    try:
        git_commit = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"]).decode().strip()
    except Exception:
        git_commit = "head"

    all_passed = True
    summary = []

    for row in target_rows:
        d_id = row["dataset_id"]
        is_neg = "1" if row["kind"] == "N" else "0"
        r1 = row["fastq_r1"]
        r2 = row["fastq_r2"]
        ref = row["reference_fasta"]
        truth_events = row["truth_events_file"]

        for tool in args.tools:
            outdir = os.path.join(args.results_dir, args.campaign, d_id, tool, git_commit, f"rep_{args.rep}")
            track = "NEG" if row["kind"] == "N" else ("MICRO" if row["kind"] == "S0" else "E2E")

            cmd = [
                "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/benchmarks/scripts/bench_compute.sh",
                "--campaign", args.campaign,
                "--dataset", d_id,
                "--track", track,
                "--tool", tool,
                "--rep", str(args.rep),
                "--threads", str(args.threads),
                "--seed", str(row["random_seed"]),
                "--outdir", outdir,
                "--r1", r1,
                "--r2", r2,
                "--ref", ref,
                "--truth-events", truth_events,
                "--is-negative", is_neg
            ]

            print(f"\n--- Running {tool} on {d_id} ---")
            p = subprocess.run(cmd)

            # Check status in run.json
            run_json = os.path.join(outdir, "run.json")
            status = "UNKNOWN"
            wall_s = "NA"
            peak_rss = "NA"
            if os.path.exists(run_json):
                with open(run_json, "r") as f:
                    r_data = json.load(f)
                    status = r_data.get("status", "FAIL")
                    wall_s = r_data.get("wall_s", "NA")
                    peak_rss = r_data.get("peak_rss_kib", "NA")

            print(f"Result for {tool} on {d_id}: status={status}, wall={wall_s}s, rss={peak_rss} KiB")
            summary.append({"dataset": d_id, "tool": tool, "status": status, "wall_s": wall_s, "rss_kib": peak_rss})
            if status != "PASS":
                all_passed = False

    print("\n================== SUITE SUMMARY ==================")
    for s in summary:
        print(f"  {s['dataset']:<15} {s['tool']:<10} -> {s['status']:<8} (time: {s['wall_s']}s, rss: {s['rss_kib']} KiB)")
    print("===================================================\n")

    sys.exit(0)

if __name__ == "__main__":
    main()
