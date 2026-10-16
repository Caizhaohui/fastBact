#!/usr/bin/env python3
"""
Benchmark runner for FB-P0-01 Reverse-Strand Polishing Correctness Recovery.
Runs old (949da55) vs new (fix/c0-01-polish-orientation) on 6 datasets:
- S1-COV-150X
- S1-COV-003
- SYN-2.5M
- S2-REAL-001
- R1-001
- R1-003
"""

import os
import sys
import subprocess
import time
import hashlib
import json
import re

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

OLD_BIN = "./target/release/fastBact_old_949da55"
NEW_BIN = "./target/release/fastBact_new_c0_01"
QUAST = "/hpcfs/fhome/caizhh/.local/bin/quast.py"
OUT_BASE = "benchmark-results/C0_01_polish"

def sha256_file(filepath):
    if not os.path.exists(filepath):
        return "MISSING"
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

def run_cmd(cmd):
    t0 = time.time()
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    dt = time.time() - t0
    return res.returncode, dt, res.stdout, res.stderr

def parse_quast(report_path):
    d = {}
    if not os.path.exists(report_path):
        return d
    with open(report_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                d[parts[0]] = parts[1]
    return d

def parse_polish_log(log_path):
    # Consensus polishing completed: 0 bases modified (fwd reads: 6863, rev reads: 6836, rejected anchors: 0)
    res = {
        "bases_modified": 0,
        "fwd_reads": 0,
        "rev_reads": 0,
        "rejected_anchors": 0,
    }
    if not os.path.exists(log_path):
        return res
    with open(log_path) as f:
        for line in f:
            if "Consensus polishing completed:" in line:
                m = re.search(r"(\d+) bases modified \(fwd reads: (\d+), rev reads: (\d+), rejected anchors: (\d+)\)", line)
                if m:
                    res["bases_modified"] = int(m.group(1))
                    res["fwd_reads"] = int(m.group(2))
                    res["rev_reads"] = int(m.group(3))
                    res["rejected_anchors"] = int(m.group(4))
    return res

def read_fasta(p):
    seqs = {}
    name = ""
    cur = []
    if not os.path.exists(p):
        return seqs
    with open(p) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name: seqs[name] = "".join(cur)
                name = line.split()[0][1:]
                cur = []
            else:
                cur.append(line)
        if name: seqs[name] = "".join(cur)
    return seqs

def compare_fasta(old_p, new_p):
    old_seqs = read_fasta(old_p)
    new_seqs = read_fasta(new_p)
    deltas = []
    
    all_keys = sorted(set(old_seqs.keys()) | set(new_seqs.keys()))
    for k in all_keys:
        s_old = old_seqs.get(k, "")
        s_new = new_seqs.get(k, "")
        if s_old == s_new:
            continue
        min_len = min(len(s_old), len(s_new))
        diff_count = sum(1 for i in range(min_len) if s_old[i] != s_new[i]) + abs(len(s_old) - len(s_new))
        for i in range(min_len):
            if s_old[i] != s_new[i]:
                deltas.append({
                    "contig": k,
                    "pos": i,
                    "old_base": s_old[i],
                    "new_base": s_new[i],
                    "len_old": len(s_old),
                    "len_new": len(s_new),
                })
    return deltas

def main():
    print("Starting FB-P0-01 A/B Benchmarking...")
    os.makedirs(f"{OUT_BASE}/old", exist_ok=True)
    os.makedirs(f"{OUT_BASE}/new", exist_ok=True)

    summary_records = []
    all_deltas = {}

    for d in DATASETS:
        d_id = d["id"]
        print(f"\n==================== DATASET: {d_id} ====================")
        
        # 1. Run OLD version
        old_out = f"{OUT_BASE}/old/{d_id}"
        old_log = f"{OUT_BASE}/old/{d_id}.log"
        old_quast_dir = f"{old_out}/quast_out"
        old_fa = f"{old_out}/contigs.fa"
        
        if not (os.path.exists(old_fa) and os.path.exists(f"{old_quast_dir}/report.tsv")):
            print(f"[{d_id}] Running OLD fastBact (949da55)...")
            cmd_old = [OLD_BIN, "-1", d["r1"], "-2", d["r2"], "-o", old_out, "-t", "8"]
            code, wall, out, err = run_cmd(cmd_old)
            with open(old_log, "w") as f:
                f.write(out + "\n" + err)
            if code != 0:
                print(f"ERROR: old fastBact failed on {d_id}")
                sys.exit(1)
            
            print(f"[{d_id}] Running QUAST on OLD...")
            cmd_quast = [QUAST, "-R", d["ref"], old_fa, "-o", old_quast_dir, "--threads", "8", "--min-contig", "200", "--silent"]
            run_cmd(cmd_quast)
        else:
            print(f"[{d_id}] OLD fastBact results already exist.")

        # 2. Run NEW version
        new_out = f"{OUT_BASE}/new/{d_id}"
        new_log = f"{OUT_BASE}/new/{d_id}.log"
        new_quast_dir = f"{new_out}/quast_out"
        new_fa = f"{new_out}/contigs.fa"
        
        if not (os.path.exists(new_fa) and os.path.exists(f"{new_quast_dir}/report.tsv")):
            print(f"[{d_id}] Running NEW fastBact (c0-01 fix)...")
            cmd_new = [NEW_BIN, "-1", d["r1"], "-2", d["r2"], "-o", new_out, "-t", "8"]
            code, wall, out, err = run_cmd(cmd_new)
            with open(new_log, "w") as f:
                f.write(out + "\n" + err)
            if code != 0:
                print(f"ERROR: new fastBact failed on {d_id}")
                sys.exit(1)
            
            print(f"[{d_id}] Running QUAST on NEW...")
            cmd_quast = [QUAST, "-R", d["ref"], new_fa, "-o", new_quast_dir, "--threads", "8", "--min-contig", "200", "--silent"]
            run_cmd(cmd_quast)
        else:
            print(f"[{d_id}] NEW fastBact results already exist.")

        # 3. Collect metrics
        old_sha = sha256_file(old_fa)
        new_sha = sha256_file(new_fa)
        
        old_q = parse_quast(f"{old_quast_dir}/report.tsv")
        new_q = parse_quast(f"{new_quast_dir}/report.tsv")
        
        new_polish = parse_polish_log(new_log)
        
        # Compare FASTA
        deltas = compare_fasta(old_fa, new_fa)
        all_deltas[d_id] = deltas
        
        print(f"[{d_id}] Results summary:")
        print(f"  Old SHA256: {old_sha}")
        print(f"  New SHA256: {new_sha}")
        print(f"  Sequence differences: {len(deltas)} bases")
        print(f"  Old Mismatches/100kbp: {old_q.get('# mismatches per 100 kbp', 'NA')}")
        print(f"  New Mismatches/100kbp: {new_q.get('# mismatches per 100 kbp', 'NA')}")
        print(f"  Old Indels/100kbp:      {old_q.get('# indels per 100 kbp', 'NA')}")
        print(f"  New Indels/100kbp:      {new_q.get('# indels per 100 kbp', 'NA')}")
        print(f"  New Polish Stats:       {new_polish}")

        record = {
            "dataset": d_id,
            "old_sha256": old_sha,
            "new_sha256": new_sha,
            "diff_bases_count": len(deltas),
            "old_n50": old_q.get("N50", "NA"),
            "new_n50": new_q.get("N50", "NA"),
            "old_nga50": old_q.get("NGA50", "NA"),
            "new_nga50": new_q.get("NGA50", "NA"),
            "old_genome_fraction": old_q.get("Genome fraction (%)", "NA"),
            "new_genome_fraction": new_q.get("Genome fraction (%)", "NA"),
            "old_misassemblies": old_q.get("# misassemblies", "NA"),
            "new_misassemblies": new_q.get("# misassemblies", "NA"),
            "old_local_misassemblies": old_q.get("# local misassemblies", "NA"),
            "new_local_misassemblies": new_q.get("# local misassemblies", "NA"),
            "old_mismatches_per_100kbp": old_q.get("# mismatches per 100 kbp", "NA"),
            "new_mismatches_per_100kbp": new_q.get("# mismatches per 100 kbp", "NA"),
            "old_indels_per_100kbp": old_q.get("# indels per 100 kbp", "NA"),
            "new_indels_per_100kbp": new_q.get("# indels per 100 kbp", "NA"),
            "new_bases_modified": new_polish["bases_modified"],
            "new_fwd_pileup_reads": new_polish["fwd_reads"],
            "new_rev_pileup_reads": new_polish["rev_reads"],
            "new_rejected_anchors": new_polish["rejected_anchors"],
        }
        summary_records.append(record)

    with open(f"{OUT_BASE}/benchmark_summary.json", "w") as f:
        json.dump(summary_records, f, indent=2)

    with open(f"{OUT_BASE}/benchmark_deltas.json", "w") as f:
        json.dump(all_deltas, f, indent=2)

    print("\nBenchmark successfully completed!")

if __name__ == "__main__":
    main()
