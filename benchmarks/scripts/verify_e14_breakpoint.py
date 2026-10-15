#!/usr/bin/env python3
"""
Verify e14 prophage P-region inversion breakpoints directly against raw FASTQ reads of R1-001 (SRR28349173).
Tests for:
1. Left Junction (J1):
   - J1-WT: NC_000913.3 left flank (..1207805) -> P-region forward start (1207806..)
   - J1-INV: NC_000913.3 left flank (..1207805) -> P-region reverse complement end
2. Right Junction (J2):
   - J2-WT: P-region forward end (..1209602) -> NC_000913.3 right flank (1209603..)
   - J2-INV: P-region reverse complement start -> NC_000913.3 right flank (1209603..)
"""

import gzip
import os
import sys

REF_PATH = "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/references/m6_identity/NC_000913.3.fasta"
R1_PATH = "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349173/SRR28349173_1.fastq.gz"
R2_PATH = "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349173/SRR28349173_2.fastq.gz"

def rev_comp(seq):
    return seq.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]

def main():
    print("Loading NC_000913.3 reference...")
    with open(REF_PATH) as f:
        ref = "".join(line.strip() for line in f if not line.startswith(">")).upper()

    p_start = 1207806 - 1  # 0-indexed 1207805
    p_end = 1209602        # 0-indexed 1209602

    p_seq = ref[p_start:p_end]
    p_rc = rev_comp(p_seq)

    flank_len = 30
    left_flank = ref[p_start - flank_len : p_start]
    right_flank = ref[p_end : p_end + flank_len]

    # Probes across breakpoints: 15bp flank + 15bp P-region (30-mer)
    k = 15
    j1_wt = left_flank[-k:] + p_seq[:k]
    j1_inv = left_flank[-k:] + p_rc[:k]

    j2_wt = p_seq[-k:] + right_flank[:k]
    j2_inv = p_rc[-k:] + right_flank[:k]

    probes = {
        "J1_WT": (j1_wt, rev_comp(j1_wt)),
        "J1_INV": (j1_inv, rev_comp(j1_inv)),
        "J2_WT": (j2_wt, rev_comp(j2_wt)),
        "J2_INV": (j2_inv, rev_comp(j2_inv)),
    }

    print("=== Junction Probes (30-mers) ===")
    for name, (p_fwd, p_rev) in probes.items():
        print(f"{name:8s} Fwd: {p_fwd}  Rev: {p_rev}")

    counts = {name: 0 for name in probes}
    total_reads = 0

    print(f"\nScanning raw reads in {R1_PATH} and {R2_PATH}...")
    for fq_path in [R1_PATH, R2_PATH]:
        with gzip.open(fq_path, "rt") as f:
            while True:
                r_id = f.readline()
                if not r_id:
                    break
                seq = f.readline().strip().upper()
                f.readline()
                f.readline()
                total_reads += 1

                for name, (p_fwd, p_rev) in probes.items():
                    if p_fwd in seq or p_rev in seq:
                        counts[name] += 1

                if total_reads % 1000000 == 0:
                    print(f"Processed {total_reads} reads... Current counts: {counts}")

    print("\n================ FINAL PROBE COUNTS ================")
    print(f"Total raw reads scanned: {total_reads}")
    for name, cnt in counts.items():
        print(f"  {name:8s}: {cnt:6d} supporting raw reads")
    print("====================================================")

    # Output JSON summary
    out_json = "benchmark-results/B1_root_cause/e14_raw_read_counts.json"
    import json
    with open(out_json, "w") as f:
        json.dump({"total_reads": total_reads, "counts": counts}, f, indent=2)
    print(f"Saved counts to {out_json}")

if __name__ == "__main__":
    main()
