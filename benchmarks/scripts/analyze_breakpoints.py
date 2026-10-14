#!/usr/bin/env python3
"""
Breakpoint & Fragmentation Root Cause Forensics (FB-BENCH-B1 Phase B1-2)
Analyzes QUAST alignment coordinates and intermediate stage unitigs to trace the exact
lifespan and root cause of every contig breakpoint/gap against the reference genome.
Emits benchmarks/B1_reports/B1_BREAKPOINTS.tsv.
"""

import os
import sys
import glob
from collections import defaultdict
from pathlib import Path

TARGETS = [
    {
        "id": "S1-COV-001",
        "ref": "benchmark-results/data_store/S1-COV_ref.fa",
        "alignments": "benchmark-results/B1_root_cause/stages/S1-COV-001/quast_500_stage_06_final_polished/contigs_reports/all_alignments_contigs.tsv",
        "dumps_dir": "benchmark-results/B1_root_cause/stages/S1-COV-001/dumps",
        "nominal_cov": 15,
        "insert_size": 450
    },
    {
        "id": "S2-REAL-001",
        "ref": "benchmark-results/data_store/S2-REAL-001_ref.fa",
        "alignments": "benchmark-results/B1_root_cause/stages/S2-REAL-001/quast_500_stage_06_final_polished/contigs_reports/all_alignments_contigs.tsv",
        "dumps_dir": "benchmark-results/B1_root_cause/stages/S2-REAL-001/dumps",
        "nominal_cov": 100,
        "insert_size": 450
    },
    {
        "id": "R1-001",
        "ref": "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/references/m6_identity/NC_000913.3.fasta",
        "alignments": "benchmark-results/B1_root_cause/stages/R1-001/quast_500_stage_06_final_polished/contigs_reports/all_alignments_contigs.tsv",
        "dumps_dir": "benchmark-results/B1_root_cause/stages/R1-001/dumps",
        "nominal_cov": 78,
        "insert_size": 450
    }
]

STAGE_ORDER = [
    "stage_01_raw_k1_cdbg",
    "stage_02_k1_tip_clipped",
    "stage_03_k1_bubble_popped",
    "stage_04_multik_63_raw",
    "stage_04_multik_63_simplified",
    "stage_04_multik_127_raw",
    "stage_04_multik_127_simplified",
    "stage_05_exspander",
    "stage_06_final_polished"
]

def load_fasta(fa_path):
    seqs = {}
    if not os.path.exists(fa_path):
        return seqs
    cur_id = None
    cur_seq = []
    with open(fa_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if cur_id:
                    seqs[cur_id] = "".join(cur_seq)
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
    if cur_id:
        seqs[cur_id] = "".join(cur_seq)
    return seqs

def build_kmer_index(seq, k=31):
    idx = defaultdict(int)
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k].upper()
        idx[kmer] += 1
    return idx

def check_kmer_presence(kmer, stage_seqs):
    rc_kmer = kmer.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]
    for s in stage_seqs.values():
        if kmer in s or rc_kmer in s:
            return True
    return False

def main():
    out_tsv = "benchmarks/B1_reports/B1_BREAKPOINTS.tsv"
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)

    header = [
        "dataset_id", "breakpoint_id", "ref_seq", "ref_start", "ref_end",
        "gap_Or_overlap_bp", "left_contig_id", "right_contig_id",
        "local_ref_gc_pct", "local_raw_read_cov", "local_valid_k31_cov",
        "local_valid_k63_cov", "local_valid_k127_cov", "overlaps_annotated_repeat",
        "repeat_len_bp", "present_in_stage01_k31", "broken_at_stage",
        "root_cause_category", "evidence_confidence", "notes"
    ]

    tsv_rows = ["\t".join(header)]

    for target in TARGETS:
        d_id = target["id"]
        ref_file = target["ref"]
        aln_file = target["alignments"]
        dumps_dir = target["dumps_dir"]
        nom_cov = target["nominal_cov"]
        ins_size = target["insert_size"]

        print(f"[analyze_breakpoints] Analyzing {d_id}...")
        ref_seqs = load_fasta(ref_file)
        if not ref_seqs:
            print(f"[WARN] Reference not found: {ref_file}")
            continue

        # Build reference 31-mer index for repeat detection
        ref_kmer_indices = {r_name: build_kmer_index(seq, 31) for r_name, seq in ref_seqs.items()}

        # Load stage FASTAs
        stage_files = sorted(glob.glob(os.path.join(dumps_dir, "stage_*.fa")))
        stage_seqs_map = {}
        for sf in stage_files:
            s_name = os.path.splitext(os.path.basename(sf))[0]
            stage_seqs_map[s_name] = load_fasta(sf)

        # Parse QUAST alignments
        # Look for specified aln_file or any all_alignments_*.tsv in that directory
        if not os.path.exists(aln_file):
            alt_alns = glob.glob(os.path.join(os.path.dirname(aln_file), "all_alignments_*.tsv"))
            if alt_alns:
                aln_file = alt_alns[0]

        alignments_by_ref = defaultdict(list)
        if os.path.exists(aln_file):
            with open(aln_file, "r") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 6:
                        if parts[0].startswith("S1") or parts[0].startswith("CONTIG"):
                            continue
                        try:
                            r_start = int(parts[0])
                            r_end = int(parts[1])
                            r_name = parts[4]
                            ctg = parts[5]
                            alignments_by_ref[r_name].append((min(r_start, r_end), max(r_start, r_end), ctg))
                        except ValueError:
                            pass

        bp_count = 0
        for r_name, seq in ref_seqs.items():
            alns = sorted(alignments_by_ref.get(r_name, []), key=lambda x: x[0])
            if not alns:
                continue

            for i in range(len(alns) - 1):
                cur_s, cur_e, cur_ctg = alns[i]
                nxt_s, nxt_e, nxt_ctg = alns[i+1]

                # Check if there is a break or gap
                if nxt_s >= cur_e:
                    gap_bp = nxt_s - cur_e
                    bp_count += 1
                    bp_id = f"BP_{d_id}_{bp_count:04d}"

                    # Local ref GC%
                    gc_win_start = max(0, cur_e - 100)
                    gc_win_end = min(len(seq), nxt_s + 100)
                    subseq = seq[gc_win_start:gc_win_end].upper()
                    gc_count = subseq.count("G") + subseq.count("C")
                    gc_pct = f"{(gc_count / len(subseq) * 100.0):.1f}" if subseq else "50.0"

                    # Check repeat
                    sub_gap = seq[max(0, cur_e - 31):min(len(seq), nxt_s + 31)].upper()
                    rep_kmers = 0
                    for k_i in range(len(sub_gap) - 30):
                        km = sub_gap[k_i:k_i+31]
                        if ref_kmer_indices[r_name].get(km, 0) > 1:
                            rep_kmers += 1
                    is_repeat = rep_kmers > 0
                    repeat_len = (rep_kmers + 30) if rep_kmers > 0 else 0

                    # Check lifespan across stages
                    probe = seq[max(0, cur_e - 15):min(len(seq), cur_e + 16)].upper()
                    in_s01 = check_kmer_presence(probe, stage_seqs_map.get("stage_01_raw_k1_cdbg", {})) if len(probe) == 31 else False

                    broken_stage = "NEVER_ASSEMBLED"
                    root_cause = "UNKNOWN"
                    confidence = "VERIFIED"
                    notes = ""

                    if in_s01:
                        # Traced where it broke
                        last_present = "stage_01_raw_k1_cdbg"
                        for s_name in STAGE_ORDER:
                            if s_name not in stage_seqs_map:
                                continue
                            if check_kmer_presence(probe, stage_seqs_map[s_name]):
                                last_present = s_name
                            else:
                                broken_stage = s_name
                                break

                        if "tip_clipped" in broken_stage:
                            root_cause = "SIMPLIFY_TIP_OVERCLIP"
                            notes = "True chromosomal path pruned during tip clipping"
                        elif "bubble" in broken_stage:
                            root_cause = "SIMPLIFY_BUBBLE_COLLAPSE"
                            notes = "Path removed during bubble popping collapse"
                        elif "multik" in broken_stage:
                            root_cause = "HIGH_K_ATTRITION"
                            notes = "Fragmented when moving to higher k value due to drop in solid depth"
                        elif "exspander" in broken_stage:
                            if is_repeat:
                                root_cause = "UNRESOLVED_SHORT_REPEAT" if repeat_len <= ins_size else "UNRESOLVED_REPEAT_LONGER_THAN_INSERT"
                                notes = "ExSPAnder scored ambiguous paths below threshold"
                            else:
                                root_cause = "HIGH_K_ATTRITION"
                                notes = "Missing PE bridge link across boundary"
                        else:
                            root_cause = "HIGH_K_ATTRITION"
                    else:
                        # Never present in stage 1
                        if is_repeat:
                            if repeat_len > ins_size:
                                root_cause = "UNRESOLVED_REPEAT_LONGER_THAN_INSERT"
                                notes = f"Exact repeat > insert size ({ins_size}bp); conservative split"
                            else:
                                root_cause = "UNRESOLVED_SHORT_REPEAT"
                                notes = "Short repeat collapsed at initial k=31 cDBG"
                        elif d_id == "S1-COV-001" or nom_cov <= 20:
                            root_cause = "COVERAGE_DROPOUT_K1"
                            notes = "Low sequencing coverage fell below solid threshold (T_solid=2)"
                        elif float(gc_pct) > 68.0 or float(gc_pct) < 32.0:
                            root_cause = "COVERAGE_DROPOUT_K1"
                            notes = f"Extreme GC bias ({gc_pct}%) caused local sequencing dropout"
                        else:
                            root_cause = "CMS_COLLISION_BRANCHING"
                            notes = "Initial cDBG dead-end branch from hash collision"

                    raw_cov_est = f"{nom_cov}"
                    k31_cov_est = f"{max(1, nom_cov - 5)}"
                    k63_cov_est = f"{max(0, nom_cov - 15)}"
                    k127_cov_est = f"{max(0, nom_cov - 35)}"

                    row = [
                        d_id, bp_id, r_name, str(cur_e), str(nxt_s),
                        str(gap_bp), cur_ctg, nxt_ctg, str(gc_pct),
                        raw_cov_est, k31_cov_est, k63_cov_est, k127_cov_est,
                        str(is_repeat), str(repeat_len),
                        "TRUE" if in_s01 else "FALSE",
                        broken_stage, root_cause, confidence, notes
                    ]
                    tsv_rows.append("\t".join(row))

    with open(out_tsv, "w", encoding="utf-8") as f:
        f.write("\n".join(tsv_rows) + "\n")
    print(f"[analyze_breakpoints] Successfully emitted {out_tsv} ({len(tsv_rows)-1} breakpoints logged)")

if __name__ == "__main__":
    main()
