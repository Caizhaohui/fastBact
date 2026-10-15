#!/usr/bin/env python3
import sys

REF_PATH = "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/references/m6_identity/NC_000913.3.fasta"

def rev_comp(seq):
    return seq.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]

with open(REF_PATH) as f:
    ref = "".join(line.strip() for line in f if not line.startswith(">")).upper()

# 1-based coordinates: 1207806 to 1209602
p_start = 1207806 - 1
p_end = 1209602

p_seq = ref[p_start:p_end]
p_rc = rev_comp(p_seq)

left_flank = ref[p_start - 60 : p_start]
right_flank = ref[p_end : p_end + 60]

print(f"P-region length: {len(p_seq)}")
print(f"Left flank (-60..0):   {left_flank}")
print(f"P-seq start (0..60):   {p_seq[:60]}")
print(f"P-seq end (-60..0):    {p_seq[-60:]}")
print(f"Right flank (0..60):   {right_flank}")

print("\n--- Inverted Repeats Analysis ---")
print(f"P-seq start:           {p_seq[:40]}")
print(f"P-seq end RC:          {rev_comp(p_seq[-40:])}")
print(f"P-seq end:             {p_seq[-40:]}")
print(f"P-seq start RC:        {rev_comp(p_seq[:40])}")

print("\n--- Checking Junction regions ---")
print(f"J1 WT junction (-20 flank + 20 p_seq):")
print(f"  {left_flank[-20:]} | {p_seq[:20]}")
print(f"J2 WT junction (-20 p_seq + 20 right_flank):")
print(f"  {p_seq[-20:]} | {right_flank[:20]}")

print("\n--- Check RC of J2 WT junction ---")
rc_j2_wt = rev_comp(p_seq[-20:] + right_flank[:20])
print(f"RC(J2 WT): {rc_j2_wt[:20]} | {rc_j2_wt[20:]}")
print(f"Left flank: {left_flank[-20:]}")

# Check similarity between left flank and rc(right flank)
print("\n--- Flank comparison ---")
print(f"left_flank[-30:]:     {left_flank[-30:]}")
print(f"rc(right_flank[:30]): {rev_comp(right_flank[:30])}")
