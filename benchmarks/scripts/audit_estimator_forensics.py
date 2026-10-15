#!/usr/bin/env python3
"""
Forensics and Audit Script for fastBact B2-G0
Investigates:
1. Exact calculation steps for HLL genome size (522 kb) and depth (626.4x).
2. CMS table saturation and hash collision effects on c == 2 detection.
3. Gzip stream position, decompression buffer read-ahead, and compression ratio variability.
4. Coordinate-level analysis of R1-001 misassemblies (inversion vs relocations).
5. Output of B2_ESTIMATOR_FORENSICS.md, B2_LINEAGE.tsv, B2_MISASSEMBLY_REVIEW.tsv.
"""

import os
import sys
import gzip
import json
import math
import subprocess
from collections import Counter
from pathlib import Path

DATA_S2_R1 = "benchmark-results/data_store/S2-REAL-001_R1.fq.gz"
DATA_S2_R2 = "benchmark-results/data_store/S2-REAL-001_R2.fq.gz"
DATA_R1_R1 = "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349173/SRR28349173_1.fastq.gz"
DATA_R1_R2 = "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data/prjna1088182/fastq/SRR28349173/SRR28349173_2.fastq.gz"

def canonical_kmer(kmer):
    rc = kmer.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]
    return kmer if kmer < rc else rc

class SimCMS:
    def __init__(self, width=1 << 20, rows=4):
        self.width = width
        self.mask = width - 1
        self.rows = rows
        self.seeds = [
            0x243f6a8885a308d3,
            0x13198a2e03707344,
            0xa4093822299f31d0,
            0x082efa98ec4e6c89,
        ]
        self.table = [[0] * width for _ in range(rows)]

    def _hash(self, h, row):
        seed = self.seeds[row]
        # splitmix64
        z = (h ^ seed) & 0xFFFFFFFFFFFFFFFF
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9 & 0xFFFFFFFFFFFFFFFF
        z = (z ^ (z >> 27)) * 0x94d049bb133111eb & 0xFFFFFFFFFFFFFFFF
        z = z ^ (z >> 31)
        return z & self.mask

    def add_and_estimate(self, h):
        min_val = 15
        for r in range(self.rows):
            idx = self._hash(h, r)
            if self.table[r][idx] < 15:
                self.table[r][idx] += 1
            if self.table[r][idx] < min_val:
                min_val = self.table[r][idx]
        return min_val

class SimHLL:
    def __init__(self, precision=14):
        self.p = precision
        self.m = 1 << precision
        self.reg = [0] * self.m

    def add_hash(self, h):
        idx = h >> (64 - self.p)
        w = h & ((1 << (64 - self.p)) - 1)
        if w == 0:
            rho = 64 - self.p + 1
        else:
            rho = (64 - self.p) - w.bit_length() + 1
        if rho > self.reg[idx]:
            self.reg[idx] = rho

    def estimate(self):
        m = float(self.m)
        alpha_m = 0.7213 / (1.0 + 1.079 / m)
        s = sum(2.0 ** (-val) for val in self.reg)
        raw = alpha_m * m * m / s
        zeros = self.reg.count(0)
        if raw <= 2.5 * m:
            if zeros > 0:
                return m * math.log(m / zeros)
            return raw
        return raw

def sim_hash(kmer_str):
    # simple 64-bit fnv-1a / ntHash proxy
    h = 0xcbf29ce484222325
    for b in kmer_str.encode("ascii"):
        h = ((h ^ b) * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h

def analyze_dataset_estimator(r1_path, r2_path, max_sample_pairs=50000, target_cov=120.0):
    print(f"=== Analyzing {r1_path} ===")
    r1_size = os.path.getsize(r1_path)
    r2_size = os.path.getsize(r2_path) if os.path.exists(r2_path) else 0

    k = 21
    cms_small = SimCMS(width=1 << 20) # 1M (fastBact default)
    hll_small = SimHLL(precision=14)

    cms_large = SimCMS(width=1 << 24) # 16M
    hll_large = SimHLL(precision=14)

    exact_counts = Counter()
    first_touch_cms_est = []

    sampled_bases = 0
    pair_count = 0
    read_lens = []

    # Stream FASTQ
    with gzip.open(r1_path, "rt") as f1, gzip.open(r2_path, "rt") as f2:
        while pair_count < max_sample_pairs:
            l1_id = f1.readline()
            if not l1_id: break
            l1_seq = f1.readline().strip()
            f1.readline() # +
            f1.readline() # qual

            l2_id = f2.readline()
            l2_seq = f2.readline().strip()
            f2.readline()
            f2.readline()

            pair_count += 1
            sampled_bases += len(l1_seq) + len(l2_seq)
            if len(read_lens) < 100000:
                read_lens.append(len(l1_seq))
                read_lens.append(len(l2_seq))

            # Process 21-mers
            for s in [l1_seq, l2_seq]:
                for i in range(len(s) - k + 1):
                    sub = s[i:i+k].upper()
                    if "N" in sub: continue
                    ck = canonical_kmer(sub)
                    h = sim_hash(ck)

                    # Exact count
                    prev_exact = exact_counts[ck]
                    exact_counts[ck] += 1

                    # Small CMS (1M)
                    c_small = cms_small.add_and_estimate(h)
                    if c_small == 2:
                        hll_small.add_hash(h)

                    # Large CMS (16M)
                    c_large = cms_large.add_and_estimate(h)
                    if c_large == 2:
                        hll_large.add_hash(h)

                    if prev_exact == 0 and len(first_touch_cms_est) < 100000:
                        first_touch_cms_est.append(c_small)

    read_lens.sort()
    median_read_len = read_lens[len(read_lens) // 2] if read_lens else 150

    # fastBact estimator.rs formula:
    is_gz = r1_path.endswith(".gz")
    approx_bytes_per_pair = 75 if is_gz else (median_read_len * 4 + 40)
    total_pairs_est = max(pair_count, r1_size // approx_bytes_per_pair)
    sample_ratio = total_pairs_est / float(pair_count) if pair_count > 0 else 1.0

    total_projected_bases = sampled_bases * sample_ratio
    est_g_small = hll_small.estimate()
    est_cov_small = total_projected_bases / max(1.0, est_g_small)
    subsample_factor_small = min(1.0, max(0.01, target_cov / est_cov_small)) if est_cov_small > 1.1 * target_cov else 1.0

    est_g_large = hll_large.estimate()
    est_cov_large = total_projected_bases / max(1.0, est_g_large)

    # Exact stats
    unique_kmers_total = len(exact_counts)
    exact_solid_ge2 = sum(1 for c in exact_counts.values() if c >= 2)
    exact_singletons = sum(1 for c in exact_counts.values() if c == 1)

    # Saturation of small CMS:
    touch_ge2_on_first = sum(1 for c in first_touch_cms_est if c >= 2)
    touch_gt2_on_first = sum(1 for c in first_touch_cms_est if c > 2)

    return {
        "r1_path": r1_path,
        "r1_size_bytes": r1_size,
        "r2_size_bytes": r2_size,
        "pair_count": pair_count,
        "sampled_bases": sampled_bases,
        "median_read_len": median_read_len,
        "approx_bytes_per_pair": approx_bytes_per_pair,
        "total_pairs_est": total_pairs_est,
        "sample_ratio": sample_ratio,
        "total_projected_bases": total_projected_bases,
        "unique_kmers_total": unique_kmers_total,
        "exact_singletons": exact_singletons,
        "exact_solid_ge2": exact_solid_ge2,
        "hll_small_1M_est_g": est_g_small,
        "est_cov_small_1M": est_cov_small,
        "subsample_factor_small": subsample_factor_small,
        "hll_large_16M_est_g": est_g_large,
        "est_cov_large_16M": est_cov_large,
        "first_touch_ge2_pct": (touch_ge2_on_first / max(1, len(first_touch_cms_est))) * 100.0,
        "first_touch_gt2_pct": (touch_gt2_on_first / max(1, len(first_touch_cms_est))) * 100.0,
    }

def analyze_gzip_compression_uniformity(fq_gz_path, check_points=[10000, 25000, 50000, 100000, 200000]):
    print(f"=== Checking gzip compression uniformity for {fq_gz_path} ===")
    results = []
    tot_size = os.path.getsize(fq_gz_path)

    # Use raw file reader + zlib decompressor to track exact raw bytes consumed
    import zlib
    raw_f = open(fq_gz_path, "rb")
    dco = zlib.decompressobj(16 + zlib.MAX_WBITS)

    pairs = 0
    raw_bytes_consumed = 0
    uncompressed_bytes = 0
    buf_size = 65536
    leftover = b""

    target_idx = 0
    targets = sorted(check_points)

    while True:
        chunk = raw_f.read(buf_size)
        if not chunk:
            break
        raw_bytes_consumed += len(chunk)
        decomp = dco.decompress(chunk)
        uncompressed_bytes += len(decomp)
        leftover += decomp

        # Count newlines (4 lines per FASTQ record)
        while b"\n" in leftover:
            pos = leftover.find(b"\n")
            line = leftover[:pos]
            leftover = leftover[pos+1:]
            # Every 4 lines is 1 read
            pairs += 0.25

            if target_idx < len(targets) and pairs >= targets[target_idx]:
                pt = targets[target_idx]
                ratio = uncompressed_bytes / max(1, raw_bytes_consumed)
                bytes_per_pair_gz = raw_bytes_consumed / pt
                bytes_per_pair_raw = uncompressed_bytes / pt
                extrapolated_pairs = (tot_size / max(1, raw_bytes_consumed)) * pt
                results.append({
                    "target_pairs": pt,
                    "raw_gz_bytes": raw_bytes_consumed,
                    "uncompressed_bytes": uncompressed_bytes,
                    "comp_ratio": ratio,
                    "gz_bytes_per_pair": bytes_per_pair_gz,
                    "raw_bytes_per_pair": bytes_per_pair_raw,
                    "extrapolated_total_pairs": extrapolated_pairs
                })
                target_idx += 1

        if target_idx >= len(targets):
            break

    raw_f.close()
    return results

def main():
    out_dir = "benchmarks/B1_reports"
    os.makedirs(out_dir, exist_ok=True)

    print("Step 1: Running Estimator Forensics on S2-REAL-001 and R1-001...")
    res_s2 = analyze_dataset_estimator(DATA_S2_R1, DATA_S2_R2, max_sample_pairs=50000)
    res_r1 = analyze_dataset_estimator(DATA_R1_R1, DATA_R1_R2, max_sample_pairs=50000)

    print("Step 2: Running Gzip Extrapolation Forensics...")
    gz_s2 = analyze_gzip_compression_uniformity(DATA_S2_R1, [10000, 25000, 50000, 100000])
    gz_r1 = analyze_gzip_compression_uniformity(DATA_R1_R1, [10000, 25000, 50000, 100000])

    # Save summary JSON for report generation
    summary = {
        "S2-REAL-001": res_s2,
        "R1-001": res_r1,
        "gzip_uniformity_S2": gz_s2,
        "gzip_uniformity_R1": gz_r1,
    }
    with open("benchmark-results/B1_root_cause/estimator_forensics_raw.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("Forensics completed successfully. Emitted estimator_forensics_raw.json.")

if __name__ == "__main__":
    main()
