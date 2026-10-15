#!/usr/bin/env python3
import os
import sys
import gzip
import json
import random
import hashlib
from pathlib import Path

DNA = ["A", "C", "G", "T"]
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

def rev_comp(seq: str) -> str:
    return "".join(COMPLEMENT.get(b, "N") for b in reversed(seq))

def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()

def make_random_seq(length: int, rng: random.Random, gc_bias: float = 0.50) -> str:
    bases = []
    for _ in range(length):
        if rng.random() < gc_bias:
            bases.append(rng.choice(["C", "G"]))
        else:
            bases.append(rng.choice(["A", "T"]))
    return "".join(bases)

def simulate_reads_pe(
    replicons: list,
    cov_nominal: float,
    read_len: int,
    frag_mean: float,
    frag_sd: float,
    out_r1_gz: str,
    out_r2_gz: str,
    rng: random.Random,
    err_rate: float = 0.005
):
    total_len = sum(len(seq) * copy for _, seq, _, copy in replicons)
    total_pairs = int((total_len * cov_nominal) / (2 * read_len))

    weights = [len(seq) * copy for _, seq, _, copy in replicons]
    total_w = sum(weights)
    norm_w = [w / total_w for w in weights]

    print(f"Generating {total_pairs} read pairs for total_len={total_len} at {cov_nominal}x...")

    CHUNK_PAIRS = 50_000
    with gzip.open(out_r1_gz, "wt", compresslevel=6) as f1, gzip.open(out_r2_gz, "wt", compresslevel=6) as f2:
        r1_buf = []
        r2_buf = []
        for pair_idx in range(1, total_pairs + 1):
            r_choice = rng.choices(replicons, weights=norm_w, k=1)[0]
            _, r_seq, is_circ, _ = r_choice
            L = len(r_seq)
            frag_len = max(read_len + 20, int(rng.gauss(frag_mean, frag_sd)))

            if is_circ:
                start = rng.randint(0, L - 1)
                if start + frag_len <= L:
                    frag = r_seq[start : start + frag_len]
                else:
                    frag = r_seq[start:] + r_seq[: (start + frag_len) % L]
            else:
                if frag_len >= L:
                    frag_len = L
                    start = 0
                else:
                    start = rng.randint(0, L - frag_len)
                frag = r_seq[start : start + frag_len]

            is_fwd = rng.random() < 0.5
            if not is_fwd:
                frag = rev_comp(frag)

            r1_seq = list(frag[:read_len])
            rc_frag = rev_comp(frag)
            r2_seq = list(rc_frag[:read_len])

            for i in range(len(r1_seq)):
                if rng.random() < err_rate:
                    orig = r1_seq[i]
                    r1_seq[i] = rng.choice([b for b in DNA if b != orig])
            for i in range(len(r2_seq)):
                if rng.random() < err_rate:
                    orig = r2_seq[i]
                    r2_seq[i] = rng.choice([b for b in DNA if b != orig])

            q_str = "I" * read_len
            r1_buf.append(f"@read_{pair_idx}/1\n{''.join(r1_seq)}\n+\n{q_str}\n")
            r2_buf.append(f"@read_{pair_idx}/2\n{''.join(r2_seq)}\n+\n{q_str}\n")

            if len(r1_buf) >= CHUNK_PAIRS:
                f1.writelines(r1_buf)
                f2.writelines(r2_buf)
                r1_buf.clear()
                r2_buf.clear()
                print(f"  Written {pair_idx}/{total_pairs} pairs...")

        if r1_buf:
            f1.writelines(r1_buf)
            f2.writelines(r2_buf)

    print(f"Done. Output files:\n  {out_r1_gz}\n  {out_r2_gz}")

def main():
    store_dir = "benchmark-results/data_store"
    os.makedirs(store_dir, exist_ok=True)

    d_id = "SYN-2.5M"
    ref_fa = f"{store_dir}/{d_id}_ref.fa"
    r1_gz = f"{store_dir}/{d_id}_R1.fq.gz"
    r2_gz = f"{store_dir}/{d_id}_R2.fq.gz"
    events_file = f"{store_dir}/{d_id}_truth_events.json"

    rng = random.Random(250)
    genome_len = 2_500_000
    print(f"Generating 2.5 Mb synthetic genome...")
    seq = make_random_seq(genome_len, rng, gc_bias=0.50)

    with open(ref_fa, "w") as f:
        f.write(f">{d_id}_chr\n{seq}\n")

    with open(events_file, "w") as f:
        json.dump({"dataset": d_id, "genome_size": genome_len, "depth": 100}, f, indent=2)

    simulate_reads_pe([("chr", seq, False, 1.0)], 100.0, 150, 450, 30, r1_gz, r2_gz, rng)

    print("Registering in benchmarks/datasets.tsv...")
    manifest_path = "benchmarks/datasets.tsv"
    r1_sha = sha256_file(r1_gz)
    r2_sha = sha256_file(r2_gz)
    ref_sha = sha256_file(ref_fa)

    # Check if already present
    rows = []
    with open(manifest_path, "r") as f:
        rows = [line for line in f]

    has_row = any(line.startswith(f"{d_id}\t") for line in rows)
    if not has_row:
        new_row = f"{d_id}\tS1\tSynthetic 2.5Mb Medium Genome\tinternal://synthetic\tNA\tNA\tNA\t{d_id}.1\tverified\t1_lin\tIllumina\tPE\t150\t450\t30\t100\t100\tbuild_fixtures_1.0\t250\t{r1_gz}\t{r2_gz}\t{ref_fa}\t{events_file}\t{r1_sha}\t{r2_sha}\t{ref_sha}\tAPPROVED\tMedium 2.5Mb bacterial genome 100x benchmark\n"
        with open(manifest_path, "a") as f:
            f.write(new_row)
        print(f"Registered {d_id} in {manifest_path}")
    else:
        print(f"{d_id} already registered in {manifest_path}")

if __name__ == "__main__":
    main()
