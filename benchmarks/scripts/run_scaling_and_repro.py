#!/usr/bin/env bash
set -uo pipefail

echo "=== Running Phase B4: Thread Scaling & Reproducibility Matrix ==="

CAMPAIGN="FB-BENCH-v0.1.0"
DATASET="S2-REAL-001"
STORE="benchmark-results/data_store"
R1="${STORE}/${DATASET}_R1.fq.gz"
R2="${STORE}/${DATASET}_R2.fq.gz"
REF="${STORE}/${DATASET}_ref.fa"

# 1. Thread Scaling: threads = 1, 2, 4 (threads = 8 already done)
for T in 1 2 4; do
  OUTDIR="benchmark-results/${CAMPAIGN}/${DATASET}_scaling/fastBact_t${T}/f328d07/rep_1"
  echo "--- Running fastBact with ${T} threads ---"
  benchmarks/scripts/bench_compute.sh \
    --campaign "${CAMPAIGN}" \
    --dataset "${DATASET}_scaling" \
    --track "E2E" \
    --tool "fastBact_t${T}" \
    --rep 1 \
    --threads "${T}" \
    --seed 301 \
    --outdir "${OUTDIR}" \
    --r1 "${R1}" \
    --r2 "${R2}" \
    --ref "${REF}" \
    --truth-events "NA" \
    --is-negative 0
done

# Copy thread 8 result into scaling structure
mkdir -p "benchmark-results/${CAMPAIGN}/${DATASET}_scaling/fastBact_t8/f328d07/rep_1"
cp -r "benchmark-results/${CAMPAIGN}/${DATASET}/fastBact/f328d07/rep_1/"* "benchmark-results/${CAMPAIGN}/${DATASET}_scaling/fastBact_t8/f328d07/rep_1/" || true

# 2. Reproducibility Repetitions (Rep 2 and Rep 3 for fastBact on S2-REAL-001)
for R in 2 3; do
  OUTDIR="benchmark-results/${CAMPAIGN}/${DATASET}/fastBact/f328d07/rep_${R}"
  echo "--- Running fastBact repetition ${R} ---"
  benchmarks/scripts/bench_compute.sh \
    --campaign "${CAMPAIGN}" \
    --dataset "${DATASET}" \
    --track "E2E" \
    --tool "fastBact" \
    --rep "${R}" \
    --threads 8 \
    --seed 301 \
    --outdir "${OUTDIR}" \
    --r1 "${R1}" \
    --r2 "${R2}" \
    --ref "${REF}" \
    --truth-events "NA" \
    --is-negative 0
done

# 3. Real dataset baseline: run Shovill on R1-001
echo "--- Running Shovill baseline on R1-001 ---"
benchmarks/scripts/run_suite.py --datasets R1-001 --tools shovill

echo "=== Phase B4 Completed Successfully ==="
