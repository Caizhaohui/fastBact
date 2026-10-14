#!/usr/bin/env bash
set -uo pipefail

echo "=== Running Thread Scaling (t=1, 2, 4) ==="

CAMPAIGN="FB-BENCH-v0.1.0"
DATASET="S2-REAL-001"
STORE="benchmark-results/data_store"
R1="${STORE}/${DATASET}_R1.fq.gz"
R2="${STORE}/${DATASET}_R2.fq.gz"
REF="${STORE}/${DATASET}_ref.fa"

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

echo "=== Thread Scaling Finished ==="
