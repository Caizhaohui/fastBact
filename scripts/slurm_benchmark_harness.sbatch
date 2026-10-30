#!/usr/bin/env bash
# F26 unified competitor benchmark — Slurm entry.
# Protocol: F26_COMPETITOR_BENCHMARK_CONTRACT.md
#
# Submit from the worktree (logs/ must exist):
#   sbatch scripts/slurm_benchmark_harness.sbatch \
#     --manifest F26_DATASET_MANIFEST.tsv \
#     --dataset R1-001 \
#     --tool all \
#     --outdir benchmarks/f26_results
#
#SBATCH --job-name=f26-competitor-bench
#SBATCH --partition=qfree
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=14G
#SBATCH --time=36:00:00
#SBATCH --output=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f26-robustness-worktree/logs/f26_benchmark_%j.out
#SBATCH --error=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f26-robustness-worktree/logs/f26_benchmark_%j.err

set -euo pipefail

ROOT="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f26-robustness-worktree"

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo "Refusing to run the benchmark harness outside Slurm. Submit this script with sbatch." >&2
  exit 2
fi

cd "${ROOT}"
mkdir -p "${ROOT}/logs"

export OMP_NUM_THREADS=8
export RAYON_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export LC_ALL=C
export TZ=UTC

echo "=== F26 competitor benchmark ==="
echo "job=${SLURM_JOB_ID} node=$(hostname) cpus=${SLURM_CPUS_PER_TASK:-unknown} mem=${SLURM_MEM_PER_NODE:-unknown}"
date -u +"%Y-%m-%dT%H:%M:%SZ"

exec python3 "${ROOT}/scripts/run_competitor_benchmark.py" "$@"
