#!/usr/bin/env bash
set -uo pipefail

# HPC Submission Wrapper for fastBact Benchmarking Suite (FB-BENCH-001)
# Submits jobs synchronously to compute node via sbatch --wait (partition: qfree, node: vnode1)
# Zero CPU / heavy memory used on login node.

mkdir -p logs

CONDA_BASE="/hpcfs/fpublic/app/miniforge3/conda"
TARGET_ENV="${HPC_ENV:-LineaBact}"
PARTITION="${HPC_PARTITION:-qfree}"
NODELIST="${HPC_NODE:-vnode1}"
CPUS="${HPC_CPUS:-8}"
MEM="${HPC_MEM:-12G}"

JOB_NAME="${1:-fastBact_bench}"
shift || true
COMMAND_TO_RUN="$*"

JOB_LOG="logs/bench_$(date +%Y%m%d_%H%M%S)_$$.log"

set +e
sbatch --wait \
  --job-name="${JOB_NAME}" \
  --partition="${PARTITION}" \
  --nodelist="${NODELIST}" \
  --nodes=1 \
  --ntasks=1 \
  --cpus-per-task="${CPUS}" \
  --mem="${MEM}" \
  --output="${JOB_LOG}" \
  --error="${JOB_LOG}" \
  --wrap="
    source \"${CONDA_BASE}/etc/profile.d/conda.sh\"
    conda activate \"${TARGET_ENV}\"
    export PATH=\"\$HOME/.cargo/bin:\$HOME/.local/bin:\$PATH\"
    export CARGO_BUILD_JOBS=${CPUS}
    export RAYON_NUM_THREADS=${CPUS}
    export OMP_NUM_THREADS=${CPUS}
    echo \"[HPC Node: \$(hostname)] [Job ID: \$SLURM_JOB_ID] Running: ${COMMAND_TO_RUN}\"
    ${COMMAND_TO_RUN}
  "
EXIT_CODE=$?
set -e

if [ -f "${JOB_LOG}" ]; then
  cat "${JOB_LOG}"
fi

exit ${EXIT_CODE}
