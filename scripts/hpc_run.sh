#!/usr/bin/env bash
set -uo pipefail
mkdir -p logs

# 1. Configure discovered Conda base & existing environment name
CONDA_BASE="/hpcfs/fpublic/app/miniforge3/conda"
TARGET_ENV="${HPC_ENV:-rust-env}"
PARTITION="${HPC_PARTITION:-qfree}"
CPUS="${HPC_CPUS:-8}"
MEM="${HPC_MEM:-10G}"

JOB_LOG="logs/job_$(date +%Y%m%d_%H%M%S)_$$.log"

COMMAND_TO_RUN="$*"

# 2. Submit synchronously to compute node via sbatch --wait (0 CPU used on login node)
set +e
sbatch --wait \
  --job-name="fastBact_ci" \
  --partition="${PARTITION}" \
  --nodes=1 \
  --ntasks=1 \
  --cpus-per-task="${CPUS}" \
  --mem="${MEM}" \
  --output="${JOB_LOG}" \
  --error="${JOB_LOG}" \
  --wrap="
    source \"${CONDA_BASE}/etc/profile.d/conda.sh\"
    conda activate \"${TARGET_ENV}\"
    export PATH=\"\$HOME/.cargo/bin:\$PATH\"
    export CARGO_BUILD_JOBS=${CPUS}
    export RAYON_NUM_THREADS=${CPUS}
    echo \"[HPC Node: \$(hostname)] Running: ${COMMAND_TO_RUN}\"
    ${COMMAND_TO_RUN}
  "
EXIT_CODE=$?
set -e

if [ -f "${JOB_LOG}" ]; then
  cat "${JOB_LOG}"
fi
exit ${EXIT_CODE}
