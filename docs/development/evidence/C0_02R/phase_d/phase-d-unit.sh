#!/bin/bash
#SBATCH --partition=qfree
#SBATCH --nodelist=vnode1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --job-name=c002r-d-unit
set -eu
cd /hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact
printf 'job=%s node=%s cpus=%s mem=%s\n' "$SLURM_JOB_ID" "$SLURMD_NODENAME" "$SLURM_CPUS_PER_TASK" "$SLURM_MEM_PER_NODE"
sha256sum benchmarks/scripts/c0_02r/verification.py tests/test_c0_02r_manifest.py
/hpcfs/fpublic/app/miniforge3/conda/bin/python3 -m pytest tests/test_c0_02r_manifest.py -v
