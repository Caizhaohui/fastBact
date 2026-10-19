#!/bin/bash
#SBATCH --partition=qfree
#SBATCH --nodelist=vnode1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=00:15:00
#SBATCH --job-name=c002r-d-verify
set -eu
root=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact
export PATH="$HOME/.cargo/bin:$HOME/.local/bin:$PATH"
cd "$root"
printf 'job=%s node=%s cpus=%s mem=%s\n' "$SLURM_JOB_ID" "$SLURMD_NODENAME" "$SLURM_CPUS_PER_TASK" "$SLURM_MEM_PER_NODE"
sha256sum benchmarks/scripts/c0_02r/verification.py benchmarks/scripts/verify_c0_02r_manifest.py tests/test_c0_02r_manifest.py
/hpcfs/fpublic/app/miniforge3/conda/bin/python3 -m pytest tests/test_c0_02r_manifest.py -v
/hpcfs/fpublic/app/miniforge3/conda/bin/python3 benchmarks/scripts/verify_c0_02r_manifest.py benchmark-results/C0_02R_takeover_20261019/results/C0_02R_RUN_MANIFEST.json
