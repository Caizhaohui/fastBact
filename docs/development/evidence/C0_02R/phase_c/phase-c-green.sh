#!/bin/bash
#SBATCH --partition=qfree
#SBATCH --nodelist=vnode1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=00:10:00
#SBATCH --job-name=c002r-c-tests
set -u
cd /hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact || exit 2
printf 'job=%s node=%s cpus=%s mem=%s\n' "$SLURM_JOB_ID" "$SLURMD_NODENAME" "$SLURM_CPUS_PER_TASK" "$SLURM_MEM_PER_NODE"
sha256sum benchmarks/scripts/c0_02r/*.py benchmarks/scripts/run_c0_02_benchmark.py tests/test_c0_02r_alignment.py tests/test_c0_02r_legacy_delta.py
/hpcfs/fhome/caizhh/.local/lib/python3.9/site-packages/quast_libs/minimap2/minimap2 --version
/hpcfs/fpublic/app/miniforge3/conda/bin/python3 -m pytest tests/test_c0_02r_alignment.py tests/test_c0_02r_legacy_delta.py -v
