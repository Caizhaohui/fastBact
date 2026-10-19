#!/bin/bash
#SBATCH --partition=qfree
#SBATCH --nodelist=vnode1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=03:00:00
#SBATCH --job-name=c002r-final-ABC
set -eu
export PATH="$HOME/.cargo/bin:$HOME/.local/bin:$PATH"
python_bin=${C0_02R_PYTHON:-/hpcfs/fpublic/app/miniforge3/conda/bin/python3}
exec "$python_bin" "$1/C/benchmarks/scripts/run_c0_02r_final_gate.py" "$1" "$2" "$3"
