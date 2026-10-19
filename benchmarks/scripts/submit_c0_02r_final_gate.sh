#!/bin/bash
set -eu
if test "$#" -ne 2; then
    echo 'Usage: bash submit_c0_02r_final_gate.sh REPOSITORY NEW_RUN_DIRECTORY' >&2
    exit 2
fi
repo=$(realpath "$1")
run_dir=$(realpath -m "$2")
mkdir "$run_dir"
mkdir "$run_dir/sources"
git -C "$repo" worktree add --detach "$run_dir/sources/A" 3950b2ab2584b24d39172b25bbbcb565232891b7
git -C "$repo" worktree add --detach "$run_dir/sources/B" b9727c3ab372f7b995e439d92b72c7a4bfd051cf
git -C "$repo" worktree add --detach "$run_dir/sources/C" HEAD
sbatch --parsable --output="$run_dir/slurm-%j.txt" "$run_dir/sources/C/benchmarks/scripts/execute_c0_02r_benchmark.sh" "$run_dir/sources" "$repo" "$run_dir/results"
