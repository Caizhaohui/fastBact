#!/bin/bash
#SBATCH -p qfree
#SBATCH -w vnode1
#SBATCH -c 8
#SBATCH --mem=12G
#SBATCH -o /hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f23-unique-anchor-worktree/.test_f23_unique_anchor_red.log
#SBATCH -e /hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f23-unique-anchor-worktree/.test_f23_unique_anchor_red.log

cd /hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f23-unique-anchor-worktree
echo "=== F23 unique-anchor gate RED run (unpatched production src, expect 2 FAIL) ==="
git rev-parse HEAD
git status --short
cargo test --test test_c0_02r_unique_anchor_gate 2>&1
echo "=== exit: $? ==="
