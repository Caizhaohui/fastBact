#!/bin/bash
#SBATCH -p qfree
#SBATCH -w vnode1
#SBATCH -c 8
#SBATCH --mem=12G
#SBATCH -o /hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f23-unique-anchor-worktree/.test_f23_horizon_repair.log
#SBATCH -e /hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f23-unique-anchor-worktree/.test_f23_horizon_repair.log

cd /hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f23-unique-anchor-worktree
echo "=== F23 horizon-guard repair verification ==="
git status --short
echo "--- [1/3] four regressed semantics tests ---"
cargo test --release --test test_c0_02r_state_semantics_gate -- test_t26 test_t27 test_t30 test_t34 2>&1 | tail -4
echo "--- [2/3] new horizon gate test ---"
cargo test --release --test test_c0_02r_unique_anchor_gate 2>&1 | tail -4
echo "--- [3/3] FULL release suite ---"
cargo test --release 2>&1 | grep -E "^(test result|running|error)|FAILED|panicked" 
echo "--- done ---"
