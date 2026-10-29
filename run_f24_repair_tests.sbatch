#!/usr/bin/env bash
#SBATCH --job-name=f24-repair-test
#SBATCH --partition=qfree
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=00:30:00
#SBATCH --output=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f24-contiguity-worktree/benchmark-results/f24_validation_20261029/logs/slurm-repair-%j.out
#SBATCH --error=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f24-contiguity-worktree/benchmark-results/f24_validation_20261029/logs/slurm-repair-%j.err

set -Eeuo pipefail
umask 027
export LC_ALL=C TZ=UTC RUST_LOG=info NO_COLOR=1
export RAYON_NUM_THREADS=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1

cd /hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f24-contiguity-worktree

echo "=== F24 FULL RELEASE TEST SUITE START ==="
date -u +"%Y-%m-%dT%H:%M:%SZ"
echo "Host: $(hostname), Slurm Job ID: ${SLURM_JOB_ID:-unknown}"

echo "--- Running full release test suite: cargo test --release ---"
cargo test --release -- --nocapture

echo "--- Building release binary: cargo build --release --bin fastBact ---"
cargo build --release --bin fastBact

echo "--- Release Binary SHA-256 ---"
sha256sum target/release/fastBact

echo "=== F24 FULL RELEASE TEST SUITE COMPLETE ==="
date -u +"%Y-%m-%dT%H:%M:%SZ"
