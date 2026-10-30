#!/usr/bin/env bash
#SBATCH --job-name=f25-release-val
#SBATCH --partition=qfree
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=14G
#SBATCH --time=00:30:00
#SBATCH --output=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f25-correctness-worktree/slurm-f25-release-%j.out
#SBATCH --error=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f25-correctness-worktree/slurm-f25-release-%j.err

set -Eeuo pipefail
umask 027
export LC_ALL=C TZ=UTC RUST_LOG=info NO_COLOR=1
export RAYON_NUM_THREADS=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1

cd /hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f25-correctness-worktree

echo "=== F25 RELEASE TEST SUITE & CLI VERIFICATION START ==="
date -u +"%Y-%m-%dT%H:%M:%SZ"
echo "Host: $(hostname), Slurm Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Git HEAD: $(git rev-parse HEAD)"
echo "Source SHA-256 (exspander.rs): $(sha256sum src/resolve/exspander.rs | awk '{print $1}')"
echo "Source SHA-256 (consensus.rs): $(sha256sum src/polish/consensus.rs | awk '{print $1}')"

# 1. Build release binary
echo "--- Step 1: Building release binary ---"
cargo build --release --bin fastBact

BIN_SHA=$(sha256sum target/release/fastBact | awk '{print $1}')
echo "Release Binary SHA-256: $BIN_SHA"

# 2. Run full workspace release test suite
echo "--- Step 2: Running full cargo test --release ---"
cargo test --release -- --nocapture

# 3. CLI smoke test
echo "--- Step 3: Running CLI Smoke Test ---"
target/release/fastBact --help > /dev/null
echo "CLI smoke test --help: OK"
target/release/fastBact --version
echo "CLI smoke test --version: OK"

echo "=== F25 RELEASE TEST SUITE & CLI VERIFICATION COMPLETE ==="
date -u +"%Y-%m-%dT%H:%M:%SZ"
