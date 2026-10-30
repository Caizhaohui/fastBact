#!/usr/bin/env bash
#SBATCH --job-name=f24-full-release
#SBATCH --partition=qfree
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=14G
#SBATCH --time=00:30:00
#SBATCH --output=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f24-contiguity-worktree/benchmark-results/f24_validation_20261029/logs/slurm-f24-full-release-%j.out
#SBATCH --error=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f24-contiguity-worktree/benchmark-results/f24_validation_20261029/logs/slurm-f24-full-release-%j.err

set -Eeuo pipefail
umask 027
export LC_ALL=C TZ=UTC RUST_LOG=info NO_COLOR=1
export RAYON_NUM_THREADS=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1

cd /hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f24-contiguity-worktree

echo "=== F24 FULL RELEASE TEST SUITE & CLI SMOKE START ==="
date -u +"%Y-%m-%dT%H:%M:%SZ"
echo "Host: $(hostname), Slurm Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Git HEAD: $(git rev-parse HEAD)"
echo "exspander.rs SHA-256: $(sha256sum src/resolve/exspander.rs | awk '{print $1}')"
echo "test_c0_02r_f24_contiguity_gate.rs SHA-256: $(sha256sum tests/test_c0_02r_f24_contiguity_gate.rs | awk '{print $1}')"

# 1. Build release binary
echo "--- Step 1: Building release binary ---"
cargo build --release

BIN_SHA=$(sha256sum target/release/fastBact | awk '{print $1}')
echo "Release Binary SHA-256: $BIN_SHA"

# 2. Run full workspace release test suite
echo "--- Step 2: Running cargo test --release ---"
cargo test --release

# 3. CLI smoke test
echo "--- Step 3: Running CLI Smoke Test ---"
target/release/fastBact --help > /dev/null
echo "CLI smoke test --help: OK"
target/release/fastBact --version
echo "CLI smoke test --version: OK"

# 4. Run targeted F24 contiguity gate test in release
echo "--- Step 4: Running F24 contiguity gate test ---"
cargo test --release --test test_c0_02r_f24_contiguity_gate -- --nocapture

# 5. Run end-to-end integration benchmark in release
echo "--- Step 5: Running synthetic bacterial assembly integration test ---"
cargo test --release --test integration_test -- --nocapture

echo "=== F24 FULL RELEASE TEST SUITE & CLI SMOKE COMPLETE ==="
date -u +"%Y-%m-%dT%H:%M:%SZ"
