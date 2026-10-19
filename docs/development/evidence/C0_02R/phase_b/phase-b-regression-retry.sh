#!/bin/bash
#SBATCH --partition=qfree
#SBATCH --nodelist=vnode1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --job-name=c002r-b-regress
set -u
export PATH="$HOME/.cargo/bin:$PATH"
evidence=/hpcfs/fhome/caizhh/.local/state/codex-with-chatgpt/fastbact-preflight-69HE09JF
export CARGO_BUILD_JOBS=4 RAYON_NUM_THREADS=8
printf 'job=%s node=%s cpus=%s mem=%s\n' "$SLURM_JOB_ID" "$SLURMD_NODENAME" "$SLURM_CPUS_PER_TASK" "$SLURM_MEM_PER_NODE"
cd "$evidence/baseline" || exit 2
export CARGO_TARGET_DIR="$evidence/target-baseline"
printf 'COMMAND (original b9727c3): cargo test --locked --offline --test integration_test -- --nocapture\n'
cargo test --locked --offline --test integration_test -- --nocapture > "$evidence/baseline-debug-integration.log" 2>&1
printf 'baseline_debug_exit=%s\n' "$?"
tail -30 "$evidence/baseline-debug-integration.log"
cd "$evidence/phase-a" || exit 2
export CARGO_TARGET_DIR="$evidence/target-phase-a"
printf 'COMMAND (fixed): cargo test --release --locked --offline --all-targets\n'
cargo test --release --locked --offline --all-targets > "$evidence/phase-b-release-regression.log" 2>&1
code=$?
printf 'fixed_release_exit=%s\n' "$code"
tail -65 "$evidence/phase-b-release-regression.log"
exit "$code"
