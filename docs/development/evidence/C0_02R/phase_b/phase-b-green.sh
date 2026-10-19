#!/bin/bash
#SBATCH --partition=qfree
#SBATCH --nodelist=vnode1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=00:25:00
#SBATCH --job-name=c002r-b-green
set -u
export PATH="$HOME/.cargo/bin:$PATH"
evidence=/hpcfs/fhome/caizhh/.local/state/codex-with-chatgpt/fastbact-preflight-69HE09JF
cd "$evidence/phase-a" || exit 2
export CARGO_BUILD_JOBS=8 CARGO_TARGET_DIR="$evidence/target-phase-a" RAYON_NUM_THREADS=8
printf 'job=%s node=%s cpus=%s mem=%s\n' "$SLURM_JOB_ID" "$SLURMD_NODENAME" "$SLURM_CPUS_PER_TASK" "$SLURM_MEM_PER_NODE"
rustfmt --edition 2021 tests/test_c0_02r_cleanup_gate.rs
sha256sum src/resolve/exspander.rs tests/test_c0_02r_cleanup_gate.rs
printf 'COMMAND: cargo test --locked --offline --test test_c0_02r_graph_atomicity --test test_c0_02r_cleanup_gate --test test_c0_02r_directed_edge_identity --test test_c0_02r_edge_gate --test test_c0_02_junction_overlap --test test_c0_01_polish_orientation -- --nocapture\n'
cargo test --locked --offline --test test_c0_02r_graph_atomicity --test test_c0_02r_cleanup_gate --test test_c0_02r_directed_edge_identity --test test_c0_02r_edge_gate --test test_c0_02_junction_overlap --test test_c0_01_polish_orientation -- --nocapture > "$evidence/phase-b-green-tests.log" 2>&1
code=$?
printf 'targeted_exit=%s\n' "$code"
tail -110 "$evidence/phase-b-green-tests.log"
if test "$code" -ne 0; then exit "$code"; fi
printf 'COMMAND: cargo test --locked --offline --all-targets\n'
cargo test --locked --offline --all-targets > "$evidence/phase-b-regression.log" 2>&1
code=$?
printf 'regression_exit=%s\n' "$code"
tail -40 "$evidence/phase-b-regression.log"
exit "$code"
