#!/bin/bash
#SBATCH --partition=qfree
#SBATCH --nodelist=vnode1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=00:10:00
#SBATCH --job-name=c002r-b-red
set -u
export PATH="$HOME/.cargo/bin:$PATH"
evidence=/hpcfs/fhome/caizhh/.local/state/codex-with-chatgpt/fastbact-preflight-69HE09JF
root=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact
export CARGO_BUILD_JOBS=8
printf 'job=%s node=%s cpus=%s mem=%s\n' "$SLURM_JOB_ID" "$SLURMD_NODENAME" "$SLURM_CPUS_PER_TASK" "$SLURM_MEM_PER_NODE"
for variant in baseline phase-a; do
    cd "$evidence/$variant" || exit 2
    cp "$root/tests/test_c0_02r_cleanup_gate.rs" "$root/tests/test_c0_02r_graph_atomicity.rs" tests/
    git rev-parse HEAD HEAD^{tree}
    sha256sum src/resolve/exspander.rs tests/test_c0_02r_cleanup_gate.rs tests/test_c0_02r_graph_atomicity.rs
    export CARGO_TARGET_DIR="$evidence/target-$variant"
    for suite in test_c0_02r_graph_atomicity test_c0_02r_cleanup_gate; do
        printf 'COMMAND: cargo test --locked --offline --test %s -- --nocapture\n' "$suite"
        cargo test --locked --offline --test "$suite" -- --nocapture > "$evidence/$variant-$suite-red.log" 2>&1
        printf 'variant=%s suite=%s cargo_exit=%s\n' "$variant" "$suite" "$?"
        tail -65 "$evidence/$variant-$suite-red.log"
    done
done
