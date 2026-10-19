#!/bin/bash
#SBATCH --partition=qfree
#SBATCH --nodelist=vnode1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=00:10:00
#SBATCH --job-name=c002r-a-base
set -u
root=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact
evidence=/hpcfs/fhome/caizhh/.local/state/codex-with-chatgpt/fastbact-preflight-69HE09JF
export PATH="$HOME/.cargo/bin:$PATH"
export CARGO_BUILD_JOBS=8 CARGO_TARGET_DIR="$evidence/target-baseline"
cp "$root/tests/test_c0_02r_directed_edge_identity.rs" "$evidence/baseline/tests/"
cd "$evidence/baseline" || exit 2
printf 'job=%s node=%s cpus=%s memory=%s\n' "$SLURM_JOB_ID" "$SLURMD_NODENAME" "$SLURM_CPUS_PER_TASK" "$SLURM_MEM_PER_NODE"
git rev-parse HEAD HEAD^{tree}
git diff --exit-code -- src Cargo.toml Cargo.lock || exit 3
sha256sum tests/test_c0_02r*.rs
for suite in test_c0_02r_directed_edge_identity test_c0_02r_edge_gate; do
    printf 'COMMAND: cargo test --locked --offline --test %s -- --nocapture\n' "$suite"
    cargo test --locked --offline --test "$suite" -- --nocapture > "$evidence/baseline-$suite-red.log" 2>&1
    code=$?
    printf '%s\t%s\n' "$suite" "$code"
    tail -70 "$evidence/baseline-$suite-red.log"
done
