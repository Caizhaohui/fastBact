#!/bin/bash
#SBATCH --partition=qfree
#SBATCH --nodelist=vnode1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=00:20:00
#SBATCH --job-name=c002r-a-red
set -u
root=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact
evidence=/hpcfs/fhome/caizhh/.local/state/codex-with-chatgpt/fastbact-preflight-69HE09JF
export PATH="$HOME/.cargo/bin:$PATH"
export CARGO_BUILD_JOBS=8
printf 'job=%s node=%s cpus=%s memory=%s\n' "$SLURM_JOB_ID" "$SLURMD_NODENAME" "$SLURM_CPUS_PER_TASK" "$SLURM_MEM_PER_NODE"
rustc --version
cargo --version
for variant in baseline inherited; do
    source_dir="$evidence/$variant"
    cp "$root/tests/test_c0_02r_directed_edge_identity.rs" "$source_dir/tests/"
    cp "$root/tests/test_c0_02r_edge_gate.rs" "$source_dir/tests/"
    cd "$source_dir" || exit 2
    git rev-parse HEAD HEAD^{tree}
    git status --porcelain=v1
    git diff --binary > "$evidence/$variant-production.patch"
    sha256sum src/resolve/exspander.rs tests/test_c0_02r*.rs
    export CARGO_TARGET_DIR="$evidence/target-$variant"
    printf 'COMMAND: cargo test --locked --offline --test test_c0_02r_directed_edge_identity --test test_c0_02r_edge_gate -- --nocapture\n'
    cargo test --locked --offline --test test_c0_02r_directed_edge_identity --test test_c0_02r_edge_gate -- --nocapture > "$evidence/$variant-a-red.log" 2>&1
    code=$?
    printf '%s\t%s\n' "$variant" "$code" | tee -a "$evidence/phase-a-red-exits.tsv"
    tail -65 "$evidence/$variant-a-red.log"
done
