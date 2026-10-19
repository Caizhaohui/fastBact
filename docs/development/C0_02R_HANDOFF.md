# C0-02R takeover

## Preflight (2026-10-19)

Actual starting HEAD: `cbf5c1a328ee71fad865da51165a7dceb353e106`, branch
`fix/c0-02r-final-gate`. Ancestors include Phase A claim `02f8236` and Phase B
claim `cbf5c1a`; these are inherited work, not independently accepted results.
`master`, `origin/main`, and `origin/master` point to `dac6a57`.
Remote URL: `git@github.com:Caizhaohui/fastBact.git`.

Starting dirty tracked paths: `src/resolve/exspander.rs` and frozen
`C0_02_JUNCTION_DELTAS.tsv`. Starting untracked paths: benchmark runner
`benchmarks/scripts/run_c0_02r_benchmark.py`, alignment events, assembly table,
manifest, and final verdict. Preserve all; the inherited final PASS is unverified.
Backup: `/hpcfs/fhome/caizhh/.local/state/codex-with-chatgpt/fastbact-preflight-69HE09JF`.
`tracked.patch` SHA-256: `bf6d89ed37b30ed47ef3c4b504bdfdeb8d5c2a0be48f1f70caf744d064c5e174`.
No reset, branch recreation, history rewrite, or frozen-artifact overwrite.

## Execution ledger

1. Phase A (in progress): run E1–E8 against b9727c3 production sources with
   test-only injection; preserve commands, source identity, exits, and logs.
   Add direct validator ambiguity/losing-edge and reciprocal-support tests,
   prove red on inherited HEAD, minimally fix, then run overlap/polish regressions.
2. Phase B (gated on A): independently reproduce G1–G6 and check complete
   snapshots, cleanup scope, retry safety, successful node conservation.
3. Phase C (gated on B): correct alignment event extraction and fail-closed
   provenance; write fresh outputs outside frozen C0-01/C0-02 artifacts.
4. Phase D (gated on C): fresh A/B/C runs across all six requested datasets,
   complete regression, review, and evidence-derived acceptance.

## Confirmed source defects / tests

At b9727c3, `validate_and_splice_path()` ignores incoming edges entirely and
backward reversal constructs `GraphEdge::new()` without graph membership checks.
Inherited HEAD checks membership but its validator bypasses edge ambiguity and
winner selection; reciprocal lookups ignore selected-edge support. New tests
in `tests/test_c0_02r_edge_gate.rs` distinguish these cases. Runtime confirmation
is pending; no correctness PASS is claimed.

## Resources and blockers

Slurm is reachable; vnode1 reports IDLE+PLANNED, 16 CPUs, 15500 MB RAM.
Use qfree/vnode1, 8 CPUs, 12 GB. All compilation/tests/assemblies/QUAST/large
data processing must run in that allocation. No production code was modified
before establishing failing tests in this takeover.

## Phase A receipt

Completed and independently reviewed: goal, code quality, security, runtime QA,
local history/caller context all returned scoped PASS. See `evidence/C0_02R/phase_a/`.
No result here accepts the inherited overall final verdict.

Slurm red job 2978126: original baseline attempt failed compilation because the
inherited test used a post-baseline `find_edge()` helper; this is NOT runtime
bug evidence. Baseline retry 2978131 uses iteration only in the test and
sequence-compatible E3/E4 fixtures: E suite 2 pass/6 fail, edge gate 1 pass/4 fail,
Cargo exits 101. Production sources/Cargo files match b9727c3 exactly.
No production diagnostic patch was needed (production patch empty).
Inherited cbf5c1a: E suite 8/8 pass, new edge gate 2 pass/3 fail, Cargo exit 101.

Green jobs 2978130 and 2978132: qfree/vnode1, 8 CPUs/12G; both completed 0:0.
Final targeted command:
`cargo test --locked --offline --test test_c0_02r_edge_gate --test test_c0_02r_directed_edge_identity --test test_c0_02_junction_overlap --test test_c0_01_polish_orientation -- --nocapture`.
52 passed: 5 edge gate + 8 E cases + 10 overlap + 29 polishing.
`cargo check --locked --offline --all-targets` exited 0.
Tests drove real resolver/validator APIs and exact sequence/graph assertions.
The separate atomicity suite and full assembler run belong to later gates.

Tested production source is cbf5c1a plus `phase_a/a-implementation.patch`;
original uncommitted short-repeat reuse edit is excluded from this phase.
Initial red test source and hashes are archived; green rustfmt changed only
new-test formatting after initial hashing. Final checksums identify archived
files; commit content identifies final tests. Baseline latest log records the
strengthened E3/E4 test hash. Raw script captures Cargo failures; the negative
Slurm wrapper exits 0 intentionally, so use recorded Cargo exit 101.

LSP hook cannot run: stable toolchain lacks rust-analyzer. This is an explicit
verification limitation; successful compute-node cargo check provides compiler
diagnostics, not an LSP result. No toolchain install or login-node build was run.

Fresh `git ls-remote --heads origin` confirms remote main/master at dac6a57,
C0-01 at 3950b2a, C0-02 at b9727c3, no remote C0-02R branch yet.
Next active phase after this commit: Phase B. Phase C/D remain gated.
