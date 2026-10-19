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

## Phase B investigation

Phase A commit: `f500d0f`. Red job 2978134: existing G1–G6 already pass on
original b9727c3 and after A (6/6 each); these are not new failure evidence.
Supplemental six cases fail 5/6 on original baseline, 3/6 after A. After A,
connected-short cleanup, incoming-only connectivity, and legal short repeat
reuse fail. The existing uncommitted short-repeat edit repairs the last case
and is adopted as inherited work after independent reproduction.

Cleanup policy: delete only unrelated isolated, noncircular nodes shorter than
500 bp that are not in rejected_nodes. Connectivity uses incident edges
between active nodes on either strand, including incoming-only edges.
Invalid/deleted targets or deleted sources do not establish active connectivity.
Candidate failure preserves graph state apart from the separately defined
cleanup of unrelated isolated nodes; a whole-resolver equality claim without
this qualification would be wrong.

Supplemental snapshots include entire CdbgGraph Debug storage (all node fields
and arena data/offsets/lengths), with explicit DNA/commit assertions. Test inputs
are deterministic; red-to-green edits are formatting only. Successful repeat
fixture expects exactly 2026 bp; retry fixture expects 1288 bp.

Job 2978136: targeted 64/64 pass; library 32/32 pass; overall job FAILED 101:0
because debug integration exceeds the unmodified 3.5s timing assertion (13.85s).
Original b9727c3 independently reproduces the same timing failure (12.58s)
in job 2978137. This timing failure is pre-existing and is not a correctness
PASS. Release all-target regression is running in 2978137; Phase C/D stay gated
until that result is assessed. No timing assertion/test is weakened or removed.

## Phase B receipt

Release regression job 2978137 COMPLETED 0:0 on qfree/vnode1, 8 CPUs/12G.
`cargo test --release --locked --offline --all-targets` passed 97 tests:
32 library + 1 full synthetic assembly + 64 targeted integration cases.
Synthetic assembly test finished in 1.58s, preserving its original timing
assertion. Original baseline debug timing failure is independently established;
debug all-target regression is still reported as failed, not relabeled PASS.
Phase B source is Phase A plus `phase_b/implementation.patch`; hashes match
targeted job 2978136. Full synthetic assembly exercised pipeline through output.
No test was removed, ignored, or weakened.

Five review lanes independently examined Phase B; static goal, code, security,
and local-context lanes PASS. Runtime lane reports the targeted gate PASS and
initial debug full-regression FAIL; the separate release receipt now supplies
the complete full-regression result. All raw outcomes are retained.
Phase B ready for independent commit; Phase C is next, Phase D follows C.

## Phase C receipt

Phase B commit: `0d13f99`. Legacy `run_c0_02_benchmark.py` now delegates to a
new modular alignment/provenance implementation, retaining compare_fasta API.
The preserved untracked `run_c0_02r_benchmark.py` is NOT the accepted runner.

Original comparator failure: job2978143 reports4540 differences for1 inserted
base, and fails reverse-complement recognition. Dotted-label receipt collision
reproduced2978159; overlapping split/merge and inversion tp:A:I holes reproduced
2978175; all repaired. Latest2978177:26/26 Python tests pass onqfree/vnode1,
8CPUs/12G. Ruff default checks and shell syntax checks pass. Static source/goal,
security, context and runtime reviewers passed their scoped checks after fixes.
New modules remain below250 pureLOC (alignment216, orchestration246; split before
future growth). No unrelated Rust refactor was introduced into the large existing
resolver module. Data boundaries have typed records; operation handling is
exhaustive, no type suppressions or broad exception swallowing were added.

Formal runner requires clean A/B/C sources, verifies original A/B SHA, builds
fresh binaries, rejects command/QUAST/parser failures, hashes source/binary/data/
FASTA/report, records exact commands/exits/versions/Slurm/resources/RSS/walltime,
and rechecks source/input identities. All18 assembly/QUAST pairs rerun; historical
FASTA/QUAST metrics are not reused. Benchmark acceptance criteria are explicitly
conservative versus bothA/B: no structural/error metric increase or genome
fraction decrease; N50 alone never fails the gate.

Builtin --dump-stages on all original versions supplies graph node/edge/base
inventories; Stage05 includes post-exSPAnder simplification. Counts are measured
stage inventories, not proof of path atomicity (the graph tests provide that).
Coordinates are0-based half-open; cs alleles follow alignment orientation.
Reference-coordinate projection does NOT establish reference-confirmed correction.
Overlapping structural relationships retain explicit ambiguity.

Formal invocation after this commit:
`bash benchmarks/scripts/submit_c0_02r_final_gate.sh "$PWD" "$PWD/benchmark-results/C0_02R_takeover_20261019"`.
This creates separate clean detached sources and submits all compute to Slurm.
Results include all requested filenames plus regenerated C0_02_JUNCTION_DELTAS.tsv
inside the NEW results directory. Starting dirty frozen/root artifacts remain
protected, outside the accepted output set. All datasets currently qualify only
as REPRODUCIBLE_ON_ORIGINAL_HPC; self-contained external packaging is not claimed.
Phase D next: formal execution, inspect acceptance, archive fresh receipts,
independent review, final documentation commit; push only after gate completion.
