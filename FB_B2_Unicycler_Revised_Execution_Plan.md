# fastBact — Revised B2 + Unicycler Integration Execution Plan

**Document ID:** FB-BENCH-B2-U-REV1  
**Status:** Proposed; requires repository-local preflight verification  
**Date:** 2026-10-15  
**Scope:** Illumina paired-end bacterial isolate de novo assembly; no long-read/hybrid development  
**Algorithm implementer:** Gemini 4  
**Review/approval:** Benchmark quality gates; no automatic multi-stage continuation

## 1. Source-of-truth and baseline identities

1. Verify the actual object referenced by tag `v0.1.0-baseline`; according to the B1 report, frozen application baseline was commit `4259fe3`.
2. Verify commit `552ff54` is the B1 telemetry/report commit, and document its diff relative to `4259fe3`. Do not silently claim these are the same baseline.
3. Never move `v0.1.0-baseline`, rewrite historical results, or overwrite prior FASTA, GFA, metrics, logs, or raw FASTQ files.
4. Record `git rev-parse HEAD`, `git status --porcelain`, tool versions, command line, dataset SHA-256, reference SHA-256, actual hostname/CPU/RAM and Slurm Job ID for every run.
5. Every production algorithm change must have one bounded causal hypothesis, its own commit, reproducible results, and a formal PASS/FAIL/INCONCLUSIVE verdict. Diagnostics-only code must be separable from algorithm changes.
6. Use the existing Slurm wrapper and compute nodes. Never compile, assemble, or run bioinformatics benchmarks on the login node. For fair CPU timings pin the benchmark to the same hardware or disclose the difference. Do not assume `vnode1` is always allocated merely from selecting `qfree`.

## 2. Upstream Unicycler reference policy

Reference upstream repository: https://github.com/rrwick/Unicycler ; inspected revision `1429301eb24df9162f3c9386e8522f942594ae42` (2026-07-29). Pin a full SHA in reports, never an unqualified moving `main` or `--depth 1` without documenting the resulting commit.

| Topic | Upstream source |
|---|---|
| Graph, topology, overlaps and cleaning | `unicycler/assembly_graph.py` |
| Graph multiplicity / copy-depth | `unicycler/assembly_graph_copy_depth.py` |
| k selection on SPAdes graph | `unicycler/spades_func.py` |
| SPAdes contig / short-read bridging | `unicycler/bridge_spades_contig.py` |
| Loop unrolling | `unicycler/bridge_loop_unroll.py` |
| Supporting bridge scoring | `unicycler/bridge_common.py` |
| Regression examples | `test/test_assembly_graph_copy_depth.py` and graph fixture files |

There is no `unicycler/bridge.py` in that checked upstream tree. `determine_copy_depth()` is not defined in `assembly_graph.py`.

Unicycler is GPL-3.0. Study concepts and published algorithms; independently design and test Rust code. Do not copy source code, line-by-line translate formulas/implementation blocks verbatim, or include GPL-protected code without a deliberate license compatibility review. Document citations and design provenance. Unicycler short-read mode is SPAdes-based: it is a **pipeline benchmark**, not a standalone graph-build comparator. Do not add SPAdes/Unicycler executables as fastBact production dependencies.

## 3. Corrected milestone sequence

**B2-G0 → B2-C1a → B2-C1b → U0 → U1 → U2a → U2b → U3a → U3b → U4a → U4b → U5a → Release Gate.**

Only B2-G0/C1a/C1b are presently authorized for implementation. U0 can be conducted as read-only research in parallel. Every subsequent stage requires separate written PASS and approval before starting implementation.

### B2-G0 — Evidence reconciliation, no production behavior change

- Inspect source and B1 records that led to HLL size 522 kb, depth 626.4x, estimated keep rate 0.1916; export intermediate counters and units.
- Verify gzip file sample-byte cursor vs decompressor buffering, both R1 and R2 gzip sizes, sampling order and variation in compression ratio across the file.
- Examine the exact multiplicity histogram of sample k=21 from a trusted exact-counter reference; test the claimed Poisson explanation before changing formulas.
- Verify which of the four new misassemblies in `R1-001 --no-subsample` are reference disagreements vs independently supported false joins.
- Validate the raw B1 metrics and the tag `4259fe3` vs telemetry `552ff54` relationship.

Deliver: `B2_ESTIMATOR_FORENSICS.md`, `B2_LINEAGE.tsv`, `B2_MISASSEMBLY_REVIEW.tsv`. Gate: evidence traceable; unknowns explicit.

### B2-C1a — Safe retention/failure fallback (one variable)

- Keep existing estimator unchanged initially; when its depth or genome-size estimate is unreliable or contradictory, **disable automatic downsampling** rather than delete reads.
- Expose `--no-subsample` / explicit `--genome-size` with clear diagnostics, preserving compatibility.
- Unknown estimate must never silently imply reliable high coverage. No forced subsampling from unreliable estimates.

A/B: frozen baseline default vs frozen baseline `--no-subsample` vs C1a default, same FASTQ/seeds/threads/evaluator. On known true <=120x data, C1a must not erroneously remove reads. If C1a retains all reads, output should be identical or nearly identical to frozen `--no-subsample`; investigate every difference.

Deliver `B2_C1A_RETAINED_READS.tsv` and `B2_C1A_ASSEMBLY.tsv`.

### B2-C1b — Accurate counts and calibrated estimation (one logical feature)

- Eliminate the fixed compressed-bytes-per-read projection. First implement an exact first-pass record/bases count for both FASTQ mates; second pass performs sampling and assembly. Zero temporary decompressed FASTQ files required.
- Keep **read totals**, **genome size**, **estimated effective depth**, and **actual retained depth** separate; precise pair count is not a precise genome-size estimate.
- Evaluate HLL and k-mer frequency models separately. A zero/one-truncated Poisson correction is an unproven experimental estimator, not an acceptance assumption; report uncertainty and test GC/coverage/repeats.
- If genome size cannot be confidently estimated, leave subsampling disabled, record reason, and permit `--genome-size` override.
- Use deterministic **pair-level** sampling with a fixed seed, never independent R1/R2 selection.
- Retain useful provenance counters and warn if effective coverage after quality filters differs materially from target.

Test 15x/30x/60x/100x/150x/300x and distinct genome sizes/GC levels; require coverage kept within predeclared tolerance when the estimate is confident. Re-evaluate 100x synthetic and confirmed pure-isolate real data; protect plasmid underrepresentation and low-depth tails.

Deliver `B2_COVERAGE_ACCURACY.tsv`, `B2_SAMPLING.tsv`, `B2_C1B_VERDICT.md`.

### U0 — Read-only Unicycler design-gap study

- Pin upstream SHA, research graph scoring, graph multiplicity, short-read bridge ranking, loop unrolling, overlaps, and their assumptions.
- Map each concept to the fastBact graph architecture and interfaces, without changing Rust production code.
- Identify unavailable data (e.g. SPAdes contig paths are not native in fastBact) and replacement evidence needed (own PE links/read paths).
- Record expected benefit, worst-case biological failure mode, implementation complexity, tests, and license provenance.

Deliver `U0_DESIGN_GAP_MATRIX.tsv`, `U0_SOURCE_MAP.md`, `U0_VERDICT.md`.

### U1 — External pipeline benchmarking (no production algorithm change)

- Add Unicycler Illumina-only to SPAdes, Shovill, and fastBact benchmark matrix.
- Pin Unicycler and external dependency versions and record nested SPAdes version.
- Distinguish environment failure from biological failure; do not score crashed or unsupported samples as accuracy losses.
- Include confirmed single-isolate truth data; keep strain-specific references and known engineered edits.
- Keep C1b default, baseline default, and baseline `--no-subsample` in result table.

Deliver `U1_PIPELINE_MATRIX.tsv`, raw logs, `U1_VERDICT.md`.

### U2a — k-by-k graph diagnostics, observation-only

Export after each graph build/upgrade: k, retained solid-kmer fraction, graph unitigs, edges, dead ends, branch nodes, assembled span, reads supporting bridge/branch, maximal graph components. Score may use Unicycler-inspired contig/dead-end heuristic as **advisory** only. Do not permit this heuristic to override biological evidence.

Deliver `U2_GRAPH_QUALITY.tsv`. Gate: graph metrics reproducible and correctly counted for circular and linear replicons.

### U2b — Evidence-guarded adaptive k selection, separately tested

- Compare existing schedule vs configurable alternatives, e.g. [31,55,79,99] and [31,55,79,103,127], one schedule per run.
- Choose based on supported k-mer spectrum, local read length after trimming, graph impact, repeat separation, and branch confidence; merged fraction alone is insufficient.
- Do **not** inject all previous unitigs as 'solid backbone' or require dead ends never increase. Preserve old graph only where independent read evidence supports it. Explicitly allow legitimate splits of unsupported joins.
- This stage does **not** modify trimming Q threshold or CMS concurrently.

Deliver `U2_K_ABLATION.tsv` and misassembly/graph diffs. Gate: no new independently verified false joins, no unacceptable genome-fraction loss, benefit on held-out isolate.

### U3a — Non-destructive multiplicity estimation, shadow mode

- Represent unknown multiplicity explicitly: `Option<u16>` or an enum with confidence and provenance; do not force 1 for every unresolved unitig.
- Distinguish vertex copy counts from directed **edge traversal counts**. Flow consistency is a conditional inference, not universal exact Kirchhoff conservation over every coverage-distorted node.
- Estimate chromosome baseline and per-replicon depth only where topology/read evidence separates them; plasmid-chromosome shared repeats can connect components and invalidate component-specific median assumptions.
- Provide synthetic truth of 1,2,3,5 repeats, low/high plasmid copies, shared insertion sequences, nonuniform coverage, and low-depth disconnected subgraphs.

Deliver `U3_COPY_NUMBER_TRUTH.tsv`, calibrated confidence errors, `U3A_VERDICT.md`. No pruning/bridging behavior change at this stage.

### U3b — Multiplicity-informed graph cleaning, single change

Only if U3a passes: use **high-confidence** multiplicity estimates for a narrowly scoped graph cleaning change. Unknown nodes remain intact. Unit tests and GFA topology invariants must pass; compare with U3a shadow-mode baseline. Do not modify PE bridging.

Deliver `U3B_GRAPH_ABLATION.tsv` and variant-level audit.

### U4a — Candidate bridge evidence and rejection (shadow mode)

Generate proposed bridges from native fastBact PE links/read paths. For each include anchors/orientation, candidate path, insert-size agreement, unique support, competing candidates, coverage consistency, multiplicity confidence, maximum repeat length, reasons to accept/reject/abstain. Quality in [0,100] is a **heuristic ranking**, not a calibrated probability. No graph edit in shadow mode.

Deliver `U4_BRIDGE_CANDIDATES.tsv`, `U4_BRIDGE_TRUTH.tsv`, `U4A_VERDICT.md`. Gate: unambiguous and impossible-to-span repeats are correctly classified; record per-class precision/false-join rates.

### U4b — Transactional conservative bridge application

- Apply one bridge at a time only when independent evidence is sufficiently strong; reject conflict/ambiguity and preserve unresolved repeat nodes.
- Check before and after every edit: orientation consistency, graph reciprocity, overlap sequence consistency, residual edge/copy capacity, circular plasmid integrity, no loss of sequence provenance.
- Permit rollback of one bridge based on local validation. Preserve rejected bridges and explicit rejection reasons.
- Compare against U4a baseline and each accepted bridge against ground truth where available.

Deliver `U4_APPLIED_BRIDGES.tsv`, `U4_ROLLBACK_LOG.tsv`, `U4B_VERDICT.md`.

### U5a — Tandem loop unrolling, separately gated

Do **not** automatically unroll a loop only from rounded multiplicity. First exclude alternate explanations (e.g., a distinct circular plasmid sharing sequence). Require path support or independently resolved copy/flow constraints, test multiple tandem copy counts, and leave ambiguous loops unexpanded. Save exact emitted graph path and evidence.

Deliver `U5_LOOP_TRUTH.tsv`, plasmid false-positive audit, `U5_VERDICT.md`.

### Release Gate — End-to-end regression and publication-ready claims

Compare: (a) immutable v0.1.0 baseline default, (b) immutable baseline `--no-subsample`, (c) B2-C1b, (d) each individually validated U-stage candidate, (e) SPAdes, Shovill, Unicycler. Run synthetic and verified true-isolate data with matched references and documented parameters. Run multi-seed/replicates, record FASTA/GFA SHA-256 and wall time/MaxRSS on same node class.

**Hard blocker:** independently confirmed **new false structural joins** relative to the last accepted quality baseline. A QUAST discrepancy against a mismatched strain reference is `REVIEW`, not automatically a false algorithm join. Do not accept an improvement on N50 alone. No claim of 'zero false joins' without coordinate-level evidence.

## 4. Mandatory benchmark structure and result fields

Base every row on dataset ID, raw FASTQ SHA, read 1/2 pairs, strain, reference SHA, whether sample is a pure isolate, genome size, plasmid truth, config, seed, thread count, commit/tag, hostname, node CPU model, RAM, external versions, elapsed seconds, MaxRSS, exit status, FASTA/GFA SHA, counts retained/discarded, stage metrics.

Quantify at least: genome fraction, NGA50, N50, contig count at >=500 bp, structural misassemblies (categorized with evidence), local misassemblies, SNP/indel errors, plasmid recall/circularity, graph dead-end count, read placement confidence, bridge abstention fraction. Include reasons for NA.

Use at minimum S1-COV-001, S1-REP-002, S1-REP-003, S1-PLAS-002, S2-REAL-001, R1-001, plus a verified pure-isolate test set independent of parameter tuning. B2-C1b must also test high depth 150x/300x and low depth 15x/30x. Do not optimize against one E. coli sample only.

**Suggested B2-C1b quality targets, conditional on validated reference/sample metadata:**
- 100x synthetic estimated genome size within documented calibrated error (initial diagnostic target <=20%); **no erroneous reduction to ~19x**.
- Default C1b <=120x must not downsample when estimate trusted; otherwise safely abstain.
- R1-001, if a pure isolate with suitable truth: target >=98.36% genome fraction, NGA50 >=58,335 bp, contigs <=165, with no new **validated** false joins relative to baseline `--no-subsample`; these are comparison targets, not guaranteed biologically achievable absolute gates.
- Preserve stable behavior on S1-COV-001 and plasmid recovery tests.
- Speed/RSS measured and reported; avoid unconditional memory cap until under full-reads workload benchmarked.

## 5. Git and handoff protocol

For each stage: (1) propose precise hypothesis and limited change set; (2) create feature branch; (3) commit code; (4) run Slurm tests; (5) record raw test manifest; (6) generate verdict; (7) stop for review. Failing experiment remains on isolated branch and is not merged. Each stage uses fresh fixed immutable benchmark comparison input. Never bundle trimming Q relaxation, CMS enlargement, high-k pruning, and bridge edits into one causal A/B commit.

Do not rename diagnostic-only stage as finished algorithm integration. Avoid `.gitignore`-ing required reproducibility artifacts; ignore cloned reference source and large intermediate datasets only, save upstream SHA and lightweight structured metrics in version control.

## 6. Required immediate Gemini instruction

Read this document and current repository. Execute **only B2-G0** first. Do not modify production algorithms, do not run commands on login node, and do not start C1a/C1b/U stages until B2-G0 audit evidence and next-step verdict are complete. Preserve frozen tags and all B1 artifacts. Output Git state, exact baseline tag target, 522-kb/626.4x provenance, sample gzip/decompressor behavior, confirmed per-sample references, raw R1-001 misassembly classifications, and a PASS/FAIL/UNKNOWN summary. Then stop.
