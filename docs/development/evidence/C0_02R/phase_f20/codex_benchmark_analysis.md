I could not save `docs/development/evidence/C0_02R/phase_f20/codex_benchmark_analysis.md`: this session permits read-only filesystem access. Local evidence retrieval also stopped when the sandbox launcher became unavailable. The report below distinguishes verified findings from causal hypotheses and identifies the analyses that remain unverified.

# fastBact Phase F20: Scientific Interpretation and Mechanical Analysis

**Baseline A:** `b9727c3ab372f7b995e439d92b72c7a4bfd051cf`  
**Candidate B:** `0606be423d187771f5925a83da2cc377c40c2a3d`  
**Disposition: C0-02R remains strictly NOT ACCEPTED.**

## 1. Principal interpretation

Candidate B changes both the evidence used to select repeat traversals and the treatment of unresolved graph material. It also changes which graph nodes participate in consensus polishing. Consequently, F20 is a comparison of bundled algorithm changes, not an isolated experiment on `classify_hop`.

The strongest supported interpretation is:

- B makes extension decisions more sensitive to authentic paired-read orientation, fragment geometry, and incompletely explored competitors.
- B preserves branches that earlier cleanup could remove, increasing residual active-node counts independently of merge withholding.
- R1-003 benefits in several final assembly metrics, consistent with better traversal selection.
- R1-001 demonstrates that these changes do not guarantee improved reference correctness.
- Both datasets retain more stage-05 graph sequence but emit slightly less final sequence. Output filtering is a plausible contributor to the genome-fraction reductions, but its exact contribution is unmeasured.
- Base-accuracy changes cannot be separated into resolver, polishing, and alignment effects from the final QUAST metrics alone.

The [metrics table](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/docs/development/evidence/C0_02R/phase_f20/assembly_metrics.tsv) and [acceptance results](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/docs/development/evidence/C0_02R/phase_f20/acceptance.tsv) establish seven regressions across fourteen comparisons.

## 2. Where the measured divergence begins

Within each dataset, A and B have identical reported input counts, trimming/merging summaries, coverage estimates, and graph summaries through the final multi-k stage.

| Dataset | Pre-exSPAnder active unitigs, A/B | Pre-exSPAnder graph bp, A/B | Effective coverage |
|---|---:|---:|---:|
| R1-001 | 5,597 / 5,597 | 5,340,838 / 5,340,838 | 149.36× |
| R1-003 | 2,133 / 2,133 | 4,853,140 / 4,853,140 | 263.13× |

Subsampling was disabled in all four runs. Thus differential read subsampling does not explain the observed outcomes.

The first divergence in these **aggregate measurements** appears at stage 05. Matching counts do not independently establish identical upstream graph sequences or topology. Furthermore, stage 05 follows repeat resolution **and subsequent graph simplification**, so this stage comparison does not isolate the resolver alone.

Source: the four requested `stage_telemetry.json` files; frozen B `src/main.rs:435`, `437`, and `443–456`.

## 3. Search admissibility and conservative withholding

### Evidence changes before classification

A uses immediate-exit support counts and best/second support comparisons. B instead evaluates routes using coordinate-aware paired-read evidence.

For a supporting pair \(p\) and route \(P\):

\[
\operatorname{span}(p,P)
=d_{\mathrm{exit}}(p)+\operatorname{sep}(P)+d_{\mathrm{entry}}(p).
\]

Only spans within the permitted insert interval contribute. Each physical read-pair ID contributes at most once. Evidence is weighted by consistency with the estimated insert distribution:

\[
w_p=\exp\left[
-\frac12\left(
\frac{\operatorname{span}(p,P)-\operatorname{median}}
{1.4826\max(\operatorname{MAD},1)}
\right)^2
\right].
\]

B also changes directional mapping, terminal seed indexing, and insert estimation. Therefore a route can lose support because its implied fragment geometry is inconsistent, even when raw linkage counts previously favored it.

Sources: frozen A `src/resolve/exspander.rs:215–264`, `293–325`; frozen B `src/resolve/exspander.rs:131–179`; A→B `src/resolve/mapper.rs` diff.

### What \(\beta(\mathrm{curr})\) means

The evidence floor is:

\[
\beta(\mathrm{curr})
=\min_{p\text{ anchored at curr}}
\left(d_{\mathrm{exit}}(p)+d_{\mathrm{entry}}(p)\right),
\]

with zero used when no anchored evidence exists.

It is a lower bound on fragment-coordinate contribution, **not a minimum-support threshold or confidence probability**.

For short repeated-node reentry, B can safely prune when:

\[
\operatorname{sep}+\beta(\mathrm{curr})>\operatorname{max\_insert}.
\]

Under this bound, recorded anchor evidence cannot support an in-window endpoint. This avoids treating an impossible continuation as an unresolved competitor.

The implementation does **not** indiscriminately apply this bound to every first visit. Ordinary first visits retain the raw separation test; the principal depth-limit check passes \(\beta=0\).

Source: frozen B `src/resolve/exspander.rs:220–236`, `250–322`, `470–485`.

### Why withholding increases

The unified classifier distinguishes sequence-invalid transitions, safely prunable transitions, live paths, repeated-node probes, unresolved circular continuation, and evidence-bearing paths that cannot be committed again.

Search is bounded: four path steps, 32 candidate paths, 256 visited states, and queue capacity 1,024. Extension also has a 64-step limit. An unresolved continuation or exhausted search budget can withhold an extension instead of accepting a best-so-far route.

The central improvement is that **an unexplored competitor is not assumed to have zero support**.

Previously consumed paths can remain visible as evidence-bearing “shadows.” This reduces the chance that earlier seed processing makes an alternative disappear from later competition.

Sources: frozen B `src/resolve/exspander.rs:181–185`, `371–624`, `668–679`, `717–723`, `823–840`, `963–978`.

Biologically, coordinate-aware evidence can distinguish alternative repeat exits when they imply different fragment lengths. Repeats lacking informative terminal anchors, possessing several compatible exits, or extending beyond the search horizon can remain unresolved. The code’s ≥500 bp “unique” category is a length heuristic, not a demonstrated guarantee of genomic uniqueness.

### Why residual unitigs are not a withholding counter

B additionally protects unresolved frontier nodes and external neighbors of merged paths from orphan cleanup. These include incoming side branches and neighbors of intermediate path nodes.

B cleanup removes short non-circular nodes only when they are unprotected and lack active connectivity. Consequently, additional residual nodes can reflect **preservation**, not merely failed merges.

Sources: frozen B `src/resolve/exspander.rs:682–699`, `1449–1492`, `1512–1539`.

A useful normalization is:

| Dataset | A residual/pre-stage-05 nodes | B residual/pre-stage-05 nodes | B/A residual count |
|---|---:|---:|---:|
| R1-001 | 5.63% | 30.87% | 5.49× |
| R1-003 | 14.06% | 33.01% | 2.35× |

B retains similar proportions on the two datasets. The larger R1-001 fold increase partly reflects its larger starting graph and A’s much smaller residual fraction. It does not establish proportionally stronger withholding in B on R1-001.

Nor do the extra 1,413 and 404 nodes represent that many independent genomic repeat loci. Their sequence identity, redundancy, and reference relevance require direct examination.

The zero ambiguous/invalid rejection counters do not contradict withholding: these counters cover particular junction-selection and validation events, not every search truncation, cycle concern, or competing-path veto. Sequence-valid overlaps also cannot establish correct repeat-copy placement.

## 4. Why R1-003 improves while R1-001 regresses

### R1-003

Observed changes include:

- Misassemblies: 2 → 1.
- Mismatches/100 kbp: 1.51 → 1.08, down 28.48%.
- Indels/100 kbp: 0.04 → 0.00.
- NGA50: 63,648 → 70,010 bp, up approximately 10%.
- Final contigs: 135 → 144.

These results are consistent with avoiding incorrect traversal and improving the distribution of correctly aligned blocks, despite fewer committed junctions and more residual graph nodes.

However, it is not established that a particular `classify_hop`, cycle, shadow, or evidence-floor decision prevented the disappearing misassembly. That requires linking QUAST breakpoints to A/B graph paths and resolver decisions.

The improvement is also metric-specific. R1-003 N90 falls from 18,556 to 16,539 bp, and genome fraction declines. NGA50 improvement therefore does not mean uniformly better contiguity or completeness.

### R1-001

Observed changes include:

- Misassemblies: 4 → 5.
- Mismatches/100 kbp: 2.14 → 2.38, up 11.21%.
- Indels/100 kbp: 0.48 → 1.03, up 114.58%.
- NGA50: 58,318 → 57,127 bp, down 2.04%.
- Final contigs: 164 → 167.

Conservative search does not make B’s joins a simple subset of A’s. Mapping, scoring, search competition, seed-order effects, cleanup, and polishing eligibility all change. B can withhold many uncertain routes while selecting a different erroneous route elsewhere.

R1-001 also enters stage 05 with more nodes and lower estimated coverage than R1-003. Those observations are compatible with differences in graph fragmentation and local evidence availability, but do not establish the causal repeat structure. Global coverage does not measure informative, uniquely anchored support at a particular repeat.

**Answer to the proposed tradeoff:** conservative withholding can prevent chimeric traversal while preserving unresolved branches. R1-003 is consistent with that benefit. R1-001 shows that this mechanism is not a universal correctness guarantee; aggregate metrics cannot identify the precise failed or avoided joins.

## 5. Genome fraction and sequence retention

### Two different length thresholds

The assembler defaults to a **500 bp output threshold**. The benchmark supplies no overriding `--min-len`.

QUAST evaluates the emitted assembly with **`--min-contig 200`**. It cannot recover 200–499 bp sequences already excluded from assembly output.

Sources: `src/cli.rs:35`; frozen B `src/polish/output.rs:65–68`; [benchmark commands](/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/docs/development/evidence/C0_02R/phase_f20/benchmark_commands.sh:251).

### More graph sequence, less emitted sequence

| Dataset | Stage-05 graph bp A → B | Final assembly bp A → B | Final aligned bp A → B |
|---|---:|---:|---:|
| R1-001 | 4,601,762 → 4,808,057 | 4,576,828 → 4,573,069 | 4,576,583 → 4,572,760 |
| R1-003 | 4,593,065 → 4,660,696 | 4,563,880 → 4,558,923 | 4,562,575 → 4,557,618 |

B retains **206,295 additional graph bp** on R1-001 and **67,631 additional graph bp** on R1-003, while emitting 3,759 and 4,957 fewer bases respectively.

This supports a distinction between graph retention and final representation. It does not establish that all retained sequence is valid, unique reference sequence.

Other output operations include terminal overlap trimming, coverage filtering, homopolymer filtering, and duplicate removal. Thus the graph-to-output difference cannot be assigned entirely to the 500 bp cutoff.

Source: frozen B `src/polish/output.rs:42–78`, `94–134`.

### Important telemetry qualification

Stage 06 sets `active_unitigs` to the curated contig count, while retaining graph-derived edge and bp fields.

Therefore **1,728 → 167** on R1-001 B is not evidence that polishing deleted 1,561 graph nodes. It switches measurement scope.

Source: frozen B `src/main.rs:480–489`.

### Is the completeness loss real?

The measured reduction in final reference coverage is real under the benchmark evaluation. Its cause is unresolved.

Using rounded genome-fraction percentages, the reductions correspond approximately to:

- R1-001: 3,853 fewer covered reference bases.
- R1-003: 4,493 fewer covered reference bases.

These are net reference-coverage differences, not exact graph-deletion counts. Genome fraction measures reference positions covered by alignments; repeat placement and alignment selection can affect it. [QUAST 5.2.0 manual](https://raw.githubusercontent.com/ablab/quast/quast_5.2.0/manual.html)

**Conclusion:** withholding valid short unitigs is plausible and architecturally supported, but not proven to explain all completeness loss. It is premature to claim either irreversible graph sequence loss or exclusively threshold-induced omission.

## 6. Base accuracy and polishing boundary effects

There are four relevant pathways:

1. Different resolver paths change sequence composition or repeat-copy selection.
2. Different contig boundaries change read recruitment and polishing context.
3. B directly changes polishing eligibility.
4. Different final sequences and alignments change the reference loci included in QUAST error measurements.

The direct polishing change is especially important. B restricts indexing, pileup allocation, and modification to nondeleted nodes meeting the output-length threshold. Its source explains that short residual repeat copies can otherwise capture an initial read seed and deprive emitted contigs of polishing coverage.

Source: frozen B `src/polish/consensus.rs:15`, `58–112`, `203`; A→B polishing diff.

Therefore the experiment contains both an **indirect boundary effect** and a **direct polishing-input change**.

The Phase E boundary-shift mechanism supplied in the task is relevant historical evidence, but does not establish the fraction of F20 errors caused by polishing. The original Phase E forensic artifacts could not be independently examined before local retrieval failed.

QUAST mismatch rates include both true reference differences and errors; indel rates count events rather than simply inserted/deleted bases. Reported R1-003 `0.00` therefore does not establish an error-free assembly. [QUAST 5.2.0 definitions](https://raw.githubusercontent.com/ablab/quast/quast_5.2.0/manual.html)

Aligned-length changes are approximately one tenth of a percent, much smaller than the reported rate changes. Denominator changes alone cannot explain those changes. Nevertheless, omission or reassignment of a small, error-enriched locus remains possible.

Historical comparisons must remain distinct:

| Dataset | Phase D mismatches/indels | F20 A | F20 B |
|---|---:|---:|---:|
| R1-001 | 2.86 / 0.96 | 2.14 / 0.48 | 2.38 / 1.03 |
| R1-003 | 1.67 / 0.04 | 1.51 / 0.04 | 1.08 / 0.00 |

R1-001 shows partial historical mismatch recovery, but worse-than-baseline mismatches and indels. R1-003 improves against both references.

## 7. Attribution limits and disposition

The following checks are necessary for definitive locus-level attribution and were **not completed**:

1. Measure stage-05 length distributions and align omitted sub-500 bp unitigs to B-only reference gaps.
2. Reconcile exact reference coverage intervals, separating gained coverage from lost coverage.
3. Compare matched A/B sequences before and after polishing at shared reference coordinates.
4. Stratify introduced, corrected, and omitted errors by proximity to changed contig boundaries.
5. Link structural breakpoints to resolver paths and withholding reasons.
6. Use controlled component ablations to separate mapper, search, preservation, and polishing changes.

Verified work comprised the supplied summary tables, four QUAST reports, four telemetry records, benchmark configuration, and relevant frozen-code inspection. No new assemblies or benchmark evaluations were run. Only one replicate exists per dataset/version.

The architectural takeaway is that stronger search admissibility and branch preservation improve the integrity of repeat-resolution decisions, while increasing the distinction between retained graph sequence and emitted assembly sequence. They do not guarantee monotonically better biological correctness, completeness, or base accuracy.

**C0-02R remains strictly NOT ACCEPTED.** R1-001 fails six acceptance comparisons, and R1-003 still fails genome fraction. R1-003’s improvements do not offset these failures.