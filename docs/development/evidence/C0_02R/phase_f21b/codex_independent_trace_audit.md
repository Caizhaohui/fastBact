**Independent audit: Phase F21-B**  
**Audit ID:** `AUDIT-20261023-017`  
**Date:** 2026-10-23  
**Overall Verdict:** **FAIL**  
**C0-02R:** **Strictly NOT ACCEPTED**

The evidence supports several observations, but material contradictions prevent certification of the claimed decision mechanism and coverage-loss attribution.

**Scope and integrity**

I directly inspected all 11 requested TSV/Markdown files, recomputed pair-ledger aggregates, and cross-checked relevant Phase F21 tables. No raw FASTAs were printed, and no files were modified.

The production `src/` tree is **100% frozen**:

```text
30967bf4940844aeb6b744732bcad6d2a72daab4
```

`HEAD:src` and `07e0e5fd927cc3c8be96b22927e33b8ae6b9e138:src` both resolve to that object. Both source comparisons returned clean, and no untracked source files were found.

**1. Configuration parity: recorded, incompletely verified**

`effective_parameter_audit.tsv` explicitly records **`BENCHMARK_CONFIGURATION_VERIFIED`**. A/B settings agree for CLI arguments, resources, k-mer schedule `[31, 63, 127]`, retained pairs, and assembly/evaluation cutoffs. The differing insert-size estimates are identified as an algorithmic difference.

However, the referenced `benchmark-results/C0_02R_f20_selective_20261022` directory is absent from this worktree. Additionally, the earlier `forensic_manifest.tsv` records different CLI syntax. Configuration parity is supported by the reconciliation table but cannot be independently certified from the underlying execution records available here.

**2. Repeat geometry: arithmetic confirmed; locus specificity contradicted**

The directed table lists seven intermediate nodes:

```text
u1926 → u3275 → u1475 → u1476 → u4349 → u4465 → u3919
```

Their raw lengths total **2,092 bp**, and every listed incoming edge overlap is **127 bp**.

The stated interflank separation reconciles as:

```text
2,092 − (8 × 127) = 1,076 bp
```

This subtracts all eight overlaps connecting the two flanks through the seven intermediate nodes. The overlap-adjusted union of the seven intermediate nodes themselves is **1,330 bp**.

The **255 bp length of `u1476` is confirmed**. Its claimed exclusive presence at 3.66 Mb is **not confirmed**: both graph tables list a 1.57 Mb occurrence and classify it as `TWO_COPY_HOMOLOGOUS_REPEAT`. This directly contradicts the report’s central asymmetric-anchor premise.

Reference coordinates also disagree between the graph tables, including those for `u1927` and `u3920`.

**3. Pair provenance and candidate decision: material reconciliation failure**

The provenance ledger contains **300 records with distinct pair IDs**.

| Relationship | Ledger pair count | Ledger weight sum | Decision-trace claim |
| :--- | ---: | ---: | :--- |
| `u1927 ↔ u169` | 0 | — | 0 spanning pairs |
| `u1476 ↔ u1926` | 16 | **14.6105** | 16 pairs; **6.3412** |
| `u1476 ↔ u4833` | 0 | — | 0 pairs |
| `u1476 ↔ u3276` | 7 | **3.2492** | 7 pairs; **0.2000** |
| `u3919 ↔ u169` | **64** | 21.9345 | **62** pairs; 21.9324 |
| `u3919 ↔ u3920` | **103** | **35.6914** | **42** pairs; **16.2510** |

Thus, the requested zero counts and the 16-pair count are confirmed **within the supplied ledger**. The claimed decision weights are not reproduced.

Candidate-specific orientation, distance, or eligibility filtering could explain different aggregates, but the supplied artifacts do not identify the contributing subsets or transformations needed to reproduce them.

The documented mechanism is backward traversal entering `u1476`, followed by selection of `u1926` over `u3276` and `u4833`. The trace’s arithmetic, `6.3412 / 0.2000 = 31.706`, is correct. Nevertheless, the contradictory locus annotation and unreconciled weights prevent independent confirmation that this false-anchor mechanism caused the jump.

**4. Authentic alternative exit: documented, not independently established**

`alternative_path_matrix.tsv` identifies **`u3920`, 2,012 bp**, as the authentic 3.67 Mb continuation and records its emission as **A `contig_142`**.

The decision trace and mechanism matrix identify **B `contig_146`**. The alternative-path matrix itself only records B’s outcome as `WITHHELD_OR_OMITTED`.

These annotations support the reported interpretation, but the bundle lacks the sequence-to-output correspondence needed to independently verify the B contig assignment. The pair-support discrepancies further qualify the resolver explanation.

**5. R1-003 Baseline A relocation: confirmed in archived event evidence**

The comparison document is corroborated by `phase_f21/misassembly_event_matrix.tsv`:

- Baseline A: **`contig_16`**, relocation **−2,094,638 bp**.
- Candidate B: **`contig_19`**, the same relocation.
- Both records share contig length **75,502 bp**, break **13066..13067**, and matching flanking reference alignments.

This confirms that the corresponding R1-003 chimeric event already existed in Baseline A. It does not independently prove the proposed search-depth explanation for the cross-isolate difference.

**6. Polishing: terminal concentration confirmed; complete causal claim overstated**

Recounting `phase_f21/base_error_loci.tsv` confirms:

```text
132 / 156 = 84.615%
```

of B’s final **error events** lie in the **101–350 bp** terminal zone. This statistic counts events, not erroneous nucleotide bases.

Of those 132 events, **130** are labeled polishing-introduced and **2** pre-existing.

`polish_boundary_evidence.tsv` contains only **25 representative records**. Its terminal examples show low estimated depth and severe strand bias, supporting the proposed mechanism, but providing no full read-vote reconstruction.

The synthesis also misstates stage arithmetic:

| R1-001 | Stage 5 mismatches/indels | Stage 6 mismatches/indels | Net increase |
| :--- | :--- | :--- | :--- |
| A | 29 / 4 | 98 / 22 | 69 / 18 |
| B | 24 / 4 | 109 / 47 | **85 / 43** |

The matrix’s statement that B introduced **61 mismatches and 39 indels** does not match these stage totals. Net count differences also cannot alone establish every individual causal event.

**7. Unique reference loss: 26,715 bp is not substantiated**

The table reports **26,715 bp**, and its percentage conversion is correct:

```text
26,715 / 4,641,652 × 100 = 0.57555% ≈ 0.576%
```

However:

- The component unique-loss rows sum to **25,574 bp**, leaving **1,141 bp unexplained**.
- Mapped span minus redundant span also gives **32,601 − 7,027 = 25,574 bp**.
- The stated QUAST decrease, **98.835% → 98.281%**, corresponds to approximately **25,715 bp**, not 26,715 bp.
- The table counts **1,561** filtered unitigs; the synthesis says **1,566**.
- The `<200` category has zero mapped span yet claims 100% redundancy. Non-alignment does not establish redundancy.
- The matrix’s “<25 kb” conclusion contradicts its own reported total.

The **234,962 bp raw filtered length must not be equated with unique reference loss**, but the replacement figure of **26,715 bp** is not auditable from this table.

**8. H1–H7 classification**

| Hypothesis | Published classification | Independent audit assessment |
| :--- | :--- | :--- |
| H1: Incorrect local PE chaining | CONFIRMED | Plausible; exact causal chain remains unreconciled. |
| H2: Repeat-internal false uniqueness | SUPPORTED | Qualified support; `u1476` specificity is contradictory. |
| H3: Authentic alternative bypassed | CONFIRMED | Documented interpretation; pair scores and output provenance need reconciliation. |
| H4: Insufficient spanning information | CONFIRMED | No spanning pair is recorded; universal impossibility without long reads is overstated. |
| H5: Terminal polishing errors | CONFIRMED | Error enrichment confirmed; detailed voting mechanism remains supported rather than fully demonstrated. |
| H6: Filtering removes unique coverage | SUPPORTED / DISSECTED | Exact unique-loss attribution fails arithmetic reconciliation. |
| H7: Configuration mismatch | REFUTED | Parity supported by summary; underlying execution verification incomplete. |

**9. Overall synthesis and completion gate**

`phase_f21b_report.md` records **`PHASE_F21B_DECISION_CAUSE_IDENTIFIED`**. I confirm that this status is present, but **do not certify the gate as satisfied**.

The report overstates certainty despite contradictory repeat annotations, unreproducible decision scores, inconsistent polishing arithmetic, and non-reconciling coverage totals. It also claims all 12 deliverables are archived although the independent audit file was absent during inspection. Its assertion that proposed interventions will permanently eliminate the failures is untested.

Certification requires reconciled graph annotations, candidate-specific pair eligibility and scoring provenance, corrected coverage accounting, and supporting output-contig and consensus-vote evidence.

**Final disposition:** **FAIL. C0-02R remains strictly NOT ACCEPTED.**

**Delivery limitation:** This is the complete stdout audit report. The requested `docs/development/evidence/C0_02R/phase_f21b/codex_independent_trace_audit.md` was **not written** because this session enforces read-only filesystem access.
