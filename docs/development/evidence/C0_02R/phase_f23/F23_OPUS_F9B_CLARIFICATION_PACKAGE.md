# F23 Opus F9(b) Clarification Package

Status: `OPUS_CLARIFICATION_REQUIRED`

Scope: resolve the cycle/re-entry boundary needed to preserve the unchanged T34
oracle. This package does not request a production change or a new paid API
call. C0-02R remains `NOT ACCEPTED`.

## 1. Opus wording

The decision says: “Re-entering a node already on the current walk marks shared
context. Cycles and tandem copies have ambiguous copy count.” Its F9(b) test
specification says: “A cycle re-entering L causes the fork to be withheld.”
The ruling does not define F9(b)'s exact graph, PE rows, declared copy truth,
whether re-entry occurs before the fork or only in a speculative candidate,
or whether attributable upstream evidence can resolve that candidate.

## 2. Unchanged T34 fixture and oracle

Source: `tests/test_c0_02r_state_semantics_gate.rs:3312-3374` in the F23 source
worktree at `f121920a9f3d76a774e8a382f2f3047f0c0c631e`.

Forward oriented nodes are S(700 bp), P(100), C(100), D(300), X(300), with
31 bp overlaps. Edges are S→P, P→C, P→X, C→C, C→D and C→X, plus the graph's
reciprocal reverse-complement edges. The unchanged oracle requires contig
S-P-C-C-X and requires D to remain unconsumed.

Insert distribution: `[150,250]`, median 200, MAD 20. PE rows:

| Directed key | Count | Endpoint distances | Relevant path separation | Span | Result |
|---|---:|---:|---:|---:|---|
| S→C | 20 | 80 + 82 | via P: 38 | 200 | Supports entry into C |
| C→X | 20 | 75 + 75 | direct: −31 | 119 | Outside window |
| C→X | 20 | 75 + 75 | via C→C: 38 | 188 | In window |
| C→D | 5 | 80 + 151 | direct: −31 | 200 | In window |
| C→D | 5 | 80 + 151 | via C→C: 38 | 269 | Outside window |

The resulting score ratio is about 3.69 in favor of C→C→X. The coordinate
arithmetic distinguishes the candidate path length. There is no S→X or other
pair row demonstrating a source-specific bridge from S across the cycle to X.
The PE keys are node-level; they do not identify which occurrence of C or
which genomic copy supplied C→X evidence.

The test helper also reuses pair IDs across unrelated directed keys:
S→C uses 34001–34020, C→D uses 34002–34006, and C→X uses 34003–34022. A real
`map_paired_reads` row writes an oriented key and its reverse-complement
reciprocal (`src/resolve/mapper.rs:208-249`), so these overlapping IDs are not
a physically consistent mapped-read ledger. T34 remains an unchanged
engineering oracle; this limits its use as biological source-attribution
evidence.

## 3. Walk and candidate-search states

At T34's second decision, `ext` contains the previously selected local prefix
S-P-C; it is not yet a graph-committed contig. Candidate search starts at C.
The first C→C hop is at search depth zero and is admitted; the later C→X
candidate is scored at separation 38. The selected steps are appended to `ext`,
and final graph mutation occurs after path construction and validation
(`src/resolve/exspander.rs:347-374`, `:1071-1095`, `:1500-1551`).

This distinction between prior extension history and speculative candidate
re-entry is observable. It does **not** by itself prove that C→X evidence
belongs to the S context. A rule that withholds every candidate self-loop
would reject T34; a rule that checks only already-selected prefix history can
accept a false source attribution when the shared C→X rows come from another
copy.

## 4. Paired counterexamples

These are declared-truth, matrix-level constructions, not executed tests or
raw-read claims. Use globally distinct pair IDs and reciprocal keys in both
worlds.

### Counterexample A: blanket cycle rejection loses an attributable bridge

Use S→P→C, C→C, C→X, C→D with the same lengths/overlaps as T34. Let the true
S context be S-P-C-C-X. Add source-specific S→X pairs with endpoint distances
75 + 75. The cycle path has separation 176 and span 326; the shorter S-P-C-X
path has separation 107 and span 257. Under insert window `[300,350]`, only the
cycle path is supported. If the candidate search remains complete and the
records are keyed at the attributable S node, a blanket “any re-entry means
withhold” rule rejects a distinguishable legal bridge. This is a falsifier for
unconditional cycle rejection, not for withholding when attributable evidence
is absent.

### Counterexample B: prefix-history-only checking can accept a chimera

Keep the same visible graph: S→C, B→C, C→C, C→D and C→X, with S and B as
distinct source flanks. In world 1, true context A is S-C-C-X and B-C-D. In
world 2, A is S-C-D and B-C-C-X. Keep the same node sequences, insert
distribution and PE coordinates/counts; allocate globally distinct pair IDs
and reciprocal keys in both worlds. S→C and B→C identify entry into the shared
node, while C→X and C→D remain pooled at C. The visible graph and matrix are
identical, but S-C-C-X is correct only in world 1. No current `PairLinkMatrix`
field distinguishes which copy supplied the C-keyed evidence. This is a
declared-truth information counterexample, not an executed raw-read
reproduction. Removing or clipping B's visible flank gives the ruling's
additional hidden-context risk.

## 5. PairLinkMatrix information boundary

`ReadPairEvidence` stores pair ID, distance from the first node's exit,
distance to the second node's entry, and read lengths. `PairLinkMatrix` keys
records by `(DirUnitigId, DirUnitigId)` and stores no occurrence index, copy
origin, placement alternatives, or path-history identity
(`src/resolve/mapper.rs:32-85`). The resolver scores the exact node key and
path separation (`src/resolve/exspander.rs:146-179`). Consequently, distance
can discriminate insert-compatible path lengths, but it cannot distinguish
two genomic origins that yield the same node key and coordinate record.

## 6. Independent review

OpenCode GLM-5.3 was attempted read-only but returned `ECONNRESET`; this is not
quota evidence. The authorized fallback analysis used Codex GPT-6.1 Sol Medium
(selected through the independent agent runner; provider-side identity is not
separately authenticated). It returned
`INDISTINGUISHABLE_WITH_CURRENT_DATA` for T34 source attribution, while
limiting that conclusion because F9(b) has no complete fixture. A separate
Codex GPT-6.1 Terra Medium read-only review independently agreed on the
observable-data limit and T34 arithmetic. It confirmed that the candidate vs
prefix distinction is observable but insufficient to prove source
attribution. Both reviews found that the same graph and matrix can accompany
different copy truth, and that T34 itself lacks an upstream source-specific
S→X record.

## 7. Required technical ruling

For a candidate path that re-enters the same oriented/physical node, under what
observable conditions may PE evidence safely resolve the path, and when must
the resolver withhold? Specifically, is an in-window pair keyed at an upstream
attributable anchor required to preserve T34-like self-loop bridges, and does
T34's current synthetic node-keyed evidence satisfy that requirement? The
ruling must define whether F9(b) concerns pre-fork selected-walk re-entry,
candidate-path re-entry, or both, while retaining T34 unchanged.
