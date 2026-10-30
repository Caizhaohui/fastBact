# Phase F24 — GPT-6 Astra Architectural Ruling

**Date**: 2026-10-29  
**Architect**: GPT-6 Astra (Codex `gpt-6-astra`, Medium reasoning, 53,107 tokens)  
**Contract**: [`F24_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md`](../../../../F24_ARCHITECTURE_AND_ACCEPTANCE_CONTRACT.md)  
**Context**: Ruling on Codex GPT-6.1 Sol Medium acceptance review findings on Candidate commit `3d44889`.

---

## 1. Attribution Lattice and Condensation DAG Propagation

The attribution lattice is strictly:
- `0 = Unknown`
- `1 = Single(lineage_id)`
- `2 = Shared`

### Accumulator State Machine ($O(1)$ per edge, $O(d)$ per junction):
```text
Empty
One(lineage_id)
Many (Shared)
```
Transitions:
- `Empty + a -> One(a)`
- `One(a) + a -> One(a)`
- `One(a) + b -> Many` (if $a \neq b$)
- `Many + x -> Many`
- Predecessor is `Many` (Shared) $\implies$ `Many`.

Processing order across condensation DAG: topological order (`(0..scc_count).rev()`).
Long nodes ($\ge 500$ bp) with in-degree 0 introduce their own Single lineage witness. Short nodes without external entries remain `Unknown (0)`.
**Crucial Requirement**: Every deciding-source path must require `Single (1)`. `Unknown (0)` and `Shared (2)` must strictly withhold decisions. `Unknown` cannot be treated as `Single`.

---

## 2. Evidence Admission Predicate ($p > b + H$)

- Centralize the strict source predicate:
  For node $U$ with entry shadow reach $b = \text{entry\_shadow\_reach}$ and horizon $H = \text{max\_insert}$, an evidence read observation starting at coordinate $p = u1\_len - dist\_u1\_exit$ must strictly satisfy:
  $$p > b + H$$
- Any read with $p \le b + H$ could have its paired mate originate inside the ambiguous upstream region $[0, b]$. It must be filtered out.
- Update positive test `test_f24_p02_interval_boundary_exact_reach` and negative test `test_f24_n02` to verify exact boundary at $b + H + 1$ (admitted) vs $b + H$ and $b + H - 1$ (rejected).

---

## 3. Scientific Gate Status & Developer Instructions

- Both R1-003 NGA50 floor (46,800 bp) and base error rates (mismatches / indels) remain mandatory under the contract.
- Developer (OpenCode GLM-5.3) must implement the technical fixes for attribution lattice and interval filtering.
- Re-run full deterministic gates and 3-replicate paired benchmark on Slurm `qfree`.
