# Mutual-commutativity checks among promoted charges (2026-07-29/30)

Pairs checked from the **stored** ℚ(t)-symbolic densities
(`qt_promotion.json`; fixed rational sections, so many-point
vanishing is strong exact evidence --- not a completed certificate:
no post-substitution degree bound is retained), plus the certified
A₂⁽²⁾ I₇ closed form and per-point normalized kernel solves for the
cylindrical I₅. Script: `code/mutual_checks.py` (fail-closed; the
skipped point for Solutions 2/3 is the t=3 decoupling locus, by
design).

## Verified — all commutators vanish at every evaluated point

| pair | Sol 1 | Sol 2 | Sol 3 |
|---|---|---|---|
| cyl [I₅, I₁₁] | 71 pts | 70 pts | 70 pts |
| cyl [I₇, I₉] | 71 pts | 71 pts | 71 pts |
| A₂⁽²⁾ [I₇, I₁₁] | 71 pts | — | — |

≈ 500 exact commutator evaluations, zero failures. Together with the
earlier direct checks [I₅,I₇] = [I₅,I₉] = 0 (two points,
`cyl_mutual_check.py`) and the defining [I₃,·] = 0 relations, the
verified commuting sets are now:
cylindrical {I₁, I₃, I₅, I₇, I₉} plus [I₅,I₁₁];
A₂⁽²⁾ {I₁, I₅, I₇, I₁₁} pairwise except [I₇,I₁₃]-class pairs.

## Attempted and stopped — out of memory

cyl [I₇,I₁₁] and [I₉,I₁₁] (all three solutions): the six workers were
OOM-killed (5 kills logged in dmesg; 125 GB host) — the per-monomial
residue lru_cache for 8×12- and 10×12-letter monomial pairs grows to
tens of GB per worker. The multiprocessing pool then hung on the dead
children; no results were lost (fail-closed), the seven completed
tasks above are recorded in the run log.

## Deferred by combinatorial size (as before)

A₂⁽²⁾ [I₇,I₁₃] and [I₁₁,I₁₃] (~10⁹–10¹⁰ Wick matchings for the worst
monomial pairs).

## Fix path (single engineering item, [road-map])

A symmetry-collapsed Wick enumerator (count matchings by
letter-multiplicity classes instead of raw positions) with bounded /
streaming accumulation would make all eight remaining pairs and the
σ = 12/14 paperclip vacancies feasible. Until then, the unchecked
pairs are labeled as such everywhere.
