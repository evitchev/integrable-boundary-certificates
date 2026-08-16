**SUPERSEDED 2026-08-11 by results/vacancy_final.md (21/21 ledger complete).**

# High-spin vacancy: partial results (batch of 2026-08-01/02, in progress)

Route: all promoted charges commute with the low charge identically in
t, so by Jacobi each pairwise commutator lies in ker ad_low at an even
charge spin; vacancy at one curve point (t = 2) is generic by
semicontinuity, and a rational section vanishing generically vanishes
identically — closing pairs with ℚ(t) proofs (`code/vacancy_high.py`).

## Verdicts (updated 2026-08-08; per-task evidence for the GF(p) verdicts in results/vacancy_a2_s24_certificate.{log,json} and results/vacancy_cyl_s18_certificate.{log,json}; the closing pass must MERGE the sigma=20 rerun output vacancy_final_modp_s20.json with the four completed verdicts into the promised vacancy_final_modp.json -- a manual step, per review)

| task | verdict |
|---|---|
| cyl σ=12, 14, 16 (Sol 1, 2, 3 each) | **VACANT** (ker dim 0) ×9 |
| A₂⁽²⁾ σ=18 | **VACANT** (ker dim 0) |
| cyl σ=18 (×3), cyl σ=20 (×3), A₂⁽²⁾ σ=24 | MemoryError in the fraction-free batch (18 GB cap trips; fail-closed); since then GF(p)-certified VACANT: cyl σ=18 ×3 (2026-08-07) and A₂⁽²⁾ σ=24 (2026-08-07); cyl σ=20 ×3 still running |
| A₂⁽²⁾ σ=20 | **VACANT** (ker dim 0; sparse+collapsed stack, 2026-08-03 run) |

## Consequences established now

- **A₂⁽²⁾ [I₇, I₁₁] = 0 identically over ℚ(t)** (σ=18 vacancy +
  Jacobi) — upgrades the earlier 71-point evidence to a proof.
- **A₂⁽²⁾ [I₇, I₁₃] = 0 identically over ℚ(t)** (σ=20 vacancy, first
  verdict of the sparse+collapsed stack).
- Cylindrical even vacancy extended to σ = 16 (all three solutions).

## Still open

Cylindrical [I₉,I₁₁] (σ=20) only.  **A₂⁽²⁾ σ=24: VACANT (GF(p)-
certified 2026-08-07, tags = ∂-rank = 5117) ⇒ [I₁₁,I₁₃] = 0
identically over ℚ(t) — the A₂⁽²⁾ pair program is COMPLETE: every
promoted pair of the family is proved, alongside the generic
infinitude theorem.**  (An earlier
revision of this list also named A₂⁽²⁾ [I₇,I₁₃], contradicting the
verdict recorded above — stale entry removed 2026-08-07, audit
catch.)  The remaining tasks run under the GF(p) vacancy certifier
(`code/vacancy_final_modp.py`; the fraction-free batch was measured
at weeks-scale by a rate probe and killed).  The three cylindrical
σ=18 sectors have already returned VACANT (certified over ℚ), closing
cylindrical [I₇,I₁₁] for all three solutions; the batch's final
artifact (`results/vacancy_final_modp.json` + write-up) will
supersede this partial file.
