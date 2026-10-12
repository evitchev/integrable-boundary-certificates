# PILLOW1: non-local maps at the exact point between Sol 1 / Sol 2 opers and the LZ pillow. Opus 5 (2), 2026-10-02.
Seal e49f9a98 (sent before computing). Exact sympy; every run capped (ulimit -v 16 GB); full logs.

## Verdict: NO transform found. Every sealed route has a stated obstruction. Two new structural facts.
KEY INVARIANT (sealed):
- Mellin multipliers, which include every Euler, Erdelyi-Kober, Riemann-Liouville and Meijer-G transform, preserve two things: the number of Mellin-recurrence terms, and the shift set together with its E-carrying subset, up to affine maps.
- So do overall factors, power maps x -> x^mu (including inversion) and algebraic gauges.
### Route (1): the pillow in Mellin form (r1_pillow_mellin.log, r1b_invariants.log)
- Only the placements w = x^(+-1) admit a finite recurrence: w^(1/2) branches at w = 0 and infinity. Then x = s^2.
  - In s the pillow has 5 terms {0, 2, 3, 4, 5}, with kappa^2 at {3, 5}.
- NEW FACT: with the momentum-dependent algebraic gauge (1-s)^h (1+s)^h, h = +-sqrt(-X/2) (the X-point exponent), (1 - s^2) divides out. The pillow becomes a THREE-TERM recurrence:
      s^0: T^2 - 1/4   (momentum-free);   s^2: -[(2T + 2 sqrt2 sqrt(-X) + 2)^2 - (4Y + 9)]/4;   s^3: -4 kappa^2.
  - Asymmetric gauges give 4-5 terms. The minimal form is unique up to the sign of h.
- Invariant: for three terms with non-E shifts {p, q} and E shift e, the pair {(q-p)/(e-p), (p-q)/(e-q)}.
  - Pillow: {-2, 2/3}.
  - Sol 1 Gamma / three-term form at t = -1/3: {-1/2, 1/3}.
  - Sol 2 at t = +1/3: {-1, 1/2}.
  - Sol 1 Schrodinger (Suzuki) form: 4 terms, unequally spaced {-2, -1, 0, 2}, E at an end. The pillow's 4-term forms ({0, 1, 2, 3}, equally spaced) do not match.
- **OBSTRUCTION: no Mellin multiplier, power map or algebraic gauge relates the pillow to Sol 1 or Sol 2 at the exact point.**
- POST-SEAL (labelled): the pillow's invariant EQUALS Sol 1's Gamma-form invariant at t' = +1/3 and at t' = -5.
  - Shift roles at t' = 1/3: the pillow's two-momentum s^2 term <-> Sol 1's base H; the momentum-free s^0 term <-> Sol 1's coefficient-1 shift; E <-> E.
  - The full cocycle test (r1c_crosspoint.py) has a free offset delta and a FREE momentum dictionary (l0, l1). It gives NONE at t' = 1/3 and t' = -5, both gauge signs.
  - Positive control (r1c_control.py): Sol 1's own three-term ODE against its Gamma form (PROOF1) gives F == 1 at t' = 1/3 and 2. The L + 1/7 tamper fails.
  - So the invariant coincidence is not an operator relation.
### Route (2): the +-alpha doublet (r2_doublet.log)
- NEW FACT: the tensor product of Suzuki(+alpha) and Suzuki(-alpha) at M = -2 is a kappa-EVEN fourth-order operator with only THREE Mellin terms. The alpha-odd parts cancel and alpha enters only as alpha^2 = -8X.
  - Its shifts (in z = 1/x) are {0, 2, 4} relative, with E at an end.
  - Its charges are the SUM over the doublet, so after monic normalisation they are the certified (doublet-average) charges.
  - This is an explicit kappa-even operator realising the doublet average.
- But its pattern (equally spaced, E at an end) does not match the pillow's {0, 2, 3}. Sym^3(pillow) has 8 terms, and pillow (x) first-order is the pillow's own pattern.
- OBSTRUCTION (same invariant): the pillow is not a multiplier image of the doublet tensor product. Its symmetric cube, or its product with a first-order factor, is not one either.
- Pillow (x) a general second-order factor: NOT tested.
### Route (3): Laplace with the pillow's X-point at p = infinity (r3_r_expgauge.log)
- L4 = ((p^2+E) chi)'''' - lam chi'' + alpha chi' - chi.
- Newton polygon at infinity: points (0,0), (1,-1), (2,-2), (3,-2), (4,-2), (4,-4). The upper hull is the single edge (0,0)-(4,-2) of slope -1/2; the other points lie strictly below it.
- So all four formal solutions are exp(c p^(1/2)) type, with NO slope-0 (regular) part.
- A right factor regular-singular at infinity (the pillow's X-point, exponent difference sqrt2 sqrt(-X)) is impossible: OBSTRUCTION.
- (The stored log's sentence "contains (1,-1) and (3,-2)?" is mis-worded: they lie BELOW the edge. The conclusion is unchanged.)
### Route (4): Euler / fractional-integral transforms
- These are Mellin multipliers, so they are covered by (1)'s obstruction.
- Structurally: B1's Erdelyi-Kober map works at M = 1 because Suzuki is then Kummer (alpha is an energy shift). At M = -2 it is a doubly-confluent Heun equation (irregular at 0 and infinity) with a 4-term recurrence, and no such parameter shift exists.
### Exponential gauges (sealed as tested, r3_r_expgauge.log)
- exp(beta s^h), h = 1, 2, 3, on the pillow's 3-term form: generic beta gives 7, 5 and 9 terms.
- For any h > 0 the new top shift 2h + 2 carries -beta^2 h^2 alone, so it cannot cancel. Then E (shift 3) is interior, while in every Sol 1 / Sol 2 form E sits at an extreme shift: OBSTRUCTION for h > 0.
- h < 0: NOT tested.
## Not covered (stated in the seal)
- kappa- or E-dependent coordinate maps.
- Exponential gauges with h < 0.
- Tensor products of the pillow with a general second-order operator.
- Transforms that are not Mellin multipliers, power maps, algebraic or exponential gauges, or Laplace.
## Disclosures
- r1c v1 hard-coded one shift role and crashed on t' = -5 / -1/3 (assertion). v2 derives the roles; the v1 logs were deleted before the v2 run.
- r2 uses z = 1/x for the doublet.
- The 'roles' line in r1c logs is printed as 'consistent role maps'.
