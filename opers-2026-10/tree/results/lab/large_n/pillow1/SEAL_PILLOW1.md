# PILLOW1 seal: a non-local map at the exact point between Suzuki + centrifugal (M = -2, Sol 1), Sol 2's n = 2 ODE, and the LZ pillow. Opus 5 (2), 2026-10-02. Before any computation.
Objects:
- Pillow (38(ad), t = -1/3): psi_zz = [kappa^2 w^(3/2)(1 - w) + V(w)] psi, with d/dz = w(1-w) d/dw and V = -X w/2 - Y w(1-w)/4 - (5/16) w(1-w) + (1/16)(1-w).
- Suzuki(M = -2): (theta^2 - theta - lam) psi = (x^-2 + alpha x^-1 - E x^2) psi, with alpha^2 = -8X and lam = 2 + Y.
- Sol 2 (k = -1/2): L3 of RED2 A3.
KEY INVARIANT (the basis for routes 1 and 4):
- A Mellin multiplier (this includes every Euler, Erdelyi-Kober, Riemann-Liouville or Meijer-G convolution transform), an overall nu-dependent factor, power maps x -> x^mu and algebraic gauges x^g all PRESERVE the number of nonzero terms of the Mellin recurrence, and the shift set up to a common scale.
- Exponential gauges exp(beta x^h) do NOT preserve it. They are treated separately: E-independent beta and h are tested for a reduction of the pillow's term count.
Routes:
(1) Put the pillow into every Moebius coordinate (6 placements of its three singular points w = 0, 1, infinity at 0 / infinity) with the minimal power map that makes its coefficients polynomial. Count the Mellin terms and compare with Sol 1 (Gamma form: 3 terms; Schrodinger form: 4) and Sol 2 (3 terms).
    - Equal counts are necessary but not sufficient. Where they match, I test the shift-set ratio, then the PROOF1-type cocycle.
    - Prediction: count >= 4 in every placement, or a ratio mismatch. Route (1) then gives an OBSTRUCTION (80%).
(2) Doublet: the tensor product of Suzuki(+alpha) and Suzuki(-alpha) (4th order). Its charges ADD, so it is charge-isospectral with the pillow after monic normalisation.
    - Compare it, via the invariant and the local exponents, with the 4th-order relatives of the pillow: the symmetric cube, and the pillow tensored with a first- or second-order factor.
    - Prediction: obstruction 75%.
(3) Laplace with the pillow's X-point sent to p = infinity: the Newton polygon of L4 at infinity.
    - Prediction: L4 has NO regular-type (slope-0) formal solution at infinity, so the pillow's regular X-point (exponent difference sqrt(2) sqrt(-X)) cannot be placed there (for algebraic identifications). Obstruction 85%.
(4) Euler / fractional-integral transforms are Mellin multipliers, hence covered by (1)'s invariant.
    - Plus a structural statement: B1's Erdelyi-Kober map exists because at M = 1 Suzuki is confluent hypergeometric (Kummer). At M = -2 it is a doubly-confluent Heun equation (irregular at 0 and infinity), and that mechanism does not apply. Verified via (1)'s counts.
Overall prediction: no transform found in these classes, 85%. A found transform is checked identically in (E, X, Y).
Limits stated in advance: kappa- or E-dependent coordinate maps, and exponential gauges beyond the tested (beta, h), are not covered.
Caps: ulimit -v 20 GB; full logs.
