# RED2: operator-level reductions at the exact point and the paperclip. Opus 5 (2), 2026-10-01. Seal 1a2b4f82 (sent before computing).
Exact sympy throughout; no solver runs (nothing to cap); full logs. Each claim is labelled OPERATOR or CHARGE level.

## (B2) Paperclip: Sol 3's k = -2 member, and the rotated pair. OPERATOR-LEVEL RELATION FOUND (sealed hypothesis partly right).
- Sol 3's member (item 298 at k = -2 = item 299's operator): S = P4(theta) - C x theta(theta-1) + C E x^2 theta(theta-1), P4 = ((theta-1/2)^2 - l0^2)((theta-1/2)^2 - l1^2).
- The rotated c = -2 pair in the family's own variable is a pair of COULOMB problems (BLZ at M = -1/2, beta^2 = 2): (theta^2 - theta - lam_pm) psi = (x - E x^2) psi.
- The oscillator form in x = y^2 does NOT work (b2_tensor.log): its x and x^2 parts D(D + 1/2) and (D + 1/2)(D + 3/2) can never both become D(D - 1).
- RESULT (b2_coulomb.log, b2_mellin.log):
  - The gauged (x^(1/2)) TENSOR (Hadamard) product T of the two Coulomb problems has EXACTLY S's theta-only part P4 iff lam_pm = ((l0 -+ l1)^2 - 1)/4, i.e. (l_pm + 1/2)^2 = (l0 -+ l1)^2/4: the ROTATED momenta.
  - Its x-parts differ from S's. But in the Mellin recurrence, sum_i p_i(nu - i) phihat(nu - i) = 0, the multiplier
        phihat_T(nu) = m(nu) psihat_S(nu),   m(nu) = rho^nu Gamma(nu + 1/2)/Gamma(nu - 1)
    turns T's recurrence into S's IDENTICALLY in nu, l0, l1 and E. This holds at C = 4 (an x-scale; item 298's normalisation C = 1 is x -> x/4) with x-scale mu = rho.
  - Controls fire: the wrong multiplier Gamma(nu + 1)/Gamma(nu - 1) has no constant solution, and unrotated momenta break P4.
- **OPERATOR LEVEL: Sol 3's k = -2 member = the Mellin conjugate, by Gamma(nu + 1/2)/Gamma(nu - 1) (a fractional, Erdelyi-Kober-type operator of order 3/2), of the tensor product of the two rotated Coulomb problems.** It is not an LCLM, as cc found.
  - This is a statement about the RECURRENCES (formal Mellin level). The contour and branch conditions under which the multiplier maps actual solutions (decaying/regular) to solutions are NOT proved.
- CHARGE level: tensor-product WKB branches add (R = R+ + R-), so the monic charges are those of (1/2) sum, consistent with cc's isospectrality (item 299). Not re-computed here.

## (B1) Paperclip: Sol 1 (Suzuki + centrifugal at M -> 1). EXACT (pencil + machine identities, b1_a1_checks.log).
- At M = 1, H_S(alpha) = -d^2 + y^2 + alpha + l(l+1)/y^2 is the radial oscillator at energy E - alpha. Its spectral determinant depends only on a = (2l + 3 - E + alpha)/4 (Kummer: proportional to 1/Gamma(a)).
- So H_S(alpha) is EXACTLY ISOSPECTRAL with the oscillator at E and l' = l + alpha/2. With the SUPER3b dictionary at M = 1 (alpha^2 = -8X, (l+1/2)^2 = -2Y), a = 1/2 + i(pi + P_X)/sqrt2: ONE member of the rotated c = -2 pair. alpha -> -alpha (Suzuki's paired +-alpha problems) gives the other.
- OPERATOR level: the regular solutions are related by the Erdelyi-Kober fractional integral M(a, b, z) -> M(a, b + alpha/2, z) in z = y^2 (term identity verified, k = 0..6). This is NON-LOCAL. No local transform exists in general, since alpha/2 is non-integer.
- CHARGE level: the single member's exact Stirling coefficients B_(n+1)(1/2 + i(pi + P_X)/sqrt2) are NOT P_X-even. Their P_X-even part equals (1/2) sum over the doublet (n+1 = 2..10 checked).
  - This is cc's formula (item 299). The cylindrical charges are the alpha-even part, i.e. the doublet average.
  - So the WKB/continued-Beta convention (alpha-odd total derivatives integrate to 0) computes the doublet average, NOT the exact log D of one member. A NEW OBSERVATION about the ODE/IM dictionary for kappa-odd opers.

## (A) Exact point
A1. Sol 1: Suzuki(M = -2) vs the pillow. NO LOCAL OPERATOR RELATION, by three independent invariants (b1_a1_checks.log, a2_laplace.log).
  - (i) kappa-parity: under x = f(z) the kappa^1 part of the normal form is exactly kappa alpha f'^2 U(f) (machine-checked), and Darboux adds a kappa-independent term. The pillow's kappa-odd part is 0.
  - (ii) Local orders of Lambda dw^2: Suzuki {-5/2, 1, -5/2} against the pillow's {-1/2, -1, -5/2} (machine-checked), so no Moebius map.
  - (iii) Suzuki(M = -2) has NO regular singular point (x = 0 and x = infinity are both irregular). The pillow has regular singular points at w = 0 (exponent difference 1/2) and at w = 1 (difference sqrt(2) sqrt(-X), momentum-dependent). Algebraic point maps preserve regular singularity and send generic irrational exponent differences to such differences, so there is no algebraic change of variables (with any gauge).
  - All of (i)-(iii) are OPERATOR-level no-go statements in the stated classes. Kappa-DEPENDENT point maps and E <-> alpha role swaps (Levi-Civita type) are not covered.
A2. Non-local: the Laplace transform L4 of Suzuki(M = -2) times x^4 (a2_laplace.log).
  - Its only finite singular points are p = +-i sqrt(E), with integer exponents {-1, 0, 1, 2}.
  - The pillow's X-point difference sqrt(2) sqrt(-X) is not an integer, so the pillow is NOT a right factor of L4 and its symmetric cube is NOT L4 (up to gauge), for every identification sending the pillow's X-point to p = +-i sqrt(E).
  - Identifications sending it to p = infinity (irregular): NOT tested.
  - A Mellin-multiplier relation (the B2 mechanism) cannot apply directly: the pillow has a (1+u)^(-5/2) branch-type irregular point at u = -1, so it has no finite Mellin recurrence.
  - VERDICT: no operator relation found. NONE in the tested classes; open beyond them.
A3. Sol 2 at t = 1/3 (third order; k = -1/2, n = 2, c = -2). NOT A SYMMETRIC SQUARE (a3_wilczynski.log).
  - The Wilczynski invariant W = q - p'/2 = (27E e^(8z) - (27 l0^2 - 9 l1^2 + 9) e^(4z) + 2) e^(-6z)/27 is not 0.
  - Controls: the symmetric square of a generic Schrodinger operator, and of Suzuki(M = -2), both give W = 0.
  - W is a projective (weight-3) invariant, so in ANY coordinate (point transformation + gauge) L3 is not the symmetric square of a second-order operator.
  - Tensor products / Mellin conjugates relating it to Suzuki(M = -2) or the pillow: NOT tested.
## Summary: operator level vs charge level
- B2: OPERATOR (Mellin-recurrence level) relation found: tensor product of the rotated Coulomb pair, then a Gamma(nu + 1/2)/Gamma(nu - 1) multiplier. CHARGE level consistent (cc).
- B1: exact OPERATOR relation (energy shift = Erdelyi-Kober fractional integral) to one member. CHARGE level = the doublet average.
- A1/A3: OPERATOR-level no-go in the stated classes.
- A2: no-go for the natural Laplace identification.
- A (all): CHARGE-level isospectrality already known.
## Deviations and disclosures
- B2's sealed hypothesis "tensor product of the two rotated OSCILLATORS" was wrong as stated (oscillator form fails). It is right with the Coulomb (M = -1/2) form plus a Mellin multiplier, which was found during the stage, not sealed.
- b2_mellin.py's printed negative-control list still shows a nu-dependent "solution"; the filter is cosmetic. Every printed entry depends on nu or is trivial, so there is no constant solution.
- b2_mellin.py first runs: C was held fixed, then not substituted into the target. Both were fixed before the reported run; the final run's log is the stored one.
- No solver runs, so no memory caps were needed.
