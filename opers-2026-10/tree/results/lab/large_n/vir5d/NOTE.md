# VIR5d (Fable seat, 2026-10-01) -- the spins in nZ: where the operator has no charge, and what else is there

EXPLORATORY stage, labelled.  Seal `SEAL_VIR5d.md` a851070a... (21:11:49Z); addendum
`SEAL_VIR5d_addendum.md` c38c786d... (21:19:35Z, a post-hoc conjecture sealed before its hold-out).
Both hashes were sent to the lead before the corresponding computations.
Notation: K = spin index (spin 2K-1), k = fibre parameter 2(t-1)/(3-t), n = k + 4, nu = (2K-1)/n;
X = PX^2, Y = pi^2 (operator variables); "drop point" = nu = i a positive integer with 2i <= K-1.

## Verdict

**Z0. The operator's vanishing locus is the record's diagonal rank-drop locus (38(dw)-(dx)).**
- (a) PROVED, top layer: the un-normalised top is (-1)^K sum_a binom(nu,a) binom(nu,K-a) X^a Y^(K-a); all
  K+1 coefficients vanish iff nu = i is an integer and 2i <= K-1; the number vanishing at nu = i is
  min(K+1, 2(K-i)), the record's multiplicity law.  (The record's remaining empirical link, grade-0
  kernel = BL weights, is not proved here.)
- (b) PROVED, loss-1 layer, via the identity
  (D - d)(phi^nu D log psi_1) = phi^nu((1/2) D log psi_1 + D^2 log psi_1) - (nu/2) phi^(nu-1) D psi_1  (d = K-1),
  which turns F1 of VIR5c into the polynomial form

      loss-1 = (nu/6) [w^(K-1)] phi^(nu-1) { 1 + (K - 3/2 - nu) e1 w + (4 nu - 2K + 2) e2 w^2 }          (F1')

  (e1 = X + Y, e2 = XY).  Machine check `z0a_f1prime.py` (exit 0): F1' == F1 (K = 1..12, five fibres);
  F1' == Codex's Sol-3 loss-1 amplitudes (K = 2..40, five fibres, 0/4095); lattice scan K = 2..12,
  i = 1..K-1: top vanishes <=> loss-1 vanishes <=> 2i <= K-1, no exception.  F1' is three terms
  phi^(nu-1) x (1, w, w^2): the reason the record needs exactly three shift amplitudes at loss 1.
- (c) CHECKED, not proved, all layers (`z0c_engine.py`): the full un-normalised R_(2K) is identically
  zero at (t, K) = (23/9, 6), (13/7, 6), (29/11, 7), (19/9, 7), (9/7, 7), (21/11, 9) -- with the four known
  from VIR4/VIR5a, ten drop points, i = 1, 2, 3.  At the lattice NON-drop point t = 1/3, K = 4, i = 2 it is
  not zero and its top is X^2 Y^2 alone.
- (d) At the boundary case K = 2i+1 (k = 1/i) the Gamma ratio is what makes the charge vanish: the
  stripped operator at k = 1/2, K = 5 gives R != 0.

**Z1. At a drop point the certified Solution-3 charge is the limit of the operator's normalised charge**
(`z1_compare.py`, exit 0).  Exact rationals at k_0 + h: max|monic - certified| = 1.2e-8, 4.3e-9, 4.8e-9,
3.7e-7 at h = 1e-8 and exactly 1e-8 times that at h = 1e-16, for (5/3, spin 5), (11/5, 7), (17/7, 9),
(7/5, 9).  A limit demonstration, not an identity: the charge there is the k-derivative of the WKB
coefficient (un-normalised leading coefficient proportional to h).

**Z2. The second kernel class.**  CFT side (`z2_kernel.py`: VIR1 kernel of ad I_3 mod d, VIR2 two-boson
conversion, VIR3 vacuum rule; gate: the Sol-3 table charge lies in every computed kernel):

| t | k | K (spin) | i | K - (2i+1) | full kernel dimension | antisymmetric class |
|---|---|---|---|---|---|---|
| 5/3 | 1 | 3 (5) | 1 | 0 | 2 | 1 |
| 7/5 | 1/2 | 5 (9) | 2 | 0 | 2 | 1 |
| 11/5 | 3 | 4 (7) | 1 | 1 | 1 | 0 |
| 17/7 | 5 | 5 (9) | 1 | 2 | 1 | 0 |
| 13/7 | 3/2 | 6 (11) | 2 | 1 | 1 | 0 |
| 23/9 | 7 | 6 (11) | 1 | 3 | 1 | 0 |
| 9/7 | 1/3 | 7 (13) | 3 | 0 | 2 | 1 |

- Sealed prediction (i) "dimension 2 at all four points" (85%) is REFUTED.  The grade-0 rank drop
  becomes a full two-dimensional kernel only at the BOUNDARY drop points K = 2i+1, k = 1/i, i.e. on the
  record's jump-point tower t_i = (2i+3)/(2i+1).  At the other drop points nothing but the Sol-3
  charge exists; the operator simply has no un-normalised charge there and the answer is Z1.
- Sealed prediction (ii) holds: the kernel is invariant under PX <-> pi and contains exactly one
  antisymmetric class.
- The lead's handle: at t = 5/3 the certified Sol-2 charge lies in the kernel and
  **Sol-2 = (1/6) Sol-3 + (5/6) E** exactly (`z4_sol2_at_5_3.py`), E the antisymmetric class below.
- **H-A is REFUTED as sealed** (`z2_compare.py` exit 2).  The derivative A of the operator under the
  antisymmetric power deformation (1 - l0^2/T^2)^(1+e)(1 - l1^2/T^2)^(1-e) has the right top (X-Y)^K and,
  at 7/5, the right loss-1 layer, but the deeper layers are wrong (at 5/3 A has no lower layer at all;
  the certified class has losses 1 and 2).  The symmetric and one-sided deformations are not in the
  kernel either.  So the second class is not a first-order response of this operator to a power
  deformation.

**Addendum, conjecture P-W (read off at i = 1, 2, then sealed).**  At t_i the antisymmetric class has
vacuum eigenvalue

      E_i = prod_{j=-i}^{i} (PX + pi - j p)(PX - pi - j p),        p^2 = 1/(i+1),

a product of two odd polynomials of degree 2i+1 in u = PX + pi and v = PX - pi with equally spaced
roots.  At t_i the VIR2 Virasoro pair has c_U = c_W = c_(1,i+1) (-2, -7, -25/2).
- T-a (t = 7/5): all 18 coefficients agree (4 were used to read the formula off; the seal said "15",
  a miscount).  T-c: no antisymmetric class at 11/5, 17/7, 13/7, 23/9.
- **T-b HOLD-OUT (i = 3, t = 9/7, K = 7, weight 14, nothing seen before the addendum seal): PASS.**
  Kernel dimension 2, one antisymmetric class, and it equals prod_{j=-3}^{3}(PX + pi - j/2)(PX - pi - j/2)
  in all 32 coefficients (`z3_pw.py`, exit 0; the weight-14 kernel took 14 min).  Gate
  (`z5_gate_K6_K7.py`): the certified Sol-3 spin-13 table charge lies in that kernel, and the spin-11
  table charge is the one-dimensional kernel at 13/7 and 23/9.
- Tamper p^2 + 1/1000: fails.
READING (not tested at the level of densities; the centre of the root string was not checked against
any representation theory; no algebra is claimed): the root spacing p in u is
1/sqrt(2(i+1)) in the momentum of U = (X + phi)/sqrt2, which is alpha_-/2 of the (1, i+1) Coulomb gas --
the spacing of the zero-mode eigenvalue of the singlet current of weight 2i+1.  The class looks like
the product W_U W_W of the two singlet currents (weight 2(2i+1) = 2K).  In operator variables
E_i = const x prod_j (l0 + l1 - j n)(l0 - l1 - j n): it vanishes exactly when two of the four dynamical
exponents differ by j n, |j| <= i, and at k = 1/i the two potential terms x^n and x^(n(i+1)) are
commensurate, so those are the Frobenius resonances of the equation.

## Answer to the question "what supplies the spin-in-nZ charges"

1. Nothing at k_0 itself: by Z0 the operator's WKB coefficient vanishes there, in every layer checked.
2. The certified Solution-3 charge is the limit of the family (Z1) -- the same operator family, one
   order in k - k_0; no other operator is needed for it.
3. A genuinely new class exists only on the jump-point tower K = 2i+1, k = 1/i.  It is odd under the
   Weyl swap, its eigenvalue is the product E_i, and at t = 5/3 the other solution's charge is a
   combination of it and the Sol-3 charge.  It is not produced by the operator's power deformations.

## Errors, deviations, disclosures

- Sealed (i) wrong (see above).  The seal's "15 coefficients" at 7/5 is 18.
- `z2_compare.py tamper` prints the string "H-A holds..." for exit 0; in tamper mode exit 0 only means
  the tampered class no longer contains the table charges (all gates flip to False, as they should).
  The wording is wrong, the control fired; log kept.
- A dead `nsimplify(...) if False else` fragment was removed from `z2_deformed.py` before it ran;
  `mellin_dual.py` is `mellin_wkb.py` with coefficients in Q[l0, l1, e]/(e^2) (diff: GENS, trunc,
  el_mul, el_scale, deformed_symbol); its e^0 part reproduces the old engine at k = 2 (R1 = 0 for the
  null deformation).
- The dual engine was not independently validated at order e beyond the top layer (the top of A is
  (X-Y)^K, as the classical formula requires).  Since H-A fails, an error there could only matter if it
  hid a pass; the certified antisymmetric class differs from A already at loss 1 at t = 5/3.
- Not done: Z1 as an exact L'Hopital identity; any density-level identification of the second class;
  Solutions 1, 2.

## Files
`SEAL_VIR5d.md`, `SEAL_VIR5d_addendum.md` (+ .sha256, .time), `z0a_f1prime.py`, `z0c_engine.py`,
`z1_point.py`, `z1_compare.py`, `z2_kernel.py`, `z2_deformed.py`, `z2_compare.py`, `z3_pw.py`,
`z4_sol2_at_5_3.py`, `mellin_wkb.py` (60795426..., unchanged), `mellin_dual.py`, `f1_lib.py`,
`vir_lib.py`, `vir2_lib.py`, `cft_lib.py` (byte-identical copies from VIR1-3), JSON outputs, logs,
`SHA256SUMS`.
