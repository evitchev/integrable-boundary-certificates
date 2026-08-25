# Lab opener: the exact multiplicity stratification of the cylindrical curve

First computation of the Mathematica discovery lab (2026-08-23).  Status
of each statement is marked: LAB (hypothesis from the Mathematica
pipeline, hashed) or CERTIFIED (re-derived in the open stack).

## Question

For each cylindrical family and each odd spin sigma, the charges
commuting with I_3 form the kernel of the linear map ad_{I_3}, whose
dimension (the multiplicity) is 1 at generic points of the curve
s^2 = (N-25)(N-1) and jumps at special points.  Three such points had
been found by sampling (results/cyl_hidden_kdv_scan.md).  Where, exactly,
does the multiplicity jump -- all of them, at each spin?

## Pipeline (wolfram/lab/)

1. `cyl_strata_build.py` (Python, open stack): ad_{I_3} in CLASS
   coordinates (non-pivot monomials of the mixed sectors modulo total
   derivatives on both sides), so the nullity is the multiplicity.  By
   linearity I_3 = sum_k c_k(N,s) mono_k, so M(N,s) = sum_k c_k M_k(N)
   with M_k the residue map of one monomial, polynomial in N through
   O(N) traces; each M_k is built by exact evaluation at 7 rational N
   and Lagrange interpolation, verified at 3 held-out N (fail-closed).
   On the curve N(t) = (t^2-25)/(t^2-1), s(t) = -24t/(t^2-1); the common
   denominator is cleared: an integer polynomial matrix in t.
   Sizes: 24x12, 71x29, 185x68, 461x160 at weights 6, 8, 10, 12.
2. `cyl_strata.wls` (Wolfram Engine 15.0): generic rank r0 at several
   rational t; a full-rank r0 x r0 submatrix chosen at a generic point;
   its determinant D(t) is a maximal minor, so EVERY rank-drop point is a
   root of D (complete candidate set); D is factored over Q; for each
   irreducible factor f the rank of M over the number field Q[t]/(f) is
   computed by elimination with entries reduced mod f (inverses from the
   extended Euclidean algorithm -- no algebraic numbers; unit-tested);
   candidates where the rank does not drop are reported as rejected.
   t = +-1 are the parametrization's poles (N -> infinity), excluded;
   N = 1 is t = infinity, checked directly in Python.
3. `code/cyl_strata_points.py` (open stack): pins every positive finding.

## Results, spins 5 / 7 / 9 / 11  (generic multiplicity 1 throughout)

| point (N, s)       | t    | Sol 1    | Sol 2    | Sol 3    | mechanism                                  |
|--------------------|------|----------|----------|----------|--------------------------------------------|
| (-2, -9)           | 3    | 1,1,1,1  | 2,2,2,2  | 2,2,2,2  | decoupling (X (+) Y at every odd spin)     |
| (-25/2, -45/2)     | 5/3  | 1,1,1,1  | 2,1,1,1  | 2,1,1,1  | self-dual embedding + one extra mixed class at spin 5 |
| (-24, -35)         | 7/5  | 1,1,1,1  | 1,1,1,1  | 1,1,2,1  | NEW: pure multiplicity jump at spin 9 (Sol 3) |
| (-143/4, -189/4)   | 9/7  | 1,1,1,1  | 1,1,1,2  | 1,1,1,1  | NEW: pure multiplicity jump at spin 11 (Sol 2) |
| (28, +-9)          | +-1/3| 1,1,1,1  | 1,1,1,1  | 1,1,1,1  | embedding points: NO multiplicity jump     |
| N = 1 endpoint     | inf  | 1,1,1    | 1,1,1    | 1,1,1    | branch point: no jump (spins 5-9 checked)  |

(dims listed at spins 5, 7, 9, 11)

COMPLETENESS (LAB): at spins 5, 7, 9, 11 these are ALL the multiplicity jumps
of the three families on the curve (finite t != +-1), Solution 1 having
none at all.  The deck-conjugate points (t -> -t) do not jump: the jumps
are sheet-specific.  Calibration: the three sampled points of the scan
are reproduced exactly, with no others at spin 5.

A PATTERN (LAB HYPOTHESIS, not a claim): the sporadic single-spin jumps
sit at t = 5/3, 7/5, 9/7 -- consecutive odd numerators and denominators
t_k = (2k+3)/(2k+1), marching toward the parametrization pole t = 1
(N -> -infinity on that branch: N = -25/2, -24, -143/4, ...) -- at
spins 5, 9, 11, the first shared by Solutions 2 and 3 (their merger),
then alternating (Sol 3 at spin 9, Sol 2 at spin 11).  Prediction to
test at spin 13: a jump at t = 11/9, i.e. (N, s) = (-238/5, -297/5),
presumably for Solution 3.  Weight-14 matrices are the next lab build.

OUTCOME (2026-08-23, preregistered in results/cyl_spin13_preregistration.md
at commit f5d6521 before any spin-13 kernel existed; certificate
code/cyl_spin13_prediction.py, open stack, exact): NOT CONFIRMED.  All
three Solutions have spin-13 kernel dimension 1 at (-238/5,-297/5),
equal to the preregistered baseline (dimension 1 on every sheet at the
t = 2 probe).  The generic spin-13 multiplicity itself is NOT yet
established: that needs the symbolic generic rank of the weight-14
matrices over Q(t), i.e. the global stratification; dimension-one
probes make it plausible, not proved.
The controls behave exactly as fixed in advance: (-2,-9) decoupling
dim 2 for Solutions 2 and 3 (so weight 14 does see jumps), Solution 1
dim 1; the single-spin points (-25/2,-45/2), (-24,-35), (-143/4,-189/4)
are back to dim 1; the off-sequence points t = 5/4 and 7/6 sit at the
baseline.
Not a refutation: the sequence already skips spin 7, so t = 11/9 may
still appear at spin 15 or beyond -- or the three points may be a
coincidence.  The hypothesis stays LAB, downgraded from "prediction" to
"open"; the global weight-14 stratification (Mathematica) decides
whether spin 13 has any new sporadic point and where.
A reading of this outcome that keeps the t-sequence alive -- the
sequence locates the points exactly but need not fix their spins, the
spin being set by a resonance level -- is written up as a FRAME in
[private-note] (Section 4), together with the
search strategy it implies for the global weight-14 run.

CERTIFIED (cyl_strata_points.py, cyl_strata_points_s11.py): every kernel
dimension in the table;
and the mechanism at (-24, -35): the spin-9 kernel of Solution 3 has no
e^{sqrt2 X}-screened direction (nor sqrt3), no pure-X or pure-Y element,
and I_3 is not X-screened at this point at all -- a jump of the
multiplicity function with no screening structure, a third mechanism
beside embedding and decoupling.  The spin-11 point (-143/4, -189/4) of
Solution 2 has the same mechanism (certified).  (N = -24 is an integer; the density
there is (dX)^4 - 7(d2X)^2 - 4(dX)^2(dY.dY) + Y-base.)

## Lessons (lab doctrine additions)

* The free Wolfram Engine licence allows few concurrent kernels: a
  wolframscript run hard-killed from outside leaves an orphaned kernel
  that blocks activation ("license-related problem").  Time limits go
  INSIDE scripts (TimeConstrained); never kill from outside.
* Position/FirstPosition default to Heads -> True (the head List matched
  Except[0] at position 0); select by index instead.
* PolynomialMod takes two arguments; reduction modulo a univariate
  polynomial is PolynomialRemainder[p, f, t].
* Rank at an algebraic candidate via Root substitution is hopeless for
  high-degree spurious factors; elimination over Q[t]/(f) is instant.
