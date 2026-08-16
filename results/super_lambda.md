# Phase 3: the background-charge family, and the naive tower identified

Script `code/super_lambda.py`, artifact `super_lambda.json`; exact
sympy in (alpha, kappa, lambda) on doubled-weight <= 9 sectors.  One
N=1 superfield with

    T_lam = 1/2 (dX)^2 + lam d2X - 1/2 psi dpsi,
    G_lam = psi dX + 2 lam dpsi,

the coefficient 2 forced by (and verified against)
Res(G_lam, G_lam) = 2 T_lam exactly, symbolic in lam.  Note
oint T_lam = oint T and oint G_lam = oint G as charges, so the whole
lambda-family lives in one free-field charge algebra; what moves is
the composite pencil P(kappa, lam) = :T_lam T_lam: + kappa
:G_lam dG_lam: and the ray selected in it.

## 1. The quantum coefficient -1/4 is universal

G-invariance (Res(G, P) in Im d) selects kappa(lam) = -1/4
independently of lam — consistent with the fixed CLO formula
I3 = int(T^2 - 1/4 G G') at general background charge.

## 2. The screening locus at general lambda

The alpha-conditions of I3(lam) factor exactly:

    (alpha^2 - 2 lam alpha - 1)(alpha^2 + 2 lam alpha - 1) = 0.

Each quadratic has root product -1 (the fermionic dual-pair invariant
a+ a- = -1, now established identically in lam on this pencil) and
root sums +-2 lam (the background charge).  At lam = 0 this
degenerates to (alpha^2 - 1)^2, the free-point locus alpha = +-1.

**Selection structure (corrected by audit round 3, and now asserted
by the script):** the single pair F- = alpha^2 - 2 lam alpha - 1
DIVIDES every free-kappa pencil condition, so one screening family is
kappa-blind at every lam — a first version of this file wrongly
claimed the opposite.  Setting kappa = -1/4 adjoins the
sign-reflected pair F+ = alpha^2 + 2 lam alpha - 1; it is the DUAL
S+- pair (CLO's psi e^{+-b Phi}) whose joint requirement selects the
super-KdV ray for lam != 0.  At lam = 0 the two pairs coincide, and
screening cannot select at all (the round-2 finding, now understood
as the degenerate case of pair-coincidence rather than an
exception).  This is certified branch-safely at lam = 3/4, where all
four roots are rational (F-: {2, -1/2}; F+: {1/2, -2}): exact
integer elimination confirms the F- roots screen at kappa = 0, 1 and
-1/4 alike, while the F+ roots screen only at kappa = -1/4 —
positive and negative controls.

## 3. The naive I3 ray is the lam* = +-sqrt(10)/4 member

The hypothesis of `super_kdv_calibration.md` section 4 is CONFIRMED
at the I3-ray level by direct solve (solution set asserted exact and
back-substituted): the naive-pencil ray (dX)^4 - 3 (dX)^2 psi dpsi
lies on the lambda-family,

    :T_lam T_lam: - 1/4 :G_lam dG_lam:  =  (1/4) I3_naive  mod Im d
    at  lam = +-sqrt(10)/4   (lam^2 = 5/8),

and at that lam the screening locus equals the naive quartic exactly:
alpha^4 - (2 + 4 lam^2) alpha^2 + 1 = alpha^4 - (9/2) alpha^2 + 1.
So the original Phase-2 computation had, unknowingly, computed the
lam^2 = 5/8 member of the standard family all along — which is why
its brackets vanished exactly and its spectral pattern matched.  The
"two towers" of the audit episode are one family seen through two
slices: the concatenation pencil happens to intersect the
lambda-family exactly at lam^2 = 5/8, the point-split pencil at
lam = 0.

Scope of the identification (audit round 3, finding 4): the solve
identifies the I3 RAY.  For the higher stored charges, note the ray
identity is a CHARGE identity — oint I3_naive = 4 oint I3(lam*) — so
the two centralizers are literally the same space, and the
multiplicity-one scans (through sigma = 11) pin each stored naive
I5, I7, I9 as the unique in-range centralizer ray at its spin.
Identifying those rays with the literature's lam*-hierarchy densities
additionally uses the literature's existence statement (a local
spin-s hierarchy charge commuting with I3(lam*)); no
coefficient-by-coefficient comparison against published deformed
densities has been performed, and none of this is an all-spin claim.

Convention-level dictionary (derived, not machine-certified): with
X = i Phi and alpha = -i b, the pair-sum invariant gives
Q = b + 1/b = 2 i lam (sign per the published currents
G = i Psi dPhi - i Q dPsi, T = -1/2 (dPhi)^2 + (Q/2) d2Phi
- 1/2 Psi dPsi; audit round 3 corrected an earlier sign slip here),
so c-hat = 1 + 2Q^2 = 1 - 8 lam^2 and c = (3/2) c-hat
= 3/2 - 12 lam^2.  The naive tower is then the
c = 3/2 - 12*(5/8) = -6 member (c-hat = -4, imaginary background
charge — a non-unitary point of the family).  The central-charge
formula itself is standard-convention arithmetic, not an
engine-verified statement (the residue engine sees no central
terms).

## 4. Scope

All statements are exact identities on weight <= 9 charge sectors
(ray selection, locus factorization, the lam* identification); no
all-spin or all-lambda centralizer theorem is claimed.  The
generic-localization caveat of the elimination applies to the locus
factorizations; the lam* identification (section 3) is an exact
membership solve, not an elimination, and needs no such caveat.

Next steps this opens: the super-paperclip scan (Phase 4) can now
fit screenings with lam as a free parameter per field — the correct
generalization of "find the quadric" — and the c(lam) dictionary can
be tested against super-Liouville spectra once vertex-operator
weights (double poles) are added to the screening module.
