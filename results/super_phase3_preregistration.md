# SUSY Phase 3 remainder: preregistration (2026-08-30)

Committed BEFORE the computations.  Engine: code/super_engine.py (N=1
super-OPE, one superfield, certified Phase 1); currents of
code/super_lambda.py: T_lam = 1/2 (dX)^2 + lam d2X - 1/2 psi dpsi,
G_lam = psi dX + 2 lam dpsi; I3(lam) = :T_lam T_lam: - 1/4 :G_lam dG_lam:
(point-split, kappa = -1/4 universal, certified super_lambda_phase3).
As charges oint G_lam = oint G and oint T_lam = oint T (the deformations
are total derivatives).  All statements are exact and range-limited
(doubled weight w2 <= 24, i.e. spin sigma <= 11); no all-spin claim.

## S1. Vacancy in the enlarged spaces at general lambda

For lam in {0 (control), 3/4, 1/2, 1, 2}: the ad_{I3(lam)} centralizer
dimensions per doubled weight w2 = 1..24 (classes mod d; fermion
bilinears make even-spin and half-integer-spin densities available).
EXPECTED (H_S1): the spectrum is lam-independent and equals the
free-point calibration -- dim 1 at sigma in {0, 1/2, 1, 3, 5, 7, 9, 11}
(oint dX, oint G, oint T and the odd tower), dim 0 at every other
half-integer and even spin in range.  Any lam at which a dimension
differs is a super jumping point: candidate, then re-checked at
neighbouring rational lam and by the direct bracket identities.

## S2. Supermultiplet pairing

The joint (ad_G, ad_{I3(lam)}) kernel at w2 = 12, 16, 20 has dim 1
(I5(lam), I7(lam), I9(lam) unique); direct brackets [I5,I7], [G,I5],
[G,I7], [I3,I5], [I3,I7], [I5,I9], [I7,I9] are total derivatives.
EXPECTED (H_S2): the tower is a tower of G-SINGLETS -- ad_{oint G}
annihilates every I_sigma(lam) -- and there are NO fermionic partners:
the half-integer sectors sigma + 1/2 (w2 = 13, 17, 21, and all odd
w2 >= 5) are vacant in ker ad_{I3(lam)}.  (The "manifestly SUSY tower
closes under ad_G" reading of the Phase-3 note is thereby refuted or
confirmed: singlets mean the closure is trivial.)

## S3. Super corank-character instance at lam = 3/4

At lam = 3/4 the screening locus (alpha^2 - 3/2 alpha - 1)(alpha^2 +
3/2 alpha - 1) has four rational roots: F- = {2, -1/2}, F+ = {1/2, -2}
(certified).  For w2 <= 16 compute the dimension of the joint kernel
{P class : Res(P V_alpha) in Im D_alpha for all alpha in the set} for
(a) the dual pair F+ = {1/2, -2}, (b) all four roots, and compare
spin by spin with the ad_{I3(3/4)} centralizer dimensions of S1.
EXPECTED (H_S3, the super analogue of the bosonic four-screening
characterization): (b) equals the centralizer spectrum at every w2 in
range; (a) is at least as large, equal at the odd-spin slots.
A strict inequality at some w2 for (b) is a finding (extra screened
classes that do not commute with I3, or vice versa).

## Outcome semantics

All three are fail-closed pins in code/super_phase3.py (registered
before its run); a failed pin is reported as a FINDING with the
observed values, not silently re-pinned.

## Pre-run addendum (2026-08-30, after a w2 <= 12 smoke run of the script)

Two adjustments, recorded before the registered run:
1. The screening-kernel helper omitted the twisted image of the
   constant (D_alpha 1 = alpha dX; the engine's weight-0 basis is
   empty), which wrongly declared oint G unscreened at w2 = 3.  Fixed
   (D_alpha 1 inserted by hand); oint G commutes with every fermionic
   screening, as the top-component argument requires.
2. S3 equality is pinned for w2 >= 3.  At w2 = 2 the exception is
   pinned explicitly: the momentum charge oint dX commutes with I3
   (dim 1) but is screened by no root -- it carries screening charge
   (Res dX(z) psi e^{alpha X}(w) = alpha psi e^{alpha X}).  This is the
   same trivial exception as the bosonic U(1) charge, not a finding.
Smoke-run observations (w2 <= 12): S1 and S2 as expected at all five
lambda; four-root kernel = centralizer for 3 <= w2 <= 12; the dual pair
alone has a strictly larger kernel (dims 1, 2, 1, 2, 4 at w2 = 7, 8, 9,
11, 12 -- including half-integer spins), so all four roots are needed.
3. H_S3(a) REFUTED by the smoke run (w2 <= 12), recorded before the
   registered run: the dual pair F+ alone has kernel dims 1,1,0,0,1,2,1,
   0,2,4 at w2 = 3..12 -- strictly larger than the centralizer at
   sigma = 3, 5.  Its commutant contains the N=1 super-Virasoro vertex algebra of
   (T_lam, G_lam), whose classes mod d number p(h) - p(h-1) from the NS
   vacuum character prod_{n>=2}(1-q^n)^{-1} prod_{r>=3/2}(1+q^r):
   1,1,0,0,1,2,1,0,2,4 at w2 = 3..12 -- EXACTLY the dual-pair kernel
   dims, at all three rational dual pairs tried, lam in {3/4, 4/3,
   5/12} (root ratios -4, -9, -9/4).  So the refutation is the expected
   structure: the pair's commutant is super-Virasoro (in range) and the
   four-root kernel selects the super-KdV charges inside it -- the
   super analogue of the bosonic four-screening characterization.
   The certificate pins (S3c, prediction) dual-pair kernel = NS vacuum
   class count for 3 <= w2 <= 16 at all three lam (w2 = 13..16 not yet
   computed when this was written).
