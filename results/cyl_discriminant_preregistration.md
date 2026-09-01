# Discriminant factorization at weights 10 and 12: preregistration (2026-08-31)

Committed BEFORE the registered run.  Prior evidence acknowledged: the
lab stratifications `wolfram/lab/cyl_strata_w10.json`, `_w12.json`
(status "LAB HYPOTHESIS -- certify in the open stack before claiming")
computed these loci; this preregisters their promotion to a
fail-closed certificate plus the beta^2/Kac reading.  Method (the
certified w14 pattern): generic rank r0 by rational probes; a
full-rank r0 x r0 submatrix at a generic point; its determinant D(t)
is a maximal minor, so EVERY rank-drop point of the weight-w ad_I3
matrix over Q-bar is a root of D (complete over-approximation); every
irreducible factor of D is then classified by exact rank -- rational
roots by direct rank, irrational factors over Q[t]/(f) with a GF(p)
screen.  This classifies the ENTIRE algebraic locus, not only rational
points: strictly stronger than the residue scans.

## Pins (H_F1, the locus over Q-bar)

    weight 10 (spin 9):  sheet 1: no jump at any point of Q-bar;
                         sheet 2: jump exactly {t = 3}, nullity 2;
                         sheet 3: jump exactly {t = 3, t = 7/5}, nullity 2 each.
    weight 12 (spin 11): sheet 1: no jump;
                         sheet 2: jump exactly {t = 3, t = 9/7}, nullity 2 each;
                         sheet 3: jump exactly {t = 3}.

Every irreducible factor of D(t) must be classified (jump, no-jump, or
pole t = +-1); an unclassified factor fails the certificate.

## Pins (H_F2, the beta^2 reading)

With x = beta^2 = (t-1)/(t+1): the jump roots map exactly to
x = 1/2 (t = 3, i.e. c_{1,2}); x = 1/6 at (w10, sheet 3) (c_{1,6}:
2q-3 = 9 at q = 6); x = 1/8 at (w12, sheet 2) (c_{1,8}:
(3q-2)/2 = 11 at q = 8).

## Pins (H_F3, the negative Kac slots -- the even-q selection, algebraic)

The Kac condition 2q-3 = 11 has the ODD solution q = 7; the
preregistered expectation is that sheet 3 at weight 12 does NOT jump
at t = 4/3 (beta^2 = 1/7) -- indeed at ANY point beyond t = 3 (H_F1).
Likewise no sheet jumps at beta^2 = 1/3 or 1/5 (t = 2, 3/2) at either
weight.  These are now ALGEBRAIC statements (all of Q-bar), not scan
silences.  A violation is a FINDING, reported, never re-pinned.

## Out of scope, recorded

The multiplicity structure of the Fitting ideal (orders of the factors
in the elementary divisors) is recorded if cheap but not pinned; the
identification of the discriminant factor with an explicit Shapovalov
determinant is the open question (mechanism note secs. 4b/5), not a
claim of this certificate.  Expensive (needs local Mathematica);
launch-attested runner record predeclared; artifacts sha-pinned.

## Domain clarification (2026-08-31, Codex round 15)

"Over Q-bar" throughout means over the parameter line MINUS the
parametrization poles t = +-1 (cleared-denominator factors there are
recorded and excluded as non-curve-points by the certificates); the
registry claims now carry this wording explicitly.
