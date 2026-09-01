# Discriminant over Q-bar at weight 14: preregistration (2026-08-31, before the run)

The third rung of the discriminant program (after the certified
weights 10, 12 -- cyl_discriminant_w1012).  Same certified method: one
full-rank maximal minor D(t) as a complete over-approximation of the
Q-bar rank-drop locus of the weight-14 ad_I3 matrices (1082 x 365,
tracked with sha-pinned sidecars); every irreducible factor decided by
exact rank (GF(p) unit-pivot screen, exact Q[t]/(f) fallback).

## Pins (H_W14)

Spin 13 (weight 14), over all of Q-bar:
    sheet 1: no jump anywhere;
    sheet 2: jump exactly {t = 3}, nullity 2;
    sheet 3: jump exactly {t = 3, t = 9/7}, nullity 2 each.
beta^2 images: {1/2} on sheet 2; {1/2, 1/8} on sheet 3 -- the sheet-3
factor at t_3 is (1 - 8x) (c_{1,8}: 2q-3 = 13 at q = 8).
Every factor classified; poles only t = +-1.

## The asymmetry prediction

(3q-2)/2 = 13 has NO integer solution (q = 28/3), so the Phi-reading
predicts sheet 2 is silent at spin 13 at every point beyond the
decoupling c_{1,2} -- certified rationally already (the spin-13 record),
now pinned over Q-bar.  A sheet-2 jump anywhere else, or any sheet-3
point beyond {c_{1,2}, c_{1,8}}, is a FINDING, reported, never
re-pinned.  Consequence if green: the discriminant table spans weights
10-14 with factors (1-2x), (1-6x), (1-8x) twice -- (1-8x) appearing at
BOTH its Kac slots, spin 11 on sheet 2 (q=8 even, (3q-2)/2 = 11) and
spin 13 on sheet 3 (2q-3 = 13).

## Domain clarification (2026-08-31, Codex round 15)

"Over Q-bar" throughout means over the parameter line MINUS the
parametrization poles t = +-1 (cleared-denominator factors there are
recorded and excluded as non-curve-points by the certificates); the
registry claims now carry this wording explicitly.
