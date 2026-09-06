# M2: the T-primariness split at t_2 -- preregistration (2026-09-01, before the run)

Tool: M1's mixed pole extraction (code/mixed_pole.py, validated).
Point: t_2 = 7/5, (N, s) = (-24, -35), Solution 3, weight 10 (spin 9),
where the certified kernel has dimension 2 (the (1 - 6x) discriminant
factor).  T = T_total = (dX)^2/2 + (dY.dY)/2, c = N + 1.

Definition.  A class [v] (v in the weight-10 kernel span K, modulo the
image D of d from weight 9) is T-PRIMARY if some representative v + dq
has vanishing poles k >= 3 in T(z)(v + dq)(w) (poles 2 and 1 are the
weight and the derivative, automatic).  The number of primary classes is
dim(ker P on K + D) - dim(ker P on D), P = the stacked pole maps
k = 3..12, exact linear algebra over Q.

EXPECTED (H_M2, the Phi_{1,3}-shadow reading): exactly ONE primary class
direction in the 2-dim kernel at t_2 -- the shadow of a singular vector
is a primary field -- and ZERO primary classes in the generic 1-dim
kernel at control points (t = 7/3, 11/5: KdV-type densities are never
primary; also the deck partner t = -7/5 as a control: generic there).
Refutation criteria: 0 primary classes at t_2 (the extra class is not a
primary shadow -- the mechanism is not a null-field image in this
sense), or 2 (both directions primary: a bigger algebra, not a
quotient), or a nonzero count at a generic control (primariness would
then be unrelated to the jump).  Exploratory, unpinned: the same count
under T_X alone and T_Y alone.  Any deviation is a FINDING, reported,
never re-pinned.

## Outcome (2026-09-01, lab run, before the registered certificate)

H_M2 REFUTED: at t_2 the 2-dim spin-9 kernel contains ZERO T_total-
primary class directions (kernel-on-(K+D) 0, on D 0); the generic
controls t = -7/5, 7/3, 11/5 give 0 as predicted.  Exploratory: under
T_X alone and T_Y alone the count is also 0 (with 20 resp. 14 primary
total-derivative directions in D, none in the kernel classes).  Reading:
the extra class is not a primary composite field -- consistent with
obstacle (i) of the mode-action note (h_{1,3} is fractional, so a
Phi_{1,3} null vector has no polynomial-field image that could be
primary).  Route 3 of that note is closed; the registered certificate
pins the observed counts in outcome form.
