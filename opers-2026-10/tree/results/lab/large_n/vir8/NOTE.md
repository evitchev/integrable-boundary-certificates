# VIR8 (Fable seat, 2026-10-01) -- direct WKB of the three-term ODEs of Solutions 3 and 2

Seal `SEAL_VIR8.md` 53429566... (22:11:03Z), sent to the lead before any run on Sols 2, 3.  Nothing here is blind.

## Verdict

**(a) The class-U equivalence holds for deg P > 2: PASS (sealed 85%).**  Exact, zero parameters:

    Sol 3 (4th order):  [ (T^2-l0^2)(T^2-l1^2) - x^(c/2) T x^(c/2) - E x^n ] phi = 0,          n = k+4
    Sol 2 (3rd order):  [ T (T^2-l1^2) - x^(c/2) (T^2-l0^2) x^(c/2) - E x^n ] phi = 0,         n = 2k+3
    T = s - theta,  M = 1/k,  c = n(1+M) = n/p^2,  l1 = n pi/p,  l0^2 : l1^2 = PX^2/a_Y : pi^2/a_X.

At k = 2 and k = 10/3, WKB orders 2..8 (spins 1..7), the charges of these ordinary differential equations are
proportional to those of the Gamma-symbol operators (items 298, 301/SOL2F): every charge is cB x B(A0, B0+1) with the
SAME Beta function as the Gamma form, the coefficient cG of the second master integral vanishes identically, and
monic(cB) == monic(R_Gamma) (`a_equiv.py`, four runs, exit 0).  Tamper (middle term x^c L(T) instead of the symmetric
x^(c/2) L x^(c/2)): disagrees already at spin 1 and produces nonzero even-spin charges (exit 0 = fires).  With VIR7's second-order case, the three certified
solutions are ordinary differential equations of orders 4, 3, 2 with a three-term potential, at every t.

**(b) My sealed answer is REFUTED (b1-b4 predicted cB != 0 at 70-75%).**  At the nZ points the three-term ODEs ALSO have
cB = cG = 0 exactly (`b_nz.py`, exit 2): Sol-3 ODE at t = 11/5 spin 7, t = 5/3 spin 5, t = 7/5 spin 9; Sol-2 ODE at t = 5/3
spin 5; and already in (a), Sol-2 at t = 2 spin 7.  The ODE form does NOT supply a charge that the Gamma form lacks.
What is actually going on (algebra + `b2_limit.py`, post-hoc, labelled):
- At every nZ point A0 = -nu k is a negative integer, so the common Beta function B(A0, B0+1) has a POLE.  The WKB
  coefficient is J = B x (polynomial): (pole) x (zero), finite along the family.  "No charge" in VIR4-VIR5d and in the
  no-charge cells of VIR7 means: the polynomial normalised by that Beta function vanishes.  Nothing is missing from J.
- Directional limits, exact rationals at k_0 + 1e-8 and 1e-16: the monic charge of each ODE converges (linearly in h) to
  the certified charge OF ITS OWN SOLUTION -- Sol-3 ODE -> Sol-3 table at 11/5 (spin 7), 5/3 (spin 5), 7/5 (spin 9);
  Sol-2 ODE -> Sol-2 table at 5/3 (spin 5).  Cross-checks (Sol-2 ODE vs Sol-3 table and vice versa at 5/3) do not converge.
- So the answer to "which one": the limit depends on the FAMILY along which the point is approached.  At t = 5/3, where
  Sols 2 and 3 share one quintic Gamma operator, their two ODE families give two different spin-5 limits, S and
  S/6 + 5E_1/6; the swapped split (P = T(T^2-l0^2), L = T^2-l1^2) gives -S/4 + 5E_1/4 (`c_splits.py`).  The two-dimensional
  kernel there is the span of the directional limits of the families meeting at that operator.
- At t = 7/5 the Sol-3 family gives S; no family is known whose limit is E_2.  Open.

**Correction to my VIR7 note.**  I wrote there that "the nZ gaps are a property of the Gamma form, not of the three-term
ODE".  That was inferred from the second-order case, where the comparison was with the SCHRODINGER form, whose Beta
function has half the argument (no pole at odd A0).  For the three-term forms in the same variable the gaps are identical
to the Gamma form's.  The correct statement: a gap is a pole of the Beta prefactor of the chosen form.

## Engine (`wkb3.py`)
Leading curve tau^dP - y^c tau^dL = y^n, uniformised by u = y^c tau^(dL-dP): tau_0 = u^(1/sigma)(1-u)^(-c/(n sigma)),
y^n = tau_0^dP (1-u), u in (0,1); theta = [n sigma (1-u)/w] u d/du, w = dP - dL u.  All WKB orders are
tau_0^(integer) x Laurent polynomial in w; J_i = (1/(n sigma)) int u^(A0-1)(1-u)^B0 w W_i(w) du is reduced exactly to the
masters M_0 = B(A0, B0+1), M_(-1) by a three-term recurrence.  Validated on (dP, dL) = (2, 1) against the Gamma form
(even and odd spins) before the seal.  Order 8 takes ~5 s; order 10 ~40 s with 16-digit rationals.

## Errors and disclosures
- Sealed (b) wrong, as above.  The seal's own remark already contained the pole; I did not draw the consequence that the
  polynomial must then vanish in the three-term form too if the two J's are proportional.
- `c_splits.py` v0 solved for (al, be) from two dependent coefficients and crashed; fixed (PX^2 pi^4 instead of PX^4 pi^2).
- The limits are demonstrations at two step sizes, not exact identities.
- Not done: excited states; a three-term form for anything beyond Sols 1-3; t = 7/5's second class.

## Exit codes
v0_validate 0; a_equiv real x4: 0; a_equiv tamper x2: 0 (fires); b_nz 2 (not as sealed); b2_limit 0; c_splits (report only).

## Files
`SEAL_VIR8.md` (+ .sha256, .time), `wkb3.py`, `v0_validate.py`, `a_equiv.py`, `b_nz.py`, `b_core.py`, `b2_limit.py`,
`c_splits.py`, `mellin_gen.py`, `mellin_wkb.py` (byte-identical copies), logs, `SHA256SUMS`.
