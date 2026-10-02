# SEAL_EXC2 -- Fable seat, 2026-10-02.  Level-1 excited-state oper for Solution 3 (fourth-order three-term ODE).
Commissioned by the lead ("GO, seal first"; class as proposed in EXC1; target = stored Sol-3 I_3 level-1 blocks at t = 2,
9/4; hold-out = Sol-3 I_5 level-1 block, uncomputed until the prediction hash is registered).  Written before any
computation of this stage.

## 0. Disclosure
Known to me: the Sol-3 I_3 level-1 singlet blocks at t = 2 and 9/4 (computed in EXC1, files d1_sol3_t*.json; I have looked
only at "irreducible quadratic" and the level-shift check, not at their coefficients).  The Sol-3 I_5 level-1 block has
never been computed by this seat.  Lead's information: Masoero-Raimondo 2312.01955 -- for twisted algebras FFH-type
trivial-monodromy conditions have no non-trivial solutions.
Hand derivations made before sealing:
(i) Vacuum, class-U normal form in xi = x^c (vartheta = xi d/dxi, b = p^2 = k/(k+1)), from (R3) of RED1 by xi -> rescaled:
        L_0 = (vartheta^2 - PX^2)(vartheta^2 - pi^2) + xi (vartheta + 1/2) - E xi^b.
(ii) TEMPLATE: Sol 1's level-1 oper of EXC1 in the same form (psi = e^(xi) phi):
        (vartheta^2 - lam~^2) + 2 xi (vartheta + 1/2 - x0) + E xi^(b_S) - 2 z xi/(xi - z)^2 - b_S xi/(xi - z),   b_S = 2/(M_S+1):
     a pure zeroth-order modification with a double pole, local exponents {-1, 2} at xi = z.
(iii) Weights: for Sol 3 xi ~ T^3, so z has weight 3 and the corrections first enter at relative order eps^2; I_1 and I_3 at
     level 1 depend only on the xi -> infinity constants of the correction (an "effective quartic" P_4), I_5 also on the
     first order in z/xi.

## 1. The class
        L = L_0 + xi [ R_2(xi) vartheta^2 + R_1(xi) vartheta + R_0(xi) ],     R_j(xi) = sum_{m=1}^{4-j} r_(j,m)/(xi - z)^m,
(no singular vartheta^3 term: that is the gauge phi -> (xi - z)^g phi; the factor xi keeps the exponents +-PX, +-pi at 0;
pole orders are the Fuchs bound).  Unknowns: r_(2,1), r_(2,2), r_(1,1..3), r_(0,1..4), z: TEN.
CONDITIONS at xi = z, for a chosen set of integer exponents e_1 < e_2 < e_3 < e_4 (sum 6 in this gauge):
   3 indicial conditions (fix r_(2,2), r_(1,3), r_(0,4));  6 no-logarithm conditions (one per pair of exponents), each
   required identically in E, i.e. coefficient by coefficient in E.
So 7 unknowns remain against at least 6 equations: a finite solution set is possible only if the E-dependence adds at
most one independent equation.  (Template count: 3 unknowns, 1 indicial, 1 no-log condition x 2 powers of E: finite.)
Exponent sets tested, in this order:  S_a = {-1, 1, 2, 4};  S_b = {-1, 0, 3, 4};  S_c = {-2, 1, 3, 4};  S_d = {-1, 0, 2, 5}.

## 2. Tests
T0 (control, known answer): the same code on the Sol-1 template class L = P_2 + 2 xi (vartheta + 1/2 - x0) + E xi^(b_S) + xi R_0,
   exponents {-1, 2}, must return exactly r_(0,2) = -2z, r_(0,1) = -b_S and EXC1's quadratic for z.
T1 (existence): for each exponent set, the solution set of the conditions at t = 2 (exact; Groebner), first at a rational
   momentum point, then symbolically if affordable.  Reported: dimension, number of solutions, and for a NEGATIVE which
   condition is inconsistent (which pair of exponents, which power of E).
T2 (target; only if T1 gives solutions): I_1 must be Delta + 1, and the I_3 values over the solutions must be the roots of
   the characteristic polynomial of the stored level-1 block, t = 2 and 9/4 (charges from wkb3 with the effective quartic).
T3 (hold-out; only if T2 passes): I_5 predicted (needs the first order in z/xi; engine extension validated on the Sol-1
   template first) -> PREDICTION_EXC2.json, hash to the lead BEFORE the Sol-3 I_5 level-1 block is computed.

## 3. Predictions
T0 passes: 90%.  T1: some sealed exponent set has solutions: 35% (Masoero-Raimondo's no-go argues against; the three-term
ODE is not an FFH connection, which is the only reason it is not lower).  T2 given T1: 60%.  T3 given T2: 85%.
If T1 is negative for all four sets: reported as a finding consistent with the twisted no-go, with the failing
condition(s) named; I will then state, labelled post-hoc, whether relaxing the class (two points; poles also in R_3;
higher pole orders) changes the count, without claiming a construction.

## 4. Controls
T0 above.  Tamper for T1/T2 (if positive): the middle term xi(vartheta + 1/2) replaced by xi(vartheta + 1) must destroy the
solutions or the I_3 match.  Exit codes: f_template.py 0 = T0 holds; f_sol3.py 0 = a sealed set has solutions, 3 = none
(clean negative), 1 = control misbehaves; g_compare.py 0 = T2 passes, 2 = fails.  Exact arithmetic; no nsimplify.
