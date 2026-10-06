# SEAL_EXC3 -- Fable (secondary) seat, 2026-10-04.  Level-one excited-state oper for Solution 2 (third-order three-term ODE).
Commissioned by the lead at the operator's request (brief EXC3_brief_2026-10-04.md, sha256 3efdb4c7...).  Written before any
computation of this stage: no data block, no engine run, no Groebner basis.  Hand derivations below are to be machine-checked
(tests T0).

## 0. Disclosure
Known to me: EXC1 (Sol 1: one apparent singularity per level in xi = x^(M+1)), EXC2 (Sol 3: level 1 = TWO apparent
singularities at xi = +-z, exponents {-1,1,2,4}, selected by invariance under the formal adjoint composed with xi -> -xi;
ONE point = half a level), the Sol 2 operator (item 303/307), its verification at t = 2, 9/4 and (my review, F29) at
t = 1/3 and other fibres.  I have NOT computed any Solution-2 excited-state block (I_3 or I_5); the lead holds the I_5
blocks as hold-out.  Code reused: frob.py (local Frobenius analysis, order-general), wkb3x.py (three-term WKB, first
order in x^(-c)), vir_lib.py / d1_data.py (data side, as in EXC1/EXC2).

## 1. Normal form (hand derivation, to be checked at T0)
Item 307's Sol 2 equation [T(T^2 - l1^2) - x^(c/2)(T^2 - l0^2)x^(c/2) - E x^n] phi = 0, T = s - theta, n = 2k+3, s = k+1,
b = p^2 = k/(k+1), c = n/b, l0^2 = n^2 P_X^2/(k p^2), l1^2 = n^2 pi^2/p^2.  With the gauge phi = x^s phi~, T = -c vartheta
(vartheta = xi d/dxi) and xi = x^c / c, dividing by -c^3:
    L_0 = vartheta (vartheta^2 - a1^2) + xi ((vartheta + 1/2)^2 - a0^2) - E' xi^b,
    a1 = l1/c = p pi,   a0 = l0/c = P_X/sqrt(k+1),   i.e.  a1^2 = b pi^2,  a0^2 = (1-b) P_X^2,   E' = c^(b-3) E (sign free).
(Check: x^(c/2) L(T) x^(c/2) = x^c L(T - c/2) and L(T - c/2) = (T - c/2)^2 - l0^2 = c^2 ((vartheta + 1/2)^2 - a0^2).)
Exponents at xi = 0: {0, +-a1} -- only the Y momentum; X enters through the middle term.  Compare Sol 3 (EXC2):
(vartheta^2 - a0^2)(vartheta^2 - a1^2) + xi (vartheta + 1/2) - E xi^b with both momenta as exponents at 0.
Dictionary to the engine (wkb3x: P(T) - x^(c/2) L(T) x^(c/2) - E x^n, c = n(1+M), M = 1/k, dP = 3, dL = 2):
    P_xi = vartheta^3 + q2 vartheta^2 + q1 vartheta + q0  ->  P_T = T^3 - c q2 T^2 + c^2 q1 T - c^3 q0;
    (1/xi) (q12 vartheta^2 + q11 vartheta + q10)  ->  x^(-c) (B2 T^2 + B1 T + B0),  B2 = -c^2 q12, B1 = c^3 q11, B0 = -c^4 q10;
    (1/xi^2) Q2(vartheta)  ->  x^(-2c) (-c^5 Q2(-T/c)).
To be checked at T0 by substitution, and by the vacuum charges (q's = 0) against the certified Sol 2 eigenvalues at t = 2.

## 2. The symmetry (hand derivation; the lead's brief expects none -- to be machine-checked at T0)
Formal adjoint for the pairing int f g dxi/xi: vartheta^+ = -vartheta, f(vartheta) xi = xi f(vartheta + 1).  Then
    P(vartheta)^+ = P(-vartheta) = -P(vartheta)   (P odd),
    [xi ((vartheta+1/2)^2 - a0^2)]^+ = ((-vartheta+1/2)^2 - a0^2) xi = xi ((vartheta + 1/2)^2 - a0^2)   (the complete square),
    (E' xi^b)^+ = E' xi^b.
So  L_0^+ = -P + xi L - E' xi^b,  and composing with xi -> -xi:  L_0^+(-xi) = -[P + xi L + (-1)^b E' xi^b] = -L_0 |_(E' -> -(-1)^b E').
The odd P only produces an overall sign; since the no-logarithm conditions are imposed identically in E', the map
    sigma: L  ->  -(L^+)(-xi)
is an anti-involution of the class, exactly as in EXC2, PROVIDED the middle term is the complete square above.  (A tamper
that breaks it: xi ((vartheta + 1)^2 - a0^2) -- the ordering tamper the brief asks for.)  On a deformation xi R_j(xi) vartheta^j
it acts as in EXC2: V_2 -> V~_2, V_1 -> 2 vartheta V~_2 - V~_1, V_0 -> vartheta^2 V~_2 - vartheta V~_1 + V~_0, V~_j(xi) = V_j(-xi),
and sends a point at z to a point at -z with exponents e -> 2 - e (the adjoint of an order-3 operator reflects exponents
about (3-1)/2 = 1).  Exponent sets invariant under sigma are those symmetric about 1.
PREDICTION S (90%): sigma is a symmetry of the Sol 2 class; sigma-invariant solutions have z2 = -z1 and the exponent set
{1-d, 1, 1+d}.  If the machine check refutes S, the stage proceeds without a symmetry and the selection is by the data
(reported as such).

## 3. The class
    L = L_0 + xi sum_points [ R_1^(k)(xi) vartheta + R_0^(k)(xi) ],   R_j = sum_(m=1)^(3-j) r_jm/(xi - z_k)^m
(no singular vartheta^2 term: that is the gauge (xi - z)^g; pole orders = Fuchs bound; the factor xi keeps the exponents
{0, +-a1} at 0).  Per point: r11, r12, r01, r02, r03, z: SIX unknowns.  Indicial polynomial at xi = z:
e(e-1)(e-2) + (r12/z) e + r03/z^2 = (e - e1)(e - e2)(e - e3), so the exponent sum is 3 and the indicial conditions fix
r12, r03.  Remaining per point: r11, r01, r02, z (FOUR) against the no-logarithm conditions at the resonances (one per
pair of exponents, each identically in E and in the free Frobenius constants).
Exponent sets, in this order: S_a = {-1, 1, 3} (the order-3 member of the pattern {-1, 1, .., m-2, m} that gave Sol 1's
{-1, 2} and Sol 3's {-1,1,2,4}; sigma-symmetric), S_b = {-2, 1, 4} (sigma-symmetric), S_c = {-1, 0, 4}, S_d = {-2, 2, 3}
(not symmetric).  Fibres t = 2 (k = 2, b = 2/3, n = 7, c = 21/2) and t = 9/4 (k = 10/3, b = 10/13, n = 29/3, c = 377/30);
momentum points (a0, a1) rational, e.g. (1/2, 1/3), (3/7, 2/9).
Two-point class: two points of type S_a (12 unknowns, 4 indicial, 8 remaining against 2 x conditions), and its
sigma-invariant subfamily z2 = -z1 (after the symmetry relations the count halves).

## 4. Level bookkeeping and selection (sealed hypotheses)
H1 (one point = ?): I_1 on a one-point solution is Delta + L1 with L1 forced by the E-coefficient condition (as r21 = -b was
in EXC2).  Priors: L1 = 1 (one point per level, as Sol 1): 45%;  L1 = 1/2 (as Sol 3): 45%;  other: 10%.
H2 (selection): if L1 = 1/2, the level-one singlet states are the sigma-invariant two-point solutions (z2 = -z1), as for
Sol 3, and on them the even-spin charges J_3, J_5 vanish: 60% given L1 = 1/2.  If L1 = 1, the level-one states are
one-point solutions, with two solutions per momentum point matching the data quadratic: 55% given L1 = 1.
H3 (even-spin charges): the VACUUM Sol 2 operator has vanishing odd-order WKB charges (even spins): 70% (to be computed
at T0; if they do not vanish, "even-spin charges vanish" cannot be a selection rule).
What I_1 and I_3 depend on: I_1 and I_3 depend on the effective cubic (the xi -> infinity constants r_j1) AND, because xi
has weight ONE for Sol 2 (xi ~ T^(dP-dL) = T), on the first order in 1/xi (weights of x^(-c) T^j relative to T^3:
4 - j).  I_5 (order 6) needs the second order in 1/xi (x^(-2c) T^j at relative order 5 - j for j = 1, 2) and the squares of
the first-order terms.  So, unlike EXC2, the engine must be extended to second order in U = x^(-c) BEFORE the hold-out.

## 5. Tests
T0 (controls, before any excited run):
  (a) normal form: the engine with q = 0 reproduces the certified Sol 2 vacuum eigenvalues (I_1, I_3, I_5; spins 1, 3, 5) at
      t = 2 from (a0, a1) above, monic;  the odd orders 3, 5 of the vacuum are reported (H3).
  (b) sigma: machine check that -(L^+)(-xi) = L for L_0 (symbolic) and that the ordering tamper breaks it.
  (c) frob.py regression: the EXC2 one-point Sol 3 system (exponents {-1,1,2,4}, t = 2, (1/2,1/3)) is reproduced (4 solutions,
      r21 = -b).
  (d) engine extension to second order in U: on the Sol 1 template at M_S = 5 the full EXC1 result R_4/c4 = -7 al^4/400 +
      3 al^2 lam^2/10 + 9 al^2/2 + 168 al z + lam^4 + 18 lam^2 + 240 z^2 + 337/5 must be reproduced INCLUDING the z^2 term
      (x0_validate reproduced everything but 240 z^2 at first order); translation covariance T -> T + d of all charges
      with the second-order term present.
T1 (existence): per exponent set, the no-log system at t = 2, (a0, a1) = (1/2, 1/3): dimension, number of solutions, and
    for a negative the inconsistent condition.  Then the two-point system and its sigma-invariant subfamily.
T2 (I_1 level; I_3 target): I_1 on the solutions -> L1 (H1).  For the class selected by H2: the characteristic polynomial
    of I_3 over the selected solutions vs the stored Sol 2 level-one singlet I_3 block (computed by d1_data for Sol 2 at
    t = 2, 9/4; the I_5 blocks are NOT computed).  Post hoc at (1/2,1/3), t = 2; then at other momentum points and at
    t = 9/4; then trace and determinant as polynomials in (P^2, Q^2) by interpolation with held-back points (as EXC2).
T3 (hold-out): I_5 trace and determinant polynomials at t = 2 and 9/4 -> PREDICTION_EXC3.json, sha256 to the lead BEFORE
    the Sol 2 I_5 level-one blocks exist; the lead compares first.
Priors: T0 all pass 80%; T1 S_a has solutions 70%; T2 I_3 match for the selected class 50% (H1 x H2); T3 given T2: 85%.

## 6. Controls and conventions
Tampers (after the hash guard): middle term xi((vartheta + 1)^2 - a0^2) (ordering) -- must destroy sigma, and the I_3 match
if one is found; exponent-set tamper ({-1,1,3} -> {-1,0,4}) on the two-point class; data tamper (+1/1000 on one matrix
entry) and cross-fibre comparison for the hold-out; sign flip of the first-order term in the I_5 map.
Exact arithmetic; Gamma atoms reduced by the functional equation before any comparison (gamma-atoms lesson); no nsimplify;
runs capped (ulimit -v 16 GB, timeout); kill by saved PID only; full logs; Singular for the eliminations.  Exit codes
reflect computed checks only.  Repo read-only; workspace ~/fable-work/exc3.
