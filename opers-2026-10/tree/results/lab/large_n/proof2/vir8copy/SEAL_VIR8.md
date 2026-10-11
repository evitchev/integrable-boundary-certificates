# SEAL_VIR8 -- Fable seat, 2026-10-01.  Direct WKB of the three-term ODEs of class U (Solutions 3 and 2).
Commissioned by the lead: (a) at a fibre where the Gamma form is exact the ODE form must reproduce the same charges;
(b) at an nZ point, does the ODE form supply the "missing" charge, and which one.  Written before any run on Sols 2, 3.

## 0. Disclosure
Done before sealing, ODE side only: (i) the exact form of the equivalence, now with all constants (same variable x):
        Gamma form  P(T) prod_(roots r of L) G_a(T - r) psi = x^n (x^(nM) + E) psi          [E -> -E done, T = s - theta]
   <=>  three-term  [ P(T) - x^(c/2) L(T) x^(c/2) - E x^n ] phi = 0,   c = n(1+M),  a = 1/M,  n = deg P + a deg L,
   by psi-hat = q phi-hat, q(nu) = sigma^(-a d nu/n) prod_r Gamma((s - nu - n - r)/sigma + (a+1)/2), sigma = nM, d = deg L.
(ii) a new engine `wkb3.py` for the three-term equation: uniformisation u = y^c tau^(dL-dP), tau_0 = u^(1/sigma)(1-u)^(-c/(n sigma)),
   everything tau_0^(integer) x Laurent polynomial in w = dP - dL u; charges reduced EXACTLY to two masters,
   M_0 = B(A0, B0+1) (the Gamma form's Beta function) and M_(-1) = int u^(A0-1)(1-u)^B0/w du:  J_i = cB M_0 + cG M_(-1).
(iii) validation on (dP, dL) = (2, 1) only (`v0_validate.py`, exit 0): for (T^2 - l1^2) - x^(c/2)(T - l0)x^(c/2) - E x^n at
   a = 5/8, 2/3, orders 2..6 (even and odd spins): cG = 0 identically and monic(cB) == the Gamma form's polynomial.
NOT done: any run with (dP, dL) = (4, 1) or (3, 2); any run at an nZ point.
Remark fixed before sealing (algebra): at every nZ point nu = spin/n integer, A0 = -nu k is a negative integer, so the common
Beta function has a POLE there.  The Gamma form's J = B x R is then (pole) x (zero); "the operator has no charge" in
VIR4-VIR5d meant R = 0 in the normalisation that divides by B.

## 1. Objects (zero parameters; l1 = n pi/p, l0^2 : l1^2 = PX^2/a_Y : pi^2/a_X)
Sol 3:  [ (T^2-l0^2)(T^2-l1^2) - x^(c/2) T x^(c/2) - E x^n ] phi = 0,       n = k+4,   M = 1/k, c = n(k+1)/k   (4th order)
Sol 2:  [ T (T^2-l1^2) - x^(c/2) (T^2-l0^2) x^(c/2) - E x^n ] phi = 0,      n = 2k+3,  M = 1/k, c = n(k+1)/k   (3rd order)

## 2. Tests
(a) EQUIVALENCE at generic fibres, exact: Sol 3 at k = 2, 10/3 and Sol 2 at k = 2, 10/3, WKB orders 2..8 (spins 1..7):
    PASS = at every order cG = 0 (or cG proportional to cB) and monic(cB) == monic(R) of the Gamma form
    (mellin_wkb for Sol 3; mellin_gen with strings of length k at +-l0 for Sol 2), odd orders vanishing in both.
    PREDICTION: PASS, 85%.
(b) nZ POINTS, exact arithmetic at k_0 itself: cB, cG of the three-term form at the order where the Gamma form's R = 0:
    b1  Sol 3 form, t = 11/5 (k = 3), spin 7 (kernel 1-dim):   cB != 0, cG = 0, monic(cB) == certified Sol-3 spin-7 charge.  75%
    b2  Sol 3 form, t = 5/3 (k = 1), spin 5 (kernel 2-dim):    monic(cB) == certified Sol-3 spin-5 charge.                    75%
    b3  Sol 2 form, t = 5/3 (k = 1), spin 5:                    monic(cB) == certified Sol-2 spin-5 charge (= S/6 + 5E_1/6).   70%
    b4  Sol 3 form, t = 7/5 (k = 1/2), spin 9 (kernel 2-dim):  monic(cB) == certified Sol-3 spin-9 charge.                    70%
    So my sealed answer to the lead's question: the three-term ODE supplies, as the residue polynomial of a resonant
    (pole) coefficient, exactly the charge of ITS OWN solution family -- the Sol-3 ODE gives the Sol-3 table charge even
    where the kernel is two-dimensional, and at t = 5/3, where Sols 2 and 3 share one Gamma operator, their two different
    three-term ODEs give the two different spin-5 charges, which together span the kernel.
    Certified inputs: vev_sol3_w6/w8/w10, vev_sol2_w6 (all opened by this seat before; nothing here is blind).
(c) exploratory, not gated: at t = 5/3 (a = 1, any subset of the five roots may be moved into L) the other symmetric splits
    P = T(T^2-l0^2), L = T^2-l1^2 and, if even, further ones: is the resulting spin-5 polynomial in span{S, E_1}?
Controls: Gamma-form engines unchanged (hash 60795426 for mellin_wkb); tamper = the middle term x^(c/2) L x^(c/2) replaced by
x^c L(T) (wrong ordering, i.e. Lc = L): must FAIL (a).  A negative for (a) would falsify the class-U derivation for
deg P > 2; I would then report the first failing order and coefficient.
Exit codes: a_equiv.py 0 = (a) passes and tamper fails; 2 = (a) fails; b_nz.py 0 = b1-b4 as predicted; 2 otherwise.
No new hold-out is created (nothing to register); exact rational arithmetic throughout; no nsimplify.
