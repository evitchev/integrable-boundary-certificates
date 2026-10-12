"""PILLOW1 positive control for the r1c cocycle conventions: Sol 1's three-term ODE (P = T^2 - l1^2, L = T - l0) vs its own Gamma form (PROOF1 Theorem 1)
must give F == 1 identically (no unknowns), at t' = 1/3 and t' = 2.  Negative control: L = T - l0 - 1/7 must fail.
Convention identical to r1c: p_h(index) = coefficient attached to psihat(index); f_c(nu) = -p_c(nu - c)/p_base(nu), f_n(nu) = p_n(nu - n)/p_base(nu)."""
import sys; sys.dont_write_bytecode = True
import sympy as sp
nu, L0, L1 = sp.symbols('nu l0 l1')
ok_all = True
for tp, shiftL, expect in ((sp.Rational(1, 3), 0, True), (sp.Integer(2), 0, True), (sp.Rational(1, 3), sp.Rational(1, 7), False)):
    a = (tp - 1) / 2; n = (tp + 3) / 2; sig = n / a; c = n + sig; s = (n - 1) / 2
    P = lambda T: T**2 - L1**2; L = lambda T: T - L0 - shiftL
    pbase = lambda v: P(s - v); pcc = lambda v: -L(s - v - c / 2); pnn = lambda v: sp.Integer(-1)
    Grat = lambda T: ((T - L0) / sig + (a + 1) / 2) / ((T - L0) / sig + (1 - a) / 2)
    def Hratio(v1, v2):
        T2 = s - v2; d = sp.simplify((s - v1) - T2); base = P(s - v1) / P(s - v2)
        return base * Grat(T2) if sp.simplify(d - sig) == 0 else base / Grat(s - v1)
    fc = lambda v: -pcc(v - c) / pbase(v); fn = lambda v: pnn(v - n) / pbase(v)
    F = fc(nu) * fn(nu - c) / (fn(nu) * fc(nu - n)) * Hratio(nu - c, nu - n)
    res = sp.simplify(F - 1) == 0
    print(f"t' = {tp}, L shift {shiftL}: F == 1 identically: {res}  (expected {expect})")
    ok_all &= (res == expect)
print('CONTROL', 'PASS' if ok_all else 'FAIL'); sys.exit(0 if ok_all else 1)
