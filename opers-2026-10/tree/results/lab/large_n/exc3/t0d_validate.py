"""EXC3 T0(d): validation of wkb3x2 (second order in x^(-c)).
 (i) regression: with q2coef = None, cB_i equals wkb3x for the Sol 3 template with a general quartic and first-order term (t = 2).
 (ii) Sol 1 template at M_S = 5 (EXC1 level-1 oper in three-term form, EXC2 seal (ii)): deformation -2 z xi/(xi-z)^2 - b xi/(xi-z)
      = -b - (2+b) z/xi - (4+b) z^2/xi^2 - ...; with the SECOND order included the engine must reproduce EXC1's Schrodinger result
      R_4/c4 = -7 al^4/400 + 3 al^2 lam^2/10 + 9 al^2/2 + 168 al z + lam^4 + 18 lam^2 + 240 z^2 + 337/5 EXACTLY (x0_validate had everything
      but the 240 z^2 at first order).
 (iii) translation covariance T -> T + d with first- and second-order terms present (dP = 3, dL = 2, the Sol 2 shape): cB_i, i = 2..6, unchanged.
 Exit 0 iff all hold."""
import hashlib, sys, time
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_EXC3.md', 'rb').read()).hexdigest() == open('SEAL_EXC3.sha256').read().split()[0]
import wkb3x as W1
import wkb3x2 as W2
ok = True
t0 = time.time()
# (i) regression against wkb3x (Sol 3 template at t = 2: k = 2, M = 1/2, general quartic + first order)
k = sp.Integer(2); M = 1 / k; n = k + 4; c = n * (1 + M)
e1 = W1.WKB3(4, 1, M, [1, 0, W1.A2, W1.A3, W1.A4], [1, -c / 2], qcoef=[W1.B0, W1.B1, W1.B2]); e1.solve(6)
e2 = W2.WKB3(4, 1, M, [1, 0, W2.A2, W2.A3, W2.A4], [1, -c / 2], qcoef=[W2.B0, W2.B1, W2.B2]); e2.solve(6)
for i in range(2, 7):
    b1, g1 = e1.J(i); b2, g2 = e2.J(i)
    same = sp.expand(b1.as_expr() - b2.as_expr()) == 0 and sp.expand(g1.as_expr() - g2.as_expr()) == 0
    ok &= same
print('(i) regression wkb3x2 == wkb3x (q2 = 0), orders 2..6: %s  (%.0fs)' % (ok, time.time() - t0), flush=True)
# (ii) Sol 1 template, second order
MS = sp.Integer(5); a = 2 / (MS - 1); n = 2 + a; M = 1 / a; c = n * (1 + M); b = 1 / (1 + M); nb = n / b
al, lam, zS = sp.symbols('alpha lam zS')
z = zS / (MS + 1)
B0_1 = 2 * nb**3 * (-(2 + b) * z)          # first order:  2 c^3 q10
B0_2 = 4 * nb**4 * (-(4 + b) * z**2)       # second order: 4 c^4 q20
eng = W2.WKB3(2, 1, M, [1, 0, -W2.l1**2 - nb**2 * b], [1, -c / 2 - W2.l0], qcoef=[W2.B0], q2coef=[W2.B1])
eng.solve(4)
cB4, cG4 = eng.J(4)
sub = {W2.l1: nb * lam / (MS + 1), W2.l0: -nb * al / (2 * (MS + 1)), W2.B0: B0_1, W2.B1: B0_2}
R4 = sp.Poly(sp.expand(cB4.as_expr().subs(sub)), al, lam, zS)
R4 = sp.Poly(sp.expand(R4.as_expr() / R4.coeff_monomial(lam**4)), al, lam, zS)
ref = sp.Poly(-7 * al**4 / 400 + 3 * al**2 * lam**2 / 10 + 9 * al**2 / 2 + 168 * al * zS + lam**4 + 18 * lam**2 + 240 * zS**2 + sp.Rational(337, 5), al, lam, zS)
diff = sp.expand(R4.as_expr() - ref.as_expr())
print('(ii) Sol 1 template, second order: wkb3x2 R_4/c4 =', R4.as_expr())
print('     EXC1 Schrodinger engine         =', ref.as_expr())
print('     difference:', diff, ' cG zero:', cG4.is_zero, ' (%.0fs)' % (time.time() - t0), flush=True)
ok &= diff == 0 and cG4.is_zero
# (iii) translation covariance, dP = 3, dL = 2 (Sol 2 shape), t = 2: k = 2, M = 1/2
k = sp.Integer(2); M = 1 / k; n = 3 + 2 * k; c = n * (1 + M); d = sp.Rational(1, 3)
T = sp.Symbol('T')
P0 = T**3 + W2.A2 * T**2 + W2.A3 * T + W2.A4
Q1 = W2.B2 * T**2 + W2.B1 * T + W2.B0
Q2 = W2.l0 * T**2 + W2.l1 * T + W2.A2 * W2.A3           # arbitrary second-order coefficients built from the available symbols
Lc = (T - c / 2)**2 - 7                                   # L(T - c/2) with l0^2 = 7 (numeric, so that l0 is free for Q2)
def coeffs(poly, var, deg, desc=True):
    p_ = sp.Poly(sp.expand(poly), var)
    return [p_.coeff_monomial(var**(deg - j)) for j in range(deg + 1)] if desc else [p_.coeff_monomial(var**j) for j in range(deg + 1)]
e0 = W2.WKB3(3, 2, M, coeffs(P0, T, 3), coeffs(Lc, T, 2), qcoef=coeffs(Q1, T, 2, False), q2coef=coeffs(Q2, T, 2, False)); e0.solve(6)
Psh, Lsh, Q1s, Q2s = [sp.expand(x.subs(T, T + d)) for x in (P0, Lc, Q1, Q2)]
es = W2.WKB3(3, 2, M, coeffs(Psh, T, 3), coeffs(Lsh, T, 2), qcoef=coeffs(Q1s, T, 2, False), q2coef=coeffs(Q2s, T, 2, False)); es.solve(6)
for i in range(2, 7):
    b0, g0 = e0.J(i); bs, gs = es.J(i)
    same = sp.expand(b0.as_expr() - bs.as_expr()) == 0 and g0.is_zero and gs.is_zero
    nterms = len(b0.terms())
    hasU2 = any(mon[5] + mon[6] + mon[7] >= 2 or (mon[0] + mon[1] >= 1) for mon, _ in b0.terms())   # B's or (l0, l1 = Q2 markers in this test) present
    print('(iii) order %d: cB unchanged under T -> T + %s: %s  (%d terms; second-order markers present: %s; cG zero: %s)' % (i, d, same, nterms, hasU2, g0.is_zero and gs.is_zero), flush=True)
    ok &= same
print('T0(d)', 'PASS' if ok else 'FAIL', '(%.0fs)' % (time.time() - t0))
sys.exit(0 if ok else 2)
