"""EXC3 T0(e): validation of wkb3x2 at MAXU = 4 (orders 1..4 in x^(-c)).
 (i) Sol 1 template (EXC1 level-1 oper, one point, M_S = 5): deformation -2 z xi/(xi-z)^2 - b xi/(xi-z) = -b - sum_m (2m + b) z^m / xi^m;
     the engine with orders m = 1..4 must reproduce EXC1's Schrodinger-form R_2, R_4, R_6 (exc1/o1_t2_n1.json, variables l0 = alpha,
     l1 = lam, z1 = z_S) EXACTLY, monic -- R_6 contains z_S^3 and z_S^4 terms, which test the orders 3 and 4.
 (ii) translation covariance T -> T + d (dP = 3, dL = 2) with all four orders present: cB_i unchanged, i = 2..6.
 (iii) regression: with q3 = q4 = 0 and MAXU = 4 the order-4 result equals the v1 engine (MAXU = 2) [order 6 differs by the U^3, U^4 products].
 Exit 0 iff (i) and (ii) hold."""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_EXC3.md', 'rb').read()).hexdigest() == open('SEAL_EXC3.sha256').read().split()[0]
import wkb3x2 as W
import wkb3x2_v1_maxu2 as W1
ok = True; t0 = time.time()
# (i)
MS = sp.Integer(5); a = 2 / (MS - 1); n = 2 + a; M = 1 / a; c = n * (1 + M); b = 1 / (1 + M); nb = n / b
al, lam, zS = sp.symbols('alpha lam zS'); z = zS / (MS + 1)
qT = [(-(2 * m + b) * z**m) * nb**2 * (2 * nb)**m for m in (1, 2, 3, 4)]       # c^2 (2c)^m q_m0
mk = [W.B0, W.C0, W.Q30, W.Q40]
eng = W.WKB3(2, 1, M, [1, 0, -W.l1**2 - nb**2 * b], [1, -c / 2 - W.l0], qcoef=[mk[0]], q2coef=[mk[1]], q3coef=[mk[2]], q4coef=[mk[3]])
eng.solve(6)
ref = json.load(open('<home>/fable-work/exc1/o1_t2_n1.json'))
l0r, l1r, z1r = sp.symbols('l0 l1 z1')
sub = {W.l1: nb * lam / (MS + 1), W.l0: -nb * al / (2 * (MS + 1))}
sub.update({mk[i]: qT[i] for i in range(4)})
for i in (2, 4, 6):
    cB, cG = eng.J(i)
    R = sp.Poly(sp.expand(cB.as_expr().subs(sub)), al, lam, zS); R = sp.Poly(sp.expand(R.as_expr() / R.coeff_monomial(lam**i)), al, lam, zS)
    Rref = sp.sympify(ref['R'][str(i)], locals={'l0': l0r, 'l1': l1r, 'z1': z1r}).subs({l0r: al, l1r: lam, z1r: zS})
    Rref = sp.Poly(sp.expand(Rref), al, lam, zS); Rref = sp.Poly(sp.expand(Rref.as_expr() / Rref.coeff_monomial(lam**i)), al, lam, zS)
    diff = sp.expand(R.as_expr() - Rref.as_expr())
    zmax = max([m[2] for m, _ in Rref.terms()])
    print('(i) order %d: wkb3x2(4 orders) == EXC1 Schrodinger engine: %s (reference has z^%d terms; %d monomials); cG zero %s  (%.0fs)' % (i, diff == 0, zmax, len(Rref.terms()), cG.is_zero, time.time() - t0), flush=True)
    if diff != 0: print('      difference:', diff)
    ok &= diff == 0 and cG.is_zero
# (ii)
k = sp.Integer(2); M = 1 / k; n = 3 + 2 * k; c = n * (1 + M); d = sp.Rational(1, 3); T = sp.Symbol('T')
P0 = T**3 + W.A2 * T**2 + W.A3 * T + W.A4; Lc = (T - c / 2)**2 - 7
Qs = [W.B1 * T + W.B0, W.C1 * T + W.C0, W.Q31 * T + W.Q30, W.Q41 * T + W.Q40]
def cf(poly, var, deg, desc=True):
    p_ = sp.Poly(sp.expand(poly), var)
    return [p_.coeff_monomial(var**(deg - j)) for j in range(deg + 1)] if desc else [p_.coeff_monomial(var**j) for j in range(deg + 1)]
e0 = W.WKB3(3, 2, M, cf(P0, T, 3), cf(Lc, T, 2), *[cf(q, T, 1, False) for q in Qs]); e0.solve(6)
es = W.WKB3(3, 2, M, cf(P0.subs(T, T + d), T, 3), cf(Lc.subs(T, T + d), T, 2), *[cf(q.subs(T, T + d), T, 1, False) for q in Qs]); es.solve(6)
for i in range(2, 7):
    b0, g0 = e0.J(i); bs, gs = es.J(i)
    same = sp.expand(b0.as_expr() - bs.as_expr()) == 0 and g0.is_zero and gs.is_zero
    syms = sorted(str(x) for x in b0.as_expr().free_symbols)
    print('(ii) order %d: cB unchanged under T -> T + %s: %s  (%d terms; symbols %s)' % (i, d, same, len(b0.terms()), syms), flush=True)
    ok &= same
# (iii)
e1 = W1.WKB3(3, 2, M, [1, W1.A2, W1.A3, W1.A4], [1, -c, c**2 / 4 - W1.l0**2], qcoef=[W1.B0, W1.B1], q2coef=[W1.C0, W1.C1]); e1.solve(4)
e2 = W.WKB3(3, 2, M, [1, W.A2, W.A3, W.A4], [1, -c, c**2 / 4 - W.l0**2], qcoef=[W.B0, W.B1], q2coef=[W.C0, W.C1]); e2.solve(4)
same4 = sp.expand(e1.J(4)[0].as_expr() - e2.J(4)[0].as_expr()) == 0
print('(iii) order 4 regression v2 == v1:', same4)
print('T0(e)', 'PASS' if ok and same4 else 'FAIL', '(%.0fs)' % (time.time() - t0)); sys.exit(0 if ok and same4 else 2)
