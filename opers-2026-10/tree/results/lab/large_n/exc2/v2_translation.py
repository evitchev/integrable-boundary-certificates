"""EXC2 control V2 (no data): translation covariance of the extended engine.  Conjugation by x^d sends T -> T + d in every
term:  P(T) -> P(T+d),  Lc(T) = T - c/2 -> T + d - c/2,  Qt(T) -> Qt(T+d).  The charges J_i, i >= 2, must not change.
This exercises the odd coefficient of P (T^3), and the j = 0, 1, 2 paths of the first-order term (x0_validate only j = 0).
Usage: v2_translation.py <t> <d> [tamper]  (tamper: Qt NOT shifted -> must fail).  Exit 0 iff cB_i equal for i = 2..6, cG = 0."""
import json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import wkb3x as W
t = sp.Rational(sys.argv[1]); d = sp.Rational(sys.argv[2]); tamper = len(sys.argv) > 3
k = 2 * (t - 1) / (3 - t); M = 1 / k; n = k + 4; c = n * (1 + M)
A2, A3, A4, B0, B1, B2 = [x for x in (W.A2, W.A3, W.A4, W.B0, W.B1, W.B2)]
T = sp.Symbol('T')
Psh = sp.Poly(sp.expand(((T + d)**4 + A2 * (T + d)**2 + A3 * (T + d) + A4)), T)
Qsh = sp.Poly(sp.expand(B2 * (T + d)**2 + B1 * (T + d) + B0), T)
pcoef = [Psh.coeff_monomial(T**(4 - j)) for j in range(5)]
qcoef = [B0, B1, B2] if tamper else [Qsh.coeff_monomial(T**j) for j in range(3)]
t0 = time.time()
ref = json.load(open('u3_engine_t%s.json' % str(t).replace('/', '_')))
loc = {'A2': A2, 'A3': A3, 'A4': A4, 'B0': B0, 'B1': B1, 'B2': B2}
eng = W.WKB3(4, 1, M, pcoef, [1, d - c / 2], qcoef=qcoef)
eng.solve(6)
ok = True
for i in (2, 3, 4, 5, 6):
    cB, cG = eng.J(i)
    if i % 2 == 0:
        same = sp.expand(cB.as_expr() - sp.sympify(ref['cB%d' % i], locals=loc)) == 0
        print('order %d: cG zero %s; cB equal to the unshifted engine: %s   (%.0fs)' % (i, cG.is_zero, same, time.time() - t0), flush=True)
        ok &= same and cG.is_zero
    else:
        print('order %d (odd): cB = %s ; cG zero %s' % (i, cB.as_expr(), cG.is_zero), flush=True)
print('V2 translation covariance (d = %s)%s: %s' % (d, ' TAMPER' if tamper else '', 'PASS' if ok else 'FAIL'))
if tamper:
    sys.exit(0 if not ok else 1)
sys.exit(0 if ok else 2)
