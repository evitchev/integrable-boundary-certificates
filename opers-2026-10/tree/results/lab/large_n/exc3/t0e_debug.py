import json, sys
import sympy as sp
sys.dont_write_bytecode = True
import wkb3x2 as W
MS = sp.Integer(5); a = 2 / (MS - 1); n = 2 + a; M = 1 / a; c = n * (1 + M); b = 1 / (1 + M); nb = n / b
al, lam, zS = sp.symbols('alpha lam zS'); z = zS / (MS + 1)
qT = [(-(2 * m + b) * z**m) * nb**2 * (2 * nb)**m for m in (1, 2, 3, 4)]
mk = [W.B0, W.C0, W.Q30, W.Q40]
eng = W.WKB3(2, 1, M, [1, 0, -W.l1**2 - nb**2 * b], [1, -c / 2 - W.l0], qcoef=[mk[0]], q2coef=[mk[1]], q3coef=[mk[2]], q4coef=[mk[3]])
eng.solve(6)
ref = json.load(open('<home>/fable-work/exc1/o1_t2_n1.json'))
l0r, l1r, z1r = sp.symbols('l0 l1 z1')
sub = {W.l1: nb * lam / (MS + 1), W.l0: -nb * al / (2 * (MS + 1))}; sub.update({mk[i]: qT[i] for i in range(4)})
cB, cG = eng.J(6)
R = sp.Poly(sp.expand(cB.as_expr().subs(sub)), al, lam, zS)
Rref = sp.Poly(sp.expand(sp.sympify(ref['R']['6'], locals={'l0': l0r, 'l1': l1r, 'z1': z1r}).subs({l0r: al, l1r: lam, z1r: zS})), al, lam, zS)
print('engine lam^6 coeff:', R.coeff_monomial(lam**6), ' reference lam^6 coeff:', Rref.coeff_monomial(lam**6))
print('engine terms:', R.as_expr())
print('reference   :', Rref.as_expr())
# proportionality test
m0 = [m for m, c_ in Rref.terms() if c_ != 0][0]
ratio = R.coeff_monomial(al**m0[0] * lam**m0[1] * zS**m0[2]) / Rref.coeff_monomial(al**m0[0] * lam**m0[1] * zS**m0[2])
print('ratio on monomial', m0, '=', ratio, '; engine - ratio*reference =', sp.expand(R.as_expr() - ratio * Rref.as_expr()))
