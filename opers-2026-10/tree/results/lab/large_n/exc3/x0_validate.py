"""EXC2: validation of the extended engine wkb3x (first-order 1/xi term) on the Sol-1 template at M_S = 5 (t = 2).
The level-1 Sol-1 oper in three-term form: P = T^2 - l1^2 - (n/b)^2 b, Lc = T - c/2 - l0, Qt = B0 = -2(b+2)(n/b)^3 z,
z = z_S/(M_S+1), l1 = (n/b) lam/(M_S+1), l0 = -(n/b) alpha/(2(M_S+1)).
Must reproduce EXC1's Schrodinger-engine result  R_4/c4 = -7 al^4/400 + 3 al^2 lam^2/10 + 9 al^2/2 + 168 al z + lam^4 + 18 lam^2 + 240 z^2 + 337/5
in every term of degree <= 1 in z (the z^2 term needs the second order in 1/xi, not implemented).  Exit 0 iff so."""
import sys
import sympy as sp
sys.dont_write_bytecode = True
import wkb3x as W
MS = sp.Integer(5)
a = 2 / (MS - 1); n = 2 + a; M = 1 / a; c = n * (1 + M); b = 1 / (1 + M)
nb = n / b
al, lam, zS = sp.symbols('alpha lam zS')
eng = W.WKB3(2, 1, M, [1, 0, -W.l1**2 - nb**2 * b], [1, -c / 2 - W.l0], qcoef=[W.B0])
eng.solve(4)
cB2, cG2 = eng.J(2); cB4, cG4 = eng.J(4)
print('cG zero at orders 2, 4:', cG2.is_zero, cG4.is_zero)
sub = {W.l1: nb * lam / (MS + 1), W.l0: -nb * al / (2 * (MS + 1)), W.B0: -2 * (b + 2) * nb**3 * zS / (MS + 1)}
R4 = sp.Poly(sp.expand(cB4.as_expr().subs(sub)), al, lam, zS)
R4 = sp.Poly(sp.expand(R4.as_expr() / R4.coeff_monomial(lam**4)), al, lam, zS)
ref = sp.Poly(-7 * al**4 / 400 + 3 * al**2 * lam**2 / 10 + 9 * al**2 / 2 + 168 * al * zS + lam**4 + 18 * lam**2 + 240 * zS**2 + sp.Rational(337, 5), al, lam, zS)
print('wkb3x :', R4.as_expr())
print('EXC1  :', ref.as_expr())
diff = sp.Poly(R4.as_expr() - ref.as_expr(), al, lam, zS)
ok = all(m[2] >= 2 for m, cc in diff.terms() if cc != 0) and cG2.is_zero and cG4.is_zero
print('difference (must contain only z^2 and higher):', diff.as_expr())
print('extended-engine validation', 'PASS' if ok else 'FAIL')
sys.exit(0 if ok else 2)
