"""Post-report check (SCAL2 screen finding): SUPER3 fibres t = 7, 4 lie on Gamma-pole sets (s/M integer at spins 5 / 7).
Are SUPER3's normalised charges there equal to the limits t -> t0 (+- 1e-9)?"""
import sys; sys.dont_write_bytecode = True
from fractions import Fraction as F
import sympy as sp
from suzuki_wkb import riccati, charge
A, L = sp.symbols('A L')
def ch(tv, k):
    M = (tv + 3) / (tv - 1); R = riccati(M, 2 * k)
    p = [p for cl, (b0, f0, p) in charge(R[2 * k], M).items() if cl[0] != 'INTEGER-f'][0]
    e = sp.expand(sum(sp.Rational(v.numerator, v.denominator) * A**(i // 2) * L**j for (i, j), v in p.items())); return sp.expand(e / sp.Poly(e, A, L).coeff_monomial(A**k))
for t0, k in ((F(7), 3), (F(4), 4)):
    e0 = ch(t0, k)
    for d in (F(1, 10**9), -F(1, 10**9), F(1, 10**12)):
        e1 = ch(t0 + d, k); diff = sp.Poly(sp.expand(e1 - e0), A, L)
        print(f't0 = {t0}, spin {2*k-1}, t0 {"+" if d > 0 else "-"} {float(abs(d)):.0e}: max |coefficient difference| = {max([abs(float(c)) for c in diff.coeffs()] + [0.0]):.3e}')
