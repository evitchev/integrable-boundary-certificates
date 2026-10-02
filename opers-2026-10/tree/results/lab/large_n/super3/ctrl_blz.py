"""SUPER3 positive control: alpha = 0 must give BLZ vacuum eigenvalues (hep-th/9812247 / hep-th/9604044 Delta-polynomials)
I3 ~ D^2 - (c+2)/12 D + c(5c+22)/2880, I5 ~ D^3 - (c+4)/8 D^2 + (c+2)(3c+20)/576 D - c(3c+14)(7c+68)/290304, c = 13 - 6(b^2 + 1/b^2), b^2 = 1/(M+1).
Delta = mu L + nu fitted on I3 (D^1, D^0 coefficients), then I3's consistency and all of I5 are checks."""
import sys; sys.dont_write_bytecode = True
from fractions import Fraction as F
import sympy as sp
from suzuki_wkb import riccati, charge
L, D, mu, nu = sp.symbols('L D mu nu')
bad = 0
for M in (F(5), F(7, 3), F(9, 5), F(-1, 3), F(11, 3), F(-2), F(-5, 2)):
    try:
        R = riccati(M, 6); [charge(R[n], M) for n in (4, 6)]
    except ZeroDivisionError:
        print(f'M = {M}: degenerate fibre (Gamma pole), screened'); continue
    b2 = 1 / (M + 1); c = 13 - 6 * (b2 + 1 / b2); c = sp.Rational(c.numerator, c.denominator)
    def ch(k):
        out = [p for cl, (b0, f0, p) in charge(R[2 * k], M).items() if cl[0] != 'INTEGER-f' and any(i == 0 and v for (i, j), v in p.items())]
        p = out[0]; e = sum(sp.Rational(v.numerator, v.denominator) * L**j for (i, j), v in p.items() if i == 0)
        e = sp.expand(e); return sp.expand(e / sp.Poly(e, L).LC())
    I3t = sp.expand((D**2 - (c + 2) / 12 * D + c * (5 * c + 22) / 2880).subs(D, mu * L + nu))
    I5t = sp.expand((D**3 - (c + 4) / 8 * D**2 + (c + 2) * (3 * c + 20) / 576 * D - c * (3 * c + 14) * (7 * c + 68) / 290304).subs(D, mu * L + nu))
    e3 = sp.Poly(sp.expand(I3t / mu**2) - ch(2), L).coeffs()
    sols = sp.solve(e3, [mu, nu], dict=True); ok = False
    for sub in sols:
        if sub.get(mu, 1) == 0: continue
        d5 = [sp.simplify(x) for x in sp.Poly(sp.expand((I5t / mu**3).subs(sub)) - ch(3), L).coeffs()]
        if all(v == 0 for v in d5): ok = True; print(f'M = {M}: c = {c}; Delta = {sp.expand((mu*L+nu).subs(sub))}; I3 fit (2 eqs, 2 unknowns), I5 (3 checks) MATCH BLZ')
        else: print(f'M = {M}: solution {sub} I5 residual {d5}')
    if not ok: bad += 1; print(f'M = {M}: NO MATCH', e3)
print('BLZ CONTROL', 'OK' if not bad else 'FAILED'); sys.exit(bad > 0)
