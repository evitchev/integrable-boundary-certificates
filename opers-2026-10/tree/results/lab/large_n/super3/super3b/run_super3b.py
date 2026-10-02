"""SUPER3b (POST-HOC DIAGNOSTIC, after the lead's observation, 2026-10-01): SUPER3's comparison with the COMMON momentum scale s free.
Dictionary alpha^2 = s X + tA, lam = s r Y + tL; 4 unknowns (s, r, tA, tL) fitted on spin 3 (XY, X, Y, 1), or on spin 5 (X^2 Y, X^2, X Y, X) where no
certified spin-3 point exists.  Everything else (spin-3 Y^2 where fitted at spin 3; spins 5, 7, 9) is checked with nothing fitted."""
import sys, json; sys.dont_write_bytecode = True
from fractions import Fraction as F
import sympy as sp
from suzuki_wkb import riccati, charge
X, Y, t = sp.symbols('X Y t'); A, L = sp.symbols('A L'); s, r, tA, tL = sp.symbols('s r tA tL')
W4 = json.load(open('vev_w4_points.json'))['1']
TAB = {W: json.load(open(f'vev_sol1_w{W}.json'))['coefficients'] for W in (6, 8, 10)}
def cert(k, tv):
    if k == 2:
        d = W4.get(str(sp.Rational(tv.numerator, tv.denominator)))
        if d is None: return None
        e = sum(sp.Rational(v) * X**int(kk.split(',')[0]) * Y**int(kk.split(',')[1]) for kk, v in d.items())
    else:
        e = 0
        for kk, v in TAB[2 * k].items():
            i, j = [int(z) for z in kk.replace('P^', '').replace('Q^', '').split()]
            e += sp.sympify(v, locals={'t': t}).subs(t, sp.Rational(tv.numerator, tv.denominator)) * X**(i // 2) * Y**(j // 2)
    e = sp.expand(e); return sp.expand(e / sp.Poly(e, X, Y).coeff_monomial(X**k))
fails = 0
for tv in [F(2), F(4), F(6), F(7), F(-4), F(5, 2), F(9, 2), F(8), F(-1, 3), F(10**6)]:
    M = (tv + 3) / (tv - 1)
    try:
        R = riccati(M, 10)
        polys = {}
        for k in (2, 3, 4, 5):
            p = [p for cl, (b0, f0, p) in charge(R[2 * k], M).items() if cl[0] != 'INTEGER-f'][0]
            polys[k] = sp.expand(sum(sp.Rational(v.numerator, v.denominator) * A**(i // 2) * L**j for (i, j), v in p.items()))
    except ZeroDivisionError as ex:
        print(f't = {tv}: degenerate ({ex}), skipped', flush=True); continue
    def mapped(k):
        e = sp.expand(polys[k].subs({A: s * X + tA, L: s * r * Y + tL}))
        return sp.expand(e / sp.Poly(e, X, Y).coeff_monomial(X**k))
    C2 = cert(2, tv)
    if C2 is not None:
        d = sp.Poly(sp.expand(sp.together(mapped(2) - C2) * 1), X, Y)
        eqs = [sp.numer(sp.together(sp.Poly(mapped(2) - C2, X, Y).coeff_monomial(m))) for m in (X * Y, X, Y, 1)]; fixed_at = 'spin 3 (XY, X, Y, 1)'
    else:
        eqs = [sp.numer(sp.together(sp.Poly(mapped(3) - cert(3, tv), X, Y).coeff_monomial(m))) for m in (X**2 * Y, X**2, X * Y, X)]; fixed_at = 'spin 5 (X^2Y, X^2, XY, X)'
    sols = [so for so in sp.solve(eqs, [s, r, tA, tL], dict=True) if so.get(s, 1) != 0]
    print(f't = {tv}, M = {M}: dictionary fixed at {fixed_at}: {len(sols)} solution(s)', flush=True)
    for so in sols:
        row = {}
        for k in (2, 3, 4, 5):
            Ck = cert(k, tv)
            if Ck is None: row[2 * k - 1] = 'no table'; continue
            diff = sp.Poly(sp.expand(mapped(k).subs(so) - Ck), X, Y)
            nz = [(m, c) for m, c in zip(diff.monoms(), diff.coeffs()) if sp.simplify(c) != 0]
            row[2 * k - 1] = 'MATCH' if not nz else f'{len(nz)} MISMATCH {[(str(m), str(c)) for m, c in nz[:3]]}'
        print('   ', {str(kk): str(v) for kk, v in so.items()}, '|', row, flush=True)
print('done')
