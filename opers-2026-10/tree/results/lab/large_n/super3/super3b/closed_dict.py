"""SUPER3b (post-hoc): ZERO-parameter check of the closed-form dictionary read off run_super3b.log:
alpha^2 = -4 M (M+1) X,  (l + 1/2)^2 = (M-1)^2/4 - (M+1) Y  (lam = l(l+1) = (l+1/2)^2 - 1/4),  M = (t+3)/(t-1).
All coefficients at spins 3, 5, 7, 9 compared with the certified Sol 1 tables, nothing fitted, at old and NEW fibres."""
import sys, json; sys.dont_write_bytecode = True
from fractions import Fraction as F
import sympy as sp
from suzuki_wkb import riccati, charge
X, Y, t = sp.symbols('X Y t'); A, L = sp.symbols('A L')
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
bad = 0; nchk = 0
for tv in [F(2), F(4), F(6), F(7), F(8), F(-4), F(5, 2), F(9, 2), F(-1, 3), F(13, 3), F(17, 5), F(25, 3), F(12), F(14), F(-7, 5), F(23, 7), F(3, 7), F(-17, 4), F(10**6)]:
    M = (tv + 3) / (tv - 1); Ms = sp.Rational(M.numerator, M.denominator)
    try:
        R = riccati(M, 10)
    except ZeroDivisionError:
        print(f't = {tv}: degenerate, skipped'); continue
    row = {}
    for k in (2, 3, 4, 5):
        Ck = cert(k, tv)
        if Ck is None: row[2 * k - 1] = '-'; continue
        try:
            p = [p for cl, (b0, f0, p) in charge(R[2 * k], M).items() if cl[0] != 'INTEGER-f'][0]
        except ZeroDivisionError:
            row[2 * k - 1] = 'degenerate'; continue
        e = sp.expand(sum(sp.Rational(v.numerator, v.denominator) * A**(i // 2) * L**j for (i, j), v in p.items()))
        e = sp.expand(e.subs({A: -4 * Ms * (Ms + 1) * X, L: (Ms - 1)**2 / 4 - (Ms + 1) * Y - sp.Rational(1, 4)}))
        e = sp.expand(e / sp.Poly(e, X, Y).coeff_monomial(X**k))
        diff = sp.Poly(sp.expand(e - Ck), X, Y); nz = [c for c in diff.coeffs() if c != 0]
        nchk += len(sp.Poly(Ck, X, Y).monoms()) - 1
        row[2 * k - 1] = 'MATCH' if not nz else f'{len(nz)} MISMATCH'; bad += len(nz) > 0
    print(f't = {tv}, M = {M}: {row}', flush=True)
print(f'ZERO-PARAMETER CHECK: {"ALL MATCH" if not bad else str(bad) + " spin-rows mismatch"}; {nchk} non-normalised coefficients compared')
sys.exit(1 if bad else 0)
