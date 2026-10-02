"""SUPER3 Step B driver: compare the extended Suzuki WKB charges with the certified Sol 1 tables. Usage: run_super3.py [tamper]."""
import sys, json; sys.dont_write_bytecode = True
from fractions import Fraction as F
import sympy as sp
from suzuki_wkb import riccati, charge, is_total_derivative
tamper = 'tamper' in sys.argv
X, Y, t = sp.symbols('X Y t')
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
def suz(k, M, R):
    """odd spin 2k-1 from R_2k: (charge poly (half-integer-f class), integer-f part is a total derivative?, alpha-parities)"""
    ch = charge(R[2 * k], M)
    half = [p for cl, (b0, f0, p) in ch.items() if cl[0] != 'INTEGER-f']
    integ = {}
    for cl, (b0, f0, p) in ch.items():
        if cl[0] == 'INTEGER-f': integ.update(p)
    assert len(half) == 1
    return half[0], is_total_derivative(integ, M), {i % 2 for (i, j) in half[0]}
def even_spin(k, M, R):
    """even spin 2k-2 from R_(2k-1): alpha-parities of the half-integer-f (charge) part; integer-f part a total derivative?"""
    ch = charge(R[2 * k - 1], M)
    half = {}; integ = {}
    for cl, (b0, f0, p) in ch.items():
        (integ if cl[0] == 'INTEGER-f' else half).update(p if cl[0] == 'INTEGER-f' else {})
        if cl[0] != 'INTEGER-f':
            for kk, v in p.items(): half[kk] = v
    return {i % 2 for (i, j), v in half.items() if v}, is_total_derivative(integ, M)
A, L = sp.symbols('A L')
def topoly(d): return sp.expand(sum(sp.Rational(v.numerator, v.denominator) * A**(i // 2) * L**j for (i, j), v in d.items()))
r, tA, tL = sp.symbols('r tA tL')
results = {}
for tv in [F(2), F(4), F(6), F(7), F(-2), F(-4), F(5, 2), F(9, 2), F(-1, 3), F(-1, 3) + F(1, 10**6), F(-1, 3) - F(1, 10**6), F(1000), F(10**6)]:
    M = (tv + 3) / (tv - 1) + (F(1, 7) if tamper else 0)
    try:
        R = riccati(M, 10)
        S = {k: suz(k, M, R) for k in (2, 3, 4, 5)}
    except ZeroDivisionError as ex:
        print(f't = {tv}: DEGENERATE ({ex}) -- skipped', flush=True); continue
    p5 = all(S[k][1] and S[k][2] == {0} for k in S)
    ev = {2 * k - 2: even_spin(k, M, R) for k in (2, 3, 4, 5)}
    ncls = {2 * k - 1: sorted(S[k][2]) for k in S}
    polys = {k: topoly(S[k][0]) for k in S}
    def mapped(k):
        e = sp.expand(polys[k].subs({A: X + tA, L: r * Y + tL}))
        return sp.expand(e / sp.Poly(e, X, Y).coeff_monomial(X**k))
    C2 = cert(2, tv)
    if C2 is None:
        # no certified spin-3 point: fix the dictionary at spin 5 top + sub-top instead (disclosed)
        C3 = cert(3, tv); m3 = mapped(3)
        eqs = [sp.Poly(m3 - C3, X, Y).coeff_monomial(mo) for mo in (X**2 * Y, X**2, X * Y)]
        fixed_at = 'spin 5 (no spin-3 point)'
    else:
        m2 = mapped(2); eqs = [sp.Poly(m2 - C2, X, Y).coeff_monomial(mo) for mo in (X * Y, X, Y)]; fixed_at = 'spin 3'
    sols = sp.solve(eqs, [r, tA, tL], dict=True)
    line = f't = {tv}, M = {M}: P5 (odd spins alpha-even, integer-f part exact) {p5}; even spins (alpha parities of charge part, int-f exact): {ev}; dictionary fixed at {fixed_at}: {len(sols)} solution(s)'
    print(line, flush=True)
    res = []
    for s_ in sols:
        row = {'dict': {str(k): str(v) for k, v in s_.items()}}
        for k in (2, 3, 4, 5):
            Ck = cert(k, tv)
            if Ck is None: row[2 * k - 1] = 'no table'; continue
            diff = sp.Poly(sp.expand(mapped(k).subs(s_) - Ck), X, Y)
            nz = [(str(mo), str(sp.nsimplify(c) if False else c)) for mo, c in zip(diff.monoms(), diff.coeffs()) if sp.simplify(c) != 0]
            top = [m for m, c in zip(diff.monoms(), diff.coeffs()) if sum(m) == k and sp.simplify(c) != 0]
            row[2 * k - 1] = 'MATCH' if not nz else f'{len(nz)} mismatched coeffs (top mismatches {len(top)}): {nz[:4]}'
        print('   dict', row['dict'], '|', {s: row[s] for s in (3, 5, 7, 9)}, flush=True)
        res.append(row)
    results[str(tv)] = {'M': str(M), 'p5': p5, 'even_spins': {str(k): [sorted(v[0]), v[1]] for k, v in ev.items()}, 'fixed_at': fixed_at, 'rows': res}
json.dump(results, open('run_super3%s.json' % ('_tamper' if tamper else ''), 'w'), indent=1)
