"""Lead's independent comparison: Fable's pre-registered VIR4d predictions (snapshot, sha256 b8446a9f...) vs the certified Sol 3
spin-11/13 tables (vev_profile_sol3_w12/w14 'vev', rational in t) at t_k = (3k+2)/(k+2); conversion P_art^2 = -PX^2, Q_art^2 = rho^2 - pi^2, monic in PX.
Also negative controls (wrong fibre).  Run from the repo root."""
import json, sympy as sp, glob, sys
t, X, Pi = sp.symbols('t X Pi')
P = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'results/lab/large_n/vir4d_prereg/PREDICTION_VIR4d.json'))['cells']
tabs = {'11': json.load(open('results/lab/anchor11/vev_profile_sol3_w12.json'))['vev'],
        '13': json.load(open(glob.glob('results/lab/anchor13/**/*sol3*w14*.json', recursive=True)[0]))['vev']}
def table(spin, tk):
    p2 = 2*(tk-1)/(tk+1); rho2 = p2/4 - 1 + 1/p2; e = 0
    for kk, c in tabs[spin].items():
        a, b = map(int, kk.split(',')); e += sp.sympify(c).subs(t, tk)*(-X)**a*(rho2 - Pi)**b
    e = sp.expand(e); n = (int(spin)+1)//2; return sp.expand(e/e.coeff(X, n).subs(Pi, 0))
def pred(cell):
    return sum(sp.Rational(c)*X**(int(k.split(',')[0])//2)*Pi**(int(k.split(',')[1])//2) for k, c in P[cell]['coefficients'].items())
def mism(d): return 0 if d == 0 else len(sp.Poly(d, X, Pi).terms())
for cell in P:
    k = int(cell.split(',')[0][2:]); spin = cell.split('=')[2]; tk = sp.Rational(3*k+2, k+2)
    print(cell, 't =', tk, 'coeffs', len(P[cell]['coefficients']), 'mismatches', mism(sp.expand(table(spin, tk) - pred(cell))))
for cell, tk in (('k=4,spin=11', sp.Rational(5, 2)), ('k=6,spin=11', sp.Rational(7, 3)), ('k=4,spin=11', sp.Rational(7, 3) + sp.Rational(1, 100))):
    print('NEG CONTROL', cell, 'vs table at t =', tk, ': mismatches', mism(sp.expand(table('11', tk) - pred(cell))))
