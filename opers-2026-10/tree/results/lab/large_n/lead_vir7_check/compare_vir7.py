"""Lead's comparison of Fable's registered VIR7 Sol 1 predictions (sha256 33af7eb7...) with the certified Sol 1 tables: spins 5, 7, 9 (vev_sol1_w6/8/10,
not blind) and spins 11, 13 (anchor11/13 vev_profile_sol1_w12/w14 -- BLIND: never opened by the seat or the lead before this run).  Conversion X = -PX^2,
Y = rho^2 - pi^2 (rho = p/2 - 1/p), monic in PX.  Negative control: t + 1/100.  Run from repo root."""
import json, sympy as sp, sys
t, X, Pi = sp.symbols('t X Pi')
P = json.load(open(sys.argv[1]))['cells']
def lw(w):
    out = {}
    for k, c in json.load(open(f'results/lab/vev/vev_sol1_w{w}.json'))['coefficients'].items():
        a, b = k.split(); out[(int(a[2:])//2, int(b[2:])//2)] = c
    return out
tabs = {'5': lw(6), '7': lw(8), '9': lw(10),
        '11': {tuple(map(int, k.split(','))): c for k, c in json.load(open('results/lab/anchor11/vev_profile_sol1_w12.json'))['vev'].items()},
        '13': {tuple(map(int, k.split(','))): c for k, c in json.load(open('results/lab/anchor13/vev_profile_sol1_w14.json'))['vev'].items()}}
def table(spin, tk):
    p2 = 2*(tk-1)/(tk+1); rho2 = p2/4 - 1 + 1/p2; e = 0
    for (a, b), c in tabs[spin].items(): e += sp.sympify(c).subs(t, tk)*(-X)**a*(rho2 - Pi)**b
    e = sp.expand(e); n = (int(spin)+1)//2; lead = e.coeff(X, n).subs(Pi, 0)
    return sp.expand(e/lead) if lead != 0 else None
def mism(d): return 0 if d == 0 else len(sp.Poly(d, X, Pi).terms())
for cell in P:
    tk = sp.Rational(cell[2:]); row = []
    for spin in ('5', '7', '9', '11', '13'):
        v = P[cell]['odd_spins'].get(spin)
        v = eval(v) if isinstance(v, str) else v
        if v.get('coefficients') is None: row.append(f'{spin}:no-charge'); continue
        pr = sum(sp.Rational(c)*X**(int(k.split(',')[0])//2)*Pi**(int(k.split(',')[1])//2) for k, c in v['coefficients'].items())
        tb = table(spin, tk)
        if tb is None: row.append(f'{spin}:table-lead0'); continue
        neg = table(spin, tk + sp.Rational(1, 100))
        row.append(f'{spin}:{mism(sp.expand(tb - pr))}/{len(sp.Poly(tb, X, Pi).terms())}(neg {mism(sp.expand(neg - pr))})')
    print(cell, ' '.join(row), flush=True)
