"""Lead's independent comparison of Fable's registered VIR5a predictions (sha256 0c12a7e9...) with the certified Sol 3 tables
(spins 5, 7, 9: results/lab/vev/vev_sol3_w6/8/10.json 'coefficients' "P^a Q^b"; spins 11, 13: anchor11/13 vev_profile 'vev' "a,b"),
evaluated at t = (3k+2)/(k+2).  Conversion P_art^2 = -PX^2, Q_art^2 = rho^2 - pi^2, rho = p/2 - 1/p, p^2 = 2(t-1)/(t+1); monic in PX.  Run from repo root."""
import json, sympy as sp, glob, sys
t, X, Pi = sp.symbols('t X Pi')
P = json.load(open(sys.argv[1]))['cells']
def load_w(w):
    d = json.load(open(f'results/lab/vev/vev_sol3_w{w}.json'))['coefficients']; out = {}
    for k, c in d.items():
        a, b = k.split(); out[(int(a[2:])//2, int(b[2:])//2)] = c
    return out
tabs = {'5': load_w(6), '7': load_w(8), '9': load_w(10)}
for spin, path in (('11', 'results/lab/anchor11/vev_profile_sol3_w12.json'), ('13', glob.glob('results/lab/anchor13/**/*sol3*w14*.json', recursive=True)[0])):
    tabs[spin] = {tuple(map(int, k.split(','))): c for k, c in json.load(open(path))['vev'].items()}
def table(spin, tk):
    p2 = 2*(tk-1)/(tk+1); rho2 = p2/4 - 1 + 1/p2; e = 0
    for (a, b), c in tabs[spin].items():
        e += sp.sympify(c).subs(t, tk)*(-X)**a*(rho2 - Pi)**b
    e = sp.expand(e); n = (int(spin)+1)//2; lead = e.coeff(X, n).subs(Pi, 0)
    return (sp.expand(e/lead) if lead != 0 else None), lead
def pred(cell, spin):
    sc = P[cell]['spins'].get(spin)
    if sc is None or sc.get('coefficients') is None: return None
    return sum(sp.Rational(c)*X**(int(k.split(',')[0])//2)*Pi**(int(k.split(',')[1])//2) for k, c in sc['coefficients'].items())
tot = 0; bad = 0
for cell in sorted(P, key=lambda c: sp.Rational(c[2:])):
    k = sp.Rational(cell[2:]); tk = (3*k+2)/(k+2); row = []
    for spin in ('5', '7', '9', '11', '13'):
        try:
            tb, lead = table(spin, tk)
        except Exception as ex:
            row.append(f'{spin}:table-error'); continue
        pr = pred(cell, spin)
        lc = P[cell]['spins'].get(spin, {}).get('leading_coefficient')
        if tb is None or pr is None:
            row.append(f'{spin}:skip(tableLead={lead},predLead={lc})'); continue
        d = sp.expand(tb - pr); m = 0 if d == 0 else len(sp.Poly(d, X, Pi).terms()); tot += 1; bad += (m > 0)
        row.append(f'{spin}:{m}/{len(sp.Poly(tb, X, Pi).terms())}' + ('' if lc is None or sp.Rational(lc) != 0 else '(predLead0)'))
    print(cell, 't =', tk, ' '.join(row), flush=True)
print('cells compared', tot, 'with mismatches', bad)
# NEGATIVE CONTROLS: each prediction against the table at a shifted fibre t_k + 1/100, and against the next fibre in the list
print('--- negative controls ---')
cells = sorted(P, key=lambda c: sp.Rational(c[2:]))
for i, cell in enumerate(cells):
    k = sp.Rational(cell[2:]); tk = (3*k+2)/(k+2)
    for tt in (tk + sp.Rational(1, 100), (3*sp.Rational(cells[(i+1) % len(cells)][2:])+2)/(sp.Rational(cells[(i+1) % len(cells)][2:])+2)):
        tb, lead = table('7', tt); pr = pred(cell, '7')
        d = sp.expand(tb - pr); print('NEG', cell, 'spin 7 vs table at t =', tt, ': mismatches', 0 if d == 0 else len(sp.Poly(d, X, Pi).terms()), '/ 15')
