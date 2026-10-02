"""Lead's BLIND spin-15 comparison of Fable's registered VIR5c predictions (sha256 81d4f9de...) with results/lab/anchor15/vev_profile_sol3_w16.json
('vev' "a,b" -> coefficient of (P_art^2)^a (Q_art^2)^b, rational in t); conversion P_art^2 = -PX^2, Q_art^2 = rho^2 - pi^2; monic in PX.
Neither the lead nor the seat had opened that table before this run.  Negative controls: shifted fibre.  Run from repo root."""
import json, sympy as sp, sys
t, X, Pi = sp.symbols('t X Pi')
P = json.load(open(sys.argv[1]))['cells']
tab = {tuple(map(int, k.split(','))): c for k, c in json.load(open('results/lab/anchor15/vev_profile_sol3_w16.json'))['vev'].items()}
def table(tk):
    p2 = 2*(tk-1)/(tk+1); rho2 = p2/4 - 1 + 1/p2; e = 0
    for (a, b), c in tab.items(): e += sp.sympify(c).subs(t, tk)*(-X)**a*(rho2 - Pi)**b
    e = sp.expand(e); lead = e.coeff(X, 8).subs(Pi, 0); return sp.expand(e/lead), lead
def pred(cell):
    sc = P[cell]['spins'].get('15')
    if sc is None or sc.get('coefficients') is None: return None
    return sum(sp.Rational(c)*X**(int(k.split(',')[0])//2)*Pi**(int(k.split(',')[1])//2) for k, c in sc['coefficients'].items())
def mism(d): return 0 if d == 0 else len(sp.Poly(d, X, Pi).terms())
for cell in sorted(P, key=lambda c: sp.Rational(c[2:])):
    k = sp.Rational(cell[2:]); tk = (3*k+2)/(k+2); pr = pred(cell)
    if pr is None: print(cell, 't =', tk, 'no spin-15 prediction'); continue
    tb, lead = table(tk)
    print(cell, 't =', tk, 'spin 15: mismatches', mism(sp.expand(tb - pr)), '/', len(sp.Poly(tb, X, Pi).terms()),
          '| NEG (t+1/100):', mism(sp.expand(table(tk + sp.Rational(1, 100))[0] - pr)), flush=True)
