"""Lead's independent spin-11 comparison of cc's registered SOL2F predictions (sha256 48f1e50b...) with results/lab/anchor11/vev_profile_sol2_w12.json
evaluated at t = 7/3, 11/5, 9/4, 5/3 (never evaluated off t = 2 before registration).  Projective comparison (raw prediction vs table, both
normalised on the first nonzero of P_X^12 / leading monomial); conversion P_art^2 = -PX^2, Q_art^2 = rho^2 - pi^2.  Run from repo root."""
import json, sympy as sp, sys
t, X, Pi = sp.symbols('t X Pi'); PX, pq = sp.symbols('P_X piq')
d = json.load(open(sys.argv[1]))
tab = {tuple(map(int, k.split(','))): c for k, c in json.load(open('results/lab/anchor11/vev_profile_sol2_w12.json'))['vev'].items()}
def table(tk):
    p2 = 2*(tk-1)/(tk+1); rho2 = p2/4 - 1 + 1/p2; e = 0
    for (a, b), c in tab.items(): e += sp.sympify(c).subs(t, tk)*(-X)**a*(rho2 - Pi)**b
    return sp.expand(e)
def norm(e):
    p = sp.Poly(e, X, Pi); lead = p.coeff_monomial(X**6)
    if lead == 0: lead = p.LC()
    return sp.expand(e/lead)
for cell in [c for c in d if c != 'note']:
    tk = sp.Rational(cell.split('t=')[1]); tb = table(tk)
    for form, v in d[cell].items():
        raw = v.get('raw') if isinstance(v, dict) else v
        if raw in (None, '0', 0):
            print(cell, form, 'predicted spin 11 = 0; certified table at t:', 'ZERO' if tb == 0 else f'NONZERO ({len(sp.Poly(tb, X, Pi).terms())} terms)'); continue
        pr = sp.expand(sp.sympify(raw, locals={'pi': pq, 'P_X': PX}).subs({PX: sp.sqrt(X), pq: sp.sqrt(Pi)}))
        dd = sp.expand(norm(tb) - norm(pr)); m = 0 if dd == 0 else len(sp.Poly(dd, X, Pi).terms())
        dneg = sp.expand(norm(table(tk + sp.Rational(1, 100))) - norm(pr))
        print(cell, form, 'spin 11 mismatches', m, '/', len(sp.Poly(tb, X, Pi).terms()), '| NEG t+1/100:', len(sp.Poly(dneg, X, Pi).terms()))
