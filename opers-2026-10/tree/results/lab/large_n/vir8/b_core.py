"""helpers shared by the VIR8 (b) scripts (copied from b_nz.py)"""
import json
import sympy as sp
import wkb3 as W3
ROOT = '<repo>/results/lab/'
l0, l1 = W3.l0, W3.l1
ts, PX, PI = sp.symbols('t PX PI')


def table(sol, s, tv):
    rho2 = sp.Rational(2) / (tv**2 - 1)
    tab = json.load(open(ROOT + 'vev/vev_sol%d_w%d.json' % (sol, s + 1)))
    expr = 0
    for key, val in tab['coefficients'].items():
        a_, b_ = [int(x.split('^')[1]) for x in key.split()]
        f = sp.sympify(val, locals={'t': ts}, rational=True)
        assert sp.fraction(sp.together(f))[1].subs(ts, tv) != 0
        expr += f.subs(ts, tv) * (-PX**2) ** (a_ // 2) * (rho2 - PI**2) ** (b_ // 2)
    P = sp.Poly(sp.expand(expr), PX, PI)
    return sp.Poly(P.as_expr() / P.coeff_monomial(PX ** (s + 1)), PX, PI)


def to_PXPI(Pp, S0, S1):
    Pp = sp.Poly(Pp, l0, l1)
    out = sp.Poly(sp.expand(sum(c * S0 ** sp.Rational(i, 2) * S1 ** sp.Rational(j, 2) * PX**i * PI**j for (i, j), c in Pp.terms())), PX, PI)
    lead = out.coeff_monomial(PX ** out.total_degree())
    return None if lead == 0 else sp.Poly(sp.expand(out.as_expr() / lead), PX, PI)


