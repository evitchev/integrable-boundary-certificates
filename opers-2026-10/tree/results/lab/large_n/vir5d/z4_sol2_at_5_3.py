"""VIR5d: at t = 5/3 express the certified Sol-2 spin-5 charge (table) in the kernel basis
{S = Sol-3 table charge (monic), E = prod_{j=-1}^{1}(PX + pi - j p)(PX - pi - j p), p^2 = 1/2}.  Exact."""
import json, sys
import sympy as sp
sys.dont_write_bytecode = True
ROOT = '<repo>/results/lab/'
t, PX, PI = sp.symbols('t PX PI')
tv = sp.Rational(5, 3)
rho2 = sp.Rational(2) / (tv**2 - 1)


def table(sol):
    tab = json.load(open(ROOT + 'vev/vev_sol%d_w6.json' % sol))
    expr = 0
    for key, val in tab['coefficients'].items():
        a_, b_ = [int(x.split('^')[1]) for x in key.split()]
        expr += sp.sympify(val, locals={'t': t}, rational=True).subs(t, tv) * (-PX**2) ** (a_ // 2) * (rho2 - PI**2) ** (b_ // 2)
    P = sp.Poly(sp.expand(expr), PX, PI)
    return sp.Poly(P.as_expr() / P.coeff_monomial(PX**6), PX, PI)


S3, S2 = table(3), table(2)
u, v = PX + PI, PX - PI
E = sp.Poly(sp.expand(u * (u**2 - sp.Rational(1, 2)) * v * (v**2 - sp.Rational(1, 2))), PX, PI)
al, be = sp.symbols('al be')
eqs = sp.Poly(S2.as_expr() - al * S3.as_expr() - be * E.as_expr(), PX, PI).coeffs()
sol = sp.solve(eqs, [al, be], dict=True)
print('Sol-2(5/3) monic =', S2.as_expr())
print('Sol-3(5/3) monic =', S3.as_expr())
print('E =', E.as_expr())
print('Sol-2 = al * Sol-3 + be * E with', sol)
assert len(sol) == 1 and sp.expand(S2.as_expr() - (al * S3.as_expr() + be * E.as_expr()).subs(sol[0])) == 0
print('exact: residual 0')
