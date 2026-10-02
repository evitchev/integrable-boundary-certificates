"""EXC1 level 2 (extension beyond the seal; labelled): TWO apparent singularities (sealed system (v)) against the level-2 data.
At t = 2, bare (P, Q) = (1/2, 1): alpha^2 = 4M(M+1) P^2 = 30, lam^2 = (M+1)(Q^2 + rho^2) = 10.
Exact elimination (Groebner, alpha algebraic) of e = Lambda_3(z1, z2) over the solutions with z1 != z2; compared with the
quintic factor of the characteristic polynomial of the 6x6 data block (d2_level2.py).  Exit 0 = equal; 2 = not."""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_EXC1.md', 'rb').read()).hexdigest() == open('SEAL_EXC1.sha256').read().split()[0]
t = sp.Integer(2); M = (t + 3) / (t - 1); u = M + 1; rho2 = 2 / (t * t - 1)
Pv, Qv = sp.Rational(1, 2), sp.Integer(1)
al2, la2 = 4 * M * (M + 1) * Pv**2, (M + 1) * (Qv**2 + rho2)
l0, l1, z1, z2, e, w = sp.symbols('l0 l1 z1 z2 e w')
t0 = time.time()
od = json.load(open('o1_t2_n2.json'))
R4 = sp.sympify(od['R']['4'], locals={'Symbol': sp.Symbol})
c4 = sp.Poly(R4, l0, l1, z1, z2).coeff_monomial(l1**4)
vac_lead = sp.Poly(sp.expand((R4 / c4).subs({z1: 0, z2: 0})), l0, l1).coeff_monomial(l0**4) * (4 * M * (M + 1))**2
L3 = sp.expand((R4 / c4 / vac_lead).subs(l1**2, la2).subs(l1**4, la2**2))
L3 = sp.expand(sp.Poly(L3, l1).as_expr().subs(l1, sp.sqrt(la2)))


def cond(zk, zo):
    r = zk / (zk - zo)
    return 2 * M * zk**2 + l0 * (M - 1) * zk - 2 * la2 + (M - 1)**2 / 2 - 2 * (u**3 * (r - r**2) * (1 - 2 * r) - 3 * u**2 * (r - r**2) + 2 * u * r)


n1 = sp.numer(sp.together(cond(z1, z2))); n2 = sp.numer(sp.together(cond(z2, z1)))
G = sp.groebner([n1, n2, sp.expand(e - L3), 1 - w * (z1 - z2), l0**2 - al2], w, z1, z2, l0, e, order='lex')
elim = [g for g in G.exprs if g.free_symbols <= {e}]
pol = sp.Poly(elim[0], e)
pol = sp.Poly(pol.as_expr() / pol.LC(), e)
print('oper: elimination polynomial for I_3 at level 2: degree %d, factor degrees %s  (%.0fs)' % (pol.degree(), [sp.Poly(f, e).degree() for f, m in sp.factor_list(pol.as_expr())[1]], time.time() - t0))
dd = json.load(open('d2_sol1_t2_P1_2_Q1.json'))
facs = [sp.Poly(sp.sympify(f, locals={'e': e}), e) for f in dd['factors']]
quint = [f for f in facs if f.degree() == 5][0]
quint = sp.Poly(quint.as_expr() / quint.LC(), e)
same = sp.expand(quint.as_expr() - pol.as_expr()) == 0
print('data: quintic factor of the 6x6 level-2 block;  oper polynomial == data quintic: %s' % same)
print('oper :', pol.as_expr())
print('data :', quint.as_expr())
sys.exit(0 if same else 2)
