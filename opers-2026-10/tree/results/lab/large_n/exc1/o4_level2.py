"""EXC1, level 2 (exploratory extension, labelled): TWO apparent singularities in the xi-plane.
(0) the coupled trivial-monodromy conditions re-derived by series at M = 5 and compared with the sealed system (v);
(1) I_1 = Delta + 2;
(2) number of solutions {z1, z2} (unordered) at a rational momentum point, and the minimal polynomial of the I_3 values
    Lambda_3(z1, z2) over the solutions (Groebner elimination, exact).
Usage: o4_level2.py <t> <alpha> <lam>   -> prints the elimination polynomial; saved to o4_t..json for comparison with data."""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_EXC1.md', 'rb').read()).hexdigest() == open('SEAL_EXC1.sha256').read().split()[0]
tstr, alv, lav = sys.argv[1], sp.Rational(sys.argv[2]), sp.Rational(sys.argv[3])
t = sp.Rational(tstr); M = (t + 3) / (t - 1); u = M + 1
l0, l1, z1, z2, e = sp.symbols('l0 l1 z1 z2 e')
t0 = time.time()


def cond(zk, zo, al, la):
    r = zk / (zk - zo)
    return 2 * M * zk**2 + al * (M - 1) * zk - 2 * la**2 + (M - 1)**2 / 2 - 2 * (u**3 * (r - r**2) * (1 - 2 * r) - 3 * u**2 * (r - r**2) + 2 * u * r)


if tstr == '2':
    # (0) series check of system (v) at M = 5
    x0, d, zz, al, la = sp.symbols('x0 d zz alpha lam_', positive=True)
    Mi = 5
    V0 = lambda xx: xx**(2 * Mi) + al * xx**(Mi - 1) + (la**2 - sp.Rational(1, 4)) / xx**2
    f1 = (x0 + d)**(Mi + 1) - x0**(Mi + 1)
    f2 = (x0 + d)**(Mi + 1) - zz
    dV = -2 * sp.diff(sp.log(f1), d, 2) - 2 * sp.diff(sp.log(f2), d, 2)
    reg = sp.series(sp.together(-2 * sp.diff(sp.log(f1), d, 2)) - 2 / d**2, d, 0, 2).removeO().coeff(d, 1) + sp.diff(-2 * sp.diff(sp.log(f2), d, 2), d).subs(d, 0)
    w1 = sp.diff(V0(x0 + d), d).subs(d, 0) + reg
    lhs = sp.simplify(w1 * x0**3)
    rhs = cond(x0**6, zz, al, la).subs(M, 5)
    print('(0) system (v) at M = 5 from the potential: %s' % (sp.simplify(lhs - rhs) == 0), flush=True)
od = json.load(open('o1_t%s_n2.json' % tstr.replace('/', '_')))
loc = {'Symbol': sp.Symbol}
R2 = sp.sympify(od['R']['2'], locals=loc); R4 = sp.sympify(od['R']['4'], locals=loc)
c2 = sp.Poly(R2, l0, l1, z1, z2).coeff_monomial(l1**2); c4 = sp.Poly(R4, l0, l1, z1, z2).coeff_monomial(l1**4)
I1c = sp.expand(R2 / c2 - (l0**2 / (4 * M) + l1**2 - (M + 1) / 6)) / (2 * (M + 1))
print('(1) two points: I_1 = Delta-part + 2 x (%s)  [level]' % I1c)
vac_lead = sp.Poly(sp.expand((R4 / c4).subs({z1: 0, z2: 0})), l0, l1).coeff_monomial(l0**4) * (4 * M * (M + 1))**2
L3 = sp.expand((R4 / c4 / vac_lead).subs({l0: alv, l1: lav}))
n1 = sp.numer(sp.together(cond(z1, z2, alv, lav))); n2 = sp.numer(sp.together(cond(z2, z1, alv, lav)))
w = sp.Symbol('w')
G = sp.groebner([n1, n2, sp.expand(e - L3), 1 - w * (z1 - z2)], w, z1, z2, e, order='lex')
elim = [g for g in G.exprs if g.free_symbols <= {e}]
pol = sp.Poly(elim[0], e)
print('(2) t = %s, alpha = %s, lam = %s: elimination polynomial for e = Lambda_3: degree %d; factor degrees %s  (%.0fs)'
      % (tstr, alv, lav, pol.degree(), [sp.Poly(f, e).degree() for f, m in sp.factor_list(pol.as_expr())[1]], time.time() - t0))
json.dump({'t': tstr, 'alpha': str(alv), 'lam': str(lav), 'elim': str(pol.as_expr())}, open('o4_t%s_a%s_l%s.json' % (tstr.replace('/', '_'), str(alv).replace('/', '_'), str(lav).replace('/', '_')), 'w'))
