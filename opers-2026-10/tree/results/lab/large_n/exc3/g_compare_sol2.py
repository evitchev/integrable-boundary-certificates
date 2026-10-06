"""EXC3 T2 (one point): I_1 (level) and I_3 of the one-point solutions of f_sol2 vs the stored Sol 2 level-one I_3 block.
Dictionary (seal sec. 1; checked at T0a): q1 = -a1^2 + r11, q0 = r01; first order q11 = r11 z + r12, q10 = r01 z + r02;
second order q21 = r11 z^2 + 2 r12 z, q20 = r01 z^2 + 2 r02 z + r03;  A3 = c^2 q1, A4 = -c^3 q0, B1 = c^3 q11, B0 = -c^4 q10,
C1 = c^4 q21, C0 = -c^5 q20, l0^2 = c^2 a0^2.  Gamma atoms: none arise (rational data).  Usage: g_compare_sol2.py <t> <exps> <a0> <a1> [tamper]
Exit 0 = I_1 = Delta + 1 (level one) on every solution and I_3 polynomial == data; 2 otherwise."""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_EXC3.md', 'rb').read()).hexdigest() == open('SEAL_EXC3.sha256').read().split()[0]
tstr, exs = sys.argv[1], sys.argv[2]; a0 = sp.Rational(sys.argv[3]); a1 = sp.Rational(sys.argv[4]); tamper = len(sys.argv) > 5
t = sp.Rational(tstr); k = 2 * (t - 1) / (3 - t); b = k / (k + 1); n = 2 * k + 3; c = n / b; rho2 = 2 / (t * t - 1)
PX2, PI2 = a0**2 / (1 - b), a1**2 / b
tt = tstr.replace('/', '_')
fs = json.load(open('f_sol2_t%s_e%s%s.json' % (tt, exs.replace(',', '_'), '_tamper' if tamper else '')))
assert sp.Rational(fs['a0']) == a0 and sp.Rational(fs['a1']) == a1
sy = {s: sp.Symbol(s) for s in ('r11', 'r12', 'r01', 'r02', 'r03', 'z', 'w')}
basis = [sp.sympify(g, locals=sy) for g in fs['basis']]
ind = {sy[k_]: sp.sympify(v, locals=sy) for k_, v in fs['indicial'].items()}
gc = json.load(open('g_engine_t%s.json' % tt))
A2, A3, A4, B0, B1, B2, C0, C1, l0 = sp.symbols('A2 A3 A4 B0 B1 B2 C0 C1 l0')
loc = {s.name: s for s in (A2, A3, A4, B0, B1, B2, C0, C1, l0)}
cB = {i: sp.sympify(gc['cB%d' % i], locals=loc) for i in (2, 3, 4, 5, 6)}
r11, r12, r01, r02, r03, z, w = [sy[s] for s in ('r11', 'r12', 'r01', 'r02', 'r03', 'z', 'w')]
q1 = -a1**2 + r11; q0 = r01
q11 = r11 * z + r12; q10 = r01 * z + r02
q21 = r11 * z**2 + 2 * r12 * z; q20 = r01 * z**2 + 2 * r02 * z + r03
val = {A2: 0, A3: c**2 * q1, A4: -c**3 * q0, B1: c**3 * q11, B0: -c**4 * q10, C1: c**4 * q21, C0: -c**5 * q20, B2: 0, l0: c * a0}
val = {k_: (v.subs(ind) if hasattr(v, 'subs') else v) for k_, v in val.items()}
P, Q = sp.symbols('P Q')
vac = {A2: 0, A3: -c**2 * b * (Q**2 + rho2), A4: 0, B1: 0, B0: 0, C1: 0, C0: 0, B2: 0, l0: c * sp.sqrt((1 - b)) * P}
lead = {i: sp.Poly(sp.expand(cB[i].subs(vac)), P, Q).coeff_monomial(P**i) for i in (2, 4, 6)}
vacv = {i: sp.expand(cB[i].subs(vac).subs({P: sp.sqrt(PX2), Q: sp.sqrt(PI2 - rho2)}) / lead[i]) for i in (2, 4, 6)}
I1 = sp.expand(cB[2].subs(val) / lead[2]); I3 = sp.expand(cB[4].subs(val) / lead[4])
J3 = sp.expand(cB[3].subs(val)); J5 = sp.expand(cB[5].subs(val))
e = sp.Symbol('e')
def elim(expr):
    G = sp.groebner(basis + [sp.expand(e - expr)], w, r11, r01, r02, z, e, order='lex')
    el = [g for g in G.exprs if g.free_symbols <= {e}]
    p_ = sp.Poly(el[0], e); return sp.Poly(p_.as_expr() / p_.LC(), e)
Gz = sp.groebner(basis, w, r11, r01, r02, z, order='lex')
elz = [g for g in Gz.exprs if g.free_symbols <= {z}][0]
print('t = %s, exponents %s, (a0, a1) = (%s, %s)%s: PX^2 = %s, pi^2 = %s' % (tstr, exs, a0, a1, ' TAMPER' if tamper else '', PX2, PI2))
print('solutions: eliminant in z: %s (degree %d)' % (sp.factor(elz), sp.Poly(elz, z).degree()))
p1 = elim(I1); print('I_1 (monic) over the solutions: roots %s; vacuum %s; level = %s' % (sp.solve(p1.as_expr(), e), vacv[2], [(v - vacv[2]) / 2 for v in sp.solve(p1.as_expr(), e)]))
p3 = elim(I3); print('I_3 (monic) minimal polynomial over the solutions: %s' % p3.as_expr())
pj3 = elim(J3); pj5 = elim(J5); print('even-spin charges J_3, J_5 over the solutions: minimal polynomials %s ; %s' % (pj3.as_expr(), pj5.as_expr()))
dd = json.load(open('d1_sol2_t%s.json' % tt))['blocks']['I3']
Mx = sp.Matrix(2, 2, [sp.sympify(dd['matrix'][i][j], locals={'P': P, 'Q': Q}) for i in range(2) for j in range(2)])
vacd = sp.Poly(sp.sympify(dd['vacuum'], locals={'P': P, 'Q': Q}), P, Q)
cp = sp.Poly(sp.expand((Mx / vacd.coeff_monomial(P**4) - e * sp.eye(2)).det()), P, Q, e)
assert all(m[0] % 2 == 0 and m[1] % 2 == 0 for m, _ in cp.terms())
cpv = sp.Poly(sum(cf * PX2 ** (m[0] // 2) * (PI2 - rho2) ** (m[1] // 2) * e ** m[2] for m, cf in cp.terms()), e)
cpv = sp.Poly(sp.expand(cpv.as_expr() / cpv.LC()), e)
print('data: level-one I_3 characteristic polynomial at this point: %s ; roots %s' % (cpv.as_expr(), [sp.nsimplify(0) if False else r_ for r_ in sp.Poly(cpv, e).nroots(n=12)]))
print('      I_3 values of the solutions (numeric): %s' % [r_ for r_ in p3.nroots(n=12)])
lvl_ok = all((v - vacv[2]) / 2 == 1 for v in sp.solve(p1.as_expr(), e))
same = sp.expand(cpv.as_expr() - p3.as_expr()) == 0
divides = p3.degree() >= cpv.degree() and sp.rem(p3, cpv, e).is_zero
print('I_1 = Delta + 1 on every solution: %s;  I_3 polynomial == data: %s;  data divides: %s' % (lvl_ok, same, divides))
ok = lvl_ok and same
sys.exit((0 if not ok else 1) if tamper else (0 if ok else 2))
