"""EXC2 T2: the solutions of the trivial-monodromy system (f_sol3.py, exponents {-1,1,2,4}) against the stored level-1
I_3 block of Sol 3.  NOTE (correction to the seal's formula (i)): the xi-exponents are a0 = p PX, a1 = p pi (not PX, pi).
Level-1 I_1, I_3 depend only on the effective quartic at xi -> infinity:
    th^4 + (-(a0^2 + a1^2) + r21) th^2 + r11 th + (a0^2 a1^2 + r01),   T = -(n/b) th.
Usage: g_compare.py <t> <a0> <a1> [tamper]     Exit 0 = I_1 = Delta + 1 and I_3 polynomial == data; 2 otherwise."""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_EXC2.md', 'rb').read()).hexdigest() == open('SEAL_EXC2.sha256').read().split()[0]
tstr = sys.argv[1]; a0 = sp.Rational(sys.argv[2]); a1 = sp.Rational(sys.argv[3]); tamper = len(sys.argv) > 4
t = sp.Rational(tstr); k = 2 * (t - 1) / (3 - t); b = k / (k + 1); n = k + 4; rho2 = 2 / (t * t - 1)
PX2, PI2 = a0**2 / b, a1**2 / b
t0 = time.time()
fs = json.load(open('f_sol3_t%s_e-1_1_2_4%s.json' % (tstr.replace('/', '_'), '_tamper' if tamper else '')))
assert sp.Rational(fs['PX']) == a0 and sp.Rational(fs['PI']) == a1
syms = {s: sp.Symbol(s) for s in ('r21', 'r11', 'r12', 'r01', 'r02', 'r03', 'z', 'w')}
basis = [sp.sympify(g, locals=syms) for g in fs['basis']]
gc = json.load(open('g_charges_t%s.json' % tstr.replace('/', '_')))
A2, A3, A4, e = sp.symbols('A2 A3 A4 e')
cB2 = sp.sympify(gc['cB2'], locals={'A2': A2, 'A3': A3, 'A4': A4}); cB4 = sp.sympify(gc['cB4'], locals={'A2': A2, 'A3': A3, 'A4': A4})
sc = n / b
def at(q2, q3, q4, expr):
    return sp.expand(expr.subs({A2: sc**2 * q2, A3: -sc**3 * q3, A4: sc**4 * q4}))
# vacuum normalisations (PX^2 = a0^2/b)
vac2 = at(-(a0**2 + a1**2), 0, a0**2 * a1**2, cB2); vac4 = at(-(a0**2 + a1**2), 0, a0**2 * a1**2, cB4)
x0s, x1s = sp.symbols('x0s x1s')            # symbolic a0^2, a1^2 to read off the PX^4 coefficient
lead4 = sp.Poly(at(-(x0s + x1s), 0, x0s * x1s, cB4), x0s, x1s).coeff_monomial(x0s**2) * b**2       # coefficient of PX^4
lead2 = sp.Poly(at(-(x0s + x1s), 0, x0s * x1s, cB2), x0s, x1s).coeff_monomial(x0s) * b            # coefficient of PX^2
r21, r11, r01, z, w = syms['r21'], syms['r11'], syms['r01'], syms['z'], syms['w']
I1 = at(-(a0**2 + a1**2) + r21, r11, a0**2 * a1**2 + r01, cB2) / lead2
L3 = at(-(a0**2 + a1**2) + r21, r11, a0**2 * a1**2 + r01, cB4) / lead4
others = [syms[s] for s in ('r12', 'r02', 'r03')]
G = sp.groebner(basis + [sp.expand(e - L3)], w, r21, r11, *others, r01, z, e, order='lex')
el = [g for g in G.exprs if g.free_symbols <= {e}]
pol = sp.Poly(el[0], e); pol = sp.Poly(pol.as_expr() / pol.LC(), e)
Gz = sp.groebner(basis, w, r21, r11, *others, r01, z, order='lex')
elz = [g for g in Gz.exprs if g.free_symbols <= {z}]
print('t = %s, xi-exponents (%s, %s) [PX^2 = %s, pi^2 = %s]%s' % (tstr, a0, a1, PX2, PI2, ' TAMPER' if tamper else ''))
print('number of solutions (degree of the eliminant in z): %d;  eliminant: %s' % (sp.Poly(elz[0], z).degree(), sp.factor(elz[0])))
# I_1: constant on the solution set?
G1 = sp.groebner(basis + [sp.expand(e - I1)], w, r21, r11, *others, r01, z, e, order='lex')
el1 = [g for g in G1.exprs if g.free_symbols <= {e}]
I1vals = sp.solve(el1[0], e)
print('I_1 (monic) over the solutions: %s;  vacuum value %s;  level = %s' % (I1vals, vac2 / lead2, [(v - vac2 / lead2) / 2 for v in I1vals]))
print('I_3: minimal polynomial over the solutions: degree %d: %s' % (pol.degree(), pol.as_expr()))
# data
P, Q = sp.symbols('P Q')
dd = json.load(open('d1_sol3_t%s.json' % tstr.replace('/', '_')))['blocks']['I3']
Mx = sp.Matrix(2, 2, [sp.sympify(dd['matrix'][i][j], locals={'P': P, 'Q': Q}) for i in range(2) for j in range(2)])
vac = sp.Poly(sp.sympify(dd['vacuum'], locals={'P': P, 'Q': Q}), P, Q)
cp = sp.expand((Mx / vac.coeff_monomial(P**4) - e * sp.eye(2)).det())
cpP = sp.Poly(cp, P, Q)
assert all(i % 2 == 0 and j % 2 == 0 for (i, j), _ in cpP.terms()), 'characteristic polynomial not even in P, Q'
cpv = sp.Poly(sum(c * PX2 ** (i // 2) * (PI2 - rho2) ** (j // 2) for (i, j), c in cpP.terms()), e)
cpv = sp.Poly(sp.expand(cpv.as_expr() / cpv.LC()), e)
print('data: level-1 characteristic polynomial at this point: %s' % cpv.as_expr())
same = sp.expand(cpv.as_expr() - pol.as_expr()) == 0
divides = pol.degree() >= 2 and sp.rem(pol, cpv, e).is_zero
lvl_ok = all((v - vac2 / lead2) / 2 == 1 for v in I1vals)
print('I_1 = Delta + 1 on every solution: %s;   I_3 polynomial == data characteristic polynomial: %s;  data polynomial divides it: %s  (%.0fs)' % (lvl_ok, same, divides, time.time() - t0))
ok = lvl_ok and same
if tamper:
    sys.exit(0 if not ok else 1)
sys.exit(0 if ok else 2)
