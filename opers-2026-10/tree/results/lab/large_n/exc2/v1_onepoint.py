"""EXC2 control V1 (no data): the I_5 map (effective quartic + first-order 1/xi term, wkb3x) on the ONE-point opers.
The one-point I_3 values are vacuum values at half-shifted exponents; if the one-point oper is a transform of the vacuum
oper at those exponents, its I_5 must be the vacuum I_5 there.  This tests the B1, B2 map (signs, scales) used in U3.
Usage: v1_onepoint.py <t> <a0> <a1> [tamper]   (tamper: B1 sign flipped -> must fail).  Exit 0 = every one-point (I_3, I_5)
pair is a vacuum pair at a half-shifted exponent pair; 2 otherwise."""
import hashlib, itertools, json, sys
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_EXC2_addendum.md', 'rb').read()).hexdigest() == open('SEAL_EXC2_addendum.sha256').read().split()[0]
tstr = sys.argv[1]; a0 = sp.Rational(sys.argv[2]); a1 = sp.Rational(sys.argv[3]); tamper = len(sys.argv) > 4
t = sp.Rational(tstr); k = 2 * (t - 1) / (3 - t); b = k / (k + 1); n = k + 4; nb = n / b
fs = json.load(open('f_sol3_t%s_e-1_1_2_4.json' % tstr.replace('/', '_')))
assert sp.Rational(fs['PX']) == a0 and sp.Rational(fs['PI']) == a1
syms = {s: sp.Symbol(s) for s in ('r21', 'r11', 'r12', 'r01', 'r02', 'r03', 'z', 'w')}
basis = [sp.sympify(g, locals=syms) for g in fs['basis']]
A2, A3, A4, B0, B1, B2, e, f = sp.symbols('A2 A3 A4 B0 B1 B2 e f')
gc = json.load(open('u3_engine_t%s.json' % tstr.replace('/', '_')))
loc = {'A2': A2, 'A3': A3, 'A4': A4, 'B0': B0, 'B1': B1, 'B2': B2}
cB4 = sp.sympify(gc['cB4'], locals=loc); cB6 = sp.sympify(gc['cB6'], locals=loc)
x0s, x1s = sp.symbols('x0s x1s')
vacmap = {A2: nb**2 * (-(x0s + x1s)), A3: 0, A4: nb**4 * x0s * x1s, B1: 0, B2: 0}
v6 = sp.expand(cB6.subs(vacmap)); v4 = sp.expand(cB4.subs(vacmap))
lead6 = sp.Poly(v6, x0s, x1s).coeff_monomial(x0s**3) * b**3; lead4 = sp.Poly(v4, x0s, x1s).coeff_monomial(x0s**2) * b**2
r21, r11, r12, r01, z = [syms[s] for s in ('r21', 'r11', 'r12', 'r01', 'z')]
sgn = -1 if tamper else 1
valmap = {A2: nb**2 * (-(a0**2 + a1**2) + r21), A3: -nb**3 * r11, A4: nb**4 * (a0**2 * a1**2 + r01),
          B2: nb**5 * (r21 * z - 4 * z), B1: -sgn * nb**6 * (r11 * z + r12)}
I3 = sp.expand(cB4.subs(valmap) / lead4); I5 = sp.expand(cB6.subs(valmap) / lead6)
others = [syms[s] for s in ('r12', 'r02', 'r03')]
gens = [syms['w'], r21, r11] + others + [r01, z]
sols = sp.solve_poly_system(basis, *gens)
print('t = %s, xi-exponents (%s, %s)%s: %d one-point solutions' % (tstr, a0, a1, ' TAMPER (B1 sign)' if tamper else '', len(sols)))

cands = {}
for d0, d1 in itertools.product((sp.Rational(1, 2), -sp.Rational(1, 2)), repeat=2):
    cands[(a0 + d0, a1 + d1)] = (sp.Rational(v4.subs({x0s: (a0 + d0)**2, x1s: (a1 + d1)**2})) / lead4, sp.Rational(v6.subs({x0s: (a0 + d0)**2, x1s: (a1 + d1)**2})) / lead6)
ok = True
for s in sols:
    sub = dict(zip(gens, s))
    i3 = sp.simplify(I3.subs(sub)); i5 = sp.simplify(I5.subs(sub))
    hit = [key for key, (c3, c5) in cands.items() if sp.simplify(c3 - i3) == 0 and sp.simplify(c5 - i5) == 0]
    hit3 = [key for key, (c3, c5) in cands.items() if sp.simplify(c3 - i3) == 0]
    print('  z = %-22s I_3 = %-10s I_5 = %-14s vacuum pair at shifted exponents: %s   (I_3 alone: %s)' % (sub[z], i3, i5, hit, hit3))
    ok &= bool(hit)
print('vacuum (I_3, I_5) at the four half-shifted exponent pairs:', cands)
print('V1 %s' % ('PASS' if ok else 'FAIL'))
if tamper:
    sys.exit(0 if not ok else 1)
sys.exit(0 if ok else 2)
