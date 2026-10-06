import json, re, sys
import sympy as sp
sys.dont_write_bytecode = True
A, B, b, z, U = sp.symbols('A B b z U'); r01, r02, s01, s02 = sp.symbols('r01 r02 s01 s02')
names = {n: sp.Symbol(n) for n in ('r11', 'r01', 'r02', 's01', 's02', 'z', 'A', 'B', 'b')}
d = json.load(open('s_sym_sol2.json')); rest = [sp.sympify(c, locals=names) for c in d['rest']]
cf = json.load(open('s_closed_sol2.json')); loc = {'A': A, 'B': B, 'b': b, 'z': z, 'r02': r02, 's02': s02}
R02 = sp.sympify(cf['r02'], locals=loc); S02 = sp.sympify(cf['s02'], locals=loc); R01 = sp.sympify(cf['r01'], locals=loc); S01 = sp.sympify(cf['s01'], locals=loc)
quad = sp.sympify(cf['c7'], locals=loc) * z**4 + sp.sympify(cf['c5'], locals=loc) * z**2 + sp.sympify(cf['c3'], locals=loc)
sub = {A: sp.Rational(1, 4), B: sp.Rational(1, 9), b: sp.Rational(2, 3)}
# exact numeric solutions of the six equations at this point
eqs = [sp.expand(c.subs(sub).subs(names['r11'], -sp.Rational(2, 3))) for c in rest]
w = sp.Symbol('w')
G = sp.groebner(eqs + [1 - w * z], w, s01, s02, r01, r02, z, order='lex')
print('numeric Groebner basis at (1/4,1/9,2/3):')
for g in G.exprs: print('   ', g)
qn = sp.Poly(quad.subs(sub), z); print('closed-form quadratic at the point:', sp.factor(qn.as_expr()))
roots = sp.solve(qn.as_expr(), z)
for zr in roots[:2]:
    vals = {k_: sp.nsimplify(0) if False else sp.simplify(v.subs(sub).subs(z, zr)) for k_, v in (('r02', R02), ('r01', R01.subs(r02, R02)), ('s02', S02), ('s01', S01.subs(s02, S02)))}
    resid = [sp.simplify(c.subs(sub).subs(names['r11'], -sp.Rational(2, 3)).subs({r02: vals['r02'], r01: vals['r01'], s02: vals['s02'], s01: vals['s01'], z: zr})) for c in rest]
    print('z =', zr, ' closed-form values:', {k_: sp.nsimplify(v) if v.is_Rational else v for k_, v in vals.items()}, ' residuals:', resid)
