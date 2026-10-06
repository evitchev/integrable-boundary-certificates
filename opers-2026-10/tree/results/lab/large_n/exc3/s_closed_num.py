"""Numeric check (30 digits) of the closed form at t = 2, (a0, a1) = (1/2, 1/3): for each root z of the quadratic, the residues
r01, r02, s01, s02 from s_closed_sol2.json, the charge I_3 through the engine mapping, and the relabelling z -> -z."""
import json, sys
import sympy as sp
sys.dont_write_bytecode = True
A, B, b, z = sp.symbols('A B b z'); r02, s02 = sp.symbols('r02 s02')
cf = json.load(open('s_closed_sol2.json')); loc = {'A': A, 'B': B, 'b': b, 'z': z, 'r02': r02, 's02': s02}
R02 = sp.sympify(cf['r02'], locals=loc); S02 = sp.sympify(cf['s02'], locals=loc); R01 = sp.sympify(cf['r01'], locals=loc).subs(r02, R02); S01 = sp.sympify(cf['s01'], locals=loc).subs(s02, S02)
c7, c5, c3 = [sp.sympify(cf[k_], locals=loc) for k_ in ('c7', 'c5', 'c3')]
t = sp.Integer(2); k = 2 * (t - 1) / (3 - t); bv = k / (k + 1); n = 2 * k + 3; c = n / bv; rho2 = 2 / (t * t - 1)
sub = {A: sp.Rational(1, 4), B: sp.Rational(1, 9), b: bv}
quad = sp.Poly((c7 * z**4 + c5 * z**2 + c3).subs(sub), z)
print('quadratic at the point:', sp.factor(quad.as_expr()))
roots = [sp.N(r_, 30) for r_ in sp.Poly(quad.as_expr(), z).nroots(n=30)]
gc = json.load(open('g_engine_t2.json'))
syms = sp.symbols('A2 A3 A4 B0 B1 B2 C0 C1 l0 Q30 Q31 Q40 Q41'); locg = {x.name: x for x in syms}
A2s, A3s, A4s, B0s, B1s, B2s, C0s, C1s, l0s, Q30s, Q31s, Q40s, Q41s = syms
cB4 = sp.sympify(gc['cB4'], locals=locg); cB2 = sp.sympify(gc['cB2'], locals=locg)
P, Q = sp.symbols('P Q')
vac = {A2s: 0, A3s: -c**2 * bv * (Q**2 + rho2), A4s: 0, B1s: 0, B0s: 0, C1s: 0, C0s: 0, B2s: 0, l0s: c * sp.sqrt(1 - bv) * P, Q30s: 0, Q31s: 0, Q40s: 0, Q41s: 0}
lead4 = sp.Poly(sp.expand(cB4.subs(vac)), P, Q).coeff_monomial(P**4); lead2 = sp.Poly(sp.expand(cB2.subs(vac)), P, Q).coeff_monomial(P**2)
def qs(r1, r2, o1, o2, o3, zz):
    d = dict(q1=r1, q0=o1)
    for m in (1, 2, 3, 4):
        d['q%d1' % m] = r1 * zz**m + m * r2 * zz**(m - 1); d['q%d0' % m] = o1 * zz**m + m * o2 * zz**(m - 1) + sp.Rational(m * (m - 1), 2) * o3 * zz**(m - 2)
    return d
for zr in roots:
    vals = {'r02': sp.N(R02.subs(sub).subs(z, zr), 30), 'r01': sp.N(R01.subs(sub).subs(z, zr), 30), 's02': sp.N(S02.subs(sub).subs(z, zr), 30), 's01': sp.N(S01.subs(sub).subs(z, zr), 30)}
    vm = {'r02': sp.N(R02.subs(sub).subs(z, -zr), 30), 'r01': sp.N(R01.subs(sub).subs(z, -zr), 30), 's02': sp.N(S02.subs(sub).subs(z, -zr), 30), 's01': sp.N(S01.subs(sub).subs(z, -zr), 30)}
    qa = qs(-bv, -3 * zr, vals['r01'], vals['r02'], 3 * zr**2, zr); qb = qs(-bv, 3 * zr, vals['s01'], vals['s02'], 3 * zr**2, -zr)
    val = {A2s: 0, A3s: c**2 * (-sp.Rational(1, 9) + qa['q1'] + qb['q1']), A4s: -c**3 * (qa['q0'] + qb['q0']), B2s: 0, l0s: c * sp.Rational(1, 2)}
    for m, (s1, s0) in zip((1, 2, 3, 4), ((B1s, B0s), (C1s, C0s), (Q31s, Q30s), (Q41s, Q40s))):
        val[s1] = c**(2 + m) * (qa['q%d1' % m] + qb['q%d1' % m]); val[s0] = -c**(3 + m) * (qa['q%d0' % m] + qb['q%d0' % m])
    I3 = sp.N(cB4.subs(val) / lead4, 20); I1 = sp.N(cB2.subs(val) / lead2, 20)
    print('z = %s: I_1 = %s, I_3 = %s' % (sp.N(zr, 12), I1, I3))
    print('    residues at z : %s' % {k_: sp.N(v, 10) for k_, v in vals.items()})
    print('    residues at -z: %s' % {k_: sp.N(v, 10) for k_, v in vm.items()})
print('data I_3 roots:', sp.Poly(sp.sympify('e**2 - 40*e/3 + 2092/9'), sp.Symbol('e')).nroots(n=12))
