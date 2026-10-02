"""EXC2 (post-hoc, after the hold-out): closed form of the level-1 oper of Sol 3 from s_sym.json.
Unknowns after the linear eliminations: r11, rho (r12 = rho z), z.  Symbolic in A = a0^2, B = a1^2, b.
Exit 0 iff (i) the final system is r11 = b(rho-b-4)/2, z = cubic(rho), F(u^2) = 0 with u = rho - (b+4), F quadratic in u^2;
(ii) at every batch point of both fibres the z-eliminant of the numeric two-point computation is reproduced."""
import json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
A, B, b, rho, z, r11, u, U2 = sp.symbols('A B b rho z r11 u U2')
d = json.load(open('s_sym.json'))
loc = {'A': A, 'B': B, 'b': b, 'rho': rho, 'z': z, 'r11': r11}
rest = [sp.sympify(c, locals=loc) for c in d['rest']]
lin = [c for c in rest if sp.Poly(c, r11, rho, z).total_degree() == 1]
assert len(lin) == 1
R11 = sp.solve(lin[0], r11)[0]
print('r11 =', sp.factor(R11))
rest2 = [sp.expand(c.subs(r11, R11)) for c in rest if c is not lin[0]]
cub = [c for c in rest2 if sp.Poly(c, rho).degree() == 3]
Z = sp.solve(cub[0], z)[0]
assert all(sp.expand(c.subs(z, Z)) == 0 for c in cub), 'the two cubic equations are not equivalent'
Zu = sp.expand(Z.subs(rho, u + b + 4))
print('z as a function of u = rho - (b+4):', sp.collect(Zu, u))
quint = [sp.expand(c.subs(z, Z).subs(rho, u + b + 4)) for c in rest2 if c not in cub]
polys = [sp.Poly(c, u) for c in quint]
print('degrees in u of the two remaining equations:', [p.degree() for p in polys])
g = sp.gcd(polys[0], polys[1])
print('gcd degree in u:', g.degree())
F = sp.Poly(sp.factor(g.as_expr()), u)
F = sp.Poly(sp.expand(F.as_expr() / F.LC()), u)
even = all(m[0] % 2 == 0 for m, _ in F.terms())
print('F(u) monic, even in u: %s; F ='% even, sp.collect(F.as_expr(), u, sp.factor))
odd_z = sp.expand(Zu + Zu.subs(u, -u)) == 0
print('z odd in u (so the relabelling u -> -u is z -> -z): %s' % odd_z)
quot = [sp.div(p, g)[0].as_expr() for p in polys]
print('cofactors of the two equations:', [sp.factor(q) for q in quot])
# z^2 eliminant: resultant of F(u) and zeta - z(u)^2
ok = even and odd_z and F.degree() == 4
json.dump({'r11': str(sp.factor(R11)), 'z_of_u': str(sp.collect(Zu, u)), 'F_of_u': str(sp.collect(F.as_expr(), u, sp.factor))}, open('s_closed.json', 'w'), indent=1)
# check against the numeric batch: z-eliminant (in z) at every point
zz = sp.Symbol('zz')
nchk = bad = 0
for tt, tstr in (('2', '2'), ('9_4', '9/4')):
    t = sp.Rational(tstr); k = 2 * (t - 1) / (3 - t); bv = k / (k + 1)
    for line in open('u3_points_t%s.txt' % tt):
        _, a0s, a1s, _ = line.split()
        p = json.load(open('u3p_t%s_%s_%s.json' % (tt, a0s.replace('/', '-'), a1s.replace('/', '-'))))
        num = sp.Poly(sp.sympify(p['z_eliminant'].replace('^', '**'), locals={'z': zz}), zz)
        sub = {A: sp.Rational(a0s)**2, B: sp.Rational(a1s)**2, b: bv}
        Fv = sp.Poly(F.as_expr().subs(sub), u); Zv = sp.Poly(Zu.subs(sub), u)
        res = sp.Poly(sp.resultant(Fv.as_expr(), zz - Zv.as_expr(), u), zz)
        nchk += 1
        same = sp.expand(res.as_expr() / res.LC() - num.as_expr() / num.LC()) == 0
        bad += not same
print('closed form against the numeric two-point batch (z-eliminant, 2 fibres): %d points, %d mismatches' % (nchk, bad))
ok &= bad == 0 and nchk == 72
print('s_closed', 'PASS' if ok else 'FAIL')
sys.exit(0 if ok else 2)
