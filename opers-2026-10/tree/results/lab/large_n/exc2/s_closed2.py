"""EXC2 (post-hoc, after the hold-out): the level-1 oper of Sol 3 in closed form -> level-1 blocks of I_3 and I_5 as EXACT
symbolic functions of the momenta (no interpolation), compared with the record's blocks (d1_sol3_t*_with5.json).
  u^2 = U root of  (3b-4) U^2 - 8(b-2)(2A+2B+b^2-2) U - 16 b ((b-2)^2-4A)((b-2)^2-4B) = 0,   A = a0^2 = b PX^2, B = a1^2 = b pi^2,
  z = (u/16)(U - 8(A+B) + 4(b-1)^2 + 4),  r11 = -s11 = b u/2,  r12 = (u+b+4) z,  s12 = (u-b-4) z, ...
Modes: real | tamper (constant of F times 11/10; must FAIL).  Exit: real 0 iff trace and det of both charges identical at both fibres."""
import json, sys
import sympy as sp
sys.dont_write_bytecode = True
mode = sys.argv[1] if len(sys.argv) > 1 else 'real'
A, B, b, rho, z, r11, u, U = sp.symbols('A B b rho z r11 u U')
d = json.load(open('s_sym.json'))
names = ['r01', 'r02', 's01', 's02']
loc = {'A': A, 'B': B, 'b': b, 'rho': rho, 'z': z, 'r11': r11}
loc.update({n: sp.Symbol(n) for n in names})
Zu = u * (u**2 - 8 * (A + B) + 4 * (b - 1)**2 + 4) / 16
cl = json.load(open('s_closed.json'))
assert sp.expand(sp.sympify(cl['z_of_u'], locals={'A': A, 'B': B, 'b': b, 'u': u}) - Zu) == 0
c1 = -8 * (b - 2) * (2 * A + 2 * B + b**2 - 2) / (3 * b - 4)
c0 = -16 * b * ((b - 2)**2 - 4 * A) * ((b - 2)**2 - 4 * B) / (3 * b - 4)
assert sp.simplify(sp.sympify(cl['F_of_u'], locals={'A': A, 'B': B, 'b': b, 'u': u}) - (u**4 + c1 * u**2 + c0)) == 0
if mode == 'tamper':
    c0 = c0 * sp.Rational(11, 10)
full = {r11: b * u / 2, rho: u + b + 4, z: Zu}
val = {n: sp.factor(sp.sympify(d['sol'][n], locals=loc).subs(full)) for n in names}
print('closed form (u = rho - b - 4):')
for n in names:
    print('   %s = %s' % (n, val[n]))
r01s01 = sp.expand(val['r01'] + val['s01'])
print('   r01 + s01 =', sp.factor(r01s01), ' [even in u: %s]' % (sp.expand(r01s01 - r01s01.subs(u, -u)) == 0))
q11 = sp.expand((b * u / 2) * Zu + (u + b + 4) * Zu + (-b * u / 2) * (-Zu) + (u - b - 4) * Zu)
print('   q11 = sum(r11 z + r12) =', sp.factor(q11))
toU = lambda ex: sp.Poly(sp.expand(ex), u)
def inU(ex):
    p = toU(ex)
    assert all(m[0] % 2 == 0 for m, _ in p.terms()), 'not even in u'
    return sum(c * U**(m[0] // 2) for m, c in p.terms())
q4U = inU(A * B + r01s01); q11U = inU(q11)
P, Q, lam = sp.symbols('P Q lam')
A2, A3, A4, B0, B1, B2 = sp.symbols('A2 A3 A4 B0 B1 B2')
ok = True
for tstr in ('2', '9/4'):
    t = sp.Rational(tstr); k = 2 * (t - 1) / (3 - t); bv = k / (k + 1); n = k + 4; nb = n / bv; rho2 = 2 / (t * t - 1)
    gc = json.load(open('u3_engine_t%s.json' % tstr.replace('/', '_')))
    locg = {'A2': A2, 'A3': A3, 'A4': A4, 'B0': B0, 'B1': B1, 'B2': B2}
    cB4 = sp.sympify(gc['cB4'], locals=locg); cB6 = sp.sympify(gc['cB6'], locals=locg)
    x0s, x1s = sp.symbols('x0s x1s')
    vacmap = {A2: nb**2 * (-(x0s + x1s)), A3: 0, A4: nb**4 * x0s * x1s, B1: 0, B2: 0}
    lead6 = sp.Poly(sp.expand(cB6.subs(vacmap)), x0s, x1s).coeff_monomial(x0s**3) * bv**3
    lead4 = sp.Poly(sp.expand(cB4.subs(vacmap)), x0s, x1s).coeff_monomial(x0s**2) * bv**2
    valmap = {A2: nb**2 * (-(A + B) - 2 * bv), A3: 0, A4: nb**4 * q4U.subs(b, bv), B2: 0, B1: -nb**6 * q11U.subs(b, bv)}
    mom = {A: bv * P**2, B: bv * (Q**2 + rho2)}
    C1 = sp.together(c1.subs(b, bv)); C0 = sp.together(c0.subs(b, bv))
    comp = sp.Matrix([[0, -C0], [1, -C1]])            # companion matrix of U^2 + C1 U + C0
    dd = json.load(open('d1_sol3_t%s_with5.json' % tstr.replace('/', '_')))
    for name, cB, lead, deg in (('I3', cB4, lead4, 4), ('I5', cB6, lead6, 6)):
        pol = sp.Poly(sp.expand(cB.subs(valmap) / lead), U)
        Mop = sp.zeros(2, 2)
        for (m,), c in pol.terms():
            Mop += c * comp**m
        tr = sp.factor(sp.together(Mop.trace()).subs(mom)); det = sp.factor(sp.together(Mop.det()).subs(mom))
        blk = dd['blocks'][name]
        Mx = sp.Matrix(2, 2, [sp.sympify(blk['matrix'][i][j], locals={'P': P, 'Q': Q}) for i in range(2) for j in range(2)])
        cN = sp.Poly(sp.sympify(blk['vacuum'], locals={'P': P, 'Q': Q}), P, Q).coeff_monomial(P**deg)
        Mn = Mx / cN
        dtr = sp.simplify(tr - Mn.trace()); ddet = sp.simplify(det - Mn.det())
        print('t = %-4s %s: degree in U %d; closed-form trace == data: %s; closed-form det == data: %s' % (tstr, name, pol.degree(), dtr == 0, ddet == 0))
        ok &= dtr == 0 and ddet == 0
print('s_closed2 [%s]: %s' % (mode, 'PASS' if ok else 'FAIL'))
if mode == 'real':
    sys.exit(0 if ok else 2)
sys.exit(0 if not ok else 1)
