"""EXC3 (post-hoc, no I_5 data): closed form of the sigma-invariant two-point system of Solution 2 from the parametric Groebner basis
(s_closed_sol2.out, ring (0,A,B,b),(s01,s02,r01,r02,z),lp):
  G[1] = z^3 (c7 z^4 + c5 z^2 + c3)   -> U = z^2 solves  c7 U^2 + c5 U + c3 = 0  (two unordered solutions = two level-one states)
  G[2]: r02 = -(poly in z)/(D z^2);  G[4]: r01 = 2 b r02/(3z) - (3-b) z/3 - b(b+3)/3;  G[5], G[8]: s02, s01;  r11 = s11 = -b.
Checks: (i) the quadratic reproduces the z-eliminant of the numeric batch at all 70 points (two fibres);
        (ii) I_3 and I_5 trace/det from the closed form (companion matrix in U) equal the interpolated polynomials (u_interp) and, for I_3,
             the data block; (iii) the exceptional points (5/8,1/6), (5/6,8/9) at t = 2: which coefficient vanishes.
Exit 0 iff (i) and (ii) hold."""
import hashlib, json, re, sys
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_EXC3.md', 'rb').read()).hexdigest() == open('SEAL_EXC3.sha256').read().split()[0]
A, B, b, z, U = sp.symbols('A B b z U')
r01, r02, s01, s02 = sp.symbols('r01 r02 s01 s02')
loc = {'A': A, 'B': B, 'b': b, 'z': z, 'r01': r01, 'r02': r02, 's01': s01, 's02': s02}
out = open('s_closed_sol2.out').read()
G = {int(m.group(1)): sp.sympify(m.group(2).replace('^', '**'), locals=loc) for m in re.finditer(r'G\[(\d+)\]=(.*)', out)}
g1 = sp.Poly(sp.expand(G[1] / z**3), z)
assert g1.degree() == 4 and all(m[0] % 2 == 0 for m, _ in g1.terms())
c7, c5, c3 = [sp.factor(g1.coeff_monomial(z**p)) for p in (4, 2, 0)]
print('quadratic in U = z^2:  c7 U^2 + c5 U + c3 = 0 with')
print('   c7 =', c7); print('   c5 =', c5); print('   c3 =', c3)
p2 = sp.Poly(G[2], r02); R02 = sp.factor(-p2.coeff_monomial(1) / p2.LC())
p5 = sp.Poly(G[5], s02); S02 = sp.factor(-p5.coeff_monomial(1) / p5.LC())
R01 = sp.factor(sp.solve(G[4], r01)[0]); S01 = sp.factor(sp.solve(G[8], s01)[0])
print('   r02 =', R02); print('   s02 = ', S02); print('   r01 =', R01); print('   s01 =', S01)
json.dump({'c7': str(c7), 'c5': str(c5), 'c3': str(c3), 'r02': str(R02), 's02': str(S02), 'r01': str(R01), 's01': str(S01), 'r11': '-b', 'r12': '-3*z', 'r03': '3*z**2'}, open('s_closed_sol2.json', 'w'), indent=1)
print('   relabelling z -> -z exchanges the two points: s02(z) == r02(-z):', sp.cancel(sp.together(S02 - R02.subs(z, -z))) == 0, ';  s01(z) == r01(-z):', sp.cancel(sp.together(S01.subs(s02, S02) - R01.subs(r02, R02).subs(z, -z))) == 0)
# residual check: all six equations of s_sym_sol2.json vanish modulo the quadratic
d = json.load(open('s_sym_sol2.json'))
names = {n: sp.Symbol(n) for n in ('r11', 'r01', 'r02', 's01', 's02', 'z', 'A', 'B', 'b')}
rest = [sp.sympify(c, locals=names) for c in d['rest']]
subs_all = {names['r11']: -b, r01: R01, r02: R02, s01: S01, s02: S02}
quad = c7 * z**4 + c5 * z**2 + c3
res_ok = True
for c in rest:
    num = sp.numer(sp.together(c.subs(subs_all)))
    rem = sp.rem(sp.Poly(sp.expand(num), z), sp.Poly(quad, z))
    res_ok &= rem.is_zero
print('all six equations vanish modulo the quadratic (symbolic in A, B, b):', res_ok)
ok = res_ok
# (i) against the numeric batch
zz = sp.Symbol('zz'); nchk = bad = 0
for tt, tstr in (('2', '2'), ('9_4', '9/4')):
    t = sp.Rational(tstr); k = 2 * (t - 1) / (3 - t); bv = k / (k + 1)
    for line in open('u_points_t%s.txt' % tt):
        _, _, a0s, a1s = line.split()
        s_ = open('u_one_sym_t%s_%s_%s.log' % (tt, a0s.replace('/', '-'), a1s.replace('/', '-'))).read()
        if int(re.search(r'vdim:\n(-?\d+)', s_).group(1)) != 4:
            continue
        num = sp.Poly(sp.sympify(re.search(r'Jz\[1\]=(.*)', s_).group(1).replace('^', '**'), locals={'z': zz}), zz)
        sub = {A: sp.Rational(a0s)**2, B: sp.Rational(a1s)**2, b: bv}
        q = sp.Poly(quad.subs(sub).subs(z, zz), zz)
        nchk += 1; bad += sp.expand(q.as_expr() / q.LC() - num.as_expr() / num.LC()) != 0
print('(i) quadratic vs numeric z-eliminants: %d points, %d mismatches' % (nchk, bad)); ok &= bad == 0 and nchk == 70
# (ii) charges from the closed form
P, Q = sp.symbols('P Q')
A2s, A3s, A4s, B0s, B1s, B2s, C0s, C1s, l0s, Q30s, Q31s, Q40s, Q41s = sp.symbols('A2 A3 A4 B0 B1 B2 C0 C1 l0 Q30 Q31 Q40 Q41')
locg = {x.name: x for x in (A2s, A3s, A4s, B0s, B1s, B2s, C0s, C1s, l0s, Q30s, Q31s, Q40s, Q41s)}
def qs(r1, r2, r3, o1, o2, o3, zz_):
    dd = dict(q1=r1, q0=o1)
    for m in (1, 2, 3, 4):
        dd['q%d1' % m] = r1 * zz_**m + m * r2 * zz_**(m - 1)
        dd['q%d0' % m] = o1 * zz_**m + m * o2 * zz_**(m - 1) + sp.Rational(m * (m - 1), 2) * o3 * zz_**(m - 2)
    return dd
for tt, tstr in (('2', '2'), ('9_4', '9/4')):
    t = sp.Rational(tstr); k = 2 * (t - 1) / (3 - t); bv = k / (k + 1); n = 2 * k + 3; c = n / bv; rho2 = 2 / (t * t - 1)
    gc = json.load(open('g_engine_t%s.json' % tt))
    cB = {i: sp.sympify(gc['cB%d' % i], locals=locg) for i in (2, 4, 6)}
    vac = {A2s: 0, A3s: -c**2 * bv * (Q**2 + rho2), A4s: 0, B1s: 0, B0s: 0, C1s: 0, C0s: 0, B2s: 0, l0s: c * sp.sqrt(1 - bv) * P, Q30s: 0, Q31s: 0, Q40s: 0, Q41s: 0}
    lead = {i: sp.Poly(sp.expand(cB[i].subs(vac)), P, Q).coeff_monomial(P**i) for i in (2, 4, 6)}
    sb = {b: bv}
    qa = qs(-bv, -3 * z, 3 * z**2, R01.subs(sb), R02.subs(sb), 3 * z**2, z)
    qb = qs(-bv, 3 * z, 3 * z**2, S01.subs(sb), S02.subs(sb), 3 * z**2, -z)
    val = {A2s: 0, A3s: c**2 * (-B + qa['q1'] + qb['q1']), A4s: -c**3 * (qa['q0'] + qb['q0']), B2s: 0, l0s: c * sp.sqrt(A)}
    for m, (s1, s0) in zip((1, 2, 3, 4), ((B1s, B0s), (C1s, C0s), (Q31s, Q30s), (Q41s, Q40s))):
        val[s1] = c**(2 + m) * (qa['q%d1' % m] + qb['q%d1' % m]); val[s0] = -c**(3 + m) * (qa['q%d0' % m] + qb['q%d0' % m])
    quadv = sp.Poly(quad.subs(sb), z)
    ui = json.load(open('u_interp_t%s.json' % tt))
    dd1 = json.load(open('d1_sol2_t%s.json' % tt))['blocks']['I3']
    Mx = sp.Matrix(2, 2, [sp.sympify(dd1['matrix'][i][j], locals={'P': P, 'Q': Q}) for i in range(2) for j in range(2)])
    Mn = Mx / sp.Poly(sp.sympify(dd1['vacuum'], locals={'P': P, 'Q': Q}), P, Q).coeff_monomial(P**4)
    for i, name in ((4, 'I3'), (6, 'I5')):
        expr = sp.together(cB[i].subs(val) / lead[i])
        num, den = sp.fraction(expr)
        # v2 method: resultant in z of the quartic with e*den - num gives, over the four roots +-z1, +-z2, the polynomial
        # prod (e - I(z_root)) = (e^2 - trace e + det)^2 up to a constant (I is the same at z and -z: relabelling r <-> s).
        ev = sp.Symbol('ev')
        res = sp.Poly(sp.resultant(quadv.as_expr(), sp.expand(ev * den - num), z), ev)
        res = sp.Poly(sp.expand(res.as_expr() / res.LC()), ev)
        tr_i = sp.sympify(ui[name + '_trace'], locals={'P': P, 'Q': Q}).subs({P**2: A / (1 - bv), Q**2: B / bv - rho2})
        det_i = sp.sympify(ui[name + '_det'], locals={'P': P, 'Q': Q}).subs({P**2: A / (1 - bv), Q**2: B / bv - rho2})
        target = sp.Poly(sp.expand((ev**2 - tr_i * ev + det_i)**2), ev)
        same = sp.expand(res.as_expr() - target.as_expr()) == 0
        extra = ''
        if name == 'I3':
            trd = sp.expand(Mn.trace().subs({P**2: A / (1 - bv), Q**2: B / bv - rho2})); detd = sp.expand(Mn.det().subs({P**2: A / (1 - bv), Q**2: B / bv - rho2}))
            same_d = sp.expand(res.as_expr() - (ev**2 - trd * ev + detd)**2) == 0
            extra = ';  == (data characteristic polynomial)^2: %s' % same_d
            ok &= same_d
        print('(ii) t = %s %s: resultant over the closed-form solutions == (e^2 - trace e + det)^2 of the interpolated block: %s (degree %d)%s' % (tstr, name, same, res.degree(), extra), flush=True)
        ok &= same
# (iii) exceptional points at t = 2
bv = sp.Rational(2, 3)
for a0s, a1s in (('5/8', '1/6'), ('5/6', '8/9'), ('1/2', '1/3')):
    sub = {A: sp.Rational(a0s)**2, B: sp.Rational(a1s)**2, b: bv}
    vals = {nm: ex.subs(sub) for nm, ex in (('c7', c7), ('c5', c5), ('c3', c3))}
    disc = sp.factor((c5**2 - 4 * c7 * c3).subs(sub))
    Dden = sp.denom(sp.together(R02)).subs(sub)
    print('(iii) t = 2 (a0, a1) = (%s, %s): c7 = %s, c5 = %s, c3 = %s; discriminant = %s; denominator of r02 = %s' % (a0s, a1s, vals['c7'], vals['c5'], vals['c3'], disc, sp.factor(Dden)))
json.dump({'c7': str(c7), 'c5': str(c5), 'c3': str(c3), 'r02': str(R02), 's02': str(S02), 'r01': str(R01), 's01': str(S01), 'r11': '-b', 'r12': '-3*z', 'r03': '3*z**2'}, open('s_closed_sol2.json', 'w'), indent=1)
print('s_closed_sol2', 'PASS' if ok else 'FAIL'); sys.exit(0 if ok else 2)
