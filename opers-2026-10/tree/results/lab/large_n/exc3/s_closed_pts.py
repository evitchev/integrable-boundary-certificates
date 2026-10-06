"""EXC3 (post-hoc): charge-level check of the closed form at rational momentum points (exact): for each point, the resultant in z of the
quartic c7 z^4 + c5 z^2 + c3 with e*den(z) - num(z), where I(z) = num/den is the engine charge on the closed-form residues, must equal
(e^2 - trace e + det)^2 with trace, det the data (I_3) / registered-and-confirmed (I_5) block invariants.  Usage: s_closed_pts.py <t> <npts>"""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_EXC3.md', 'rb').read()).hexdigest() == open('SEAL_EXC3.sha256').read().split()[0]
tstr = sys.argv[1]; npts = int(sys.argv[2]); tt = tstr.replace('/', '_')
t = sp.Rational(tstr); k = 2 * (t - 1) / (3 - t); bv = k / (k + 1); n = 2 * k + 3; c = n / bv; rho2 = 2 / (t * t - 1)
A, B, b, z, ev = sp.symbols('A B b z ev'); r02, s02 = sp.symbols('r02 s02')
cf = json.load(open('s_closed_sol2.json')); loc = {'A': A, 'B': B, 'b': b, 'z': z, 'r02': r02, 's02': s02}
R02 = sp.sympify(cf['r02'], locals=loc); S02 = sp.sympify(cf['s02'], locals=loc)
R01 = sp.sympify(cf['r01'], locals=loc).subs(r02, R02); S01 = sp.sympify(cf['s01'], locals=loc).subs(s02, S02)
quad = sp.sympify(cf['c7'], locals=loc) * z**4 + sp.sympify(cf['c5'], locals=loc) * z**2 + sp.sympify(cf['c3'], locals=loc)
gc = json.load(open('g_engine_t%s.json' % tt))
syms = sp.symbols('A2 A3 A4 B0 B1 B2 C0 C1 l0 Q30 Q31 Q40 Q41'); locg = {x.name: x for x in syms}
A2s, A3s, A4s, B0s, B1s, B2s, C0s, C1s, l0s, Q30s, Q31s, Q40s, Q41s = syms
cB = {i: sp.sympify(gc['cB%d' % i], locals=locg) for i in (4, 6)}
P, Q = sp.symbols('P Q')
vac = {A2s: 0, A3s: -c**2 * bv * (Q**2 + rho2), A4s: 0, B1s: 0, B0s: 0, C1s: 0, C0s: 0, B2s: 0, l0s: c * sp.sqrt(1 - bv) * P, Q30s: 0, Q31s: 0, Q40s: 0, Q41s: 0}
lead = {i: sp.Poly(sp.expand(cB[i].subs(vac)), P, Q).coeff_monomial(P**i) for i in (4, 6)}
ui = json.load(open('u_interp_t%s.json' % tt))
dd = json.load(open('d1_sol2_t%s_with5.json' % tt))['blocks']
def qs(r1, r2, o1, o2, o3, zz):
    d = dict(q1=r1, q0=o1)
    for m in (1, 2, 3, 4):
        d['q%d1' % m] = r1 * zz**m + m * r2 * zz**(m - 1); d['q%d0' % m] = o1 * zz**m + m * o2 * zz**(m - 1) + sp.Rational(m * (m - 1), 2) * o3 * zz**(m - 2)
    return d
pts = [l.split()[2:4] for l in open('u_points_t%s.txt' % tt)][:npts]
ok = True; t0 = time.time()
for a0s, a1s in pts:
    a0 = sp.Rational(a0s); a1 = sp.Rational(a1s); sub = {A: a0**2, B: a1**2, b: bv}
    qv = sp.Poly(quad.subs(sub), z)
    if qv.degree() < 4 or qv.coeff_monomial(1) == 0:
        print('(%s, %s): exceptional, skipped' % (a0s, a1s)); continue
    qa = qs(-bv, -3 * z, R01.subs(sub), R02.subs(sub), 3 * z**2, z); qb = qs(-bv, 3 * z, S01.subs(sub), S02.subs(sub), 3 * z**2, -z)
    val = {A2s: 0, A3s: c**2 * (-a1**2 + qa['q1'] + qb['q1']), A4s: -c**3 * (qa['q0'] + qb['q0']), B2s: 0, l0s: c * a0}
    for m, (s1, s0) in zip((1, 2, 3, 4), ((B1s, B0s), (C1s, C0s), (Q31s, Q30s), (Q41s, Q40s))):
        val[s1] = c**(2 + m) * (qa['q%d1' % m] + qb['q%d1' % m]); val[s0] = -c**(3 + m) * (qa['q%d0' % m] + qb['q%d0' % m])
    X = a0**2 / (1 - bv); Y = a1**2 / bv - rho2
    line = '(%s, %s):' % (a0s, a1s)
    for i, name in ((4, 'I3'), (6, 'I5')):
        num, den = sp.fraction(sp.together(cB[i].subs(val) / lead[i]))
        res = sp.Poly(sp.resultant(qv.as_expr(), sp.expand(ev * den - num), z), ev); res = sp.Poly(sp.expand(res.as_expr() / res.LC()), ev)
        blk = dd[name]
        Mx = sp.Matrix(2, 2, [sp.sympify(blk['matrix'][r][cc], locals={'P': P, 'Q': Q}) for r in range(2) for cc in range(2)])
        Mn = Mx / sp.Poly(sp.sympify(blk['vacuum'], locals={'P': P, 'Q': Q}), P, Q).coeff_monomial(P**i)
        evl = lambda ex: sum(cq * X**(mm[0] // 2) * Y**(mm[1] // 2) for mm, cq in sp.Poly(sp.expand(ex), P, Q).terms())
        target = sp.Poly(sp.expand((ev**2 - evl(Mn.trace()) * ev + evl(Mn.det()))**2), ev)
        same = sp.expand(res.as_expr() - target.as_expr()) == 0
        line += '  %s resultant == (data charpoly)^2: %s' % (name, same); ok &= same
    print(line + '  (%.0fs)' % (time.time() - t0), flush=True)
print('s_closed_pts t = %s: %s' % (tstr, 'PASS' if ok else 'FAIL')); sys.exit(0 if ok else 2)
