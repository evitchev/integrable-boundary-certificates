"""JOINT1: state-by-state (joint) check of the excited-state operators -- sealed method (SEAL_JOINT1.md).
Cells: C1/C2 EXC1 Sol 1 (t = 2, 9/4); C3/C4 EXC2 Sol 3; C5/C6 EXC3 Sol 2; all level one, 2 x 2 singlet blocks of I_3 and I_5.
Data: archived blocks (inputs/d1_sol<s>_t<t>_with5.json), each divided by the P^(2K) coefficient of its vacuum eigenvalue.
Operator: archived closed forms (inputs/).  Joint invariant: tr(M_3 M_5) = a1 b1 + a2 b2 (and tr(M_3^2 M_5), tr(M_3 M_5^2)) vs the operator's
sum over its two selected solutions; explicit pairs over the quadratic splitting field at rational points.
Usage: joint1.py <cell C1..C6> [crossed | shift]   (tampers: crossed pairing on the operator side; P -> P + 1/10 on the operator side).
Exit 0 iff commute, non-degenerate at the points, mixed traces and pairs match (real); tampers: 0 iff they fire."""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_JOINT1.md', 'rb').read()).hexdigest() == open('SEAL_JOINT1.sha256').read().split()[0]
cell = sys.argv[1]; mode = sys.argv[2] if len(sys.argv) > 2 else 'real'
CELLS = {'C1': (1, '2'), 'C2': (1, '9/4'), 'C3': (3, '2'), 'C4': (3, '9/4'), 'C5': (2, '2'), 'C6': (2, '9/4')}
sol, tstr = CELLS[cell]
t = sp.Rational(tstr); k = 2 * (t - 1) / (3 - t); b = k / (k + 1); rho2 = 2 / (t * t - 1); tt = tstr.replace('/', '_')
P, Q = sp.symbols('P Q'); U = sp.Symbol('U'); z = sp.Symbol('z')
POINTS = [(sp.Rational(1, 2), sp.Rational(1, 3)), (sp.Rational(2, 5), sp.Rational(1, 4))]
SHIFT = sp.Rational(1, 10) if mode == 'shift' else 0
t0 = time.time()
# ---------------------------------------------------------------- data side
dd = json.load(open('inputs/d1_sol%d_t%s_with5.json' % (sol, tt)))
M = {}
for name, K in (('I3', 2), ('I5', 3)):
    blk = dd['blocks'][name]
    Mx = sp.Matrix(2, 2, [sp.sympify(blk['matrix'][i][j], locals={'P': P, 'Q': Q}) for i in range(2) for j in range(2)])
    lead = sp.Poly(sp.sympify(blk['vacuum'], locals={'P': P, 'Q': Q}), P, Q).coeff_monomial(P**(2 * K))
    M[name] = (Mx / lead).applyfunc(sp.expand)
M3, M5 = M['I3'], M['I5']
comm = (M3 * M5 - M5 * M3).applyfunc(sp.expand)
commute = comm == sp.zeros(2, 2)
e = sp.Symbol('e')
cp3 = sp.Poly(sp.expand((M3 - e * sp.eye(2)).det()), e)
disc3 = sp.factor(sp.discriminant(cp3.as_expr(), e))
data = {'T11': sp.expand((M3 * M5).trace()), 'T21': sp.expand((M3 * M3 * M5).trace()), 'T12': sp.expand((M3 * M5 * M5).trace()),
        's3': sp.expand(M3.trace()), 'p3': sp.expand(M3.det()), 's5': sp.expand(M5.trace()), 'p5': sp.expand(M5.det())}
data['crossed_T11'] = sp.expand(data['s3'] * data['s5'] - data['T11'])      # a1 b2 + a2 b1
print('%s: Sol %d, t = %s, block 2 x 2.  [M_3, M_5] = 0 symbolically: %s.  disc(I_3 charpoly) = %s' % (cell, sol, tstr, commute, disc3), flush=True)


# ---------------------------------------------------------------- operator side: symmetric functions over the two solutions
def sym_traces(f3, f5, modpoly, var):
    """f3, f5: expressions (polynomial or rational) in var; modpoly: monic polynomial in var (degree d) whose roots are the solutions
    (each unordered solution counted mult times).  Returns traces of f3 f5, f3^2 f5, f3 f5^2 and of f3, f5 (sum over the roots) via the
    multiplication matrices in K[var]/(modpoly); exact."""
    q = sp.Poly(modpoly, var); d = q.degree()
    def rep(f):
        num, den = sp.fraction(sp.together(f))
        num = sp.Poly(num, var); den = sp.Poly(den, var)
        if den.degree() == 0:
            r = sp.rem(num, q); return sp.Poly(sp.expand(r.as_expr() / den.as_expr()), var)
        inv = sp.invert(den.as_expr(), q.as_expr(), var)
        return sp.rem(sp.Poly(sp.expand(num.as_expr() * inv), var), q)
    def mult(fpoly):
        cols = []
        for j in range(d):
            colp = sp.rem(sp.Poly(sp.expand(fpoly.as_expr() * var**j), var), q)
            cols.append([colp.coeff_monomial(var**i) for i in range(d)])
        return sp.Matrix(d, d, lambda i, j: cols[j][i])
    A3, A5 = mult(rep(f3)), mult(rep(f5))
    return {'T11': sp.expand((A3 * A5).trace()), 'T21': sp.expand((A3 * A3 * A5).trace()), 'T12': sp.expand((A3 * A5 * A5).trace()),
            's3': sp.expand(A3.trace()), 's5': sp.expand(A5.trace()), 'T20': sp.expand((A3 * A3).trace())}


opsym = None        # symbolic (P, Q) operator traces, if available
oppts = {}          # per point: operator traces
if sol == 1:
    Mv = (t + 3) / (t - 1)
    od = json.load(open('inputs/o1_t%s_n1.json' % tt))
    l0, l1, z1 = sp.symbols('l0 l1 z1')
    R = {i: sp.sympify(od['R'][str(i)], locals={'Symbol': sp.Symbol}) for i in (4, 6)}
    L = {}
    for i, K in ((4, 2), (6, 3)):
        c = sp.Poly(R[i], l0, l1, z1).coeff_monomial(l1**i)
        vac_lead = sp.Poly(sp.expand((R[i] / c).subs(z1, 0)), l0, l1).coeff_monomial(l0**i) * (4 * Mv * (Mv + 1))**K
        L[i] = sp.expand(R[i] / c / vac_lead)
    cond = sp.expand(2 * Mv * z1**2 + l0 * (Mv - 1) * z1 - 2 * l1**2 + (Mv - 1)**2 / 2)
    cond = sp.expand(cond / sp.Poly(cond, z1).LC())
    tr = sym_traces(L[4], L[6], cond, z1)
    PX, PI = sp.symbols('PX PI')
    def to_bare(expr):
        p_ = sp.Poly(sp.expand(expr), l0)
        assert all(m[0] % 2 == 0 for m, _ in p_.terms()), 'odd power of alpha in a symmetric function'
        ex = sum(cf * (4 * Mv * (Mv + 1) * (P + SHIFT)**2)**(m[0] // 2) for m, cf in p_.terms())
        return sp.expand(ex.subs(l1**2, (Mv + 1) * (Q**2 + rho2)))
    opsym = {k_: to_bare(v) for k_, v in tr.items()}
elif sol == 3:
    A, B, bb, rho, zz, r11, u = sp.symbols('A B b rho z r11 u')
    d = json.load(open('inputs/s_sym.json')); names = ['r01', 'r02', 's01', 's02']
    loc = {'A': A, 'B': B, 'b': bb, 'rho': rho, 'z': zz, 'r11': r11}; loc.update({n: sp.Symbol(n) for n in names})
    Zu = u * (u**2 - 8 * (A + B) + 4 * (bb - 1)**2 + 4) / 16
    c1 = -8 * (bb - 2) * (2 * A + 2 * B + bb**2 - 2) / (3 * bb - 4); c0 = -16 * bb * ((bb - 2)**2 - 4 * A) * ((bb - 2)**2 - 4 * B) / (3 * bb - 4)
    full = {r11: bb * u / 2, rho: u + bb + 4, zz: Zu}
    val = {n: sp.sympify(d['sol'][n], locals=loc).subs(full) for n in names}
    r01s01 = sp.expand(val['r01'] + val['s01'])
    q11 = sp.expand((bb * u / 2) * Zu + (u + bb + 4) * Zu + (-bb * u / 2) * (-Zu) + (u - bb - 4) * Zu)
    def inU(ex):
        p_ = sp.Poly(sp.expand(ex), u); assert all(m[0] % 2 == 0 for m, _ in p_.terms())
        return sum(cf * U**(m[0] // 2) for m, cf in p_.terms())
    q4U = inU(A * B + r01s01); q11U = inU(q11)
    n = k + 4; nb = n / b
    gc = json.load(open('inputs/u3_engine_t%s.json' % tt))
    A2, A3, A4, B0, B1, B2 = sp.symbols('A2 A3 A4 B0 B1 B2'); locg = {s.name: s for s in (A2, A3, A4, B0, B1, B2)}
    cB4 = sp.sympify(gc['cB4'], locals=locg); cB6 = sp.sympify(gc['cB6'], locals=locg)
    x0s, x1s = sp.symbols('x0s x1s')
    vacmap = {A2: nb**2 * (-(x0s + x1s)), A3: 0, A4: nb**4 * x0s * x1s, B1: 0, B2: 0}
    lead6 = sp.Poly(sp.expand(cB6.subs(vacmap)), x0s, x1s).coeff_monomial(x0s**3) * b**3
    lead4 = sp.Poly(sp.expand(cB4.subs(vacmap)), x0s, x1s).coeff_monomial(x0s**2) * b**2
    valmap = {A2: nb**2 * (-(A + B) - 2 * b), A3: 0, A4: nb**4 * q4U.subs(bb, b), B2: 0, B1: -nb**6 * q11U.subs(bb, b)}
    I3U = sp.expand(cB4.subs(valmap) / lead4); I5U = sp.expand(cB6.subs(valmap) / lead6)
    modU = sp.expand(U**2 + c1.subs(bb, b) * U + c0.subs(bb, b))
    tr = sym_traces(I3U, I5U, modU, U)
    mom = {A: b * (P + SHIFT)**2, B: b * (Q**2 + rho2)}
    opsym = {k_: sp.expand(sp.cancel(sp.together(v)).subs(mom)) for k_, v in tr.items()}
else:   # sol 2 (EXC3): rational points only
    A, B, bb = sp.symbols('A B b'); r02, s02 = sp.symbols('r02 s02')
    cf = json.load(open('inputs/s_closed_sol2.json')); loc = {'A': A, 'B': B, 'b': bb, 'z': z, 'r02': r02, 's02': s02}
    R02 = sp.sympify(cf['r02'], locals=loc); S02 = sp.sympify(cf['s02'], locals=loc)
    R01 = sp.sympify(cf['r01'], locals=loc).subs(r02, R02); S01 = sp.sympify(cf['s01'], locals=loc).subs(s02, S02)
    quad = sp.sympify(cf['c7'], locals=loc) * z**4 + sp.sympify(cf['c5'], locals=loc) * z**2 + sp.sympify(cf['c3'], locals=loc)
    n = 2 * k + 3; c = n / b
    gc = json.load(open('inputs/g_engine_t%s.json' % tt))
    syms = sp.symbols('A2 A3 A4 B0 B1 B2 C0 C1 l0 Q30 Q31 Q40 Q41'); locg = {x.name: x for x in syms}
    A2s, A3s, A4s, B0s, B1s, B2s, C0s, C1s, l0s, Q30s, Q31s, Q40s, Q41s = syms
    cB = {i: sp.sympify(gc['cB%d' % i], locals=locg) for i in (4, 6)}
    vac = {A2s: 0, A3s: -c**2 * b * (Q**2 + rho2), A4s: 0, B1s: 0, B0s: 0, C1s: 0, C0s: 0, B2s: 0, l0s: c * sp.sqrt(1 - b) * P, Q30s: 0, Q31s: 0, Q40s: 0, Q41s: 0}
    lead = {i: sp.Poly(sp.expand(cB[i].subs(vac)), P, Q).coeff_monomial(P**i) for i in (4, 6)}
    def qs(r1, r2, o1, o2, o3, zz_):
        dq = dict(q1=r1, q0=o1)
        for m in (1, 2, 3, 4):
            dq['q%d1' % m] = r1 * zz_**m + m * r2 * zz_**(m - 1); dq['q%d0' % m] = o1 * zz_**m + m * o2 * zz_**(m - 1) + sp.Rational(m * (m - 1), 2) * o3 * zz_**(m - 2)
        return dq
    def op_traces_at(Pv, Qv):
        Pv = Pv + SHIFT
        sub = {A: (1 - b) * Pv**2, B: b * (Qv**2 + rho2), bb: b}
        qv = sp.Poly(sp.expand(quad.subs(sub)), z); qv = sp.Poly(sp.expand(qv.as_expr() / qv.LC()), z)
        qa = qs(-b, -3 * z, R01.subs(sub), R02.subs(sub), 3 * z**2, z); qb = qs(-b, 3 * z, S01.subs(sub), S02.subs(sub), 3 * z**2, -z)
        val = {A2s: 0, A3s: c**2 * (-sub[B] + qa['q1'] + qb['q1']), A4s: -c**3 * (qa['q0'] + qb['q0']), B2s: 0, l0s: c * sp.sqrt(sub[A])}
        for m, (s1, s0) in zip((1, 2, 3, 4), ((B1s, B0s), (C1s, C0s), (Q31s, Q30s), (Q41s, Q40s))):
            val[s1] = c**(2 + m) * (qa['q%d1' % m] + qb['q%d1' % m]); val[s0] = -c**(3 + m) * (qa['q%d0' % m] + qb['q%d0' % m])
        I3z = sp.together(cB[4].subs(val) / lead[4]); I5z = sp.together(cB[6].subs(val) / lead[6])
        trz = sym_traces(I3z, I5z, qv.as_expr(), z)              # over the 4 roots +-z1, +-z2: each unordered solution twice
        return {k_: sp.Rational(v) / 2 for k_, v in trz.items()}


def report(label, dat, op):
    keys = ('s3', 's5', 'T11', 'T21', 'T12')
    res = {k_: sp.expand(dat[k_] - op[k_]) == 0 for k_ in keys}
    crossed_vs_op = sp.expand(dat['crossed_T11'] - op['T11']) == 0
    sep = sp.expand(dat['T11'] - dat['crossed_T11'])              # (a1 - a2)(b1 - b2)
    print('   %s: marginals tr/det-level (s3, s5) equal: %s %s; MIXED tr(M3 M5): %s; tr(M3^2 M5): %s; tr(M3 M5^2): %s; data crossed invariant equals operator T11: %s; (a1-a2)(b1-b2) = %s'
          % (label, res['s3'], res['s5'], res['T11'], res['T21'], res['T12'], crossed_vs_op, sep if len(str(sep)) < 80 else '(nonzero polynomial, %d terms)' % len(sp.Poly(sep, P, Q).terms())), flush=True)
    return res, crossed_vs_op, sep


ok = commute
if opsym is not None:
    res, crossed, sep = report('SYMBOLIC in (P, Q)', data, opsym)
    if mode == 'crossed':
        ok &= crossed and all(res[k_] for k_ in ('s3', 's5'))          # the crossed pairing would match only if T11 == crossed
    elif mode == 'shift':
        ok &= all(res.values())
    else:
        ok &= all(res.values()) and sep != 0
# explicit pairs at the points
for (Pv, Qv) in POINTS:
    sub = {P: Pv, Q: Qv}
    M3v = M3.subs(sub); M5v = M5.subs(sub)
    dat_pt = {k_: sp.Rational(v.subs(sub)) for k_, v in data.items()}
    if opsym is not None:
        op_pt = {k_: sp.Rational(v.subs(sub)) for k_, v in opsym.items()}
    else:
        op_pt = op_traces_at(Pv, Qv)
    disc_pt = sp.Rational(disc3.subs(sub))
    # data pairs: eigenvectors of M3 over Q(sqrt(disc))
    a1, a2 = sp.solve(cp3.as_expr().subs(sub), e)
    pairs_data = []
    for a in (a1, a2):
        v = (M3v - a * sp.eye(2)).nullspace()[0]
        w = M5v * v
        bval = sp.radsimp(sp.cancel(w[0] / v[0])) if v[0] != 0 else sp.radsimp(sp.cancel(w[1] / v[1]))
        pairs_data.append((sp.radsimp(a), sp.radsimp(bval)))
    p3_op = (op_pt['s3']**2 - op_pt['T20']) / 2        # a1 a2 from the operator's own traces
    print('   point (P,Q) = (%s, %s): disc(I_3) = %s%s; DATA pairs (I_3, I_5): %s' % (Pv, Qv, disc_pt, ' DEGENERATE' if disc_pt == 0 else '', pairs_data), flush=True)
    bsum = op_pt['s5']; T11 = op_pt['T11']
    pairs_op = None
    a_op = sorted([sp.radsimp(x) for x in sp.solve(e**2 - op_pt['s3'] * e + p3_op, e)], key=lambda x: sp.N(x, 30).as_real_imag())
    if a_op[0] != a_op[1]:
        b1 = sp.radsimp((T11 - a_op[1] * bsum) / (a_op[0] - a_op[1])); b2 = sp.radsimp(bsum - b1)
        pairs_op = [(a_op[0], b1), (a_op[1], b2)]
    same_pairs = pairs_op is not None and all(any(sp.simplify(pa[0] - po[0]) == 0 and sp.simplify(pa[1] - po[1]) == 0 for po in pairs_op) for pa in pairs_data)
    crossed_pairs = pairs_op is not None and all(any(sp.simplify(pa[0] - po[0]) == 0 and sp.simplify(pa[1] - (bsum - po[1])) == 0 for po in pairs_op) for pa in pairs_data)
    print('      OPER pairs %s: %s;  sets of pairs EQUAL: %s;  crossed pairing would match: %s' % ('(shifted momentum)' if mode == 'shift' else '', pairs_op, same_pairs, crossed_pairs), flush=True)
    if opsym is None:
        res, crossed, sep = report('at this point', dat_pt, op_pt)
        if mode == 'real':
            ok &= all(res.values()) and sep != 0
        elif mode == 'crossed':
            ok &= crossed
        else:
            ok &= all(res.values())
    ok &= (same_pairs if mode == 'real' else (crossed_pairs if mode == 'crossed' else same_pairs)) and disc_pt != 0
print('%s [%s]: %s  (%.0fs)' % (cell, mode, ('PASS' if ok else 'FAIL') if mode == 'real' else ('tamper DOES NOT FIRE' if ok else 'tamper fires'), time.time() - t0))
sys.exit((0 if ok else 2) if mode == 'real' else (0 if not ok else 1))
