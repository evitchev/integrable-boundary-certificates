"""(v2: q = 1/(z1 - z2) kept as a variable with the relation q (z1 - z2) = 1; v1 cleared denominators and reached degree 35.)
EXC2 addendum: TWO apparent singularities of type {-1,1,2,4} for the Sol-3 three-term operator.
Generates the 14 no-logarithm equations (7 at each point; the other point enters through its regular expansion) and
writes a Singular script for the elimination.   Usage: h_two.py <t> <a0> <a1> [onepoint]
'onepoint' = control: second point's residues set to zero (must reproduce the four one-point solutions)."""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import frob as Fb
assert hashlib.sha256(open('SEAL_EXC2_addendum.md', 'rb').read()).hexdigest() == open('SEAL_EXC2_addendum.sha256').read().split()[0]
t = sp.Rational(sys.argv[1]); a0 = sp.Rational(sys.argv[2]); a1 = sp.Rational(sys.argv[3]); onepoint = len(sys.argv) > 4
k = 2 * (t - 1) / (3 - t); b = k / (k + 1); n = k + 4
z, E, e = Fb.z, Fb.E, Fb.e
z2 = sp.Symbol('z2')
QQ = sp.Symbol('qq')
th = sp.Symbol('th')
names = ['21', '22', '11', '12', '13', '01', '02', '03', '04']
r = {nm: sp.Symbol('r' + nm) for nm in names}
s = {nm: sp.Symbol('s' + nm) for nm in names}
t0 = time.time()
P = sp.Poly(sp.expand((th**2 - a0**2) * (th**2 - a1**2)), th)
Pcoef = [P.coeff_monomial(th**j) for j in range(5)]
Lcoef = [sp.Rational(1, 2), 1]
K = 2          # span 5: need l_k up to k = 5 - 4 = 1; take 2


def local_two(rr, ss, zz, zo):
    """operator around xi = zz: point-1 residues rr (singular), point-2 residues ss at zo (regular expansion)."""
    order = 4
    powers = [{0: sp.Integer(1)}]
    def theta_on(series):
        out = {}
        for kk, g in series.items():
            f = sp.expand(g * (e + kk))
            out[kk - 1] = sp.expand(out.get(kk - 1, 0) + zz * f)
            out[kk] = sp.expand(out.get(kk, 0) + f)
        return {kk: v for kk, v in out.items() if v != 0}
    for j in range(order):
        powers.append(Fb.trunc(theta_on(powers[-1]), K))
    xi = {0: zz, 1: sp.Integer(1)}
    out = {}
    for j, c in enumerate(Pcoef):
        out = Fb.add_series(out, {kk: c * v for kk, v in powers[j].items()})
    for j, c in enumerate(Lcoef):
        out = Fb.add_series(out, Fb.mul_series({kk: c * v for kk, v in powers[j].items()}, xi))
    xib = {i: sp.binomial(b, i) / zz**i for i in range(0, K + order + 1)}
    out = Fb.add_series(out, {kk: -E * v for kk, v in xib.items()})
    for nm, c in rr.items():
        j, m = int(nm[0]), int(nm[1])
        out = Fb.add_series(out, Fb.mul_series(Fb.mul_series(powers[j], xi), {-m: c}))
    qq = QQ if zz == z else -QQ          # 1/(zz - zo)
    for nm, c in ss.items():
        j, m = int(nm[0]), int(nm[1])
        ser = {i: sp.binomial(-m, i) * qq ** (m + i) for i in range(0, K + order + 1)}
        out = Fb.add_series(out, Fb.mul_series(Fb.mul_series(powers[j], xi), {i: c * v for i, v in ser.items()}))
    return Fb.trunc(out, K)


exps = [-1, 1, 2, 4]
eq_all = []
for (rr, ss, zz, zo, tag) in ((r, s, z, z2, 'A'), (s, r, z2, z, 'B')):
    ssub = {nm: (sp.Integer(0) if onepoint and tag == 'A' else v) for nm, v in ss.items()}
    if onepoint and tag == 'B':
        continue
    lk = local_two(rr, ssub, zz, zo)
    lk = {kk: v.subs(Fb.z, zz) if False else v for kk, v in lk.items()}
    lead, conds, kaps = Fb.nolog_conditions(lk, exps, 4)
    ind = {rr['22']: -4 * zz, rr['13']: 8 * zz**2, rr['04']: -8 * zz**3}
    assert sp.expand((lead - zz**4 * sp.prod([(e - x) for x in exps])).subs(ind)) == 0
    eq_all.append((tag, conds, kaps, ind))
    print('point %s: conditions generated (%.0fs)' % (tag, time.time() - t0), flush=True)
ind_all = {}
for _, _, _, ind in eq_all:
    ind_all.update(ind)
if onepoint:
    ind_all.update({v: 0 for v in s.values()})
polys = []
for tag, conds, kaps, ind in eq_all:
    for pos, expr in conds:
        num = sp.numer(sp.together(sp.expand(expr.subs(ind_all))))
        Pp = sp.Poly(sp.expand(num), E, *kaps)
        for mon, cf in Pp.terms():
            if cf != 0:
                polys.append((tag, pos, mon, sp.expand(cf)))
unk1 = [r[nm] for nm in ('21', '11', '12', '01', '02', '03')] + [z]
unk2 = [] if onepoint else [s[nm] for nm in ('21', '11', '12', '01', '02', '03')] + [z2]
unknowns = unk1 + unk2
print('%d equations in %d unknowns (+ q); total degrees %s  (%.0fs)' % (len(polys), len(unknowns), [sp.Poly(p_, *(unknowns + [QQ])).total_degree() for _, _, _, p_ in polys], time.time() - t0))
# charges: effective quartic
gc = json.load(open('g_charges_t%s.json' % str(t).replace('/', '_')))
A2, A3, A4 = sp.symbols('A2 A3 A4')
cB2 = sp.sympify(gc['cB2'], locals={'A2': A2, 'A3': A3, 'A4': A4}); cB4 = sp.sympify(gc['cB4'], locals={'A2': A2, 'A3': A3, 'A4': A4})
sc = n / b
x0s, x1s = sp.symbols('x0s x1s')
def at(q2, q3, q4, expr):
    return sp.expand(expr.subs({A2: sc**2 * q2, A3: -sc**3 * q3, A4: sc**4 * q4}))
lead4 = sp.Poly(at(-(x0s + x1s), 0, x0s * x1s, cB4), x0s, x1s).coeff_monomial(x0s**2) * b**2
lead2 = sp.Poly(at(-(x0s + x1s), 0, x0s * x1s, cB2), x0s, x1s).coeff_monomial(x0s) * b
tot = lambda nm: r[nm] + (0 if onepoint else s[nm])
I1 = at(-(a0**2 + a1**2) + tot('21'), tot('11'), a0**2 * a1**2 + tot('01'), cB2) / lead2
L3 = at(-(a0**2 + a1**2) + tot('21'), tot('11'), a0**2 * a1**2 + tot('01'), cB4) / lead4
vacI1 = at(-(a0**2 + a1**2), 0, a0**2 * a1**2, cB2) / lead2
# Singular script
ev, iv, wv = sp.symbols('ee ii ww')
gens = [wv] + unknowns + ([] if onepoint else [QQ]) + [iv, ev]
sat = wv * z * (1 if onepoint else z2) - 1
def sing(p_):
    return str(sp.expand(sp.numer(sp.together(p_)))).replace('**', '^')
lines = ['ring R = 0,(%s),dp;' % ','.join(str(g) for g in gens), 'option(redSB);', 'ideal I =']
allp = [p_ for _, _, _, p_ in polys] + [sat, sp.expand(iv - I1), sp.expand(ev - L3)] + ([] if onepoint else [sp.expand(QQ * (z - z2) - 1)])
lines.append(',\n'.join(sing(p_) for p_ in allp) + ';')
elimvars = '*'.join(str(g) for g in [wv] + unknowns + ([] if onepoint else [QQ]))
lines += ['ideal G = std(I);', 'print("dim:"); dim(G);', 'print("vdim:"); vdim(G);',
          'ideal Je = eliminate(G, %s*ii);' % elimvars, 'print("eliminant in ee:"); Je;',
          'ideal Ji = eliminate(G, %s*ee);' % elimvars, 'print("eliminant in ii:"); Ji;',
          'ideal Jz = eliminate(G, %s*ii*ee);' % '*'.join(str(g) for g in [wv] + [u for u in unknowns if u != z] + ([] if onepoint else [QQ])), 'print("eliminant in z:"); Jz;', 'quit;']
fn = 'h_two_v2_t%s_%s_%s%s.sing' % (str(t).replace('/', '_'), str(a0).replace('/', '-'), str(a1).replace('/', '-'), '_onepoint' if onepoint else '')
open(fn, 'w').write('\n'.join(lines) + '\n')
json.dump({'vacI1': str(vacI1), 'onepoint': onepoint}, open(fn.replace('.sing', '.json'), 'w'))
print('written %s; vacuum monic I_1 = %s  (%.0fs)' % (fn, vacI1, time.time() - t0))
