"""EXC3 (seal sec. 3, two-point class): TWO apparent singularities of type S_a = {-1,1,3} for the Sol 2 normal form, at xi = z
(residues r) and xi = z2 (residues s); the other point enters each local expansion through its regular Taylor series.
Modes:  free2 -> z2 free, r11 = s11 = -b substituted (kk := z + z2, the sigma-pair test);
        sym  -> z2 = -z imposed (the sigma-invariant candidates; 7 unknowns r11 r01 r02 s11 s01 s02 z);
        free -> z2 free, with qq = 1/(z - z2) as a variable (8 unknowns + qq).
Indicial conditions fix r12 = -3 z, r03 = 3 z^2 (and the same at z2).  Writes h_two_<mode>_<TAG>.sing (Singular: dim, vdim,
eliminants in ee (I_3 monic), ii (I_1 monic), jj (J_3), kk (sigma-check s11 - r11), zz) and runs it.
Usage: h_two_sol2.py <mode> <t> <a0> <a1> [tamper]   (tamper: middle term xi((vartheta+1)^2 - a0^2))."""
import hashlib, json, re, subprocess, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import frob as Fb
assert hashlib.sha256(open('SEAL_EXC3.md', 'rb').read()).hexdigest() == open('SEAL_EXC3.sha256').read().split()[0]
mode, tstr = sys.argv[1], sys.argv[2]; a0 = sp.Rational(sys.argv[3]); a1 = sp.Rational(sys.argv[4]); tamper = len(sys.argv) > 5
t = sp.Rational(tstr); k = 2 * (t - 1) / (3 - t); b = k / (k + 1); n = 2 * k + 3; c = n / b; rho2 = 2 / (t * t - 1)
TAG = 't%s_%s_%s%s' % (tstr.replace('/', '_'), sys.argv[3].replace('/', '-'), sys.argv[4].replace('/', '-'), '_tamper' if tamper else '')
z, E, e = Fb.z, Fb.E, Fb.e
z2 = sp.Symbol('z2'); QQ = sp.Symbol('qq')
names = ['11', '12', '01', '02', '03']
r = {nm: sp.Symbol('r' + nm) for nm in names}; s = {nm: sp.Symbol('s' + nm) for nm in names}
t0 = time.time()
Pcoef = [0, -a1**2, 0, 1]
shift = 1 if tamper else sp.Rational(1, 2)
Lcoef = [shift**2 - a0**2, 2 * shift, 1]
exps = [-1, 1, 3]; K = exps[-1] - exps[0] - 3 + 1


def local_two(rr, ss, zz, zo):
    order = 3
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
    for j, cc in enumerate(Pcoef):
        out = Fb.add_series(out, {kk: cc * v for kk, v in powers[j].items()})
    for j, cc in enumerate(Lcoef):
        out = Fb.add_series(out, Fb.mul_series({kk: cc * v for kk, v in powers[j].items()}, xi))
    xib = {i: sp.binomial(b, i) / zz**i for i in range(0, K + order + 1)}
    out = Fb.add_series(out, {kk: -E * sp.expand_func(v) for kk, v in xib.items()})
    for nm, cc in rr.items():
        j, m = int(nm[0]), int(nm[1])
        out = Fb.add_series(out, Fb.mul_series(Fb.mul_series(powers[j], xi), {-m: cc}))
    qq = (1 / (zz - zo)) if mode == 'sym' else (QQ if zz == z else -QQ)   # free, free2: qq variable
    for nm, cc in ss.items():
        j, m = int(nm[0]), int(nm[1])
        ser = {i: sp.binomial(-m, i) * qq ** (m + i) for i in range(0, K + order + 1)}
        out = Fb.add_series(out, Fb.mul_series(Fb.mul_series(powers[j], xi), {i: cc * v for i, v in ser.items()}))
    return Fb.trunc(out, K)


zB = -z if mode == 'sym' else z2
ind = {r['12']: -3 * z, r['03']: 3 * z**2, s['12']: -3 * zB, s['03']: 3 * zB**2}
if mode == 'free2':      # free points with r11 = s11 = -b substituted (the E-linear conditions, local to each point)
    ind.update({r['11']: -b, s['11']: -b})
polys = []
for (rr, ss, zz, zo, tag) in ((r, s, z, zB, 'A'), (s, r, zB, z, 'B')):
    lk = local_two(rr, ss, zz, zo)
    lead, conds, kaps = Fb.nolog_conditions(lk, exps, 3)
    assert sp.expand((lead - zz**3 * sp.prod([(e - x) for x in exps])).subs(ind)) == 0, 'indicial'
    for pos, expr in conds:
        num = sp.numer(sp.together(sp.expand(expr.subs(ind))))
        Pp = sp.Poly(sp.expand(num), E, *kaps)
        for mon, cf in Pp.terms():
            if cf != 0:
                core = sp.Integer(1)
                for fac, mult in sp.factor_list(sp.expand(cf))[1]:
                    if fac not in (z, z2):
                        core *= fac**mult
                polys.append(sp.expand(core))
    print('point %s: conditions generated (%.0fs)' % (tag, time.time() - t0), flush=True)
unk = ([r['11'], s['11']] if mode != 'free2' else []) + [r['01'], r['02'], s['01'], s['02'], z] + ([] if mode == 'sym' else [z2, QQ])
print('%d equations in %d unknowns; total degrees %s' % (len(polys), len(unk), [sp.Poly(p_, *unk).total_degree() for p_ in polys]), flush=True)
# charges (engine json) and observables
gc = json.load(open('g_engine_t%s.json' % tstr.replace('/', '_')))
A2, A3, A4, B0, B1, B2, C0, C1, l0, Q30, Q31, Q40, Q41 = sp.symbols('A2 A3 A4 B0 B1 B2 C0 C1 l0 Q30 Q31 Q40 Q41')
loc = {x.name: x for x in (A2, A3, A4, B0, B1, B2, C0, C1, l0, Q30, Q31, Q40, Q41)}
cB = {i: sp.sympify(gc['cB%d' % i], locals=loc) for i in (2, 3, 4, 5, 6)}
P, Q = sp.symbols('P Q')
vac = {A2: 0, A3: -c**2 * b * (Q**2 + rho2), A4: 0, B1: 0, B0: 0, C1: 0, C0: 0, B2: 0, l0: c * sp.sqrt(1 - b) * P, Q30: 0, Q31: 0, Q40: 0, Q41: 0}
lead = {i: sp.Poly(sp.expand(cB[i].subs(vac)), P, Q).coeff_monomial(P**i) for i in (2, 4, 6)}
def qs(rr, zz):
    # coefficient of xi^(-m) in xi [r_j1/(xi-z) + r_j2/(xi-z)^2 + r_j3/(xi-z)^3]:  r_j1 z^m + m r_j2 z^(m-1) + m(m-1)/2 r_j3 z^(m-2)
    d = dict(q1=rr['11'], q0=rr['01'])
    for m in (1, 2, 3, 4):
        d['q%d1' % m] = rr['11'] * zz**m + m * rr['12'] * zz**(m - 1)
        d['q%d0' % m] = rr['01'] * zz**m + m * rr['02'] * zz**(m - 1) + sp.Rational(m * (m - 1), 2) * rr['03'] * zz**(m - 2)
    return d
qa, qb = qs(r, z), qs(s, zB)
val = {A2: 0, A3: c**2 * (-a1**2 + qa['q1'] + qb['q1']), A4: -c**3 * (qa['q0'] + qb['q0']), B2: 0, l0: c * a0}
# (1/xi^m) (q_m1 vartheta + q_m0)  ->  x^(-mc) (c^(2+m) q_m1 T - c^(3+m) q_m0)
for m, (s1, s0) in zip((1, 2, 3, 4), ((B1, B0), (C1, C0), (Q31, Q30), (Q41, Q40))):
    val[s1] = c**(2 + m) * (qa['q%d1' % m] + qb['q%d1' % m]); val[s0] = -c**(3 + m) * (qa['q%d0' % m] + qb['q%d0' % m])
val = {k_: sp.expand(sp.sympify(v).subs(ind)) for k_, v in val.items()}
obs = {'ii': sp.expand(cB[2].subs(val) / lead[2]), 'ee': sp.expand(cB[4].subs(val) / lead[4]), 'ff': sp.expand(cB[6].subs(val) / lead[6]),
       'jj': sp.expand(cB[3].subs(val)), 'kk': (s['11'] - r['11']) if mode != 'free2' else (z + z2)}
syms = {nm: sp.Symbol(nm) for nm in obs}
def rel(nm):
    num, den = sp.fraction(sp.together(syms[nm] - obs[nm]))
    core = sp.Integer(1)
    for fac, mult in sp.factor_list(sp.expand(num))[1]:
        if fac not in (z, z2):
            core *= fac**mult
    return sp.expand(core)
txt = lambda p_: str(sp.expand(p_)).replace('**', '^')
vars_ = [str(u) for u in unk] + ['ww'] + list(obs)
allp = polys + [sp.expand(sp.Symbol('ww') * z * (1 if mode == 'sym' else z2) - 1)] + ([] if mode == 'sym' else [sp.expand(QQ * (z - z2) - 1)]) + [rel(nm) for nm in obs]
elim_core = '*'.join(str(u) for u in unk) + '*ww'
lines = ['ring R = 0,(%s),dp;' % ','.join(vars_), 'option(redSB);', 'ideal I =', ',\n'.join(txt(p_) for p_ in allp) + ';',
         'ideal G = std(I);', 'print("dim:"); dim(G);', 'print("vdim:"); vdim(G);']
for nm in obs:
    others = '*'.join(o for o in obs if o != nm)
    lines += ['ideal J%s = eliminate(G, %s*%s);' % (nm, elim_core, others), 'print("eliminant %s:"); J%s;' % (nm, nm)]
zel = '*'.join([str(u) for u in unk if u != z] + ['ww'] + list(obs))
lines += ['ideal Jz = eliminate(G, %s);' % zel, 'print("eliminant z:"); Jz;', 'quit;']
fn = 'h_two_%s_%s.sing' % (mode, TAG)
open(fn, 'w').write('\n'.join(lines) + '\n')
print('written %s (%.0fs); running Singular' % (fn, time.time() - t0), flush=True)
res = subprocess.run(['<home>/ib-lab/solvers/env/bin/Singular', '-q'], stdin=open(fn), capture_output=True, text=True, timeout=3000).stdout
open(fn.replace('.sing', '.out'), 'w').write(res)
print(res[:3000])
print('done (%.0fs)' % (time.time() - t0))
