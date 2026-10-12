"""EXC4 B: level two for Solution 3 -- FOUR apparent singularities at xi = +-z1, +-z3, type {-1,1,2,4}, sigma-invariant pairs
(EXC2 normal form (th^2 - A)(th^2 - B) + xi(th + 1/2) - E xi^b; sigma relations within each pair from EXC2: s21 = r21 = -b, s22 = +4z,
s13 = 8z^2, s04 = 8z^3, s11 = -r11, s12 = r12 - 2(b+4)z; s01, s02, s03 free and eliminated by the point's own conditions).
Each point's local expansion includes the other three points' regular expansions (frob.py machinery, as EXC2's h_two_v2).
Linear eliminations per point (r03 = -z(r12 + 4z) from the kap2 condition at the resonance 2; r02 and r01 from the conditions linear in
them; the same for the s's) leave (r11, r12, z) per pair; the rest goes to Singular with the observables I_1 (ii) and I_3 (ee) from the
effective quartic (q2 = -(A+B) + sum r21, q3 = sum r11, q4 = AB + sum r01; cB2, cB4 of EXC2's g_charges_t2.json at t = 2).
Usage: h_four_sol3.py <t> <P> <Q> [mode]   mode: four (default) | pair (ONE sigma pair: EXC2 regression, must give the EXC2 quadratic) | tamper (xi(th+1)).
Writes h_four_<TAG>.sing/.out and prints dim, vdim, eliminants."""
import hashlib, json, re, subprocess, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import frob as Fb
assert hashlib.sha256(open('SEAL_EXC4.md', 'rb').read()).hexdigest() == open('SEAL_EXC4.sha256').read().split()[0]
tstr, Ps, Qs = sys.argv[1], sys.argv[2], sys.argv[3]; mode = sys.argv[4] if len(sys.argv) > 4 else 'four'
t = sp.Rational(tstr); k = 2 * (t - 1) / (3 - t); b = k / (k + 1); n = k + 4; rho2 = 2 / (t * t - 1)
Pv, Qv = sp.Rational(Ps), sp.Rational(Qs); A = b * Pv**2; B = b * (Qv**2 + rho2)
TAG = 't%s_P%s_Q%s_%s' % (tstr.replace('/', '_'), Ps.replace('/', '-'), Qs.replace('/', '-'), mode)
z, E, e = Fb.z, Fb.E, Fb.e
z3 = sp.Symbol('z3')
names = ['21', '22', '11', '12', '13', '01', '02', '03', '04']
PTS = []     # (position symbol/expr, residue dict, tag)
r = {nm: sp.Symbol('r' + nm) for nm in names}; s = {nm: sp.Symbol('s' + nm) for nm in names}
u = {nm: sp.Symbol('u' + nm) for nm in names}; v = {nm: sp.Symbol('v' + nm) for nm in names}
PTS = [(z, r, 'A'), (-z, s, 'B')] + ([] if mode == 'pair' else [(z3, u, 'C'), (-z3, v, 'D')])
t0 = time.time()
th = sp.Symbol('th')
Pcoef = [A * B, 0, -(A + B), 0, 1]
Lcoef = [sp.Integer(1) if mode == 'tamper' else sp.Rational(1, 2), 1]
exps = [-1, 1, 2, 4]; K = 2


def local_at(idx):
    zz, rr, tag = PTS[idx]
    order = 4
    powers = [{0: sp.Integer(1)}]
    def theta_on(series):
        out = {}
        for kk, g in series.items():
            f = sp.expand(g * (e + kk))
            out[kk - 1] = sp.expand(out.get(kk - 1, 0) + zz * f)
            out[kk] = sp.expand(out.get(kk, 0) + f)
        return {kk: vv for kk, vv in out.items() if vv != 0}
    for j in range(order):
        powers.append(Fb.trunc(theta_on(powers[-1]), K))
    xi = {0: zz, 1: sp.Integer(1)}
    out = {}
    for j, cc in enumerate(Pcoef):
        out = Fb.add_series(out, {kk: cc * vv for kk, vv in powers[j].items()})
    for j, cc in enumerate(Lcoef):
        out = Fb.add_series(out, Fb.mul_series({kk: cc * vv for kk, vv in powers[j].items()}, xi))
    xib = {i: sp.binomial(b, i) / zz**i for i in range(0, K + order + 1)}
    out = Fb.add_series(out, {kk: -E * sp.expand_func(vv) for kk, vv in xib.items()})
    for nm, cc in rr.items():
        j, m = int(nm[0]), int(nm[1])
        out = Fb.add_series(out, Fb.mul_series(Fb.mul_series(powers[j], xi), {-m: cc}))
    for jdx, (zo, ss, tago) in enumerate(PTS):
        if jdx == idx:
            continue
        qq = 1 / (zz - zo)
        for nm, cc in ss.items():
            j, m = int(nm[0]), int(nm[1])
            ser = {i: sp.binomial(-m, i) * qq ** (m + i) for i in range(0, K + order + 1)}
            out = Fb.add_series(out, Fb.mul_series(Fb.mul_series(powers[j], xi), {i: cc * vv for i, vv in ser.items()}))
    return Fb.trunc(out, K)


# sigma relations and indicial values
def fixed(rr, ss, zz):
    d = {rr['21']: -b, rr['22']: -4 * zz, rr['13']: 8 * zz**2, rr['04']: -8 * zz**3,
         ss['21']: -b, ss['22']: 4 * zz, ss['13']: 8 * zz**2, ss['04']: 8 * zz**3,
         ss['11']: -rr['11'], ss['12']: rr['12'] - 2 * (b + 4) * zz}
    d[rr['03']] = -zz * (d.get(rr['12'], rr['12']) + 4 * zz)
    d[ss['03']] = zz * (d[ss['12']] - 4 * zz)       # the kap2 condition at the mirrored point (z -> -z)
    return d
sym = fixed(r, s, z)
if mode != 'pair':
    sym.update(fixed(u, v, z3))
polys = []
for idx in range(len(PTS)):
    lk = local_at(idx)
    lead, conds, kaps = Fb.nolog_conditions(lk, exps, 4)
    zz = PTS[idx][0]
    assert sp.expand((lead - zz**4 * sp.prod([(e - x) for x in exps])).subs(sym)) == 0, 'indicial at point %d' % idx
    cnt = 0
    for pos, expr in conds:
        num = sp.numer(sp.together(sp.expand(expr.subs(sym).subs(sym))))
        Pp = sp.Poly(sp.expand(num), E, *kaps)
        for mon, cf in Pp.terms():
            core = sp.Integer(1)
            for fac, mult in sp.factor_list(sp.expand(cf))[1]:
                if fac not in (z, z3):
                    core *= fac**mult
            core = sp.expand(core)
            if core != 0 and core not in polys and sp.expand(-core) not in polys:
                polys.append(core); cnt += 1
    print('point %s: %d new conditions (%.0fs)' % (PTS[idx][2], cnt, time.time() - t0), flush=True)
unk = [r['11'], r['12'], r['01'], r['02'], s['01'], s['02'], z] + ([] if mode == 'pair' else [u['11'], u['12'], u['01'], u['02'], v['01'], v['02'], z3])
# linear eliminations of r01, r02, s01, s02 (u, v): pick equations linear in the variable with a coefficient free of the other unknowns
sol = {}
def try_lin(var):
    for c in polys:
        cc = sp.expand(c.subs(sol))
        if cc == 0 or not cc.has(var):
            continue
        p_ = sp.Poly(cc, var)
        if p_.degree() == 1 and not (p_.LC().free_symbols & (set(unk) - {z, z3})):
            vv = sp.factor(-p_.coeff_monomial(1) / p_.LC())
            for k_ in list(sol):
                sol[k_] = sp.factor(sol[k_].subs(var, vv))
            sol[var] = vv
            return True
    return False
for var in [r['02'], s['02'], r['01'], s['01']] + ([] if mode == 'pair' else [u['02'], v['02'], u['01'], v['01']]):
    ok_ = try_lin(var)
    print('   %s %s' % (var, 'eliminated' if ok_ else 'NOT linear -> kept'), flush=True)
rest = []
for c in polys:
    cc = sp.factor(sp.numer(sp.together(c.subs(sol))))
    if cc == 0:
        continue
    core = sp.Integer(1)
    for fac, mult in sp.factor_list(cc)[1]:
        if fac not in (z, z3) and (fac.free_symbols & set(unk)):
            core *= fac
    core = sp.expand(core)
    if core != 1 and core not in rest and sp.expand(-core) not in rest:
        rest.append(core)
remaining = [x for x in unk if x not in sol]
print('%d equations in %d unknowns %s (%.0fs)' % (len(rest), len(remaining), remaining, time.time() - t0), flush=True)
# observables from the effective quartic (EXC2 engine at t = 2)
gc = json.load(open('<home>/fable-work/exc2/g_charges_t%s.json' % tstr.replace('/', '_')))
A2, A3, A4 = sp.symbols('A2 A3 A4'); loc = {'A2': A2, 'A3': A3, 'A4': A4}
cB2 = sp.sympify(gc['cB2'], locals=loc); cB4 = sp.sympify(gc['cB4'], locals=loc)
nb = n / b; x0s, x1s = sp.symbols('x0s x1s')
lead4 = sp.Poly(sp.expand(cB4.subs({A2: nb**2 * (-(x0s + x1s)), A3: 0, A4: nb**4 * x0s * x1s})), x0s, x1s).coeff_monomial(x0s**2) * b**2
lead2 = sp.Poly(sp.expand(cB2.subs({A2: nb**2 * (-(x0s + x1s)), A3: 0, A4: nb**4 * x0s * x1s})), x0s, x1s).coeff_monomial(x0s) * b
res_all = [rr for _, rr, _ in PTS]
q2 = -(A + B) + sum(rr['21'] for rr in res_all); q3 = sum(rr['11'] for rr in res_all); q4 = A * B + sum(rr['01'] for rr in res_all)
val = {A2: nb**2 * q2, A3: -nb**3 * q3, A4: nb**4 * q4}
I1 = sp.together(cB2.subs(val).subs(sym).subs(sol) / lead2); I3 = sp.together(cB4.subs(val).subs(sym).subs(sol) / lead4)
vac2 = sp.expand(cB2.subs({A2: nb**2 * (-(A + B)), A3: 0, A4: nb**4 * A * B}) / lead2)
ii, ee = sp.symbols('ii ee')
def rel(symb, expr):
    num, den = sp.fraction(sp.together(symb - expr)); core = sp.Integer(1)
    for fac, mult in sp.factor_list(sp.expand(num))[1]:
        if fac not in (z, z3):
            core *= fac**mult
    return sp.expand(core)
ww = sp.Symbol('ww')
allp = rest + [sp.expand(ww * z * (1 if mode == 'pair' else z3) - 1), rel(ii, I1), rel(ee, I3)]
if mode != 'pair':
    allp.append(sp.expand(sp.Symbol('vv') * (z**2 - z3**2) - 1))       # distinct pairs
vars_ = [str(x) for x in remaining] + ['ww'] + ([] if mode == 'pair' else ['vv']) + ['ii', 'ee']
txt = lambda p_: str(sp.expand(p_)).replace('**', '^')
elim_core = '*'.join(str(x) for x in remaining) + '*ww' + ('' if mode == 'pair' else '*vv')
lines = ['ring R = 0,(%s),dp;' % ','.join(vars_), 'option(redSB);', 'ideal I =', ',\n'.join(txt(p_) for p_ in allp) + ';',
         'ideal G = std(I);', 'print("dim:"); dim(G);', 'print("vdim:"); vdim(G);',
         'ideal Ji = eliminate(G, %s*ee);' % elim_core, 'print("eliminant ii:"); Ji;',
         'ideal Je = eliminate(G, %s*ii);' % elim_core, 'print("eliminant ee:"); Je;',
         'ideal Jz = eliminate(G, %s*ii*ee);' % '*'.join([str(x) for x in remaining if x != z] + ['ww'] + ([] if mode == 'pair' else ['vv'])), 'print("eliminant z:"); Jz;', 'quit;']
fn = 'h_four_%s.sing' % TAG
open(fn, 'w').write('\n'.join(lines) + '\n')
print('written %s; vacuum monic I_1 = %s; running Singular (%.0fs)' % (fn, vac2, time.time() - t0), flush=True)
res = subprocess.run(['<home>/ib-lab/solvers/env/bin/Singular', '-q'], stdin=open(fn), capture_output=True, text=True, timeout=21000).stdout
open(fn.replace('.sing', '.out'), 'w').write(res)
print(res[:4000]); print('done (%.0fs)' % (time.time() - t0))
