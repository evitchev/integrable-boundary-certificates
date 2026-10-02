"""EXC2 (after the hold-out; labelled post-hoc): the SELF-ADJOINT two-point system in closed form, symbolic in
A = a0^2, B = a1^2 and b.  Points at xi = +z and xi = -z, exponents {-1,1,2,4}; adjoint symmetry (derived by hand, checked
here by the equations themselves): s21 = r21 = -b, s11 = -r11, s12 = r12 - 2(b+4) z, s22 = +4z, s13 = 8z^2, s04 = 8 z^3.
All 14 no-log equations are generated; the linear ones are solved; the rest are reduced to polynomials in (rho, zeta = z^2),
r12 = rho z.  Output s_sym.json."""
import json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import frob as Fb
A, B, b = sp.symbols('A B b')
z, E, e = Fb.z, Fb.E, Fb.e
th = sp.Symbol('th')
names = ['21', '22', '11', '12', '13', '01', '02', '03', '04']
r = {nm: sp.Symbol('r' + nm) for nm in names}
s = {nm: sp.Symbol('s' + nm) for nm in names}
t0 = time.time()
Pcoef = [A * B, 0, -(A + B), 0, 1]
Lcoef = [sp.Rational(1, 2), 1]
K = 2


def local_two(rr, ss, zz, zo):
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
    out = Fb.add_series(out, {kk: -E * sp.expand_func(v) for kk, v in xib.items()})
    for nm, c in rr.items():
        j, m = int(nm[0]), int(nm[1])
        out = Fb.add_series(out, Fb.mul_series(Fb.mul_series(powers[j], xi), {-m: c}))
    qq = 1 / (zz - zo)
    for nm, c in ss.items():
        j, m = int(nm[0]), int(nm[1])
        ser = {i: sp.binomial(-m, i) * qq ** (m + i) for i in range(0, K + order + 1)}
        out = Fb.add_series(out, Fb.mul_series(Fb.mul_series(powers[j], xi), {i: c * v for i, v in ser.items()}))
    return Fb.trunc(out, K)


exps = [-1, 1, 2, 4]
rho = sp.Symbol('rho'); zeta = sp.Symbol('zeta')
sym = {r['21']: -b, s['21']: -b, r['22']: -4 * z, s['22']: 4 * z, r['13']: 8 * z**2, s['13']: 8 * z**2, r['04']: -8 * z**3, s['04']: 8 * z**3,
       s['11']: -r['11'], r['12']: rho * z, s['12']: rho * z - 2 * (b + 4) * z}
sym[r['03']] = -z * (sym[r['12']] + 4 * z)                      # kap2 coefficient at the resonance 2 (local)
sym[s['03']] = z * (sym[s['12']] - 4 * z)                       # the same at xi = -z
eqs = {}
for (rr, ss, zz, zo, tag) in ((r, s, z, -z, 'A'), (s, r, -z, z, 'B')):
    lk = local_two(rr, ss, zz, zo)
    lead, conds, kaps = Fb.nolog_conditions(lk, exps, 4)
    out = []
    for pos, expr in conds:
        num = sp.numer(sp.together(sp.expand(expr.subs(sym).subs(sym))))
        Pp = sp.Poly(sp.expand(num), E, *kaps)
        for mon, cf in Pp.terms():
            out.append((pos, mon, sp.factor(cf)))
    eqs[tag] = out
    print('point %s: %d conditions (%.0fs)' % (tag, len(out), time.time() - t0), flush=True)
    for pos, mon, cf in out:
        print('   res %d %s: free symbols %s, length %d' % (pos, mon, sorted(str(x) for x in cf.free_symbols), len(str(cf))))
U = {'r11': r['11'], 'r01': r['01'], 'r02': r['02'], 's01': s['01'], 's02': s['02']}
EA = [sp.numer(sp.together(cf)) for _, _, cf in eqs['A']]; EB = [sp.numer(sp.together(cf)) for _, _, cf in eqs['B']]
for tag, L in (('A', EA), ('B', EB)):
    for i_, c in enumerate(L):
        print('%s%d: %s' % (tag, i_, sp.factor(c)))
sol = {}
def solve_lin(c, sym_):
    p_ = sp.Poly(sp.expand(sp.numer(sp.together(c.subs(sol)))), sym_)
    assert p_.degree() == 1, (sym_, p_.degree())
    v = sp.factor(-p_.coeff_monomial(1) / p_.LC())
    for k_ in list(sol):
        sol[k_] = sp.factor(sol[k_].subs(sym_, v))
    sol[sym_] = v
    print('%s = %s' % (sym_, v), flush=True)
solve_lin(EA[0], r['02']); solve_lin(EB[0], s['02']); solve_lin(EA[1], r['01']); solve_lin(EB[1], s['01'])
rest = []
for c in EA[2:] + EB[2:]:
    cc = sp.factor(sp.numer(sp.together(c.subs(sol))))
    if cc == 0:
        print('   (identically satisfied)'); continue
    core = sp.Integer(1)
    for fac, mult in sp.factor_list(cc)[1]:
        if fac != z and (fac.free_symbols & {rho, z, r['11']}):
            core *= fac
    core = sp.expand(core)
    if core not in rest and sp.expand(-core) not in rest:
        rest.append(core)
print('remaining distinct equations in (r11, rho, z): %d (%.0fs)' % (len(rest), time.time() - t0))
for c in rest:
    pz = sp.Poly(c, r['11'], rho, z)
    print('   degrees (r11, rho, z) = %s; %d terms: %s' % ((pz.degree(r['11']), pz.degree(rho), pz.degree(z)), len(sp.Poly(c, r['11'], rho, z, A, B, b).terms()), sp.collect(c, [r['11'], rho, z]) if len(str(c)) < 1500 else '...'))
json.dump({'sol': {str(k_): str(v) for k_, v in sol.items()}, 'sym': {str(k_): str(v) for k_, v in sym.items()}, 'rest': [str(c) for c in rest]}, open('s_sym.json', 'w'), indent=1)
