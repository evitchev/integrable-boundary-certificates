"""EXC3 (post-hoc, no data): the sigma-invariant two-point system of Solution 2 in closed form, symbolic in A = a0^2, B = a1^2, b.
Points at xi = +-z, exponents {-1,1,3}; indicial: r12 = -3z, r03 = 3z^2, s12 = 3z, s03 = 3z^2; sigma: s11 = r11 (checked on the
numeric solutions: kk = 0).  All no-log equations of both points generated with frob.py; linear eliminations; the rest printed.
Output s_sym_sol2.json."""
import json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import frob as Fb
A, B, b = sp.symbols('A B b')
z, E, e = Fb.z, Fb.E, Fb.e
names = ['11', '12', '01', '02', '03']
r = {nm: sp.Symbol('r' + nm) for nm in names}; s = {nm: sp.Symbol('s' + nm) for nm in names}
t0 = time.time()
Pcoef = [0, -B, 0, 1]; Lcoef = [sp.Rational(1, 4) - A, 1, 1]
exps = [-1, 1, 3]; K = 2


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
    qq = 1 / (zz - zo)
    for nm, cc in ss.items():
        j, m = int(nm[0]), int(nm[1])
        ser = {i: sp.binomial(-m, i) * qq ** (m + i) for i in range(0, K + order + 1)}
        out = Fb.add_series(out, Fb.mul_series(Fb.mul_series(powers[j], xi), {i: cc * v for i, v in ser.items()}))
    return Fb.trunc(out, K)


sym = {r['12']: -3 * z, r['03']: 3 * z**2, s['12']: 3 * z, s['03']: 3 * z**2, s['11']: r['11']}
eqs = {}
for (rr, ss, zz, zo, tag) in ((r, s, z, -z, 'A'), (s, r, -z, z, 'B')):
    lk = local_two(rr, ss, zz, zo)
    lead, conds, kaps = Fb.nolog_conditions(lk, exps, 3)
    out = []
    for pos, expr in conds:
        num = sp.numer(sp.together(sp.expand(expr.subs(sym).subs(sym))))
        Pp = sp.Poly(sp.expand(num), E, *kaps)
        for mon, cf in Pp.terms():
            core = sp.Integer(1)
            for fac, mult in sp.factor_list(sp.expand(cf))[1]:
                if fac != z:
                    core *= fac**mult
            out.append((pos, mon, sp.expand(core)))
    eqs[tag] = out
    print('point %s: %d conditions (%.0fs)' % (tag, len(out), time.time() - t0), flush=True)
    for pos, mon, cf in out:
        print('   res %d %s: %s' % (pos, mon, sp.factor(cf) if len(str(cf)) < 300 else '%d terms' % len(sp.Poly(cf).terms())))
unk = [r['11'], r['01'], r['02'], s['01'], s['02'], z]
allp = [cf for tag in 'AB' for _, _, cf in eqs[tag] if cf != 0]
sol = {}
def solve_lin(sym_):
    for c in allp:
        cc = sp.expand(c.subs(sol))
        if cc == 0:
            continue
        p_ = sp.Poly(cc, sym_)
        if p_.degree() == 1 and not (p_.LC().free_symbols & set(unk)):
            v = sp.factor(-p_.coeff_monomial(1) / p_.LC())
            for k_ in list(sol):
                sol[k_] = sp.factor(sol[k_].subs(sym_, v))
            sol[sym_] = v
            print('%s = %s' % (sym_, v), flush=True)
            return True
    return False
for sym_ in (r['02'], s['02'], r['01'], s['01'], r['11']):
    solve_lin(sym_)
rest = []
for c in allp:
    cc = sp.factor(sp.numer(sp.together(c.subs(sol))))
    if cc == 0:
        continue
    core = sp.Integer(1)
    for fac, mult in sp.factor_list(cc)[1]:
        if fac != z and (fac.free_symbols & set(unk)):
            core *= fac
    core = sp.expand(core)
    if core not in rest and sp.expand(-core) not in rest and core != 1:
        rest.append(core)
print('remaining distinct equations: %d (%.0fs)' % (len(rest), time.time() - t0))
for c in rest:
    print('   ', sp.collect(c, [r['11'], z], sp.factor) if len(str(c)) < 1500 else '%d terms' % len(sp.Poly(c).terms()))
json.dump({'sol': {str(k_): str(v) for k_, v in sol.items()}, 'sym': {str(k_): str(v) for k_, v in sym.items()}, 'rest': [str(c) for c in rest]}, open('s_sym_sol2.json', 'w'), indent=1)
