"""EXC4 A1 (second part): TWO apparent singularities in the sixth-order form of Solution 3 at t = 2, at free positions z1, z2, each of type
S_6 = {-1,1,2,3,4,6}; residues free (no symmetry imposed; the adjoint invariance is tested afterwards through kk = sum of the odd parts).
Unknowns: 2 x 16 + 2 = 34 (after the indicial conditions); equations 2 x 16.  Singular (no short timeout; PID saved by the caller).
Observables: ii (monic I_1), ee (monic I_3) from the effective sextic (sum over both points of r_j1).  The sealed prediction: solutions with
z2 != -z1 exist and their I_3 values contain the level-one data quadratic's roots (I_1 = vacuum + 2 on them).
Usage: f_six_two.py <a0> <a1>"""
import hashlib, json, re, subprocess, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import frob as Fb
import mellin_gen6 as MG
assert hashlib.sha256(open('SEAL_EXC4.md', 'rb').read()).hexdigest() == open('SEAL_EXC4.sha256').read().split()[0]
a0s, a1s = sys.argv[1], sys.argv[2]; a0 = sp.Rational(a0s); a1 = sp.Rational(a1s); A, B = a0**2, a1**2
t = sp.Integer(2); k = sp.Integer(2); b = k / (k + 1); n = k + 4; M = 1 / k; c = n * (1 + M); f = sp.Rational(3, 2) / c; rho2 = 2 / (t * t - 1)
t0 = time.time()
th, T = sp.symbols('th T')
z, E, e = Fb.z, Fb.E, Fb.e
z2 = sp.Symbol('z2'); QQ = sp.Symbol('qq')
names = [(j, m) for j in range(0, 5) for m in range(1, 7 - j)]
r = {(j, m): sp.Symbol('r%d%d' % (j, m)) for (j, m) in names}; s = {(j, m): sp.Symbol('s%d%d' % (j, m)) for (j, m) in names}
Q6 = sp.Poly(sp.expand((th**2 - A) * (th**2 - B) * (th**2 - f**2)), th)
Pcoef = [Q6.coeff_monomial(th**j) for j in range(7)]
Lcoef = [sp.Integer(-1)]
exps = [-1, 1, 2, 3, 4, 6]; K = exps[-1] - exps[0] - 6 + 1


def local_two(rr, ss, zz, zo, sign):
    order = 6
    powers = [{0: sp.Integer(1)}]
    def theta_on(series):
        out = {}
        for kk, g in series.items():
            fx = sp.expand(g * (e + kk))
            out[kk - 1] = sp.expand(out.get(kk - 1, 0) + zz * fx)
            out[kk] = sp.expand(out.get(kk, 0) + fx)
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
    for (j, m), cc in rr.items():
        out = Fb.add_series(out, Fb.mul_series(Fb.mul_series(powers[j], xi), {-m: cc}))
    qq = sign * QQ          # 1/(zz - zo)
    for (j, m), cc in ss.items():
        ser = {i: sp.binomial(-m, i) * qq ** (m + i) for i in range(0, K + order + 1)}
        out = Fb.add_series(out, Fb.mul_series(Fb.mul_series(powers[j], xi), {i: cc * v for i, v in ser.items()}))
    return Fb.trunc(out, K)


ind = {}
for rr, zz in ((r, z), (s, z2)):
    ind.update({rr[4, 2]: -6 * zz, rr[3, 3]: 24 * zz**2, rr[2, 4]: -72 * zz**3, rr[1, 5]: 144 * zz**4, rr[0, 6]: -144 * zz**5})
polys = []
for (rr, ss, zz, zo, sign, tag) in ((r, s, z, z2, 1, 'A'), (s, r, z2, z, -1, 'B')):
    lk = local_two(rr, ss, zz, zo, sign)
    lead, conds, kaps = Fb.nolog_conditions(lk, exps, 6)
    assert sp.expand((lead - zz**6 * sp.prod([(e - x) for x in exps])).subs(ind)) == 0
    for pos, expr in conds:
        num = sp.numer(sp.together(sp.expand(expr.subs(ind))))
        Pp = sp.Poly(sp.expand(num), E, *kaps)
        for mon, cf in Pp.terms():
            core = sp.Integer(1)
            for fac, mult in sp.factor_list(sp.expand(cf))[1]:
                if fac not in (z, z2):
                    core *= fac**mult
            polys.append(sp.expand(core))
    print('point %s done (%.0fs)' % (tag, time.time() - t0), flush=True)
unk = [x for x in r.values() if x not in ind] + [z] + [x for x in s.values() if x not in ind] + [z2, QQ]
print('%d equations in %d unknowns (%.0fs)' % (len(polys), len(unk), time.time() - t0), flush=True)
# observables
hs = [MG.h1, MG.h2, MG.h3, MG.h4, MG.h5, MG.h6]
eng = MG.MellinWKB(n, M, [MG.ONE] + [MG.P_(h) for h in hs][:4]); eng.solve(4)
cB = {i: sp.expand(eng.R(i)[0]) for i in (2, 4)}
P, Q = sp.symbols('P Q')
Hv = sp.Poly(sp.expand((T**2 - n**2 * P**2 / b) * (T**2 - n**2 * (Q**2 + rho2) / b) * (T**2 - sp.Rational(9, 4))), T)
hv = {hs[j - 1]: Hv.coeff_monomial(T**(6 - j)) for j in range(1, 7)}
lead = {i: sp.Poly(sp.expand(cB[i].subs(hv)), P, Q).coeff_monomial(P**i) for i in (2, 4)}
Qeff = Q6.as_expr() + sum((r[j, 1] + s[j, 1]) * th**j for j in range(0, 5))
He = sp.Poly(sp.expand(c**6 * Qeff.subs(th, -T / c)), T)
he = {hs[j - 1]: He.coeff_monomial(T**(6 - j)) for j in range(1, 7)}
I1 = sp.expand(cB[2].subs(he) / lead[2]); I3 = sp.expand(cB[4].subs(he) / lead[4])
ii, ee, kk, ww = sp.symbols('ii ee kk ww')
def txt(p_):
    # integer coefficients (Singular parses 'x^2/4' as x^(2/4)): clear denominators
    pp = sp.Poly(sp.expand(p_), *sorted(sp.expand(p_).free_symbols, key=str)) if sp.expand(p_).free_symbols else sp.Poly(sp.expand(p_), sp.Symbol('dummy'))
    _, pc = pp.clear_denoms()
    return str(pc.as_expr()).replace('**', '^')
vars_ = [str(x) for x in unk] + ['ww', 'ii', 'ee', 'kk']
allp = polys + [sp.expand(QQ * (z - z2) - 1), sp.expand(ww * z * z2 - 1), sp.expand(ii - I1), sp.expand(ee - I3), sp.expand(kk - (z + z2))]
core = '*'.join(str(x) for x in unk) + '*ww'
lines = ['ring R = 0,(%s),dp;' % ','.join(vars_), 'option(redSB);', 'ideal I =', ',\n'.join(txt(p_) for p_ in allp) + ';',
         'ideal G = std(I);', 'print("dim:"); dim(G);', 'print("vdim:"); vdim(G);',
         'ideal Ji = eliminate(G, %s*ee*kk);' % core, 'print("eliminant ii:"); Ji;',
         'ideal Je = eliminate(G, %s*ii*kk);' % core, 'print("eliminant ee:"); Je;',
         'ideal Jk = eliminate(G, %s*ii*ee);' % core, 'print("eliminant kk:"); Jk;', 'quit;']
fn = 'f_six_two_%s_%s.sing' % (a0s.replace('/', '-'), a1s.replace('/', '-'))
open(fn, 'w').write('\n'.join(lines) + '\n')
print('written %s (%.0fs); running Singular' % (fn, time.time() - t0), flush=True)
res = subprocess.run(['<home>/ib-lab/solvers/env/bin/Singular', '-q'], stdin=open(fn), capture_output=True, text=True).stdout
open(fn.replace('.sing', '.out'), 'w').write(res)
print(res[:3000]); print('done (%.0fs)' % (time.time() - t0))
