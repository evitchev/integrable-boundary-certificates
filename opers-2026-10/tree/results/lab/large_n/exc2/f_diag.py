"""EXC2 T1 diagnosis for the negative exponent sets: print the E-dependent no-logarithm conditions explicitly."""
import sys
import sympy as sp
sys.dont_write_bytecode = True
import frob as Fb
t = sp.Integer(2); k = 2 * (t - 1) / (3 - t); b = k / (k + 1)
a0, a1 = sp.Rational(1, 2), sp.Rational(1, 3)
z, E, e = Fb.z, Fb.E, Fb.e
th = sp.Symbol('th')
r = {(j, m): sp.Symbol('r%d%d' % (j, m)) for j in (2, 1, 0) for m in range(1, 5 - j)}
R = {j: {m: r[j, m] for m in range(1, 5 - j)} for j in (2, 1, 0)}
P = sp.Poly(sp.expand((th**2 - a0**2) * (th**2 - a1**2)), th)
Pcoef = [P.coeff_monomial(th**j) for j in range(5)]
for exps in ([-1, 1, 2, 4], [-1, 0, 3, 4], [-1, 0, 2, 5]):
    lk = Fb.local_operator(Pcoef, [sp.Rational(1, 2), 1], b, -1, R, exps[-1] - exps[0] - 3)
    lead, conds, kaps = Fb.nolog_conditions(lk, exps, 4)
    sol_ind = sp.solve([c for c in sp.Poly(sp.expand(lead - z**4 * sp.prod([(e - x) for x in exps])), e).coeffs() if c != 0], [r[2, 2], r[1, 3], r[0, 4]], dict=True)[0]
    eqs = Fb.split_conditions([(p_, sp.expand(c_.subs(sol_ind))) for p_, c_ in conds], kaps)
    print('exponents %s: indicial %s' % (exps, sol_ind))
    for pos, mon, cf in eqs:
        if mon[0] >= 1:
            kn = ' '.join('%s^%d' % (s, p) for s, p in zip(['E'] + [str(x) for x in kaps], mon) if p)
            print('   resonance at exponent %d, coefficient of [%s]:  %s = 0' % (pos, kn, cf if len(str(cf)) < 300 else str(cf)[:300] + ' ...'))
