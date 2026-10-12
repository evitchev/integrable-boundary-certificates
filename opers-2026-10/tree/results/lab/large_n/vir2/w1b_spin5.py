"""VIR2 step b: the I_5 system at a few fibres -- which I_3 solutions survive."""
import sys, time
import sympy as sp
sys.dont_write_bytecode = True
from w1_lib import *
t0 = time.time()
CAND = {1: [('sealed S1', lambda u_: {q: 1 - u_**2 / 2, g: -u_ / 2})],
        2: [('found at spin 3: a^2 = 1 - 2u^2, beta = -sqrt2 u', lambda u_: {q: 1 - 2 * u_**2, g: -u_})],
        3: [('sealed S3', lambda u_: {q: 1 / (2 * u_**2), g: 1 / (2 * u_)}),
            ('second I_3 solution: a^2 = 2u^2, beta = -sqrt2 u', lambda u_: {q: 2 * u_**2, g: -u_})]}
for uval in (sp.Rational(3), sp.Rational(2, 5), sp.Rational(-7, 3)):
    t, N, s, r = fibre(uval)
    print('\n=== u = %s: t = %s, N = %s' % (uval, t, N), flush=True)
    for sol in (1, 2, 3):
        d3 = I3_two_boson(sol, uval)
        d5 = I5_two_boson(sol, uval)
        sys3 = system(d3, r)
        sys5 = system(d5, r)
        print('Sol %d: I_5 two-boson image has %d monomials; %d conditions  (%.1fs)' % (sol, len(d5), len(sys5), time.time() - t0), flush=True)
        for name, f in CAND[sol]:
            sub = f(uval)
            ok3 = all(sp.expand(e.subs(sub)) == 0 for e in sys3)
            ok5 = all(sp.expand(e.subs(sub)) == 0 for e in sys5)
            print('   %-55s I_3: %s   I_5: %s' % (name, ok3, ok5), flush=True)
        G = sp.groebner(sys3 + sys5, q, g, order='lex')
        print('   Groebner basis of the joint {I_3, I_5} system (a != 0, beta != 0):', [str(sp.factor(e)) for e in G.exprs], flush=True)
        # a = 0 control
        k3 = reduce_conditions(W.conditions(d3), r, keep_a=True)
        k5 = reduce_conditions(W.conditions(d5), r, keep_a=True)
        g0 = sp.factor(sp.gcd_list([sp.expand(e.subs(a, 0)) for e in k3 + k5]))
        print('   a = 0: gcd over I_3 and I_5 conditions: %s   (g = u: %s, g = -1/u: %s)' % (g0, uval, -1 / uval), flush=True)
print('done %.1fs' % (time.time() - t0))
