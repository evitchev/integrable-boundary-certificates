"""VIR2 step a: the I_3 system at a few fibres -- size, sealed pairs, full solution sets."""
import sys, time
import sympy as sp
sys.dont_write_bytecode = True
from w1_lib import *
t0 = time.time()
for uval in (sp.Rational(3), sp.Rational(2, 5), sp.Rational(-7, 3)):
    t, N, s, r = fibre(uval)
    print('\n=== u = %s: t = %s, N = %s, s = %s' % (uval, t, N, s))
    for sol in (1, 2, 3):
        d3 = I3_two_boson(sol, uval)
        sys3 = system(d3, r)
        print('Sol %d: I_3 two-boson image has %d monomials; %d conditions  (%.1fs)' % (sol, len(d3), len(sys3), time.time() - t0))
        sealed = SEALED(sol, uval)
        if sealed:
            print('   sealed pair satisfies all I_3 conditions: %s' % all(sp.expand(e.subs(sealed)) == 0 for e in sys3))
        # a = 0 solutions
        a0 = [sp.factor(sp.gcd_list([sp.expand(e.subs(q, 0)) for e in sys3]))]
        print('   a = 0: gcd of conditions in g: %s   (Virasoro: g = u = %s, g = -1/u = %s)' % (a0[0], uval, -1 / uval))
        G = sp.groebner(sys3 + [], q, g, order='lex')
        print('   Groebner basis (lex q > g) of the I_3 system:')
        for e in G.exprs:
            print('      ', sp.factor(e))
print('done %.1fs' % (time.time() - t0))
