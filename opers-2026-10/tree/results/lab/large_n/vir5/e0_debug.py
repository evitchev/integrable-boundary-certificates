import sys
import sympy as sp
sys.dont_write_bytecode = True
import mellin_wkb as MW
import wkb_lib as WL
l0 = MW.l0
# BLZ: (theta - 1/2 - l)(theta - 1/2 + l) psi = x^2 (x^{2M} + E) psi  ->  Ht(T) = T^2 - l^2, s = 1/2
for Mv in (sp.Integer(1), sp.Rational(1, 3)):
    eng = MW.MellinWKB(2, Mv, [MW.P_(1), MW.P_(-l0**2)])
    eng.solve(4)
    for i in (1, 2, 3, 4):
        print('M = %s: W_%d keys %s' % (Mv, i, sorted(eng.W[i].keys())))
    for i in (2, 3, 4):
        Ri, ok = eng.R(i)
        print('   Mellin R_%d = %s (homog %s)' % (i, sp.factor(Ri), ok))
    pr = WL.Problem([WL.g0, WL.g1 - 1], 'P', -1, 2, (WL.F(1, 2), -1)); pr.solve(4)
    for i in (2, 3, 4):
        Rk = pr.R(i)[0].subs({WL.g0: sp.Rational(1, 2) - WL.l0, WL.g1: sp.Rational(1, 2) + WL.l0}, simultaneous=True).subs(WL.M, Mv)
        print('   chain  R_%d = %s' % (i, sp.factor(sp.cancel(Rk))))
