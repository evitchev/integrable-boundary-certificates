"""EXC2 T2, step 1: I_1 and I_3 of the three-term operator with a GENERAL quartic
    [ T^4 + A2 T^2 + A3 T + A4 - x^(c/2) T x^(c/2) - E x^n ] phi = 0     (n = k+4, M = 1/k)
from wkb3 (exact): cB_2, cB_4 as polynomials in (A2, A3, A4).  ODE side only.   Usage: g_charges.py <t>"""
import json, sys
import sympy as sp
sys.dont_write_bytecode = True
import wkb3g as W
t = sp.Rational(sys.argv[1]); k = 2 * (t - 1) / (3 - t); M = 1 / k; n = k + 4; c = n * (1 + M)
eng = W.WKB3(4, 1, M, [1, 0, W.A2, W.A3, W.A4], [1, -c / 2])
eng.solve(4)
out = {'t': str(t), 'k': str(k), 'n': str(n)}
for i in (2, 3, 4):
    cB, cG = eng.J(i)
    assert cG.is_zero
    out['cB%d' % i] = str(cB.as_expr())
    print('order %d: cB = %s' % (i, cB.as_expr()))
json.dump(out, open('g_charges_t%s.json' % str(t).replace('/', '_'), 'w'), indent=1)
