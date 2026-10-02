"""EXC2 U3, step 1 (ODE side only): I_5 (order 6) of the three-term operator with a general quartic and the first-order
1/xi term:  [T^4 + A2 T^2 + A3 T + A4 - x^(c/2) T x^(c/2) - E x^n + x^(-c)(B2 T^2 + B1 T + B0)] phi = 0.   Usage: u3_engine.py <t>"""
import json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import wkb3x as W
t = sp.Rational(sys.argv[1]); k = 2 * (t - 1) / (3 - t); M = 1 / k; n = k + 4; c = n * (1 + M)
t0 = time.time()
eng = W.WKB3(4, 1, M, [1, 0, W.A2, W.A3, W.A4], [1, -c / 2], qcoef=[W.B0, W.B1, W.B2])
eng.solve(6)
out = {'t': str(t)}
for i in (2, 4, 6):
    cB, cG = eng.J(i)
    assert cG.is_zero, 'second master at order %d' % i
    out['cB%d' % i] = str(cB.as_expr())
    print('order %d: %d terms  (%.0fs)' % (i, len(cB.terms()), time.time() - t0), flush=True)
print('cB6 =', out['cB6'])
json.dump(out, open('u3_engine_t%s.json' % str(t).replace('/', '_'), 'w'), indent=1)
