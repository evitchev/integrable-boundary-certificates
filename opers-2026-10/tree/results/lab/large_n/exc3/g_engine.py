"""EXC3: WKB charges (orders 2, 4, 6 and the odd 3, 5) of the deformed Sol 2 three-term operator, symbolic in the deformation:
  P_T = T^3 + A2 T^2 + A3 T + A4,  Lc = (T - c/2)^2 - l0^2,  x^(-c)(B1 T + B0),  x^(-2c)(C1 T + C0),  x^(-3c)(Q31 T + Q30),  x^(-4c)(Q41 T + Q40)   (no vartheta^2 deformation).
v2: four orders in x^(-c) (v1 had two; its outputs g_engine_t*_v1.json are kept).
Fibre t -> g_engine_t<t>.json.  Usage: g_engine.py <t>"""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_EXC3.md', 'rb').read()).hexdigest() == open('SEAL_EXC3.sha256').read().split()[0]
import wkb3x2 as W
t = sp.Rational(sys.argv[1]); k = 2 * (t - 1) / (3 - t); M = 1 / k; n = 2 * k + 3; c = n * (1 + M)
t0 = time.time()
eng = W.WKB3(3, 2, M, [1, W.A2, W.A3, W.A4], [1, -c, c**2 / 4 - W.l0**2], qcoef=[W.B0, W.B1], q2coef=[W.C0, W.C1], q3coef=[W.Q30, W.Q31], q4coef=[W.Q40, W.Q41])
eng.solve(6)
out = {'t': str(t), 'gens': 'l0 (Lc), A2 A3 A4 (P_T), B0 B1 (x^-c), C0 C1 (x^-2c), Q30 Q31 (x^-3c), Q40 Q41 (x^-4c)', 'engine': 'wkb3x2 v2 MAXU=4'}
for i in (2, 3, 4, 5, 6):
    cB, cG = eng.J(i)
    assert cG.is_zero, 'second master at order %d' % i
    out['cB%d' % i] = str(cB.as_expr())
    print('order %d: %d terms  (%.0fs)' % (i, len(cB.terms()), time.time() - t0), flush=True)
json.dump(out, open('g_engine_t%s.json' % str(t).replace('/', '_'), 'w'), indent=1)
print('cB2 =', out['cB2']); print('cB3 =', out['cB3']); print('cB4 =', out['cB4'])
