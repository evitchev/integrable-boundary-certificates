"""VIR5a step 1: predictions of the continued operator at one fibre -> prediction_<tag>.json.
Opens no data file.  Usage: e1_predict.py <k as p/q> [tamper_gamma | tamper_M]"""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import mellin_wkb as MW
SEAL = open('SEAL_VIR5a.sha256').read().split()[0]
assert hashlib.sha256(open('SEAL_VIR5a.md', 'rb').read()).hexdigest() == SEAL
k = sp.Rational(sys.argv[1])
mode = sys.argv[2] if len(sys.argv) > 2 else 'real'
PX, PI = sp.symbols('PX PI')
t0 = time.time()
K = 14
n = k + 4
tv = (3 * k + 2) / (k + 2)
h = MW.family_symbol(k, K // 2 + 1)
Mv = 1 / k
if mode == 'tamper_gamma':
    T2 = sp.Symbol('T2')
    expr = sp.expand((1 - MW.l0**2 * T2) * (1 - MW.l1**2 * T2))
    h = [MW.P_(expr.coeff(T2, j)) for j in range(K // 2 + 2)]
if mode == 'tamper_M':
    Mv = 1 / k + sp.Rational(1, 10)
eng = MW.MellinWKB(n, Mv, h)
eng.solve(K)
c2 = n**2 * (k + 1) / k
cell = {'k': str(k), 't': str(tv), 'n': str(n), 'M': str(Mv), 'scale_c2': str(c2), 'mode': mode, 'seal': SEAL, 'spins': {}, 'odd_orders_vanish': {}}
for i in range(2, K + 1):
    Ri, ok = eng.R(i)
    assert ok
    if i % 2 == 1:
        cell['odd_orders_vanish'][str(i)] = bool(sp.simplify(Ri) == 0)
        continue
    Pl = sp.Poly(Ri, MW.l0, MW.l1)
    assert all(a % 2 == 0 and b % 2 == 0 for (a, b), _ in Pl.terms())
    Rm = sp.Poly(sp.expand(sum(cf * c2 ** ((a + b) // 2) * PX**a * PI**b for (a, b), cf in Pl.terms())), PX, PI)
    lead = Rm.coeff_monomial(PX**i)
    if lead == 0:
        cell['spins'][str(i - 1)] = {'leading_coefficient': '0', 'coefficients': None}
    else:
        Rn = sp.Poly(sp.expand(Rm.as_expr() / lead), PX, PI)
        cell['spins'][str(i - 1)] = {'leading_coefficient': str(lead),
                                     'coefficients': {'%d,%d' % (a, b): str(c) for (a, b), c in sorted(Rn.terms())}}
    print('k = %s spin %d: leading %s  (%.0fs)' % (k, i - 1, lead, time.time() - t0), flush=True)
tag = str(k).replace('/', '_').replace('-', 'm') + ('' if mode == 'real' else '_' + mode)
fn = 'prediction_k%s.json' % tag
json.dump(cell, open(fn, 'w'), indent=1, sort_keys=True)
print('written', fn, hashlib.sha256(open(fn, 'rb').read()).hexdigest(), 'odd orders vanish:', cell['odd_orders_vanish'])
