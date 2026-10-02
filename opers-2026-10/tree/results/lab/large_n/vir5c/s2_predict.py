"""VIR5c S2: engine predictions at one fibre up to a given engine order -> pred_<tag>.json.
Opens no data file.  Usage: s2_predict.py <k as p/q> <order> [stripped]"""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import mellin_wkb as MW
SEAL = open('SEAL_VIR5c.sha256').read().split()[0]
assert hashlib.sha256(open('SEAL_VIR5c.md', 'rb').read()).hexdigest() == SEAL
assert hashlib.sha256(open('mellin_wkb.py', 'rb').read()).hexdigest().startswith('60795426')
k = sp.Rational(sys.argv[1])
K = int(sys.argv[2])
mode = sys.argv[3] if len(sys.argv) > 3 else 'real'
PX, PI = sp.symbols('PX PI')
t0 = time.time()
n = k + 4
tv = (3 * k + 2) / (k + 2)
h = MW.family_symbol(k, K // 2 + 1)
if mode == 'stripped':
    T2 = sp.Symbol('T2')
    expr = sp.expand((1 - MW.l0**2 * T2) * (1 - MW.l1**2 * T2))
    h = [MW.P_(expr.coeff(T2, j)) for j in range(K // 2 + 2)]
eng = MW.MellinWKB(n, 1 / k, h)
c2 = n**2 * (k + 1) / k
cell = {'k': str(k), 't': str(tv), 'n': str(n), 'mode': mode, 'order': K, 'seal': SEAL, 'spins': {}, 'odd_orders_vanish': {}}
tag = str(k).replace('/', '_').replace('-', 'm') + ('' if mode == 'real' else '_' + mode)
fn = 'pred_k%s_o%d.json' % (tag, K)
eng.W = {}
for N in range(1, K + 1):
    rest = eng.stage(N)
    eng.W[N] = MW.el_scale(rest, sp.Rational(-1) / MW.R_(eng.n))
    i = N
    if i < 2:
        continue
    Ri, ok = eng.R(i)
    assert ok
    if i % 2 == 1:
        cell['odd_orders_vanish'][str(i)] = bool(sp.expand(Ri) == 0)
        print('k = %s order %d (odd): vanishes %s  (%.0fs)' % (k, i, cell['odd_orders_vanish'][str(i)], time.time() - t0), flush=True)
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
    print('k = %s spin %d done  (%.0fs)' % (k, i - 1, time.time() - t0), flush=True)
    json.dump(cell, open(fn, 'w'), indent=1, sort_keys=True)      # checkpoint after every spin
print('written', fn, hashlib.sha256(open(fn, 'rb').read()).hexdigest())
