"""VIR5d Z0(c),(d): does the operator's un-normalised WKB coefficient R_(2K) vanish identically (all layers)?
Usage: z0c_engine.py <k as p/q> <engine order> [stripped]     (ODE side only; opens no data file)"""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import mellin_wkb as MW
SEAL = open('SEAL_VIR5d.sha256').read().split()[0]
assert hashlib.sha256(open('SEAL_VIR5d.md', 'rb').read()).hexdigest() == SEAL
assert hashlib.sha256(open('mellin_wkb.py', 'rb').read()).hexdigest().startswith('60795426')
k = sp.Rational(sys.argv[1]); order = int(sys.argv[2]); mode = sys.argv[3] if len(sys.argv) > 3 else 'real'
n = k + 4
tv = (3 * k + 2) / (k + 2)
h = MW.family_symbol(k, order // 2 + 1)
if mode == 'stripped':
    T2 = sp.Symbol('T2')
    expr = sp.expand((1 - MW.l0**2 * T2) * (1 - MW.l1**2 * T2))
    h = [MW.P_(expr.coeff(T2, j)) for j in range(order // 2 + 2)]
eng = MW.MellinWKB(n, 1 / k, h)
t0 = time.time()
eng.solve(order)
out = {'k': str(k), 't': str(tv), 'n': str(n), 'mode': mode, 'rows': {}}
for i in range(2, order + 1, 2):
    Ri, ok = eng.R(i)
    assert ok
    Pl = sp.Poly(Ri, MW.l0, MW.l1)
    K = i // 2
    nu = sp.Rational(2 * K - 1) / n
    terms = {m: c for m, c in Pl.terms() if c != 0}
    by_loss = {}
    for (a, b), c in terms.items():
        by_loss.setdefault(K - (a + b) // 2, []).append((a // 2, b // 2))
    row = {'K': K, 'spin': 2 * K - 1, 'nu': str(nu), 'nu_integer': bool(nu.is_integer), 'identically_zero': len(terms) == 0,
           'nonzero_monomials_by_loss': {str(l_): sorted(v) for l_, v in sorted(by_loss.items())}}
    out['rows'][str(K)] = row
    print('k = %s (t = %s, n = %s) %s  K = %d spin %2d  nu = %-6s  R identically zero: %s   top monomials (X^a Y^b): %s'
          % (k, tv, n, mode, K, 2 * K - 1, nu, row['identically_zero'], by_loss.get(0, [])), flush=True)
tag = str(k).replace('/', '_').replace('-', 'm') + ('' if mode == 'real' else '_' + mode)
json.dump(out, open('z0c_k%s.json' % tag, 'w'), indent=1, sort_keys=True)
print('done in %.0fs' % (time.time() - t0))
