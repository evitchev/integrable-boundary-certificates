"""VIR5d Z1: monic operator polynomial at k = k0 + 10^-h, one spin (exact rationals; ODE side only).
Usage: z1_point.py <k0> <spin> <h>"""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import mellin_wkb as MW
SEAL = open('SEAL_VIR5d.sha256').read().split()[0]
assert hashlib.sha256(open('SEAL_VIR5d.md', 'rb').read()).hexdigest() == SEAL
k0 = sp.Rational(sys.argv[1]); spin = int(sys.argv[2]); hexp = int(sys.argv[3])
k = k0 + sp.Rational(1, 10**hexp)
PX, PI = sp.symbols('PX PI')
t0 = time.time()
polys, odd = MW.family_monic(k, spin + 1, PX, PI)
P, lead = polys[spin]
assert P is not None
out = {'k0': str(k0), 'h': '1e-%d' % hexp, 'spin': spin, 'leading_coefficient': str(lead),
       'coefficients': {'%d,%d' % (a, b): str(c) for (a, b), c in sorted(P.terms())}}
tag = str(k0).replace('/', '_')
json.dump(out, open('z1_k%s_s%d_h%d.json' % (tag, spin, hexp), 'w'), indent=1, sort_keys=True)
print('k0 = %s + 1e-%d spin %d: |leading coefficient| = %.3e  (%.0fs)' % (k0, hexp, spin, abs(float(lead)), time.time() - t0))
