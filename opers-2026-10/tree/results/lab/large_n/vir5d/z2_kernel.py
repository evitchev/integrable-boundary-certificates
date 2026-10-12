"""VIR5d Z2, CFT side: the kernel of ad I_3 (Solution 3) modulo total derivatives at a drop point, and
the cylinder vacuum eigenvalue of every kernel class as a polynomial in (PX, pi).
Pipeline: VIR1 engine (kernel in Heis_X (x) Y-currents) -> VIR2 conversion to two bosons (requires the
class to lie in Heis (x) Vir_N + d; reported if not) -> VIR3 vacuum rule; rho reduced by rho^2 = 2/(t^2-1).
Usage: z2_kernel.py <t as p/q> <K>    -> z2_kernel_t<..>_K<..>.json"""
import hashlib, json, sys, time
from fractions import Fraction as F
import sympy as sp
sys.dont_write_bytecode = True
import vir_lib as V
import vir2_lib as W
import cft_lib as C
SEAL = open('SEAL_VIR5d.sha256').read().split()[0]
assert hashlib.sha256(open('SEAL_VIR5d.md', 'rb').read()).hexdigest() == SEAL
tq = F(sys.argv[1]); K = int(sys.argv[2]); w = 2 * K
t0 = time.time()
N, s = V.curve_point(tq)
P = V.P4(3, N, s)
com = V.Commutant(w)
ker = com.kernel(P, N)
print('t = %s (N = %s), weight %d: kernel dimension %d  (%.0fs)' % (tq, N, w, len(ker), time.time() - t0), flush=True)
rho2 = sp.Rational(2) / (sp.Rational(tq.numerator, tq.denominator) ** 2 - 1)
rho = W.rho
out = {'t': str(tq), 'K': K, 'kernel_dimension': len(ker), 'classes': []}
for idx, vec in enumerate(ker):
    try:
        img, hv = W.to_two_boson([sp.Rational(x.numerator, x.denominator) for x in vec], w)
    except AssertionError as e:
        print('  class %d: NOT in Heis (x) Vir_N + d (%s)' % (idx, e), flush=True)
        out['classes'].append({'in_HV': False})
        continue
    ev = sp.expand(C.vacuum_eigenvalue(img))
    Pr = sp.Poly(ev, C.PX, C.PI, rho)
    red = {}
    for (a, b, r), c in Pr.terms():
        key = (a, b, r % 2)
        red[key] = red.get(key, 0) + c * rho2 ** (r // 2)
    odd = {k_: v for k_, v in red.items() if k_[2] == 1 and v != 0}
    even = {(a, b): v for (a, b, r), v in red.items() if r == 0 and v != 0}
    par = all(a % 2 == 0 and b % 2 == 0 for a, b in even)
    print('  class %d: in HV + d; eigenvalue monomials %d; odd-in-rho part zero: %s; even in PX and pi: %s  (%.0fs)'
          % (idx, len(even), not odd, par, time.time() - t0), flush=True)
    out['classes'].append({'in_HV': True, 'odd_rho_zero': not odd, 'even': par,
                           'eigenvalue': {'%d,%d' % k_: str(v) for k_, v in sorted(even.items())}})
tag = str(tq).replace('/', '_')
json.dump(out, open('z2_kernel_t%s_K%d.json' % (tag, K), 'w'), indent=1, sort_keys=True)
print('done (%.0fs)' % (time.time() - t0))
