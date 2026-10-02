"""VIR8 (b), POST-HOC (labelled): after b_nz.py showed cB = cG = 0 at the nZ points for the three-term forms too,
the directional limits.  For each case: monic(cB) of the three-term ODE at k0 + h, h = 1e-8 and 1e-16 (exact rationals),
against the certified charge of the named solution at t0; and the size of the un-normalised leading coefficient.
PASS(case) = max|difference| is O(h) (ratio of the two values within a factor 2 of 1e8)."""
import hashlib, json, sys, time
import sympy as sp
sys.dont_write_bytecode = True
import wkb3 as W3
from b_core import table, to_PXPI
l0, l1 = W3.l0, W3.l1
from b_core import PX, PI  # same symbols as the helpers
_unused = sp.symbols('PX PI')
cases = [('b1', 3, sp.Integer(3), 7, 3), ('b2', 3, sp.Integer(1), 5, 3), ('b3', 2, sp.Integer(1), 5, 2), ('b4', 3, sp.Rational(1, 2), 9, 3),
         ('x1 (cross-check: Sol-2 ODE limit vs the SOL-3 table)', 2, sp.Integer(1), 5, 3), ('x2 (cross-check: Sol-3 ODE limit vs the SOL-2 table)', 3, sp.Integer(1), 5, 2)]
ok = True
t0 = time.time()
for name, sol, k0, spin, cert_sol in cases:
    tv0 = (3 * k0 + 2) / (k0 + 2)
    cert = table(cert_sol, spin, tv0)
    errs, leads = {}, {}
    for hexp in (8, 16):
        k = k0 + sp.Rational(1, 10**hexp)
        M = 1 / k
        tv = (3 * k + 2) / (k + 2)
        p2 = 2 * (tv - 1) / (tv + 1)
        if sol == 3:
            n = k + 4; c = n * (1 + M)
            eng = W3.WKB3(4, 1, M, [1, 0, -(l0**2 + l1**2), 0, l0**2 * l1**2], [1, -c / 2])
            S1 = n * n / p2; S0 = S1
        else:
            n = 2 * k + 3; c = n * (1 + M)
            eng = W3.WKB3(3, 2, M, [1, 0, -l1**2, 0], [1, -c, c * c / 4 - l0**2])
            S1 = n * n / p2; S0 = S1 / k
        i = spin + 1
        eng.solve(i)
        cB, cG = eng.J(i)
        assert cG.is_zero
        mB = to_PXPI(cB, S0, S1)
        mons = set(mB.monoms()) | set(cert.monoms())
        errs[hexp] = max(abs(mB.coeff_monomial(m) - cert.coeff_monomial(m)) for m in mons)
        leads[hexp] = abs(sp.Poly(cB, l0, l1).coeff_monomial(l0 ** i))
    r = errs[8] / errs[16] if errs[16] != 0 else sp.oo
    good = errs[8] < sp.Rational(1, 10**4) and sp.Rational(10**8, 2) < r < 2 * 10**8
    crossed = name.startswith('x')
    if not crossed:
        ok &= bool(good)
    print('%s: Sol-%d ODE, k0 = %s (t0 = %s), spin %d vs certified Sol-%d: max|monic - certified| = %.4e (h=1e-8), %.4e (h=1e-16); un-normalised lead %.3e, %.3e  -> %s  (%.0fs)'
          % (name, sol, k0, tv0, spin, cert_sol, float(errs[8]), float(errs[16]), float(leads[8]), float(leads[16]),
             'converges O(h)' if good else 'does NOT converge to this charge', time.time() - t0), flush=True)
print('directional limits', 'PASS' if ok else 'FAIL')
sys.exit(0 if ok else 2)
