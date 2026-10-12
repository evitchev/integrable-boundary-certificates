"""NZERO1: trace where the zero comes from. Exact cells per brief:
  t=4: k=-6, n=-9, M=1/k=-1/6; spin 9 (nu=spin/n=-9/-9=1, k*nu=-6 -> nu=-1... wait)
  Actually: nu = spin/n. At t=4, n=-9: spin 9 -> nu = 9/(-9) = -1. k*nu = (-6)*(-1) = 6 (positive integer). CHECK.
  At t=7, n=-3: spin 3 -> nu = 3/(-3) = -1, k*nu = (-3/5... wait k at t=7)
  k at t=7: 2(7-1)/(3-7) = 12/(-4) = -3. n = 2k+3 = -3. spin 3: nu=-1, k*nu=3. spin 9: nu=-3, k*nu=9.
  
  So B(-k*nu, nu*(1+k)) = B(k*nu, nu + nu*k):
  t=4 sp9: B(6, -1+(-6)) = B(6,-7)?? That's Gamma(6)*Gamma(-7)/Gamma(-1) -> both poles.
  Wait, the brief says B(-k*nu, nu*(1+k)). Let me re-read.
  nu = spin/n. At t=4: n=-9, spin=9, nu=9/(-9)=-1.
  -k*nu = -(-6)(-1) = -6. nu*(1+k) = -1*(1+(-6)) = -1*(-5) = 5.
  B(-6, 5) = Gamma(-6)*Gamma(5)/Gamma(-1). Both Gamma(-6) and Gamma(-1) have poles -> infinity/infinity."""
import sys, os
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
os.chdir(os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
sys.path.insert(0, '.')
import common, gwkb, thwkb
from common import *

def trace_cell(tval, K_val=10):
    tv = sp.Rational(tval)
    kq = 2 * (tv - 1) / (3 - tv)
    k = F(kq.p, kq.q)
    g = classG(2, tv)
    n = g['D']
    Mv = 1 / k
    bb = n * Mv
    print(f"\n=== t={tval}, k={k}, n={n}, M={Mv}, aX={g['aX']}, aY={g['aY']} ===")
    bl = blocks_of(g, Mv)
    sv = gwkb.wkb(bl, g['s'], Mv, K_val)
    results = {}
    for kk in range(2, K_val + 1, 2):
        spin = kk - 1
        sk = sv[kk]
        has_pole = any(a.denominator == 1 and a <= 0 for (a, c) in sk)
        if not has_pole:
            r = thwkb.integrate(sk, bb)
            poly = {k_: sp.Rational(v.numerator, v.denominator) for k_, v in r[2].items()} if r else {}
            Cv = 2 * n * (1 + Mv)
            ch = sp.expand(sum(v * (sp.Rational(Cv.numerator, Cv.denominator) * PX**2 / g['aY'])**(i // 2) *
                                (sp.Rational(Cv.numerator, Cv.denominator) * PI**2 / g['aX'])**(j // 2)
                            for (i, j), v in poly.items()))
            is_zero = (ch == 0)
            print(f"  spin {spin}: REGULAR path, charge {'ZERO' if is_zero else 'non-zero'}")
            results[spin] = ('regular', is_zero)
            continue
        pw, poly, nonpole = thwkb.integrate_reg(sk, bb)
        poly = dict(poly)
        Cv = 2 * n * (1 + Mv)
        ch = sp.expand(sum(v * (sp.Rational(Cv.numerator, Cv.denominator) * PX**2 / g['aY'])**(i // 2) *
                            (sp.Rational(Cv.numerator, Cv.denominator) * PI**2 / g['aX'])**(j // 2)
                        for (i, j), v in poly.items()))
        is_zero = (ch == 0)
        nu = F(spin) / F(int(n))
        knu = k * nu
        n1k = nu * (1 + k)
        print(f"  spin {spin}: POLE path, charge {'*** ZERO ***' if is_zero else 'non-zero'}")
        print(f"    nu={nu}, -k*nu={-knu}, nu*(1+k)={n1k}")
        print(f"    B(-k*nu, nu*(1+k)) = Gamma({-knu})*Gamma({n1k})/Gamma({nu}) -> " + 
              ("POLE/POLE" if (-knu).is_Integer and -knu <= 0 and nu.is_Integer and nu <= 0 else "evaluable"))
        # Show which terms die
        for (a, c), p in sorted(sk.items()):
            a_r = sp.Rational(a.numerator, a.denominator)
            c_r = sp.Rational(c.numerator, c.denominator)
            if a_r.is_Integer and a_r <= 0:
                m = int(-a_r)
                neg_c = -c_r
                killed = (neg_c.is_Integer and neg_c <= 0)
                print(f"      v^{a_r}(v+E')^{c_r}: res(Gamma)={F((-1)**m, sp.factorial(m))}, Gamma(-c)=Gamma({neg_c}) {'KILLED' if killed else '= '+str(sp.gamma(neg_c))}")
        results[spin] = ('pole', is_zero)
    return results

# The three exceptional cells
print("Exceptional cells (operator coefficient = 0, certified charge non-zero):")
trace_cell('4')
trace_cell('7')
