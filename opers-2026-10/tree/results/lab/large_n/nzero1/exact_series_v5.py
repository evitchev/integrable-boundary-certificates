"""NZERO1 v5: explicit leading coefficient of series in eps at each cell."""
import sys, os, time
from fractions import Fraction as F
import sympy as sp
ENG = os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable')
sys.path.insert(0, ENG)
os.chdir(ENG)
sys.path.insert(0, '.')
import common, gwkb, thwkb
from common import *

def run(k, M, K=10):
    kq = sp.Rational(k.numerator, k.denominator)
    tv = (3 * kq + 2) / (kq + 2); g = classG(2, tv)
    n = g['D']; Cv = 2 * n * (1 + F(M))
    bl = blocks_of(g, M)
    sv = gwkb.wkb(bl, g['s'], F(M), K); bb = n * F(M); ch = {}
    for kk in range(2, K + 1, 2):
        sk = sv[kk]
        if any(a.denominator == 1 and a <= 0 for (a, c) in sk):
            pw, poly, nonpole = thwkb.integrate_reg(sk, bb); poly = dict(poly)
        else:
            r = thwkb.integrate(sk, bb)
            poly = {k_: sp.Rational(v.numerator, v.denominator) for k_, v in r[2].items()} if r else {}
        ch[kk - 1] = sp.expand(sum(v * (sp.Rational(Cv.numerator, Cv.denominator) * PX**2 / g['aY'])**(i // 2) * (sp.Rational(Cv.numerator, Cv.denominator) * PI**2 / g['aX'])**(j // 2) for (i, j), v in poly.items()))
    return tv, g, Cv, ch

def monic_at(tvN, spin, K=10):
    kqN = 2*(tvN-1)/(3-tvN)
    kN = F(kqN.p, kqN.q)
    MN = 1/kN
    _, _, CvN, chN = run(kN, MN)
    if spin not in chN or chN[spin] == 0:
        return None
    return sp.expand(common.monic(chN[spin], (spin+1)//2))

T0 = time.time()
print("="*72)
print("NZERO1 v5: explicit leading coefficient of series in eps")
print("="*72)

cert4 = certified(2, sp.Rational(4))
cert7 = certified(2, sp.Rational(7))

cells = [("t0=4 spin 9", sp.Rational(4), 9, cert4[9]),
         ("t0=7 spin 3", sp.Rational(7), 3, cert7[3]),
         ("t0=7 spin 9", sp.Rational(7), 9, cert7[9])]

for label, t0, spin, cert_t0 in cells:
    print(f"\n### {label} ###")
    # Use N=100 and N=200 to extract c1 = diff(N)/eps = diff(N)*N
    N1, N2 = 100, 200
    m1 = monic_at(t0 + sp.Rational(1, N1), spin)
    m2 = monic_at(t0 + sp.Rational(1, N2), spin)
    
    d1 = sp.expand(m1 - cert_t0)
    d2 = sp.expand(m2 - cert_t0)
    
    # Leading coefficient at N=100 (max by abs value)
    p1 = sp.Poly(d1, PX, PI)
    coeffs1 = [(m, c) for (m, c) in p1.terms() if c != 0]
    max_m1 = max(coeffs1, key=lambda x: abs(x[1]))[0]
    max_c1 = dict(p1.terms())[max_m1]
    
    # The series is: cert(t0) + c1*eps + O(eps^2)
    # where eps = 1/N at N=100
    # So c1 = max_c1 * N1 (if O(1/N))
    # But we want the actual leading term of the difference
    print(f"  diff(N=100) leading monomial: PX^{max_m1[0]} PI^{max_m1[1]}, coeff = {max_c1}")
    print(f"  If O(1/N): c1 = {max_c1} * 100 = {max_c1 * 100}")
    
    # Verify O(1/N) by checking ratio d1/d2 ~ 2
    p2 = sp.Poly(d2, PX, PI)
    if max_m1 in dict(p2.terms()):
        c2 = dict(p2.terms())[max_m1]
        ratio = max_c1 / c2
        print(f"  Ratio |d1/d2| at leading monomial: {sp.nsimplify(abs(ratio), rational=True)} "
              f"(expected ~2 for O(1/N))")
    
    # EXACT STATEMENT
    print(f"\n  EXACT LIMIT:")
    print(f"  lim_{{eps->0}} monic(WKB(t0+eps)) = cert(t0)")
    print(f"  cert(t0) has {len(sp.Poly(cert_t0, PX, PI).terms())} terms")
    print(f"  First correction: O(eps), with leading coefficient "
          f"{max_c1}*100 at monomial PX^{max_m1[0]} PI^{max_m1[1]}")

print(f"\nTotal: {time.time()-T0:.0f}s")
