import sys, os
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
os.chdir(os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
sys.path.insert(0, '.')
import common, gwkb, thwkb
from common import *

# Replicate compute_limit_v3.py's run() EXACTLY
def run(k, M, K=10):
    kq = sp.Rational(k.numerator, k.denominator)
    tv = (3 * kq + 2) / (kq + 2); g = classG(2, tv)
    n = g['D']; Cv = 2 * n * (1 + F(M))
    bl = blocks_of(g, M)
    sv = gwkb.wkb(bl, g['s'], F(M), K); bb = n * F(M); ch = {}
    for kk in range(2, K + 1, 2):
        sk = sv[kk]
        if any(a.denominator == 1 and a <= 0 for (a, c) in sk): pw, poly, nonpole = thwkb.integrate_reg(sk, bb); poly = dict(poly)
        else:
            r = thwkb.integrate(sk, bb); poly = {k_: sp.Rational(v.numerator, v.denominator) for k_, v in r[2].items()} if r else {}
        ch[kk - 1] = sp.expand(sum(v * (sp.Rational(Cv.numerator, Cv.denominator) * PX**2 / g['aY'])**(i // 2) * (sp.Rational(Cv.numerator, Cv.denominator) * PI**2 / g['aX'])**(j // 2) for (i, j), v in poly.items()))
    return tv, g, Cv, ch

# At t=401/100, spin 9
tv = sp.Rational('401/100')
kq = 2*(tv-1)/(3-tv)
print(f"t={tv}, kq={kq}, k={F(kq.p, kq.q)}")
# compute_limit_v3 used: run(kN, MN) where kN=F(kqN.p,kqN.q), MN=1/kN
kN = F(kq.p, kq.q)
MN = 1/kN
_, _, CvN, chN = run(kN, MN)
P = sp.Poly(chN[9], PX, PI)
mono = P.coeff_monomial(PX**5)  # (spin+1)//2 = 5
print(f"ch[9] has {len(P.terms())} terms, PX powers = {sorted(m[0] for m in P.terms())}")
print(f"monic coeff (PX^5) = {mono}")
if mono != 0:
    m = sp.expand(chN[9]/mono)
    print(f"monic OK: {len(sp.Poly(m,PX,PI).terms())} terms")
else:
    print("MONIC ZERO - bug!")
    # What power IS the leading term?
    maxpx = max(m[0] for m in P.terms())
    print(f"Max PX power is {maxpx}, not 5")
    mono2 = P.coeff_monomial(PX**maxpx)
    print(f"Leading coeff (PX^{maxpx}) = {mono2}")

# Now try my exact_series_v2 approach
def run_at_tv(tv2, K=10):
    kq2 = 2*(tv2-1)/(3-tv2)
    k2 = F(kq2.p, kq2.q)
    g2 = classG(2, tv2)
    n2 = g2['D']
    Mv2 = 1/k2
    Cv2 = 2*n2*(1+Mv2)
    bl2 = blocks_of(g2, Mv2)
    sv2 = gwkb.wkb(bl2, g2['s'], Mv2, K)
    bb2 = n2*Mv2
    ch2 = {}
    for kk in range(2, K+1, 2):
        sk2 = sv2[kk]
        if any(a.denominator == 1 and a <= 0 for (a,c) in sk2):
            pw2, poly2, nonpole2 = thwkb.integrate_reg(sk2, bb2); poly2 = dict(poly2)
        else:
            r2 = thwkb.integrate(sk2, bb2)
            poly2 = {k_: sp.Rational(v.numerator, v.denominator) for k_, v in r2[2].items()} if r2 else {}
        ch2[kk-1] = sp.expand(sum(v * (sp.Rational(Cv2.numerator, Cv2.denominator)*PX**2/g2['aY'])**(i//2) *
                                   (sp.Rational(Cv2.numerator, Cv2.denominator)*PI**2/g2['aX'])**(j//2)
                               for (i,j), v in poly2.items()))
    return g2, Cv2, ch2

g3, Cv3, ch3 = run_at_tv(tv, 10)
P3 = sp.Poly(ch3[9], PX, PI)
print(f"\nrun_at_tv: ch[9] has {len(P3.terms())} terms, PX powers = {sorted(m[0] for m in P3.terms())}")
mono3 = P3.coeff_monomial(PX**5)
print(f"monic coeff (PX^5) = {mono3}")
# Compare the two charges
d = sp.expand(chN[9] - ch3[9])
print(f"ch from run() vs run_at_tv(): diff = {'ZERO' if d==0 else 'DIFFER'}")
