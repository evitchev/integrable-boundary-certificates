"""Debug: replicate ext_fibres.score() EXACTLY and find where my comparison diverges."""
import sys, os
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
os.chdir(os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
sys.path.insert(0, '.')
import common, gwkb, thwkb
from common import *

# EXACT copy of ext_fibres functions
def run(k, M, K=10, delta_shift=0):
    kq = sp.Rational(k.numerator, k.denominator)
    tv = (3 * kq + 2) / (kq + 2); g = classG(2, tv)
    assert g['alpha'] == F(k) and g['gamma'] == 1 and g['D'] == 2 * F(k) + 3, g
    n = g['D']; Cv = 2 * n * (1 + F(M))
    bl = blocks_of(g, M)
    if delta_shift: bl = [(a, sh, d + delta_shift) if a != 1 else (a, sh, d) for (a, sh, d) in bl]
    sv = gwkb.wkb(bl, g['s'], F(M), K); bb = n * F(M); ch = {}
    for kk in range(2, K + 1, 2):
        sk = sv[kk]
        if any(a.denominator == 1 and a <= 0 for (a, c) in sk): pw, poly, nonpole = thwkb.integrate_reg(sk, bb); poly = dict(poly)
        else:
            r = thwkb.integrate(sk, bb); poly = {k_: sp.Rational(v.numerator, v.denominator) for k_, v in r[2].items()} if r else {}
        ch[kk - 1] = sp.expand(sum(v * (sp.Rational(Cv.numerator, Cv.denominator) * PX**2 / g['aY'])**(i // 2) * (sp.Rational(Cv.numerator, Cv.denominator) * PI**2 / g['aX'])**(j // 2) for (i, j), v in poly.items()))
    return tv, g, Cv, ch

def gnorm(expr):
    rep = {}
    for g_ in expr.atoms(sp.gamma):
        q = g_.args[0]
        if not q.is_Rational or q.is_Integer: continue
        fac = sp.Integer(1); x = q
        while x > 1: x -= 1; fac *= x
        while x < 0: fac /= x; x += 1
        rep[g_] = fac * sp.gamma(x)
    return sp.cancel(sp.together(expr.xreplace(rep)))

# Test at t=9/4 (known good)
print("=== t=9/4, spin 3: step-by-step ===")
tv = sp.Rational('9/4'); kq = 2*(tv-1)/(3-tv); k = F(kq.p, kq.q)
print(f"k = {k}")
_, g, Cv, ch = run(k, 1/k)
cert = certified(2, tv)
s = 3
ch_s = ch[s]
print(f"WKB charge: {sp.expand(ch_s)}")
print(f"Certified:  {sp.expand(cert[s])}")
# EXACT score() logic:
gn_ch = gnorm(ch_s)
print(f"\ngnorm(WKB): {gn_ch}")
if gn_ch == 0:
    print("  -> WKB charge identically 0")
else:
    d = sp.expand(gnorm(sp.expand(monic(ch_s, (s+1)//2) - cert[s])))
    print(f"score diff: {d}")
    print(f"  {'MATCH' if d==0 else 'DIFFER'}")
    
# My version:
my_monic = sp.expand(ch_s / sp.Poly(ch_s, PX, PI).coeff_monomial(PX**2))
my_diff = sp.expand(gnorm(my_monic - cert[s]))
print(f"\nMy monic: {my_monic}")
print(f"My diff: {'ZERO' if my_diff==0 else f'{len(sp.Poly(my_diff,PX,PI).terms())} terms'}")

# Check: is monic() the same as dividing by coeff_monomial?
print(f"\ncommon.monic(ch): {monic(ch_s, 2)}")
print(f"ch/coeff:          {sp.expand(ch_s / sp.Poly(ch_s,PX,PI).coeff_monomial(PX**2))}")
print(f"Same: {sp.simplify(monic(ch_s,2) - sp.expand(ch_s / sp.Poly(ch_s,PX,PI).coeff_monomial(PX**2))) == 0}")

# The issue might be that certified() uses a different normalisation
print(f"\ncert[3] leading PX^2 coeff: {sp.Poly(cert[3],PX,PI).coeff_monomial(PX**2)}")
print(f"WKB  leading PX^2 coeff:     {sp.Poly(ch_s,PX,PI).coeff_monomial(PX**2)}")
