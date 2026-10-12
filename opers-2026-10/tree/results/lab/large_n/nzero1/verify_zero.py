"""Verify: does the sum of the 9 surviving terms cancel exactly?
And compute the proper limit via symbolic expansion in eps = t - t0."""
import sys, os
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
os.chdir(os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
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
        if any(a.denominator == 1 and a <= 0 for (a, c) in sk): pw, poly, nonpole = thwkb.integrate_reg(sk, bb); poly = dict(poly)
        else:
            r = thwkb.integrate(sk, bb); poly = {k_: sp.Rational(v.numerator, v.denominator) for k_, v in r[2].items()} if r else {}
        ch[kk - 1] = sp.expand(sum(v * (sp.Rational(Cv.numerator, Cv.denominator) * PX**2 / g['aY'])**(i // 2) * (sp.Rational(Cv.numerator, Cv.denominator) * PI**2 / g['aX'])**(j // 2) for (i, j), v in poly.items()))
    return tv, g, Cv, ch, ch.get(9)

# At t=4, look at the individual terms before summing
print("=== t=4 spin 9: individual term values ===")
tv = sp.Rational(4); k = F(-6); Mv = F(-1, 6)
g = classG(2, tv)
Cv = 2 * g['D'] * (1 + Mv)
bl = blocks_of(g, Mv)
sv = gwkb.wkb(bl, g['s'], Mv, 10)
sk = sv[10]  # spin 9
bb = g['D'] * Mv
pw, poly, nonpole = thwkb.integrate_reg(sk, bb)
poly = dict(poly)
print(f"integrate_reg returned {len(poly)} surviving terms:")
total = sp.Integer(0)
for (i, j), v in sorted(poly.items()):
    val = sp.Rational(v.numerator, v.denominator) * (sp.Rational(Cv.numerator, Cv.denominator) * PX**2 / g['aY'])**(i // 2) * (sp.Rational(Cv.numerator, Cv.denominator) * PI**2 / g['aX'])**(j // 2)
    total += val
    print(f"  (i,j)=({i},{j}): coeff={v}, value = {sp.simplify(val)}")
print(f"\nSUM = {sp.simplify(total)}")
print(f"Is zero: {sp.simplify(total) == 0}")

# Now the LIMIT: compute at t = 4 + eps with eps symbolic
print("\n\n=== LIMIT via symbolic eps ===")
# k(t) = 2(t-1)/(3-t); at t=4+eps: k = 2(3+eps)/(-1-eps) = -2(3+eps)/(1+eps)
eps = sp.Symbol('eps')
t_eps = 4 + eps
k_eps = 2*(t_eps-1)/(3-t_eps)
# M = 1/k = (3-t)/2(t-1) = (-1-eps)/(6+2*eps)
M_eps = 1/k_eps
# n = 2k+3
n_eps = 2*k_eps + 3
print(f"k(eps) = {sp.simplify(k_eps)}")
print(f"M(eps) = {sp.simplify(M_eps)}")
print(f"n(eps) = {sp.simplify(n_eps)}")
# aX = -2(t-3)/(t+5) = -2(1+eps)/(9+eps)
aX_eps = -2*(t_eps-3)/(t_eps+5)
# aY = 4(t-1)/(t+5) = 4(3+eps)/(9+eps)
aY_eps = 4*(t_eps-1)/(t_eps+5)
print(f"aX(eps) = {sp.simplify(aX_eps)}")
print(f"aY(eps) = {sp.simplify(aY_eps)}")

# The key: at generic t, the charge has a PX^(spin-1) term. At t=4 it vanishes.
# The LIMIT of the monic normalised charge should be finite.
# Strategy: compute the charge at t = 4 + 1/N for large N (already done), 
# and check if it converges to the certified value.
# From the earlier run: at t=801/200 (= 4.005), the monic had 21 terms but MISMATCHED.
# This suggests the approach is NOT simply continuous -- there's a branch issue.

# Let me check: what does the certified charge look like as t -> 4?
# cert[9] at t=4 has 21 terms. Does cert[9](t) -> cert[9](4) as t->4?
cert4 = certified(2, sp.Rational(4))
c9_4 = sp.expand(cert4[9])
print(f"\nCertified spin 9 at t=4: {len(sp.Poly(c9_4, PX, PI).terms())} terms")
print(f"Leading: {sp.Poly(c9_4, PX, PI).coeff_monomial(PX**8)} * PX^8 * ...")

# The WKB at t=4+eps gives a charge whose monic tends to something.
# Let me compute at several eps values and check convergence of the monic.
print("\n=== Convergence check: monic at t=4+1/N for N=10..50 ===")
for N in [10, 20, 50, 100]:
    tvN = sp.Rational(4) + sp.Rational(1, N)
    kqN = 2*(tvN-1)/(3-tvN)
    kN = F(kqN.p, kqN.q)
    try:
        _, _, CvN, chN = run(kN, 1/kN)
        P = sp.Poly(chN[9], PX, PI)
        mono = P.coeff_monomial(PX**8)
        if mono != 0:
            mN = sp.expand(chN[9] / mono)
            # Compare with certified at same t
            certN = certified(2, tvN)
            dN = sp.expand(sp.simplify(mN - certN[9]))
            print(f"  N={N}: t={tvN}, monic has {len(sp.Poly(mN,PX,PI).terms())} terms, diff vs cert = {'ZERO' if dN==0 else f'{len(sp.Poly(dN,PX,PI).terms())} terms'}")
        else:
            print(f"  N={N}: monic still zero")
    except Exception as ex:
        print(f"  N={N}: ERROR {ex}")
