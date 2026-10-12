"""Debug: compare my charge computation vs ext_fibres.run() using direct function calls."""
import sys, os
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
os.chdir(os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
sys.path.insert(0, '.')
import common, gwkb, thwkb
from common import *

# Replicate ext_fibres.run() exactly (copy-paste from ext_fibres.py)
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

# Test at t=4 spin 9
print("=== ext_fibres.run() at t=4 (k=-6, M=-1/6) ===")
tv, g, Cv, ch = run(F(-6), F(-1, 6))
print(f"t={tv}, n={g['D']}, C={Cv}")
for s in [3, 5, 7, 9]:
    e = ch.get(s)
    if e is None: print(f"  spin {s}: MISSING"); continue
    P = sp.Poly(e, PX, PI)
    pxpows = sorted(m[0] for m in P.terms())
    mono = P.coeff_monomial(PX**(s-1))
    print(f"  spin {s}: {len(P.terms())} terms, PX powers={pxpows}, PX^{s-1} coeff={'ZERO' if mono==0 else str(sp.simplify(mono))[:40]}")

# Now check: what does score() do?
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

print("\n=== score() path at t=4 ===")
cert = certified(2, sp.Rational(4))
for s in [3, 5, 7, 9]:
    e = ch.get(s)
    if e is None: continue
    ge = gnorm(e)
    if ge == 0:
        print(f"  spin {s}: gnorm(ch) == 0 -> 'WKB charge identically 0 (certified: %s)' % {'ZERO' if cert[s]==0 else 'non-zero'}")
    else:
        d = sp.expand(gnorm(sp.expand(monic(e, (s+1)//2) - cert[s])))
        print(f"  spin {s}: monic diff = {'ZERO (MATCH)' if d==0 else f'{len(sp.Poly(d,PX,PI).terms())} terms DIFFER'}")

# The key question: at t=4 spin 9, the raw charge has PX powers [0] only.
# But the CERTIFIED charge has PX powers up to 8. So the WKB is giving a DIFFERENT
# object than the certified charge. The issue is that the WKB at this point gives
# the "subleading" part (the PI^2-only terms), not the full top-layer charge.
print("\n=== Key insight check ===")
e9 = ch[9]
P9 = sp.Poly(e9, PX, PI)
print(f"WKB spin 9 at t=4: terms = {list(P9.terms())}")
c9 = sp.expand(cert[9])
Pc9 = sp.Poly(c9, PX, PI)
print(f"Certified spin 9 at t=4 (first 5): {list(Pc9.terms())[:5]}")
print(f"\nThe WKB gives ONLY the PI^2 term (PX^0), while certified has full polynomial.")
print("This means the WKB 'charge' at the exceptional point is NOT the full charge,")
print("but only the subleading (in PX) component. The monic normalisation then fails")
print("because it expects the leading PX^(s-1) term which has vanished.")
