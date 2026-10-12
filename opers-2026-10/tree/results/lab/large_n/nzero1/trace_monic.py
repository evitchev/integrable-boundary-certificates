"""Check whether the monic coefficient (PX^(spin-1)) vanishes in the raw charge,
causing the gnorm(ch[s]) == 0 path in score()."""
import sys, os
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
os.chdir(os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
sys.path.insert(0, '.')
import common, gwkb, thwkb
from common import *

def check_monic(tval, K_val=10):
    tv = sp.Rational(tval)
    kq = 2 * (tv - 1) / (3 - tv)
    k = F(kq.p, kq.q)
    g = classG(2, tv)
    n = g['D']
    Mv = 1 / k
    bb = n * Mv
    Cv = 2 * n * (1 + Mv)
    bl = blocks_of(g, Mv)
    sv = gwkb.wkb(bl, g['s'], Mv, K_val)
    print(f"\n=== t={tval}, k={k}, n={n} ===")
    for kk in range(2, K_val + 1, 2):
        spin = kk - 1
        sk = sv[kk]
        has_pole = any(a.denominator == 1 and a <= 0 for (a, c) in sk)
        if has_pole:
            pw, poly, nonpole = thwkb.integrate_reg(sk, bb)
            poly = dict(poly)
        else:
            r = thwkb.integrate(sk, bb)
            poly = {k_: sp.Rational(v.numerator, v.denominator) for k_, v in r[2].items()} if r else {}
        ch_raw = sp.expand(sum(v * (sp.Rational(Cv.numerator, Cv.denominator) * PX**2 / g['aY'])**(i // 2) *
                                (sp.Rational(Cv.numerator, Cv.denominator) * PI**2 / g['aX'])**(j // 2)
                            for (i, j), v in poly.items()))
        # Check the PX^(spin-1) coefficient
        P = sp.Poly(ch_raw, PX, PI)
        mono_coeff = P.coeff_monomial(PX**(spin-1))
        # What monomials are present?
        exponents = set()
        for (mon, coeff) in P.terms():
            exponents.add(mon[0] // 2)  # power of PX
        print(f"  spin {spin}: raw charge has PX powers {sorted(exponents)}")
        print(f"    PX^{spin-1} coeff = {sp.simplify(mono_coeff)}")
        if mono_coeff == 0:
            print(f"    *** MONIC COEFFICIENT IS ZERO -> score() reports 'WKB charge identically 0' ***")
        else:
            # Try the full monic
            ch_monic = sp.expand(ch_raw / mono_coeff)
            print(f"    monic OK, leading terms: {len(sp.Poly(ch_monic, PX, PI).terms())} terms")

check_monic('4')
check_monic('7')
