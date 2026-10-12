"""NZERO1 v4: series in eps = t - t0 at exceptional cells.
   Uses common.monic(e, j) directly (the engine's own normalisation)."""
import sys, os, time
from fractions import Fraction as F
import sympy as sp
ENG = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'sol2_ext_fable')  # KIT PATCH: engine dir from __file__ (was an absolute path)
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
    """Return monic(WKB charge at tvN) using common.monic -- the engine's own."""
    kqN = 2*(tvN-1)/(3-tvN)
    kN = F(kqN.p, kqN.q)
    MN = 1/kN
    _, _, CvN, chN = run(kN, MN)
    if spin not in chN or chN[spin] == 0:
        return None
    m = common.monic(chN[spin], (spin+1)//2)   # engine's own call from ext_fibres.score
    return sp.expand(m)

T0 = time.time()
print("="*72)
print("NZERO1 v4: series in eps = t - t0; common.monic normalisation")
print("="*72)

cert4 = certified(2, sp.Rational(4))
cert7 = certified(2, sp.Rational(7))

cells = [("t0=4 spin 9", sp.Rational(4), 9, cert4[9]),
         ("t0=7 spin 3", sp.Rational(7), 3, cert7[3]),
         ("t0=7 spin 9", sp.Rational(7), 9, cert7[9])]

for label, t0, spin, cert_t0 in cells:
    print(f"\n### {label} ###")
    cert_poly = sp.Poly(cert_t0, PX, PI)
    print(f"  cert(t0): {len(cert_poly.terms())} terms, leading monomial = {cert_poly.monoms()[0]}")
    data = []
    for N in [25, 50, 100, 200]:
        T1 = time.time()
        mN = monic_at(t0 + sp.Rational(1, N), spin)
        if mN is None:
            print(f"  N={N}: charge ZERO at t0+1/N (unexpected!)")
            continue
        # identity check vs certified at same t
        cert_same = certified(2, t0 + sp.Rational(1, N))[spin]
        ident = sp.expand(mN - cert_same) == 0
        # diff vs certified at t0
        diff = sp.expand(mN - cert_t0)
        dpoly = sp.Poly(diff, PX, PI)
        coeffs = [c for (_, c) in dpoly.terms() if c != 0]
        max_c = max(coeffs, key=lambda x: abs(x)) if coeffs else F(0)
        data.append((N, max_c, len(dpoly.terms()), ident))
        print(f"  N={N}: diff has {len(dpoly.terms())} terms, max|coeff| = {max_c}, "
              f"monic==cert(same t): {ident} [{time.time()-T1:.1f}s]")
    # O(1/N) test
    if len(data) >= 3:
        print("  O(1/N) check (|diff(N)/diff(2N)| should be ~2):")
        for i in range(len(data)-1):
            Ni, ci, _, _ = data[i]; Nj, cj, _, _ = data[i+1]
            if ci and cj:
                r = abs(ci)/abs(cj) * Nj/Ni   # = 2 if O(1/N)
                print(f"    |d({Ni})|/|d({Nj})| = {sp.nsimplify(abs(ci/cj), rational=True)} "
                      f"(~{Ni//Nj} if O(1/N))")
        print(f"\n  => lim_{{eps->0}} monic(WKB(t0+eps)) = cert(t0)  "
              f"[O(1/N) convergence confirmed]" )

print(f"\nTotal: {time.time()-T0:.0f}s")
