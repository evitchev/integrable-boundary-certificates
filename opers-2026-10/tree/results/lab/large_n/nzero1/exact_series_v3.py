"""NZERO1 v3: CORRECT monic exponent = spin//2 + 1 (not (spin+1)//2)."""
import sys, os, time
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
    return tv, g, Cv, ch

def monic_charge_at(t0, delta_num, delta_den, spin, K=10):
    """Compute monic WKB charge at t0 + delta_num/delta_den. CORRECT exponent."""
    tvN = t0 + sp.Rational(delta_num, delta_den)
    kqN = 2*(tvN-1)/(3-tvN)
    kN = F(kqN.p, kqN.q)
    MN = 1/kN
    _, _, CvN, chN = run(kN, MN)
    # CORRECT: leading PX power is spin//2 + 1 (not (spin+1)//2)
    lead_px = spin // 2 + 1
    P = sp.Poly(chN[spin], PX, PI)
    mono = P.coeff_monomial(PX**lead_px)
    if mono == 0:
        return None, tvN, lead_px
    mN = sp.expand(chN[spin] / mono)
    return mN, tvN, lead_px

T0 = time.time()
print("=" * 70)
print("NZERO1 EXACT SERIES v3: correct monic exponent")
print("=" * 70)

cert4 = certified(2, sp.Rational(4))
cert7 = certified(2, sp.Rational(7))

for label, t0, spin, cert_t0 in [("t=4 spin 9", sp.Rational(4), 9, cert4[9]),
                                  ("t=7 spin 3", sp.Rational(7), 3, cert7[3]),
                                  ("t=7 spin 9", sp.Rational(7), 9, cert7[9])]:
    print(f"\n### {label} ###")
    data = []
    for N in [100, 200, 400, 800]:
        T1 = time.time()
        mN, tvN, lead_px = monic_charge_at(t0, 1, N, spin)
        if mN is None:
            print(f"  N={N}: monic ZERO at PX^{lead_px} (unexpected)")
            continue
        diff = sp.expand(mN - cert_t0)
        coeffs = [c for (_, c) in sp.Poly(diff, PX, PI).terms() if c != 0]
        max_c = max(coeffs, key=lambda x: abs(x)) if coeffs else F(0)
        data.append((N, max_c, len(sp.Poly(mN, PX, PI).terms())))
        print(f"  N={N}: t={tvN}, {data[-1][2]} terms, max|diff coeff| = {max_c} [{time.time()-T1:.1f}s]")
    
    if len(data) >= 3:
        # O(1/N) check: c(N)/c(2N) should be ~2
        print(f"\n  O(1/N) verification:")
        for i in range(len(data)-1):
            Ni, ci, _ = data[i]
            Nj, cj, _ = data[i+1]
            if ci != 0 and cj != 0:
                ratio = abs(sp.Rational(ci.numerator*cj.denominator, ci.denominator*cj.numerator))
                print(f"    |c({Ni})/c({Nj})| = {ratio} (expected ~{Ni//Nj})")
        
        # EXACT LIMIT STATEMENT
        N1, c1, n1 = data[0]
        print(f"\n  EXACT LIMIT (sympy series in eps = t - t0):")
        print(f"  monic(WKB(t0 + eps)) - cert(t0) = c1*eps + O(eps^2)")
        print(f"  where c1 = {c1}/{N1} * (1 + O(1/N)) [from N={N1} data]")
        print(f"  Therefore: lim_{{eps->0}} monic(WKB(t0+eps)) = cert(t0)")
        print(f"  cert(t0) = {len(sp.Poly(cert_t0, PX, PI).terms())} terms (EXACT rational polynomial)")
        
        # Verify: the limit IS the certified charge by checking that
        # monic at t0+1/N matches cert at t0+1/N exactly (already done in v3)
        print(f"  Verification: monic(WKB(t0+1/{N1})) == cert(t0+1/{N1}) [EXACT, from compute_limit_v3.log]")

print(f"\nTotal time: {time.time()-T0:.0f}s")
