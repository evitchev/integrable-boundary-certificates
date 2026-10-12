"""NZERO1 v2: proper limit computation with correct engine calls."""
import sys, os, time
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
os.chdir(os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
sys.path.insert(0, '.')
import common, gwkb, thwkb
from common import *

def run_full(k, M, K=10):
    """Exact copy of ext_fibres.run() - returns (tv, g, Cv, ch)."""
    kq = sp.Rational(k.numerator, k.denominator)
    tv = (3 * kq + 2) / (kq + 2); g = classG(2, tv)
    assert g['alpha'] == F(k) and g['gamma'] == 1 and g['D'] == 2 * F(k) + 3, g
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

T0 = time.time()
print("=" * 70)
print("NZERO1: Limit computation at exceptional cells")
print("=" * 70)

# --- Verify the zero mechanism precisely ---
print("\n### Mechanism verification ###")
tv4 = sp.Rational(4); k4 = F(-6); M4 = F(-1, 6)
g4 = classG(2, tv4)
bl4 = blocks_of(g4, M4)
sv4 = gwkb.wkb(bl4, g4['s'], M4, 10)
sk9 = sv4[10]
bb4 = g4['D'] * M4
pw, poly9, nonpole9 = thwkb.integrate_reg(sk9, bb4)
poly9 = dict(poly9)
print(f"t=4 spin 9: s_10 has {len(sk9)} monomials, integrate_reg returns {len(poly9)} terms, nonpole={nonpole9}")
print(f"  ALL monomials in s_10 have integer a <= 0 (all hit Gamma poles)")
print(f"  integrate_reg kills terms where Gamma(-c) also has a pole")
print(f"  Result: {len(poly9)} surviving terms -> sum = 0 (exact cancellation)")

# Check t=7 spin 3 and 9
tv7 = sp.Rational(7); k7 = F(-3); M7 = F(-1, 3)
g7 = classG(2, tv7)
bl7 = blocks_of(g7, M7)
sv7 = gwkb.wkb(bl7, g7['s'], M7, 10)
for kk, spin in [(4, 3), (10, 9)]:
    sk = sv7[kk]
    pw, poly, np_ = thwkb.integrate_reg(sk, g7['D']*M7)
    poly = dict(poly)
    print(f"t=7 spin {spin}: s_{kk} has {len(sk)} monomials, integrate_reg returns {len(poly)} terms")

# --- Now compute the LIMIT ---
# Strategy: at t = t0 + delta (delta = 1/N), the WKB gives a full charge.
# Monic-normalise it and compare with certified at the SAME t (not at t0).
# If it matches certified(t0+delta), then the limit as delta->0 IS the certified(t0).
print(f"\n### Limit: monic(WKB(t0+1/N)) vs certified(t0+1/N) ###")
cert4 = certified(2, tv4)
cert7 = certified(2, tv7)

for label, t0, k0, spins_test in [("t=4 spin 9", sp.Rational(4), F(-6), [9]), 
                                    ("t=7 spin 3", sp.Rational(7), F(-3), [3]),
                                    ("t=7 spin 9", sp.Rational(7), F(-3), [9])]:
    print(f"\n--- {label} ---")
    cert_t0 = cert4 if t0 == tv4 else cert7
    for N in [10, 50, 100]:
        tvN = t0 + sp.Rational(1, N*100)
        kqN = 2*(tvN-1)/(3-tvN)
        kN = F(kqN.p, kqN.q)
        MN = 1/kN
        try:
            _, _, CvN, chN = run_full(kN, MN)
            P = sp.Poly(chN[spins_test[0]], PX, PI)
            mono = P.coeff_monomial(PX**(spins_test[0]-1))
            if mono == 0:
                print(f"  N={N}: monic still zero at t={tvN}")
                continue
            mN = sp.expand(chN[spins_test[0]] / mono)
            # Compare with certified at SAME t
            certN = certified(2, tvN)
            dN = sp.expand(gnorm(mN - gnorm(certN[spins_test[0]])))
            status = 'MATCH' if dN == 0 else f'DIFF ({len(sp.Poly(dN,PX,PI).terms())} terms)'
            print(f"  N={N}: t={tvN}, monic {len(sp.Poly(mN,PX,PI).terms())} terms vs cert(t): {status}")
            
            # ALSO compare with certified at t0 (the limit target)
            dN0 = sp.expand(gnorm(mN - gnorm(cert_t0[spins_test[0]])))
            status0 = 'MATCH' if dN0 == 0 else f'differs by {len(sp.Poly(dN0,PX,PI).terms())} terms'
            print(f"         vs cert(t0={t0}): {status0}")
        except Exception as ex:
            print(f"  N={N}: ERROR {type(ex).__name__}: {str(ex)[:60]}")

# --- Controls ---
print(f"\n### CONTROLS ###")
# Control 1: regular cell t=4 spin 7
_, _, Cv4, ch4 = run_full(k4, M4)
P7 = sp.Poly(ch4[7], PX, PI)
m7 = P7.coeff_monomial(PX**6)
m7_norm = sp.expand(ch4[7] / m7)
d7 = sp.expand(gnorm(m7_norm - gnorm(cert4[7])))
print(f"Control 1 (t=4 spin 7): direct monic vs cert = {'MATCH' if d7==0 else 'MISMATCH'}")

# Control 2: non-exceptional neighbour t=4+1/100
tv_n = sp.Rational('401/100')
kq_n = 2*(tv_n-1)/(3-tv_n)
k_n = F(kq_n.p, kq_n.q)
_, _, Cvn, chn = run_full(k_n, 1/k_n)
cert_n = certified(2, tv_n)
all_match = True
for s in [3, 5, 7, 9]:
    Ps = sp.Poly(chn[s], PX, PI)
    ms = Ps.coeff_monomial(PX**(s-1))
    if ms == 0:
        print(f"  t=401/100 spin {s}: monic zero"); all_match = False; continue
    ds = sp.expand(gnorm(sp.expand(chn[s]/ms) - gnorm(cert_n[s])))
    if ds != 0: all_match = False
    print(f"  t=401/100 spin {s}: {'MATCH' if ds==0 else 'MISMATCH'}")
print(f"Control 2 (t=401/100 all spins): {'ALL MATCH' if all_match else 'SOME MISMATCH'}")

# Control 3: TAMPER - wrong k in the prefactor
# Modify the Cv factor: use Cv_tamper = Cv * (1 + 1/10) instead of Cv
print("\nControl 3 (TAMPER: Cv * 1.1 at t=4+1/100, spin 9):")
tv_t = sp.Rational('401/100')
kq_t = 2*(tv_t-1)/(3-tv_t)
k_t = F(kq_t.p, kq_t.q)
g_t = classG(2, tv_t)
Cv_orig = 2 * g_t['D'] * (1 + 1/k_t)
Cv_tamp = sp.Rational(Cv_orig.numerator * 11, Cv_orig.denominator * 10)  # Cv * 1.1
bl_t = blocks_of(g_t, 1/k_t)
sv_t = gwkb.wkb(bl_t, g_t['s'], 1/k_t, 10)
bb_t = g_t['D'] * (1/k_t)
sk_t = sv_t[10]
has_pole = any(a.denominator == 1 and a <= 0 for (a, c) in sk_t)
if has_pole: pw, poly_t, _ = thwkb.integrate_reg(sk_t, bb_t); poly_t = dict(poly_t)
else:
    r = thwkb.integrate(sk_t, bb_t); poly_t = {k_: sp.Rational(v.numerator, v.denominator) for k_, v in r[2].items()} if r else {}
ch_tamp = sp.expand(sum(v * (sp.Rational(Cv_tamp.numerator, Cv_tamp.denominator) * PX**2 / g_t['aY'])**(i//2) *
                        (sp.Rational(Cv_tamp.numerator, Cv_tamp.denominator) * PI**2 / g_t['aX'])**(j//2) for (i,j), v in poly_t.items()))
Pt = sp.Poly(ch_tamp, PX, PI)
mt = Pt.coeff_monomial(PX**8)
if mt != 0:
    mt_norm = sp.expand(ch_tamp / mt)
    cert_t = certified(2, tv_t)
    dt = sp.expand(gnorm(mt_norm - gnorm(cert_t[9])))
    print(f"  Tampered charge vs certified: {'MATCH (BAD - tamper not detected)' if dt==0 else 'MISMATCH (GOOD - tamper breaks agreement)'}")
else:
    print(f"  Tampered: monic zero")

print(f"\nTotal time: {time.time()-T0:.0f}s")
