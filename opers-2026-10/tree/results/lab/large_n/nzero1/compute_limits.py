"""NZERO1: compute the WKB charge as an exact limit t -> t0 at the three exceptional cells,
compare with certified charges. Also run controls."""
import sys, os, time
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
os.chdir(os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
sys.path.insert(0, '.')
import common, gwkb, thwkb
from common import *

def gnorm(expr):
    """Reduce Gamma(q) for rational non-integer q to Gamma(frac), 0<frac<1."""
    rep = {}
    for g_ in expr.atoms(sp.gamma):
        q = g_.args[0]
        if not q.is_Rational or q.is_Integer: continue
        fac = sp.Integer(1); x = q
        while x > 1: x -= 1; fac *= x
        while x < 0: fac /= x; x += 1
        rep[g_] = fac * sp.gamma(x)
    return sp.cancel(sp.together(expr.xreplace(rep)))

def wkb_charge(tval, K_val=10):
    """Compute the WKB charge dict at a generic t (sympy symbol)."""
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
    ch = {}
    for kk in range(2, K_val + 1, 2):
        sk = sv[kk]
        has_pole = any(a.denominator == 1 and a <= 0 for (a, c) in sk)
        if has_pole:
            pw, poly, nonpole = thwkb.integrate_reg(sk, bb)
            poly = dict(poly)
        else:
            r = thwkb.integrate(sk, bb)
            poly = {k_: sp.Rational(v.numerator, v.denominator) for k_, v in r[2].items()} if r else {}
        ch[kk-1] = sp.expand(sum(v * (sp.Rational(Cv.numerator, Cv.denominator) * PX**2 / g['aY'])**(i // 2) *
                                   (sp.Rational(Cv.numerator, Cv.denominator) * PI**2 / g['aX'])**(j // 2)
                               for (i, j), v in poly.items()))
    return g, ch

def monic_normalize(ch, spin):
    """Divide by PX^(spin-1) coefficient; return None if that coefficient is zero."""
    P = sp.Poly(sp.expand(ch), PX, PI)
    mono = P.coeff_monomial(PX**(spin-1))
    if mono == 0:
        return None
    return sp.expand(ch / mono)

def limit_charge(t0, spin, K_val=10, eps_steps=3):
    """Compute the charge at t0+1/eps for several eps, find the pattern, extrapolate.
    For exactness: compute at t0 + 1/N for N = 2,3,4 and check convergence."""
    results = []
    for N in range(2, 5):
        tv = sp.Rational(t0) + sp.Rational(1, N*100)  # t0 + 1/(100N)
        try:
            g, ch = wkb_charge(tv, K_val)
            m = monic_normalize(ch[spin], spin)
            if m is not None:
                results.append((tv, m))
                print(f"    t={tv} ({tv.numerator}/{tv.denominator}): monic OK, {len(sp.Poly(m, PX, PI).terms())} terms", flush=True)
            else:
                results.append((tv, None))
                print(f"    t={tv}: monic coefficient still zero", flush=True)
        except Exception as ex:
            print(f"    t={tv}: ERROR {type(ex).__name__}: {str(ex)[:80]}", flush=True)
            results.append((tv, None))
    return results

def certified_and_compare(t0, spins):
    """Get certified charges at t0 and compare with the limit values."""
    tv0 = sp.Rational(t0)
    cert = certified(2, tv0)
    print(f"\n  Certified charges at t={t0}:")
    for s in spins:
        if s in cert:
            c = gnorm(cert[s])
            print(f"    spin {s}: {'ZERO' if c == 0 else f'{len(sp.Poly(c, PX, PI).terms())} terms'}")
    return cert

print("=" * 70)
print("NZERO1: Exact limits at exceptional cells")
print("=" * 70)

# --- Cell 1: t=4, spin 9 ---
print("\n### Cell 1: t=4 (k=-6, n=-9), spin 9 ###")
T0 = time.time()
cert = certified_and_compare('4', [3, 5, 7, 9])
print(f"  (certified computed in {time.time()-T0:.1f}s)")

print("\n  Computing WKB charges near t=4:")
res = limit_charge('4', 9)
if res and res[0][1] is not None:
    # Compare the first successful monic with certified
    tv_test, m_test = res[0]
    diff = sp.expand(gnorm(m_test - cert[9]))
    if diff == 0:
        print(f"\n  *** LIMIT MATCHES CERTIFIED at t={tv_test} ***")
    else:
        print(f"\n  MISMATCH vs certified at t={tv_test}: {len(sp.Poly(diff, PX, PI).terms())} terms differ")

# --- Cell 2: t=7, spin 3 ---
print("\n### Cell 2: t=7 (k=-3, n=-3), spin 3 ###")
cert7 = certified_and_compare('7', [1, 3, 5, 7, 9])
print("\n  Computing WKB charges near t=7:")
res7 = limit_charge('7', 3)
if res7 and res7[0][1] is not None:
    tv_test, m_test = res7[0]
    diff = sp.expand(gnorm(m_test - cert7[3]))
    if diff == 0:
        print(f"\n  *** LIMIT MATCHES CERTIFIED at t={tv_test} ***")
    else:
        print(f"\n  MISMATCH vs certified at t={tv_test}: {len(sp.Poly(diff, PX, PI).terms())} terms differ")

# --- Cell 3: t=7, spin 9 ---
print("\n### Cell 3: t=7 (k=-3, n=-3), spin 9 ###")
print("\n  Computing WKB charges near t=7:")
res7b = limit_charge('7', 9)
if res7b and res7b[0][1] is not None:
    tv_test, m_test = res7b[0]
    diff = sp.expand(gnorm(m_test - cert7[9]))
    if diff == 0:
        print(f"\n  *** LIMIT MATCHES CERTIFIED at t={tv_test} ***")
    else:
        print(f"\n  MISMATCH vs certified at t={tv_test}: {len(sp.Poly(diff, PX, PI).terms())} terms differ")

# --- Control 1: regular cell t=4, spin 7 ---
print("\n### CONTROL 1: t=4, spin 7 (regular, should match directly) ###")
g4, ch4 = wkb_charge('4')
m7 = monic_normalize(ch4[7], 7)
diff = sp.expand(gnorm(m7 - cert[7]))
print(f"  Direct value matches certified: {'YES' if diff == 0 else 'NO'}")

# --- Control 2: non-exceptional neighbour t=4+1/100 ---
print("\n### CONTROL 2: t=4+1/100 (non-exceptional neighbour) ###")
tv_n = sp.Rational('401/100')
kq_n = 2 * (tv_n - 1) / (3 - tv_n)
k_n = F(kq_n.p, kq_n.q)
print(f"  t=401/100, k={k_n}")
gn, chn = wkb_charge(tv_n)
for s in [3, 5, 7, 9]:
    ms = monic_normalize(chn[s], s)
    if ms is not None:
        cert_n = certified(2, tv_n)
        d = sp.expand(gnorm(ms - gnorm(cert_n[s])))
        print(f"  spin {s}: matches certified at t=401/100: {'YES' if d == 0 else 'NO'}")
    else:
        print(f"  spin {s}: monic coefficient zero (unexpected for non-exceptional)")

print(f"\nTotal time: {time.time()-T0:.0f}s")
