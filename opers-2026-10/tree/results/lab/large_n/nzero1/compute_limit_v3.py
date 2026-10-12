"""NZERO1 v3: CORRECT monic normalisation (divide by PX^(spin+1), not PX^(spin-1))."""
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

def score_check(ch_s, cert_s, s):
    """Exact ext_fibres.score() logic."""
    gn = gnorm(ch_s)
    if gn == 0:
        return 'WKB_ZERO', None
    d = sp.expand(gnorm(sp.expand(monic(ch_s, (s+1)//2) - cert_s)))
    return ('MATCH' if d == 0 else f'DIFF({len(sp.Poly(d,PX,PI).terms())})'), d

T0 = time.time()
print("=" * 70)
print("NZERO1 v3: Limits with CORRECT monic (PX^(spin+1))")
print("=" * 70)

# --- Mechanism (confirmed from earlier runs) ---
print("""
### Mechanism (from trace_zero.log + verify_zero.log) ###
At t=4 (k=-6, n=-9) spin 9 and t=7 (k=-3, n=-3) spins 3,9:
- The WKB s_k contains ONLY monomials v^a(v+E')^c with a a non-positive integer.
- thwkb.integrate_reg() replaces Gamma(a) by its residue (-1)^m/m! for a=-m.
- For terms where Gamma(-c) ALSO has a pole (i.e., -c is a non-positive integer),
  the product res(Gamma(a)) * Gamma(-c) = finite * infinity -> killed (set to 0).
- At these exceptional cells, ALL surviving terms cancel exactly in the sum.
  integrate_reg returns an EMPTY dict (0 terms).
- File:line: thwkb.py ~L218-225 (integrate_reg), the `if g == 0: continue` line.
- The zero is then detected in ext_fibres.py score() via gnorm(ch[s]) == 0.
""")

# --- Limit computation ---
print("### Limit: monic(WKB(t0+delta)) vs certified(t0+delta) ###")
cert4 = certified(2, sp.Rational(4))
cert7 = certified(2, sp.Rational(7))

cells = [
    ("t=4 spin 9", sp.Rational(4), F(-6), 9),
    ("t=7 spin 3", sp.Rational(7), F(-3), 3),
    ("t=7 spin 9", sp.Rational(7), F(-3), 9),
]

for label, t0, k0, spin in cells:
    print(f"\n--- {label} ---")
    cert_t0 = cert4 if t0 == sp.Rational(4) else cert7
    for delta_denom in [100, 500, 1000]:
        tvN = t0 + sp.Rational(1, delta_denom)
        kqN = 2*(tvN-1)/(3-tvN)
        kN = F(kqN.p, kqN.q)
        MN = 1/kN
        T1 = time.time()
        try:
            _, _, CvN, chN = run(kN, MN)
            status, diff = score_check(chN[spin], certified(2, tvN)[spin], spin)
            # Also vs cert at t0
            d0, _ = score_check(chN[spin], cert_t0[spin], spin) if gnorm(chN[spin]) != 0 else ('WKB_ZERO', None)
            print(f"  t={tvN} (+1/{delta_denom}): vs cert(same t): {status} | vs cert(t0): {d0} [{time.time()-T1:.1f}s]")
        except Exception as ex:
            print(f"  t={tvN}: ERROR {type(ex).__name__}: {str(ex)[:60]}")

# --- Controls ---
print(f"\n### CONTROLS ###")
# C1: regular cell t=4 spin 7
_, _, Cv4, ch4 = run(F(-6), F(-1, 6))
s, d = score_check(ch4[7], cert4[7], 7)
print(f"C1 (t=4 spin 7, regular): {s}")

# C2: non-exceptional neighbour t=4+1/100, all spins
tv_n = sp.Rational('401/100')
kq_n = 2*(tv_n-1)/(3-tv_n); k_n = F(kq_n.p, kq_n.q)
_, _, Cvn, chn = run(k_n, 1/k_n)
cert_n = certified(2, tv_n)
all_ok = True
for s_ in [1, 3, 5, 7, 9]:
    st, _ = score_check(chn[s_], cert_n[s_], s_)
    if st != 'MATCH': all_ok = False
    print(f"  C2 t=401/100 spin {s_}: {st}")
print(f"C2 summary: {'ALL MATCH' if all_ok else 'SOME MISMATCH'}")

# C3: TAMPER - shift delta by 1/10 in one of the blocks (as ext_fibres does)
# Actually, the brief says "wrong k in the Beta prefactor". The Beta prefactor is
# implicit in the integrate_reg step. Tamper: use k+1/10 instead of k for M=1/k.
print("\nC3 (TAMPER: M -> M/(1+eps) i.e. k -> k*(1+eps), eps=1/10):")
tv_t = sp.Rational('401/100')
kq_t = 2*(tv_t-1)/(3-tv_t); k_t = F(kq_t.p, kq_t.q)
M_tamp = 1/(k_t * sp.Rational(11, 10))  # k_tamper = k*1.1, so M_tamper = 1/(k*1.1)
try:
    _, _, Cv_t, ch_t = run(k_t, M_tamp)
    cert_t = certified(2, tv_t)
    st, _ = score_check(ch_t[9], cert_t[9], 9)
    print(f"  Tampered at t=401/100 spin 9: {st} ({'GOOD - breaks' if st != 'MATCH' else 'BAD - still matches'})")
except Exception as ex:
    print(f"  Tamper raised: {type(ex).__name__}: {str(ex)[:60]} (GOOD - engine inconsistent)")

# C4: Tamper at the exceptional cell itself
print("\nC4 (TAMPER at t=4 spin 9: perturb M by 1%):")
M_tamp2 = F(-1, 6) + F(1, 600)
try:
    _, _, Cv_t2, ch_t2 = run(F(-6), M_tamp2)
    st2, d2 = score_check(ch_t2[9], cert4[9], 9)
    print(f"  Tampered at t=4 spin 9: {st2}")
except Exception as ex:
    print(f"  Tamper raised: {type(ex).__name__}: {str(ex)[:60]}")

print(f"\nTotal time: {time.time()-T0:.0f}s")
