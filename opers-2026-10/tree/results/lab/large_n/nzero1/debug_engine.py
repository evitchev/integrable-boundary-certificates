"""Debug: compare my wkb_charge against ext_fibres.run() directly at a known-good point."""
import sys, os
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
os.chdir(os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
sys.path.insert(0, '.')
import common, gwkb, thwkb
from common import *
# Import ext_fibres run() 
sys.argv = ['ext_fibres.py']  # prevent __main__ block
exec(open('ext_fibres.py').read().split("if __name__ != '__main__':")[0])

# Test at t=9/4 (known good fibre)
print("=== Comparing my impl vs ext_fibres.run() at t=9/4 ===")
tv_test = sp.Rational('9/4')
kq_test = 2 * (tv_test - 1) / (3 - tv_test)
k_test = F(kq_test.p, kq_test.q)
print(f"t=9/4, k={k_test}")

# Method A: ext_fibres.run()
tv_A, g_A, Cv_A, ch_A = run(k_test, 1/k_test)
print(f"A (run): n={g_A['D']}, C={Cv_A}")

# Method B: my wkb_charge equivalent
g_B = classG(2, tv_test)
n_B = g_B['D']
M_B = 1/k_test
bb_B = n_B * M_B
Cv_B = 2 * n_B * (1 + M_B)
bl_B = blocks_of(g_B, M_B)
sv_B = gwkb.wkb(bl_B, g_B['s'], M_B, 10)
ch_B = {}
for kk in range(2, 11, 2):
    sk = sv_B[kk]
    has_pole = any(a.denominator == 1 and a <= 0 for (a, c) in sk)
    if has_pole:
        pw, poly, nonpole = thwkb.integrate_reg(sk, bb_B)
        poly = dict(poly)
    else:
        r = thwkb.integrate(sk, bb_B)
        poly = {k_: sp.Rational(v.numerator, v.denominator) for k_, v in r[2].items()} if r else {}
    ch_B[kk-1] = sp.expand(sum(v * (sp.Rational(Cv_B.numerator, Cv_B.denominator) * PX**2 / g_B['aY'])**(i // 2) *
                                 (sp.Rational(Cv_B.numerator, Cv_B.denominator) * PI**2 / g_B['aX'])**(j // 2)
                             for (i, j), v in poly.items()))

print(f"B (manual): n={n_B}, C={Cv_B}")
print(f"Same g? {g_A == g_B}  Same Cv? {Cv_A == Cv_B}")
for s in [1, 3, 5]:
    dA = sp.Poly(ch_A[s], PX, PI).terms()
    dB = sp.Poly(ch_B[s], PX, PI).terms()
    print(f"  spin {s}: A has {len(dA)} terms, B has {len(dB)} terms")
    # Check if they match up to monic
    mA = sp.Poly(ch_A[s], PX, PI).coeff_monomial(PX**(s-1))
    mB = sp.Poly(ch_B[s], PX, PI).coeff_monomial(PX**(s-1))
    if mA != 0 and mB != 0:
        diff = sp.expand(sp.simplify(ch_A[s]/mA - ch_B[s]/mB))
        print(f"    monic diff: {'ZERO' if diff == 0 else f'{len(sp.Poly(diff,PX,PI).terms())} terms'}")
    else:
        print(f"    monic coeff A={mA}, B={mB}")

# Now test at t=4 (exceptional)
print("\n=== At t=4 ===")
tv_4 = sp.Rational(4)
k_4 = F(-6)
tv_A2, g_A2, Cv_A2, ch_A2 = run(k_4, F(-1,6))
print(f"run(): spin 9 charge = {ch_A2[9] if ch_A2.get(9) is not None else 'MISSING'}")
P9 = sp.Poly(ch_A2[9], PX, PI) if ch_A2.get(9) else sp.Poly(0, PX, PI)
print(f"  PX powers in spin 9: {sorted(m[0] for m in P9.terms())}")
print(f"  PX^8 coeff: {P9.coeff_monomial(PX**8)}")

# Check what the certified side looks like
cert = certified(2, tv_4)
c9 = sp.expand(cert[9])
Pc9 = sp.Poly(c9, PX, PI)
print(f"\n  Certified spin 9: {len(Pc9.terms())} terms, PX powers = {sorted(m[0] for m in Pc9.terms())}")
print(f"  Certified spin 9 (first 3 terms): {list(Pc9.terms())[:3]}")
