import sys, os
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
os.chdir(os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
sys.path.insert(0, '.')
import common, gwkb, thwkb
from common import *

def is_neg_int(x): return x.denominator == 1 and x <= 0

def trace_cell(tval, K_val=10):
    tv = sp.Rational(tval)
    kq = 2 * (tv - 1) / (3 - tv)
    k = F(kq.p, kq.q)
    g = classG(2, tv)
    n = g['D']
    Mv = 1 / k
    bb = n * Mv
    print(f"\n=== t={tval}, k={k}, n={n}, M={Mv}, aX={g['aX']}, aY={g['aY']} ===")
    bl = blocks_of(g, Mv)
    sv = gwkb.wkb(bl, g['s'], Mv, K_val)
    for kk in range(2, K_val + 1, 2):
        spin = kk - 1
        sk = sv[kk]
        has_pole = any(a.denominator == 1 and a <= 0 for (a, c) in sk)
        if not has_pole:
            r = thwkb.integrate(sk, bb)
            poly = {k_: sp.Rational(v.numerator, v.denominator) for k_, v in r[2].items()} if r else {}
            Cv = 2 * n * (1 + Mv)
            ch = sp.expand(sum(v * (sp.Rational(Cv.numerator, Cv.denominator) * PX**2 / g['aY'])**(i // 2) *
                                (sp.Rational(Cv.numerator, Cv.denominator) * PI**2 / g['aX'])**(j // 2)
                            for (i, j), v in poly.items()))
            print(f"  spin {spin}: REGULAR, charge {'ZERO' if ch==0 else 'non-zero'}")
            continue
        pw, poly, nonpole = thwkb.integrate_reg(sk, bb)
        poly = dict(poly)
        Cv = 2 * n * (1 + Mv)
        ch = sp.expand(sum(v * (sp.Rational(Cv.numerator, Cv.denominator) * PX**2 / g['aY'])**(i // 2) *
                            (sp.Rational(Cv.numerator, Cv.denominator) * PI**2 / g['aX'])**(j // 2)
                        for (i, j), v in poly.items()))
        nu = F(spin) / F(int(n))
        knu = k * nu
        n1k = nu * (1 + k)
        both_pole = is_neg_int(-knu) and is_neg_int(nu)
        label = "POLE/POLE (inf/inf)" if both_pole else "evaluable"
        print(f"  spin {spin}: POLE path, charge {'*** ZERO ***' if ch==0 else 'non-zero'}")
        print(f"    nu={nu}, -k*nu={-knu}, nu*(1+k)={n1k}  B args -> {label}")
        killed_count = 0
        alive_count = 0
        for (a, c), p in sorted(sk.items()):
            a_r = F(a.numerator, a.denominator)
            c_r = F(c.numerator, c.denominator)
            neg_c = -c_r
            if a_r.denominator == 1 and a_r <= 0:
                m = int(-a_r)
                if is_neg_int(neg_c):
                    killed_count += 1
                    print(f"      v^{a_r}(v+E')^{c_r}: res={F((-1)**m, sp.factorial(m))} x Gamma({neg_c})=KILLED (pole)")
                else:
                    alive_count += 1
                    gv = sp.gamma(sp.Rational(neg_c.numerator, neg_c.denominator))
                    print(f"      v^{a_r}(v+E')^{c_r}: res={F((-1)**m, sp.factorial(m))} x Gamma({neg_c})={gv}")
            else:
                alive_count += 1
        print(f"    total terms: {len(sk)}, killed: {killed_count}, alive: {alive_count}")

trace_cell('4')
trace_cell('7')
# Controls: regular cells
print("\n--- CONTROLS ---")
trace_cell('4')  # spin 7 should be regular
trace_cell('9/4')  # nearby non-exceptional fibre
