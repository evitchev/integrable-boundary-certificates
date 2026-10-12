"""NZERO1: trace where the zero comes from at t=4 spin 9, t=7 spins 3,9.
The engine path: classG -> blocks_of -> gwkb.wkb -> thwkb.integrate_reg -> charges.
The zero is in the Gamma/Beta prefactor of integrate_reg."""
import sys, os
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
os.chdir(os.path.expanduser('<repo>/results/lab/large_n/sol2_ext_fable'))
sys.path.insert(0, '.')
import common, gwkb, thwkb
from common import *

def trace_cell(tval, M_val, K_val=10):
    tv = sp.Rational(tval)
    kq = 2 * (tv - 1) / (3 - tv)
    k = F(kq.p, kq.q)
    g = classG(2, tv)
    n = g['D']
    bb = n * F(M_val)
    print(f"\n=== t={tval}, k={k}, M={F(M_val)}, n={n}, aX={g['aX']}, aY={g['aY']} ===")
    bl = blocks_of(g, F(M_val))
    sv = gwkb.wkb(bl, g['s'], F(M_val), K_val)
    for kk in range(2, K_val + 1, 2):
        spin = kk - 1
        sk = sv[kk]
        has_pole = any(a.denominator == 1 and a <= 0 for (a, c) in sk)
        if not has_pole:
            print(f"  spin {spin}: no pole terms -> regular integrate() path")
            continue
        # This is the interesting path
        pw, poly, nonpole = thwkb.integrate_reg(sk, bb)
        # Now build the charge
        Cv = 2 * n * (1 + F(M_val))
        ch = sp.expand(sum(v * (sp.Rational(Cv.numerator, Cv.denominator) * PX**2 / g['aY'])**(i // 2) *
                            (sp.Rational(Cv.numerator, Cv.denominator) * PI**2 / g['aX'])**(j // 2)
                        for (i, j), v in poly.items()))
        # Check if zero
        is_zero = (ch == 0)
        # Now trace WHY: look at each term in sk
        print(f"  spin {spin}: pole path, integrate_reg result = {'ZERO' if is_zero else 'non-zero'}")
        print(f"    common power (a+c) = {pw}, nonpole = {nonpole}")
        print(f"    s_k has {len(sk)} monomials:")
        for (a, c), p in sorted(sk.items()):
            a_r = sp.Rational(a.numerator, a.denominator)
            c_r = sp.Rational(c.numerator, c.denominator)
            # The Gamma prefactor in integrate() would be gamma(a) * gamma(-a-c) / gamma(-c)
            # In integrate_reg: residue of gamma(a) = (-1)^m / m! where a = -m
            if a_r.is_Integer and a_r <= 0:
                m = int(-a_r)
                res = F((-1)**m, sp.factorial(m))
                g_c = sp.gamma(sp.Rational(-c_r.numerator, -c_r.denominator)) if not ((-c_r).is_Integer and -c_r <= 0) else 0
                status = "POLE" if g_c != 0 else "ZERO (Gamma(-c) also poles)"
                print(f"      v^({a_r}) (v+E')^({c_r}): Gamma(a) pole, res={res}, Gamma(-c)={g_c} [{status}]")
            else:
                print(f"      v^({a_r}) (v+E')^({c_r}): non-pole term")
        if is_zero:
            print("    *** CHARGE IS ZERO ***")
            # Dig into integrate_reg internals
            expr = sk
            tot = {a + c for (a, c) in expr}
            print(f"    common a+c = {tot.pop()}")
            for (a, c), p in expr.items():
                if a.denominator == 1 and a <= 0:
                    m = int(-a)
                    res = sp.Rational((-1)**m, sp.factorial(m))
                    c_val = sp.Rational(c.numerator, c.denominator)
                    neg_c = -c_val
                    if neg_c.is_Integer and neg_c <= 0:
                        print(f"      TERM KILLED: a={a}, c={c}, -c={neg_c}: Gamma(-c)=Gamma({neg_c}) has a pole -> g=0")
                    else:
                        g_val = sp.gamma(neg_c)
                        print(f"      term: a={a}, c={c}, res*Gamma(-c) = {res}*{g_val} = {sp.simplify(res*g_val)}")

# Test the exceptional cells
trace_cell('4', F(1, 9))    # t=4, k=-9/7... wait, let me compute k correctly
