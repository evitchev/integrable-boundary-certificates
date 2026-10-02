"""VIR8 engine validation (ODE side only; no certified data):
(i)  BLZ/Suzuki family (dP, dL) = (2, 1): the three-term form  [ (T^2 - l1^2) - x^(c/2)(T - l0)x^(c/2) - E x^n ] phi = 0
     against the Gamma form (T^2 - l1^2) G_a(T - l0) (mellin_gen), which VIR7 showed equal to Suzuki + centrifugal and to
     the Sol-1 tables.  Printed: the two master coefficients cB, cG of the three-term charge and the Gamma-form R_i.
Exit 0 iff at every even order the three-term charge is proportional to the Gamma-form charge (cG either zero or
proportional to cB, and monic(cB) == monic(R))."""
import sys, time
import sympy as sp
sys.dont_write_bytecode = True
import wkb3 as W3
import mellin_gen as MG

def monic(P, var):
    P = sp.Poly(P, W3.l0, W3.l1)
    if P.is_zero:
        return None
    lead = P.coeff_monomial(var)
    return None if lead == 0 else sp.Poly(sp.expand(P.as_expr() / lead), W3.l0, W3.l1)

ok = True
t0 = time.time()
for a in (sp.Rational(5, 8), sp.Rational(2, 3)):
    M = 1 / a
    n = 2 + a
    c = n * (1 + M)
    K = 6
    # three-term: P = T^2 - l1^2 ; Lc(T) = T - c/2 - l0
    eng = W3.WKB3(2, 1, M, [1, 0, -W3.l1**2], [1, -c / 2 - W3.l0])
    eng.solve(K)
    # Gamma form
    blocks = [[1, 0, -MG.l1**2] + [0] * (K - 1), MG.string_series(a, MG.l0, n / a, K + 1)]
    g = MG.MellinWKB(n, M, MG.symbol_from_blocks(blocks, K + 1)); g.solve(K)
    for i in range(2, K + 1):
        cB, cG = eng.J(i)
        Rg, _ = g.R(i)
        Rg = sp.Poly(sp.expand(Rg.subs({MG.l0: W3.l0, MG.l1: W3.l1})), W3.l0, W3.l1)
        var = W3.l0**i if i % 2 == 0 else W3.l0**i
        mB, mG, mR = monic(cB, W3.l0**i), monic(cG, W3.l0**i), monic(Rg, W3.l0**i)
        sameB = (mB is None and mR is None) or (mB is not None and mR is not None and mB == mR)
        sameG = cG.is_zero or (mG is not None and mR is not None and mG == mR)
        ok &= sameB and sameG
        print('a = %s order %d (spin %d): cG zero: %s; monic(cB) == monic(Gamma form): %s; monic(cG) == monic(Gamma form) or cG = 0: %s   (%.0fs)'
              % (a, i, i - 1, cG.is_zero, sameB, sameG, time.time() - t0), flush=True)
print('three-term engine validation', 'PASS' if ok else 'FAIL')
sys.exit(0 if ok else 2)
