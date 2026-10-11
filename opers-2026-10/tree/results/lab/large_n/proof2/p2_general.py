"""PROOF2 (2) VERIFIED extension beyond Sols 2, 3: general three-term operators, using Fable's engines (copied, read-only).
P = (T^2 - l1^2) prod (T - r) [extra rational roots], L = (T - l0) prod (T - q) [extra rational roots]; non-integer k.
ODE engine W3.WKB3(dP, dL, M, P coeffs, Lc = L(T - c/2) coeffs) vs Gamma form MG (blocks: P factors and one string G_k(T - root) per root of L).
Usage: p2_general.py k 'r1,r2' 'q1' order [tamper]"""
import sys, time; sys.dont_write_bytecode = True
sys.path.insert(0, '.')
import sympy as sp
import wkb3 as W3, mellin_gen as MG
k = sp.Rational(sys.argv[1]); rs = [sp.Rational(x) for x in sys.argv[2].split(',') if x not in ('', '-')]; qs = [sp.Rational(x) for x in sys.argv[3].split(',') if x not in ('', '-')]
Kmax = int(sys.argv[4]); tamper = len(sys.argv) > 5
M = 1 / k; l0, l1 = W3.l0, W3.l1; T = sp.Symbol('T')
P = sp.expand((T**2 - l1**2) * sp.Mul(*[(T - r) for r in rs])); Lp = sp.expand((T - l0) * sp.Mul(*[(T - q) for q in qs]))
dP, dL = sp.degree(P, T), sp.degree(Lp, T); n = dP + k * dL; c = n * (1 + M)
sh = 0 if tamper else c / 2
Lc = sp.expand(Lp.subs(T, T - sh))
pco = sp.Poly(P, T).all_coeffs(); lco = sp.Poly(Lc, T).all_coeffs()
eng = W3.WKB3(dP, dL, M, pco, lco); assert eng.n == n, (eng.n, n)
eng.solve(Kmax)
blocks = [[1, 0, -MG.l1**2] + [0] * (Kmax - 1)] + [[1, -r] + [0] * Kmax for r in rs]
blocks += [MG.string_series(k, MG.l0, n / k, Kmax + 1)] + [MG.string_series(k, q, n / k, Kmax + 1) for q in qs]
g = MG.MellinWKB(n, M, MG.symbol_from_blocks(blocks, Kmax + 1)); g.solve(Kmax)
def monic(Pp):
    Pp = sp.Poly(Pp, l0, l1)
    if Pp.is_zero: return None
    lead = Pp.coeff_monomial(l0 ** Pp.total_degree())
    return None if lead == 0 else sp.Poly(sp.expand(Pp.as_expr() / lead), l0, l1)
ok = True
for i in range(2, Kmax + 1):
    cB, cG = eng.J(i); Rg, _ = g.R(i)
    Rg = sp.Poly(sp.expand(sp.sympify(Rg).subs({MG.l0: l0, MG.l1: l1})), l0, l1)
    mB, mG, mR = monic(cB), monic(cG), monic(Rg)
    sameB = (cB.is_zero and Rg.is_zero) or (mB is not None and mR is not None and mB == mR)
    sameG = cG.is_zero or (mG is not None and mR is not None and mG == mR)
    A0 = (1 - i) / (n / k); degen = A0.is_integer and A0 <= 0
    print(f'k = {k}, P extra roots {rs}, L extra roots {qs}{" TAMPER" if tamper else ""}: order {i}: A0 = {A0}{" (BETA-POLE CELL)" if degen else ""}; cB zero {cB.is_zero}; R zero {Rg.is_zero}; cG zero {cG.is_zero}; monic(cB) == monic(R): {sameB and sameG}', flush=True)
    if degen: continue
    ok &= sameB and sameG
print('GENERAL', 'AGREE' if ok else 'DISAGREE'); sys.exit((0 if ok else 2) if not tamper else (0 if not ok else 2))
