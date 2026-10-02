"""VIR8 (c), exploratory (sealed as such): at t = 5/3 the quintic Gamma operator T (T^2-l0^2)(T^2-l1^2) psi = x^5 (x^5 + E) psi
is the k = 1 member of several class-U families (a = 1: any subset of roots may sit in L).  For the symmetric splits,
the limit k -> 1 of the monic spin-5 charge of the family, expressed in the kernel basis {S, E_1}:
   S = certified Sol-3 charge, E_1 = prod_{j=-1..1}(PX + pi - j p)(PX - pi - j p), p^2 = 1/2.
Family conventions: strings of length k in L, unit roots in P, n = deg P + k deg L, M = 1/k, l1 = n pi/p and
l0^2 : l1^2 from the powers (power-alpha factor <-> PX^2/alpha etc.)."""
import sys, time
import sympy as sp
sys.dont_write_bytecode = True
import wkb3 as W3
from b_core import table, PX, PI
l0, l1 = W3.l0, W3.l1
tv0 = sp.Rational(5, 3)
S = table(3, 5, tv0)
u, v = PX + PI, PX - PI
E1 = sp.Poly(sp.expand(u * (u**2 - sp.Rational(1, 2)) * v * (v**2 - sp.Rational(1, 2))), PX, PI)
S2 = table(2, 5, tv0)
al, be = sp.symbols('al be')
t0 = time.time()
splits = {
    'Sol 3:  P = (T^2-l0^2)(T^2-l1^2), L = T':        (4, 1, lambda c: ([1, 0, -(l0**2 + l1**2), 0, l0**2 * l1**2], [1, -c / 2]), (1, 1)),
    'Sol 2:  P = T(T^2-l1^2), L = T^2-l0^2':           (3, 2, lambda c: ([1, 0, -l1**2, 0], [1, -c, c * c / 4 - l0**2]), ('k', 1)),
    'swap:   P = T(T^2-l0^2), L = T^2-l1^2':           (3, 2, lambda c: ([1, 0, -l0**2, 0], [1, -c, c * c / 4 - l1**2]), (1, 'k')),
    'both:   P = T, L = (T^2-l0^2)(T^2-l1^2)':         None,
}
for name, spec in splits.items():
    if spec is None:
        print('%-45s not in the class (deg P must exceed deg L)' % name); continue
    dP, dL, mk, powers = spec
    res = {}
    for hexp in (8, 16):
        k = 1 + sp.Rational(1, 10**hexp)
        M = 1 / k
        n = dP + dL * k
        c = n * (1 + M)
        pc, lc = mk(c)
        eng = W3.WKB3(dP, dL, M, pc, lc)
        eng.solve(6)
        cB, cG = eng.J(6)
        assert cG.is_zero
        tv = (3 * k + 2) / (k + 2); p2 = 2 * (tv - 1) / (tv + 1)
        c2 = n * n / p2
        aX = k if powers[0] == 'k' else 1
        aYp = k if powers[1] == 'k' else 1
        Pp = sp.Poly(cB, l0, l1)
        out = sp.Poly(sp.expand(sum(cf * (c2 / aX) ** sp.Rational(i, 2) * (c2 / aYp) ** sp.Rational(j, 2) * PX**i * PI**j for (i, j), cf in Pp.terms())), PX, PI)
        out = sp.Poly(sp.expand(out.as_expr() / out.coeff_monomial(PX**6)), PX, PI)
        # least-squares-free: solve al, be from two coefficients, report the max residual over all
        eqs = sp.Poly(out.as_expr() - al * S.as_expr() - be * E1.as_expr(), PX, PI)
        sol = sp.solve([eqs.coeff_monomial(PX**6), eqs.coeff_monomial(PX**2 * PI**4)], [al, be], dict=True)[0]
        resid = max(abs(cf) for cf in sp.Poly(eqs.as_expr().subs(sol), PX, PI).coeffs())
        res[hexp] = (sol[al], sol[be], resid)
    print('%-45s limit = al S + be E_1 with al = %.10f, be = %.10f (h = 1e-16); residual %.2e (h=1e-8), %.2e (h=1e-16)  (%.0fs)'
          % (name, float(res[16][0]), float(res[16][1]), float(res[8][2]), float(res[16][2]), time.time() - t0), flush=True)
