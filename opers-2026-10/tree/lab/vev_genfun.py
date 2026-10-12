"""Vacuum expectation values by generating functions (2026-09-05): the same quantity as code/cyl_vev.py (the certified
Wick-pairing recipe, validated on hep-th/0404195 eqs. (59)-(60)), computed without enumerating pairings.

  * X-sector: letters d^k J (J = dX) with contraction C(a, b) = (-1)^b R^{(a+b)}(0), R(eps) = 1/(4 sinh^2(eps/2)) - 1/eps^2,
    and zero mode P for the letter J itself: the sum over partial pairings is the coefficient of prod t_a^{n_a}/n_a! in
    exp(P t_0 + (1/2) sum_{ab} t_a t_b C(a, b))  -- computed here by a memoized recursion over the multiset of orders.
  * Y-sector: pair letters (d^a Y . d^b Y) of an O(N) vector field with the same contraction C and zero mode Q e for a = 0:
    <exp(sum_j s_j Y_{a_j}.Y_{b_j})> = det(1 - 2 S C)^{-N/2} exp( Q^2 [ 2 S (1 - 2 C S)^{-1} ]_{00} ) with S the symmetric
    matrix of the quadratic form; closed cycles give the powers of N, zero-mode chains the powers of Q^2 -- the recipe's
    combinatorics.  The value of a monomial is the coefficient of prod s_j (with multiplicities) in the series.
  * Contractions pair only same-species letters, so the value of a mixed monomial is the product of the two sectors' values.
Article momenta: P_art^2 = -P^2, Q_art^2 = -Q^2 (as in cyl_vev.vev_polynomial).  Exact (sympy rationals); N, P, Q symbolic or numeric.
Validated against cyl_vev.vev_polynomial on the certified weight-4..10 densities of the three families (all coefficients)."""
import sys, os
from functools import lru_cache
import sympy as sp
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(ROOT, 'code'))
from cyl_vev import _Rd, C as C_recipe
P, Q, Nsym = sp.symbols('P Q N')
def Cc(a, b): return sp.Rational(C_recipe(a, b))
@lru_cache(maxsize=None)
def x_value(orders):
    """sum over partial pairings of X-letters with J-derivative orders `orders` (sorted tuple), zero mode P for order 0"""
    if not orders: return sp.Integer(1)
    first, rest = orders[0], orders[1:]
    total = (P if first == 0 else 0) * x_value(rest)
    for i, o in enumerate(rest):
        total += Cc(first, o) * x_value(rest[:i] + rest[i + 1:])
    return sp.expand(total)
def y_value(pairs, N):
    """pairs: tuple of (a, b) J-derivative orders of the Y pair-letters; returns the O(N) vacuum value (polynomial in Q^2, N)"""
    if not pairs: return sp.Integer(1)
    orders = sorted({o for pr in pairs for o in pr}); idx = {o: i for i, o in enumerate(orders)}; n = len(orders)
    from collections import Counter
    mult = Counter(tuple(sorted(pr)) for pr in pairs); kinds = sorted(mult)
    s = sp.symbols(f's0:{len(kinds)}'); eps = sp.Symbol('eps')
    S = sp.zeros(n, n)
    for sk, (a, b) in zip(s, kinds):
        if a == b: S[idx[a], idx[a]] += sk
        else: S[idx[a], idx[b]] += sk / 2; S[idx[b], idx[a]] += sk / 2
    Cm = sp.Matrix(n, n, lambda i, j: Cc(orders[i], orders[j]))
    M = 2 * eps * S * Cm; order = len(pairs)
    # log det(1 - M) = - sum_k tr(M^k)/k ; series to order `order` in eps
    Mk = sp.eye(n); logdet = sp.Integer(0)
    for k in range(1, order + 1):
        Mk = Mk * M; logdet -= Mk.trace() / k
    expo = -sp.Rational(1, 2) * N * logdet
    if 0 in idx:
        # Q^2 [ 2S (1 - 2 C S)^{-1} ]_{00} = Q^2 sum_k [2S (2CS)^k]_{00}
        A = 2 * eps * Cm * S; Ak = sp.eye(n); chain = sp.Integer(0)
        for k in range(0, order + 1):
            chain += (eps * S * Ak)[idx[0], idx[0]]; Ak = Ak * A      # (1/2) mu^T (2S)(1 - C 2S)^{-1} mu = Q^2 [S (1 - 2CS)^{-1}]_{00}
        expo += Q**2 * chain
    ser = sp.series(sp.exp(expo), eps, 0, order + 1).removeO()
    coeff = sp.expand(ser.coeff(eps, order))
    for sk, kind in zip(s, kinds):
        m = mult[kind]; coeff = sp.expand(sp.factorial(m) * sp.Poly(coeff, sk).coeff_monomial(sk**m))      # (Y_a.Y_b)^m = m! x coefficient of s^m in exp(s Y_a.Y_b)
    return sp.expand(coeff)
def vev_monomial(mono, N):
    xs, yp = mono
    xo = tuple(sorted(k - 1 for k in xs)); yo = tuple((k - 1, l - 1) for k, l in yp)
    return sp.expand(x_value(xo) * y_value(yo, N))
def vev_polynomial(dens, N):
    """eigenvalue in the ARTICLE momenta (P_art = i P, Q_art = i Q), as cyl_vev.vev_polynomial"""
    tot = sp.Integer(0)
    for mono, c in dens.items():
        cc = sp.sympify(c) if isinstance(c, sp.Basic) else sp.Rational(c.numerator, c.denominator)
        tot += cc * vev_monomial(mono, N)
    Pa, Qa = sp.symbols('P_art Q_art')
    out = sp.expand(tot.subs({P: sp.I * Pa, Q: sp.I * Qa}, simultaneous=True))
    return sp.expand(out.subs({Pa: P, Qa: Q}, simultaneous=True))
if __name__ == "__main__":
    from fractions import Fraction as F
    from mixed_engine import gen_mixed_basis, genuine_kernel_mixed, cyl_P4
    from cyl_vev import vev_polynomial as recipe
    import time
    def point(t):
        t = F(t); return (t * t - 25) / (t * t - 1), F(-24) * t / (t * t - 1)
    ok = True
    for sol in (1, 2, 3):
        tv = F(2); Nv, sv = point(tv); P4 = cyl_P4(sol, Nv, sv)
        for W in (4, 6, 8, 10):
            dens = P4 if W == 4 else genuine_kernel_mixed(P4, gen_mixed_basis(W), Nv)[0]
            t0 = time.time(); a = vev_polynomial(dens, sp.Rational(Nv.numerator, Nv.denominator)); t1 = time.time()
            b = recipe(dens, Nv, W); t2 = time.time()
            same = sp.expand(a - b) == 0; ok &= same
            print(f"Solution {sol} t = {tv} weight {W}: generating-function == recipe: {same}  ({t1-t0:.1f}s vs {t2-t1:.1f}s)", flush=True)
    print("VALIDATION", "PASSED" if ok else "FAILED")
