"""Generating-function VEV in exact rational arithmetic with per-variable truncation (2026-09-06).

The same quantity as lab/vev_genfun.py (which was validated exactly against the certified Wick-pairing recipe
code/cyl_vev.py on the twelve certified weight-4..10 densities), reorganised so that memory stays bounded and speed
is polynomial in practice:

  * X-sector (letters d^k J, contraction C(a, b), zero mode P for k = 0): the memoized recursion over the sorted multiset
    of orders, with the value held as a tuple of Fractions (coefficients of P^d) instead of a sympy expression.
  * Y-sector (pair letters (d^a Y . d^b Y) of the O(N) vector, zero mode Q e for a = 0):
        <exp(sum_j s_j Y_{a_j}.Y_{b_j})> = det(1 - 2 S C)^{-N/2} exp( Q^2 [ S (1 - 2 C S)^{-1} ]_{00} ),
    the coefficient of prod_j s_j^{m_j} (times prod m_j!) extracted from truncated multivariate polynomials in the
    kind variables s_j, every variable capped at its multiplicity m_j (so the total degree is capped at the number of
    pairs), the two exponentials expanded separately: exp((N/2) sum_k tr(M^k)/k) with M = 2 S C, and
    exp(Q^2 Ch) = sum_j Q^{2j} Ch^j / j! with Ch = sum_k [S (2 C S)^k]_{00}.  N numeric (Fraction); nothing symbolic.
  * Mixed monomials: product of the two sectors' values (contractions pair only same-species letters).
Article momenta at the end: P -> i P, Q -> i Q (as cyl_vev.vev_polynomial / vev_genfun.vev_polynomial).
Output: a sympy polynomial in P, Q identical to vev_genfun.vev_polynomial(dens, N)."""
import sys, os
from fractions import Fraction as F
from functools import lru_cache
from collections import Counter
from math import factorial
import sympy as sp
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'code'))
import cyl_vev
# contraction R(eps) = 1/(4 sinh^2(eps/2)) - 1/eps^2 and C(a, b) = (-1)^b R^{(a+b)}(0), as cyl_vev; the derivative table is
# extended beyond cyl_vev's 22 entries (weight 18 needs formal diagonal entries up to order 32) and checked on the overlap
_eps = sp.Symbol('eps'); _NR = 40
_R = sp.series(1 / (4 * sp.sinh(_eps / 2) ** 2) - 1 / _eps ** 2, _eps, 0, _NR + 2).removeO()
_RD = [F(str(sp.diff(_R, _eps, k).subs(_eps, 0))) for k in range(_NR)]
if _RD[:len(cyl_vev._Rd)] != [F(v) for v in cyl_vev._Rd]:
    raise RuntimeError("contraction table disagrees with cyl_vev on the overlap")
@lru_cache(maxsize=None)
def Cc(a, b):
    return (-1) ** b * _RD[a + b]

# ---- X-sector: univariate polynomials in P as tuples of Fractions ----
def padd(a, b):
    n = max(len(a), len(b))
    return tuple((a[i] if i < len(a) else F(0)) + (b[i] if i < len(b) else F(0)) for i in range(n))
@lru_cache(maxsize=None)
def x_value(orders):
    if not orders: return (F(1),)
    first, rest = orders[0], orders[1:]
    total = ((F(0),) + x_value(rest)) if first == 0 else (F(0),)
    for i, o in enumerate(rest):
        c = Cc(first, o)
        if c: total = padd(total, tuple(c * z for z in x_value(rest[:i] + rest[i + 1:])))
    return total

# ---- Y-sector: truncated multivariate polynomials in the kind variables (dict exponent-tuple -> Fraction) ----
def tmul(a, b, caps):
    out = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            e = tuple(x + y for x, y in zip(ea, eb))
            if all(x <= c for x, c in zip(e, caps)):
                out[e] = out.get(e, F(0)) + ca * cb
    return {e: c for e, c in out.items() if c}
def tadd(a, b, scale=F(1)):
    out = dict(a)
    for e, c in b.items():
        out[e] = out.get(e, F(0)) + scale * c
    return {e: c for e, c in out.items() if c}
def tscale(a, s): return {e: s * c for e, c in a.items()} if s else {}
def matmul(A, B, caps):
    n = len(A); out = [[{} for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for l in range(n):
            if not A[i][l]: continue
            for j in range(n):
                if B[l][j]: out[i][j] = tadd(out[i][j], tmul(A[i][l], B[l][j], caps))
    return out
def texp(a, caps, order):
    """exp(a) truncated, a without constant term"""
    k = len(caps); one = {(0,) * k: F(1)}
    out, term = dict(one), dict(one)
    for i in range(1, order + 1):
        term = tscale(tmul(term, a, caps), F(1, i))
        if not term: break
        out = tadd(out, term)
    return out
def y_value(pairs, N):
    """returns {j: coefficient of Q^(2j)} (Fractions) for the O(N) vacuum value of prod (Y_a.Y_b) over `pairs`"""
    if not pairs: return {0: F(1)}
    orders = sorted({o for pr in pairs for o in pr}); idx = {o: i for i, o in enumerate(orders)}; n = len(orders)
    mult = Counter(tuple(sorted(pr)) for pr in pairs); kinds = sorted(mult); k = len(kinds)
    caps = tuple(mult[kd] for kd in kinds); order = len(pairs)
    def unit(j): return {tuple(1 if i == j else 0 for i in range(k)): F(1)}
    S = [[{} for _ in range(n)] for _ in range(n)]
    for j, (a, b) in enumerate(kinds):
        if a == b: S[idx[a]][idx[a]] = tadd(S[idx[a]][idx[a]], unit(j))
        else:
            S[idx[a]][idx[b]] = tadd(S[idx[a]][idx[b]], unit(j), F(1, 2)); S[idx[b]][idx[a]] = tadd(S[idx[b]][idx[a]], unit(j), F(1, 2))
    Cm = [[Cc(orders[i], orders[j]) for j in range(n)] for i in range(n)]
    M = [[{} for _ in range(n)] for _ in range(n)]          # M = 2 S C
    A = [[{} for _ in range(n)] for _ in range(n)]          # A = 2 C S
    for i in range(n):
        for j in range(n):
            for l in range(n):
                if S[i][l] and Cm[l][j]: M[i][j] = tadd(M[i][j], S[i][l], 2 * Cm[l][j])
                if Cm[i][l] and S[l][j]: A[i][j] = tadd(A[i][j], S[l][j], 2 * Cm[i][l])
    I = [[({(0,) * k: F(1)} if i == j else {}) for j in range(n)] for i in range(n)]
    Mk = I; L = {}
    for kk in range(1, order + 1):
        Mk = matmul(Mk, M, caps)
        tr = {}
        for i in range(n): tr = tadd(tr, Mk[i][i])
        if tr: L = tadd(L, tr, F(1, kk))          # odd traces can vanish while higher ones do not: no early exit
    EL = texp(tscale(L, F(N) / 2), caps, order)
    target = caps; mfac = F(1)
    for c in caps: mfac *= factorial(c)
    out = {}
    if 0 in idx:
        i0 = idx[0]; Ak = I; Ch = {}
        for kk in range(order + 1):
            Ch = tadd(Ch, matmul(S, Ak, caps)[i0][i0]); Ak = matmul(Ak, A, caps)
        Chj = {(0,) * k: F(1)}
        for j in range(order + 1):
            if j: Chj = tscale(tmul(Chj, Ch, caps), F(1, j))
            if not Chj: break
            c = tmul(EL, Chj, caps).get(target, F(0))
            if c: out[j] = mfac * c
    else:
        c = EL.get(target, F(0))
        if c: out[0] = mfac * c
    return out
def vev_monomial(mono, N):
    """{(d, j): Fraction} = coefficient of P^d Q^(2j) (bare momenta)"""
    xs, yp = mono
    xo = tuple(sorted(k - 1 for k in xs)); yo = tuple((k - 1, l - 1) for k, l in yp)
    xp = x_value(xo); yv = y_value(yo, N); out = {}
    for d, cx in enumerate(xp):
        if not cx: continue
        for j, cy in yv.items(): out[(d, j)] = out.get((d, j), F(0)) + cx * cy
    return out
def vev_table(dens, N):
    tot = {}
    for mono, c in dens.items():
        cc = c if isinstance(c, F) else F(int(sp.Rational(c).p), int(sp.Rational(c).q))
        if not cc: continue
        for key, v in vev_monomial(mono, N).items(): tot[key] = tot.get(key, F(0)) + cc * v
    return {key: v for key, v in tot.items() if v}
def vev_polynomial(dens, N):
    """eigenvalue in the ARTICLE momenta (P_art = i P, Q_art = i Q), as vev_genfun.vev_polynomial"""
    P, Q = sp.symbols('P Q'); tot = sp.Integer(0)
    for (d, j), v in vev_table(dens, N).items():
        tot += sp.Rational(v.numerator, v.denominator) * sp.I**d * (-1)**j * P**d * Q**(2 * j)
    return sp.expand(tot)
