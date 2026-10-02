"""SUPER3 Step B engine: WKB local IMs of the EXTENDED Suzuki ODE  psi'' = [kappa^2 (x^(2M) - 1) + kappa alpha x^(M-1) + lam x^-2] psi
(Suzuki quant-ph/0003066 eq. (1) rescaled; the lam term is our extension).  Riccati R' + R^2 = Q, R = sum_n kappa^(1-n) R_n.
Monomial x^e Lambda^f (Lambda = x^(2M) - 1); coefficient = polynomial in (alpha, lam): {(i, j): Fraction}.
int x^e Lambda^f dx := (-1)^(f-1/2) B((e+1)/(2M), f+1) / (2M)   (continued Beta on (0, 1), the record's convention)."""
import sys; sys.dont_write_bytecode = True
from fractions import Fraction as F
from collections import defaultdict
def cadd(A, B, c=1):
    for k, v in B.items():
        A[k] = A.get(k, 0) + c * v
        if A[k] == 0: del A[k]
    return A
def cmul(A, B):
    o = {}
    for (i, j), u in A.items():
        for (k, l), v in B.items():
            o[(i + k, j + l)] = o.get((i + k, j + l), 0) + u * v
    return {k: v for k, v in o.items() if v}
def fadd(P, Q, c=1):
    for m, cf in Q.items():
        cur = cadd(dict(P.get(m, {})), cf, c)
        if cur: P[m] = cur
        elif m in P: del P[m]
    return P
def fmul(P, Q):
    o = {}
    for (e1, f1), c1 in P.items():
        for (e2, f2), c2 in Q.items():
            fadd(o, {(e1 + e2, f1 + f2): cmul(c1, c2)})
    return o
def fder(P, M):
    o = {}
    for (e, f), c in P.items():
        if e + 2 * M * f: fadd(o, {(e - 1, f): c}, e + 2 * M * f)
        if f: fadd(o, {(e - 1, f - 1): c}, 2 * M * f)
    return o
def riccati(M, nmax):
    half = F(1, 2)
    R = [{(F(0), half): {(0, 0): F(1)}}]
    inv = {(F(0), -half): {(0, 0): half}}
    S1 = {(M - 1, F(0)): {(1, 0): F(1)}}; fadd(S1, fder(R[0], M), -1)     # order kappa^1: R_0' + 2 R_0 R_1 = alpha U
    R.append(fmul(inv, S1))
    for n in range(2, nmax + 1):
        S = {}
        if n == 2: fadd(S, {(F(-2), F(0)): {(0, 1): F(1)}})
        fadd(S, fder(R[n - 1], M), -1)
        for i in range(1, n):
            fadd(S, fmul(R[i], R[n - i]), -1)
        R.append(fmul(inv, S))
    return R
def rfa(x, n):
    r = F(1)
    if n >= 0:
        for j in range(n): r *= x + j
    else:
        for j in range(-n):
            d = x - 1 - j
            if d == 0: raise ZeroDivisionError('degenerate Pochhammer')
            r /= d
    return r
def charge(P, M):
    """sum over monomials of the continued Beta, grouped by the class of b mod 1; each class relative to its first b0, f0.
    Returns {class: (b0, f0, poly)}: value = Gamma(b0)Gamma(f0+1)/Gamma(b0+f0+1)/(2M) * poly."""
    out = {}
    for (e, f), c in P.items():
        b = (e + 1) / (2 * M); cls = (b - (b.numerator // b.denominator), f - (f.numerator // f.denominator))
        if cls[1] == 0: out.setdefault(('INTEGER-f', cls[0]), [None, None, {}]); fadd(out[('INTEGER-f', cls[0])][2], {(e, f): c}); continue
        if cls not in out: out[cls] = [b, f, {}]
        b0, f0, poly = out[cls]
        db, df = b - b0, f - f0
        assert db.denominator == 1 and df.denominator == 1
        db, df = int(db), int(df)
        # Gamma(b)/Gamma(b0) * Gamma(f+1)/Gamma(f0+1) * Gamma(b0+f0+1)/Gamma(b+f+1)
        ratio = rfa(b0, db) * rfa(f0 + 1, df) / rfa(b0 + f0 + 1, db + df) * (1 if df % 2 == 0 else -1)
        cadd(poly, c, ratio)
    return out

def is_total_derivative(P, M):
    """Exact: is P (a sum of monomials x^e Lambda^f) = d/dx G for G in the span of monomials x^(e+1) Lambda^f', f' in {f, f+1}?"""
    from eng_null import nullspace_solve
    targets = {}
    for (e, f), c in P.items():
        for k, v in c.items(): targets.setdefault(k, {})[(e, f)] = v
    for k, T in targets.items():
        cands = sorted({(e + 1, f + d) for (e, f) in T for d in (0, 1, 2, 3)})
        cols = [fder({m: {(0, 0): F(1)}}, M) for m in cands]
        keys = sorted({m for col in cols for m in col} | set(T))
        A = [[col.get(m, {}).get((0, 0), F(0)) for col in cols] for m in keys]
        b = [T.get(m, F(0)) for m in keys]
        if nullspace_solve(A, b) is None: return False
    return True
