"""NUM1 WKB engine (own, new): theta-form Riccati for -psi'' + (x^(2M) + alpha x^(M-1) + lam/x^2 - E) psi = 0, E = -e:
S^2 + S' - S - lam = x^(2M+2) + alpha x^(M+1) + e x^2,   S = theta log psi, ' = d/du, u = log x.
Grading: alpha, d/du, the -S term: grade 1; lam: grade 2.  S_0 = -x W^(1/2) (decaying branch), W = x^(2M) + e.
Terms x^a W^b with coefficient polynomials in (alpha, lam) (exact Fractions).  int_0^inf du x^a W^b = (1/2M) e^(b + a/2M) B(a/2M, -b - a/2M) (continued).
Grade k integrates to C_k e^(mu (1-k)), mu = (M+1)/(2M)."""
from fractions import Fraction as F
import mpmath as mp
def padd(p, q, s=1):
    out = dict(p)
    for k, v in q.items():
        out[k] = out.get(k, 0) + s * v
        if out[k] == 0: del out[k]
    return out
def pmul(p, q):
    out = {}
    for k1, a in p.items():
        for k2, b in q.items():
            k = (k1[0] + k2[0], k1[1] + k2[1]); out[k] = out.get(k, 0) + a * b
    return {k: v for k, v in out.items() if v != 0}
def eadd(E1, E2, s=1):
    out = {k: dict(v) for k, v in E1.items()}
    for k, v in E2.items():
        out[k] = padd(out.get(k, {}), v, s)
        if not out[k]: del out[k]
    return out
def emul(E1, E2):
    out = {}
    for (a1, b1), p1 in E1.items():
        for (a2, b2), p2 in E2.items():
            k = (a1 + a2, b1 + b2); out[k] = padd(out.get(k, {}), pmul(p1, p2))
            if not out[k]: del out[k]
    return out
def escale(E1, c): return {k: {kk: vv * c for kk, vv in v.items()} for k, v in E1.items()}
def eD(E1, M):
    out = {}
    for (a, b), p in E1.items():
        if a != 0: out = eadd(out, {(a, b): {k: v * a for k, v in p.items()}})
        if b != 0: out = eadd(out, {(a + 2 * M, b - 1): {k: v * 2 * M * b for k, v in p.items()}})
    return out
def wkb(M, K):
    """S_0..S_K as expressions; polys keyed (i, j) = alpha^i lam^j"""
    M = F(M)
    S = [{(F(1), F(1, 2)): {(0, 0): F(-1)}}]
    inv2S0 = {(F(-1), F(-1, 2)): {(0, 0): F(-1, 2)}}          # 1/(2 S_0)
    for k in range(1, K + 1):
        rhs = {}
        if k == 1: rhs = eadd(rhs, {(M + 1, F(0)): {(1, 0): F(1)}})
        if k == 2: rhs = eadd(rhs, {(F(0), F(0)): {(0, 1): F(1)}})
        rhs = eadd(rhs, eD(S[k - 1], M), -1); rhs = eadd(rhs, S[k - 1])
        for i in range(1, k): rhs = eadd(rhs, emul(S[i], S[k - i]), -1)
        S.append(emul(rhs, inv2S0))
    return S
def integral_coeff(expr, M, alpha, lam):
    """C with int du expr = C e^(mu(1-k)); returns (C, poles) with poles = list of offending terms"""
    M = F(M); tot = mp.mpf(0); poles = []
    for (a, b), p in expr.items():
        A_ = a / (2 * M); B_ = -b - a / (2 * M)
        val = sum(mp.mpf(v.numerator) / v.denominator * alpha ** i * lam ** j for (i, j), v in p.items())
        if val == 0: continue
        if (A_.denominator == 1 and A_ <= 0) or (B_.denominator == 1 and B_ <= 0): poles.append((a, b)); continue
        g = mp.gamma(mp.mpf(A_.numerator) / A_.denominator) * mp.gamma(mp.mpf(B_.numerator) / B_.denominator)
        g /= mp.gamma(mp.mpf((-b).numerator) / (-b).denominator) if not ((-b).denominator == 1 and -b <= 0) else mp.inf
        tot += val * g / (2 * M)
    return tot, poles
def epower(expr, M):
    s = {b + a / (2 * F(M)) for (a, b) in expr}; assert len(s) <= 1, s; return s.pop() if s else None
