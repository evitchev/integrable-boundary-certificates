"""Modular linear algebra and exact helpers (certificate cyl_stratification_w1014).

numpy int64 arithmetic modulo primes p < 2^31: every entry is kept in [0, p), so a product of two entries is < 2^62 and is
reduced mod p before any addition (the overflow rule).  Exact helpers (CRT, rational reconstruction, exact sparse
products) use Python integers.
"""
import math
import random
from fractions import Fraction as F

import numpy as np
from sympy import isprime

P31 = 2147483647          # 2^31 - 1 (prime)


def primes_below(bound, k, skip=()):
    """The k largest primes below `bound` (descending), excluding `skip`; sympy.isprime is deterministic below 2^64."""
    out, q = [], bound - 1
    while len(out) < k:
        if q not in skip and isprime(q):
            out.append(q)
        q -= 1
    return out


# ---------------------------------------------------------------- matrices of polynomials mod p

def coeff_arrays(d, p):
    """Sparse M_poly (JSON dict) -> list of int64 arrays C_k = coefficient of t^k mod p, k = 0..dmax."""
    m, n = d["m"], d["n"]
    dmax = max(len(e) - 1 for _, _, e in d["entries"])
    C = [np.zeros((m, n), dtype=np.int64) for _ in range(dmax + 1)]
    for ii, jj, e in d["entries"]:
        for k, c in enumerate(e):
            if c:
                C[k][ii, jj] = c % p
    return C


def tmod(t0, p):
    """A Fraction or int t0 as an element of F_p (asserts the denominator is a unit)."""
    t0 = F(t0)
    assert t0.denominator % p != 0
    return (t0.numerator % p) * pow(t0.denominator % p, p - 2, p) % p


def eval_mod(C, t0, p):
    """M(t0) mod p by Horner, t0 an integer in [0, p)."""
    A = C[-1].copy()
    for k in range(len(C) - 2, -1, -1):
        A = (A * t0) % p
        A = (A + C[k]) % p
    return A


def rref_mod(A, p):
    """Reduced row echelon form of A mod p.  Returns (R, pivot columns, original indices of the pivot rows); the original rows
    at those indices restricted to the pivot columns form a nonsingular submatrix mod p."""
    A = np.array(A, dtype=np.int64) % p
    m, n = A.shape
    perm = np.arange(m)
    piv = []
    r = 0
    for c in range(n):
        if r == m:
            break
        nz = np.nonzero(A[r:, c])[0]
        if nz.size == 0:
            continue
        i = r + int(nz[0])
        if i != r:
            A[[r, i]] = A[[i, r]]
            perm[[r, i]] = perm[[i, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        col = A[:, c].copy()
        col[r] = 0
        nzr = np.nonzero(col)[0]
        if nzr.size:
            A[nzr] = (A[nzr] - (col[nzr, None] * A[r][None, :]) % p) % p
        piv.append(c)
        r += 1
    return A, piv, [int(x) for x in perm[:r]]


def rank_mod(A, p):
    return len(rref_mod(A, p)[1])


def det_mod(A, p):
    """Determinant of a square matrix mod p (Gaussian elimination)."""
    A = np.array(A, dtype=np.int64) % p
    n = A.shape[0]
    det = 1
    for c in range(n):
        nz = np.nonzero(A[c:, c])[0]
        if nz.size == 0:
            return 0
        i = c + int(nz[0])
        if i != c:
            A[[c, i]] = A[[i, c]]
            det = (-det) % p
        pv = int(A[c, c])
        det = det * pv % p
        inv = pow(pv, p - 2, p)
        if c + 1 < n:
            f = (A[c + 1:, c] * inv) % p
            A[c + 1:, c:] = (A[c + 1:, c:] - (f[:, None] * A[c, c:][None, :]) % p) % p
    return det


def nullspace_mod(A, p):
    """Null space of A mod p in RREF normalization: (pivot cols, free cols, list of vectors (int arrays))."""
    R, piv, _ = rref_mod(A, p)
    n = A.shape[1]
    pset = set(piv)
    free = [c for c in range(n) if c not in pset]
    vecs = []
    for f in free:
        v = np.zeros(n, dtype=np.int64)
        v[f] = 1
        for i, c in enumerate(piv):
            v[c] = (-int(R[i, f])) % p
        vecs.append(v)
    return piv, free, vecs


# ---------------------------------------------------------------- univariate polynomials mod p (lists, low -> high)

def ptrim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def pmod_eval(a, x, p):
    acc = 0
    for c in reversed(a):
        acc = (acc * x + c) % p
    return acc


def pmulp(a, b, p):
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] = (out[i + j] + x * y) % p
    return ptrim(out)


def psubp(a, b, p):
    n = max(len(a), len(b))
    return ptrim([((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % p for i in range(n)])


def pdivmodp(a, b, p):
    a = [x % p for x in a]
    b = ptrim([x % p for x in b])
    assert b, "division by zero polynomial"
    inv = pow(b[-1], p - 2, p)
    q = [0] * max(len(a) - len(b) + 1, 0)
    for i in range(len(a) - len(b), -1, -1):
        c = a[i + len(b) - 1] * inv % p
        q[i] = c
        if c:
            for j, y in enumerate(b):
                a[i + j] = (a[i + j] - c * y) % p
    return ptrim(q), ptrim(a)


def pmonic(a, p):
    a = ptrim([x % p for x in a])
    inv = pow(a[-1], p - 2, p)
    return [x * inv % p for x in a]


def pgcdp(a, b, p):
    a, b = ptrim([x % p for x in a]), ptrim([x % p for x in b])
    while b:
        a, b = b, pdivmodp(a, b, p)[1]
    return pmonic(a, p) if a else []


def interp_mod(xs, ys, p):
    """Newton interpolation mod p through the points (xs distinct mod p); coefficients low -> high."""
    n = len(xs)
    c = [y % p for y in ys]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            c[i] = (c[i] - c[i - 1]) * pow((xs[i] - xs[i - j]) % p, p - 2, p) % p
    poly = []
    for i in range(n - 1, -1, -1):
        poly = pmulp(poly, [(-xs[i]) % p, 1], p) if poly else []
        poly = ptrim([(poly[0] + c[i]) % p] + poly[1:]) if poly else ptrim([c[i] % p])
    return poly


def ratrecon_poly(xs, ys, p, D):
    """Rational function a/q with deg a <= D, deg q <= D, q monic, a(x) = y q(x) at all points (needs len(xs) >= 2D + 1):
    extended Euclid on (prod (t - x_i), interpolant), stopped at the first remainder of degree <= D.  None if it fails."""
    Z = interp_mod(xs, ys, p)
    Mo = [1]
    for x in xs:
        Mo = pmulp(Mo, [(-x) % p, 1], p)
    r0, r1 = Mo, Z
    s0, s1 = [], [1]          # cofactors of Z
    while r1 and len(r1) - 1 > D:
        q, r = pdivmodp(r0, r1, p)
        r0, r1 = r1, r
        s0, s1 = s1, psubp(s0, pmulp(q, s1, p), p)
    if not s1 or len(s1) - 1 > D:
        return None
    inv = pow(s1[-1], p - 2, p)
    a = [x * inv % p for x in r1]
    q = [x * inv % p for x in s1]
    # a = Z*q mod Mo; it must agree with y*q at every point (true by construction) and q must not vanish there
    if any(pmod_eval(q, x, p) == 0 for x in xs):
        return None
    return a, q


# ---------------------------------------------------------------- CRT and rational reconstruction

def crt_pair(r1, m1, r2, m2):
    """x = r1 mod m1, x = r2 mod m2 (coprime) -> x mod m1*m2."""
    k = (r2 - r1) * pow(m1, -1, m2) % m2
    return r1 + m1 * k, m1 * m2


def symmetric(x, m):
    x %= m
    return x - m if x > m // 2 else x


def ratrecon(u, m):
    """Wang's rational reconstruction: a/b with |a|, b <= sqrt(m/2), a = b u mod m; None if none exists."""
    u %= m
    bound = math.isqrt(m // 2)
    r0, r1 = m, u
    s0, s1 = 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    if s1 < 0:
        r1, s1 = -r1, -s1
    if math.gcd(r1, s1) != 1:
        return None
    return F(r1, s1)


# ---------------------------------------------------------------- exact evaluation and products

def exact_int_matrix(d, t0):
    """Sparse integer matrix b^dmax * M_poly(a/b) as {(ii, jj): int} for t0 = a/b (a positive multiple of M_poly(t0))."""
    t0 = F(t0)
    a, b = t0.numerator, t0.denominator
    dmax = max(len(e) - 1 for _, _, e in d["entries"])
    pw_a = [a**k for k in range(dmax + 1)]
    pw_b = [b**k for k in range(dmax + 1)]
    out = {}
    for ii, jj, e in d["entries"]:
        v = sum(c * pw_a[k] * pw_b[dmax - k] for k, c in enumerate(e))
        if v:
            out[(ii, jj)] = v
    return out


def int_vector(v):
    """A list of Fractions -> primitive integer vector (common denominator cleared, content divided out)."""
    den = 1
    for x in v:
        den = math.lcm(den, F(x).denominator)
    iv = [int(F(x) * den) for x in v]
    g = 0
    for x in iv:
        g = math.gcd(g, x)
    return [x // g for x in iv] if g else iv


def sparse_matvec_is_zero(Mint, m, v):
    acc = [0] * m
    for (ii, jj), c in Mint.items():
        x = v[jj]
        if x:
            acc[ii] += c * x
    return not any(acc)


def exact_nullity(d, t0, rng, nprimes_rank=2, max_primes=4000, verbose=False):
    """Exact nullity of M_poly(t0), t0 rational (t0 != +-1), by my own computation:
      rank lower bound: rank of M_poly(t0) mod two 31-bit primes (a nonzero minor mod p is nonzero over Q);
      nullity lower bound: exact rational null vectors (RREF normalization, CRT over 31-bit primes + rational reconstruction),
      verified by exact integer products.  Returns a dict with the proof data."""
    n, m = d["n"], d["m"]
    ps = primes_below(P31 + 1, nprimes_rank + max_primes)
    rank_ps = ps[:nprimes_rank]
    ranks = []
    for p in rank_ps:
        C = coeff_arrays(d, p)
        ranks.append(rank_mod(eval_mod(C, tmod(t0, p), p), p))
    rho = max(ranks)
    nul_upper = n - rho
    # reconstruction primes: reduce the exact integer matrix b^dmax M_poly(t0) mod p (same kernel as M_poly(t0))
    Mint = exact_int_matrix(d, t0)
    acc = None
    used = 0
    piv0 = None
    vecs_exact = None
    for p in ps[nprimes_rank:]:
        A = np.zeros((m, n), dtype=np.int64)
        for (ii, jj), c in Mint.items():
            A[ii, jj] = c % p
        piv, free, vecs = nullspace_mod(A, p)
        if len(free) != nul_upper:
            continue                       # this prime sees a smaller rank: skip it
        if piv0 is None:
            piv0, free0 = piv, free
            acc = [[(int(x), p) for x in v] for v in vecs]
        else:
            if piv != piv0:
                continue
            acc = [[crt_pair(r, mm, int(x), p) for (r, mm), x in zip(av, v)] for av, v in zip(acc, vecs)]
        used += 1
        if used % 4 == 0 or used == 1:
            cand = []
            ok = True
            for av in acc:
                vr = [ratrecon(r, mm) for r, mm in av]
                if any(x is None for x in vr):
                    ok = False
                    break
                cand.append(vr)
            if ok:
                ivs = [int_vector(vr) for vr in cand]
                if all(sparse_matvec_is_zero(Mint, m, iv) for iv in ivs):
                    vecs_exact = ivs
                    break
        if used >= max_primes:
            break
    res = {"t0": str(F(t0)), "rank_primes": rank_ps, "ranks_mod_p": ranks, "nullity_upper": nul_upper,
           "recon_primes_used": used}
    if vecs_exact is not None:
        # independence: in RREF normalization each vector is nonzero exactly at its own free column among free0
        indep = all((iv[f] != 0) == (k == kk) for k, iv in enumerate(vecs_exact) for kk, f in enumerate(free0))
        res.update({"null_vectors_verified": len(vecs_exact), "independent": indep, "free_cols": free0,
                    "max_bits": max(abs(x).bit_length() for iv in vecs_exact for x in iv),
                    "nullity": nul_upper if indep else None})
        res["vectors"] = vecs_exact
    else:
        res.update({"null_vectors_verified": 0, "nullity": None})
    return res
