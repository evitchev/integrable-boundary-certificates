"""stratification_w6_w8: exact rank of an integer-polynomial matrix over the number field Q[t]/(f), f irreducible over Q.
Dense univariate arithmetic over Fractions (my own code). Elements: coefficient lists low -> high, reduced mod f."""
from fractions import Fraction as F


def trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def pmul(a, b):
    if not a or not b:
        return []
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            if y != 0:
                out[i + j] += x * y
    return trim(out)


def psub(a, b):
    n = max(len(a), len(b))
    return trim([(a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0) for i in range(n)])


def pdivmod(a, f):
    a = trim([F(x) for x in a])
    f = trim([F(x) for x in f])
    q = [F(0)] * max(len(a) - len(f) + 1, 0)
    lc = f[-1]
    while len(a) >= len(f):
        c = a[-1] / lc
        sh = len(a) - len(f)
        q[sh] = c
        for i, y in enumerate(f):
            a[sh + i] -= c * y
        a = trim(a)
    return trim(q), a


def pmod(a, f):
    return pdivmod(a, f)[1]


def pinv(a, f):
    """Inverse of a mod f (f irreducible, a not 0 mod f), by the extended Euclidean algorithm over Q."""
    r0, r1 = trim([F(x) for x in f]), pmod(a, f)
    s0, s1 = [], [F(1)]
    while r1:
        q, r = pdivmod(r0, r1)
        r0, r1 = r1, r
        s0, s1 = s1, psub(s0, pmul(q, s1))
    # r0 is a nonzero constant (gcd), s0 * a = r0 mod f
    assert len(r0) == 1, "not invertible: f is not irreducible or a = 0 mod f"
    return pmod([c / r0[0] for c in s0], f)


def rank_mod(Mpoly, f):
    """Exact rank of Mpoly (integer coefficient lists) over Q[t]/(f)."""
    f = trim([F(x) for x in f])
    rows = [[pmod([F(c) for c in e], f) for e in row] for row in Mpoly]
    nr, nc = len(rows), len(rows[0])
    rank = 0
    for c in range(nc):
        piv = next((r for r in range(rank, nr) if rows[r][c]), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        pr = rows[rank]
        inv = pinv(pr[c], f)
        pr = [pmod(pmul(x, inv), f) for x in pr]
        rows[rank] = pr
        for r in range(rank + 1, nr):
            if rows[r][c]:
                fac = rows[r][c]
                rows[r] = [psub(x, pmod(pmul(fac, y), f)) if y else x for x, y in zip(rows[r], pr)]
        rank += 1
        if rank == nr:
            break
    return rank
