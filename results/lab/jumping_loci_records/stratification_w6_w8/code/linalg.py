"""stratification_w6_w8: exact linear algebra helpers (my own; Fractions / integers only)."""
from fractions import Fraction as F


def rank_and_pivots(A):
    """Exact rank of a matrix of Fractions/ints by Gaussian elimination; returns (rank, pivot_rows, pivot_cols)
    where pivot_rows are indices into the ORIGINAL row order (rows are scanned in order, so the pivot rows
    are linearly independent and the pivot_rows x pivot_cols submatrix is nonsingular)."""
    rows = [[F(x) for x in r] for r in A]
    nr = len(rows)
    nc = len(rows[0]) if nr else 0
    orig = list(range(nr))
    rank = 0
    pcols = []
    for c in range(nc):
        piv = next((r for r in range(rank, nr) if rows[r][c] != 0), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        orig[rank], orig[piv] = orig[piv], orig[rank]
        pr = rows[rank]
        inv = 1 / pr[c]
        for r in range(rank + 1, nr):
            if rows[r][c] != 0:
                f = rows[r][c] * inv
                rr = rows[r]
                for cc in range(c, nc):
                    if pr[cc] != 0:
                        rr[cc] -= f * pr[cc]
        pcols.append(c)
        rank += 1
        if rank == nr:
            break
    return rank, sorted(orig[:rank]), pcols


def rank(A):
    return rank_and_pivots(A)[0]


def det(A):
    """Exact determinant (Fractions)."""
    n = len(A)
    rows = [[F(x) for x in r] for r in A]
    d = F(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if rows[r][c] != 0), None)
        if piv is None:
            return F(0)
        if piv != c:
            rows[c], rows[piv] = rows[piv], rows[c]
            d = -d
        pr = rows[c]
        d *= pr[c]
        inv = 1 / pr[c]
        for r in range(c + 1, n):
            if rows[r][c] != 0:
                f = rows[r][c] * inv
                rr = rows[r]
                for cc in range(c, n):
                    if pr[cc] != 0:
                        rr[cc] -= f * pr[cc]
    return d


def solve(A, b):
    """Solve A x = b exactly for nonsingular square A; returns (x, det A)."""
    n = len(A)
    rows = [[F(x) for x in r] + [F(bb)] for r, bb in zip(A, b)]
    d = F(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if rows[r][c] != 0), None)
        if piv is None:
            raise ZeroDivisionError("singular")
        if piv != c:
            rows[c], rows[piv] = rows[piv], rows[c]
            d = -d
        pr = rows[c]
        d *= pr[c]
        inv = 1 / pr[c]
        for r in range(c + 1, n):
            if rows[r][c] != 0:
                f = rows[r][c] * inv
                rr = rows[r]
                for cc in range(c, n + 1):
                    if pr[cc] != 0:
                        rr[cc] -= f * pr[cc]
    x = [F(0)] * n
    for c in range(n - 1, -1, -1):
        s = rows[c][n] - sum(rows[c][k] * x[k] for k in range(c + 1, n))
        x[c] = s / rows[c][c]
    return x, d


def interpolate(xs, ys):
    """Exact polynomial (coefficient list low -> high, Fractions) of degree < len(xs) through the points (Newton form)."""
    n = len(xs)
    coef = [F(y) for y in ys]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) / (xs[i] - xs[i - j])
    # expand Newton form
    poly = [F(0)]
    for i in range(n - 1, -1, -1):
        # poly = poly * (x - xs[i]) + coef[i]
        new = [F(0)] * (len(poly) + 1)
        for k, c in enumerate(poly):
            new[k + 1] += c
            new[k] -= c * xs[i]
        new[0] += coef[i]
        poly = new
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly
