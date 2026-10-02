"""Exact Fraction linear solve: returns a solution x of A x = b or None."""
from fractions import Fraction as F
def nullspace_solve(A, b):
    n = len(A[0]) if A else 0
    M = [list(r) + [bb] for r, bb in zip(A, b)]
    piv = []; r = 0
    for c in range(n):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]; pv = M[r][c]; M[r] = [x / pv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]; M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        piv.append(c); r += 1
    for i in range(r, len(M)):
        if M[i][n] != 0: return None
    x = [F(0)] * n
    for i, c in enumerate(piv): x[c] = M[i][n]
    return x
