"""VIR2-R -- O(N) pair-language Wick engine with FORMAL N (numeric rational N), built on wick.py (my WINF1/VIR1-R explicit Wick code).
Field = {(xs, yp): Fraction}: xs = sorted tuple of X derivative orders (>= 1), yp = sorted tuple of Y pairs (a, b), a <= b, (a, b) = d^aY . d^bY.
Propagator d^m phi(z) d^n phi(w) ~ (-1)^(m-1) (m+n-1)! (z-w)^-(m+n) for X-X and for Y_i-Y_j (delta_ij).  Y letters carry dummy labels (one per pair);
a Y-Y contraction unifies labels; after contraction each label class keeps 2 letters (-> new pair) or 0 (closed cycle -> factor N)."""
from fractions import Fraction as F
from math import factorial
def _letters(key, side):
    xs, yp = key; L = [('x', a, None) for a in xs]
    for k, (a, b) in enumerate(yp): L += [('y', a, (side, k)), ('y', b, (side, k))]
    return L
def _matchings(A, B):
    def rec(ia, used, pairs):
        if ia == len(A): yield pairs; return
        yield from rec(ia + 1, used, pairs)
        for jb in range(len(B)):
            if jb not in used and B[jb][0] == A[ia][0]:
                yield from rec(ia + 1, used | {jb}, pairs + [(ia, jb)])
    yield from rec(0, frozenset(), [])
def _taylor(n, T):
    if n == 0:
        if T == 0: yield ()
        return
    for k in range(T + 1):
        for rest in _taylor(n - 1, T - k): yield (k,) + rest
def ope_pairs(A, B, Nval, poles):
    out = {}
    for kA, cA in A.items():
        LA = _letters(kA, 'A')
        for kB, cB in B.items():
            LB = _letters(kB, 'B')
            for pairs in _matchings(LA, LB):
                fac = F(cA) * F(cB); P = 0
                parent = {}
                def find(x):
                    while parent.get(x, x) != x: x = parent[x]
                    return x
                for ia, jb in pairs:
                    m, n = LA[ia][1], LB[jb][1]
                    fac *= (-1) ** (m - 1) * factorial(m + n - 1); P += m + n
                    if LA[ia][0] == 'y':
                        ra, rb = find(LA[ia][2]), find(LB[jb][2])
                        if ra != rb: parent[ra] = rb
                usedA = {ia for ia, _ in pairs}; usedB = {jb for _, jb in pairs}
                remA = [LA[i] for i in range(len(LA)) if i not in usedA]; remB = [LB[j] for j in range(len(LB)) if j not in usedB]
                # closed cycles: label classes that lost all their letters
                labels = set(l for (_, _, l) in LA + LB if l is not None)
                roots_all = {find(l) for l in labels}
                for T in range(P + 1):
                    r = P - T
                    if r not in poles: continue
                    for sh in _taylor(len(remA), T):
                        cf = fac
                        for k in sh: cf /= factorial(k)
                        letters = [(t_, d + k, l) for (t_, d, l), k in zip(remA, sh)] + remB
                        xs = tuple(sorted(d for (t_, d, l) in letters if t_ == 'x'))
                        groups = {}
                        for (t_, d, l) in letters:
                            if t_ == 'y': groups.setdefault(find(l), []).append(d)
                        ok = True; yp = []
                        for g in groups.values():
                            if len(g) != 2: ok = False; break
                            yp.append(tuple(sorted(g)))
                        assert ok, 'label class with != 2 letters'
                        ncyc = len(roots_all - set(groups.keys()))
                        key = (xs, tuple(sorted(yp)))
                        d_ = out.setdefault(r, {}); d_[key] = d_.get(key, 0) + cf * F(Nval) ** ncyc
    return {r: {k: v for k, v in f.items() if v} for r, f in out.items()}
def dfield_pairs(A):
    out = {}
    for (xs, yp), c in A.items():
        for i in range(len(xs)):
            k = (tuple(sorted(xs[:i] + (xs[i] + 1,) + xs[i + 1:])), yp); out[k] = out.get(k, 0) + c
        for j, (a, b) in enumerate(yp):
            for na, nb in ((a + 1, b), (a, b + 1)):
                k = (xs, tuple(sorted(yp[:j] + (tuple(sorted((na, nb))),) + yp[j + 1:]))); out[k] = out.get(k, 0) + c
    return {k: v for k, v in out.items() if v}
def mulx(xs_extra, A):
    return {(tuple(sorted(xs + tuple(xs_extra))), yp): c for (xs, yp), c in A.items()}
def fadd(A, B, s=1):
    out = dict(A)
    for k, v in B.items():
        out[k] = out.get(k, 0) + s * v
        if out[k] == 0: del out[k]
    return out
def to_explicit(A, Nc):
    """pair field -> explicit letters (comp, p) with p = derivative order of a = d(boson) (i.e. order - 1); X comp -1"""
    out = {}
    for (xs, yp), c in A.items():
        poly = {tuple(sorted((-1, a - 1) for a in xs)): F(c)}
        for (a, b) in yp:
            new = {}
            for m, cc in poly.items():
                for i in range(Nc):
                    nm = tuple(sorted(m + ((i, a - 1), (i, b - 1)), key=str)); new[nm] = new.get(nm, 0) + cc
            poly = new
        for m, cc in poly.items():
            m2 = tuple(sorted(m, key=str)); out[m2] = out.get(m2, 0) + cc
    return {m: v for m, v in out.items() if v}
