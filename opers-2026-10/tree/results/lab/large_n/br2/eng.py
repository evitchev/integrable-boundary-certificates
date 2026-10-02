"""SUPER2 engine: exact two-boson (X, phi) vertex algebra, dPhi_a dPhi_b ~ delta_ab/(z-w)^2, e^{u Phi}e^{v Phi} ~ (z-w)^{u.v}.
Field = polynomial in x_(a,k) = d^k Phi_a (a = 0: X, 1: phi), times e^{m Phi}. Monomial = sorted tuple of (a, k)."""
import sys; sys.dont_write_bytecode = True
from fractions import Fraction as F
from math import factorial
from itertools import combinations
from functools import lru_cache
def parts(w, kmax=None):
    if w == 0: yield (); return
    kmax = kmax or w
    for k in range(min(w, kmax), 0, -1):
        for rest in parts(w - k, k): yield (k,) + rest
@lru_cache(None)
def basis(w):
    out = set()
    for wx in range(w + 1):
        for px in parts(wx):
            for pf in parts(w - wx): out.add(tuple(sorted([(0, k) for k in px] + [(1, k) for k in pf])))
    return sorted(out)
def padd(A, B, c=1):
    for m, v in B.items():
        A[m] = A.get(m, 0) + c * v
        if A[m] == 0: del A[m]
    return A
def pmul(A, B):
    out = {}
    for m1, c1 in A.items():
        for m2, c2 in B.items():
            m = tuple(sorted(m1 + m2)); out[m] = out.get(m, 0) + c1 * c2
    return {m: c for m, c in out.items() if c}
def Eseries(s, jmax):
    """E_j, coefficients of zeta^j in exp(sum_n s.x_n zeta^n / n!)."""
    E = [{(): F(1)}]
    for j in range(1, jmax + 1):
        acc = {}
        for n in range(1, j + 1):
            bn = {((a, n),): F(s[a]) / factorial(n) for a in (0, 1) if s[a] != 0}
            padd(acc, pmul(bn, E[j - n]), F(n, j))
        E.append(acc)
    return E
def residue(mono, s, e=0):
    """Coefficient of zeta^{-1} in S(z) A(w)/e^{(s+m)Phi(w)}, S = e^{s Phi}, A = mono e^{m Phi}, e = s.m (integer)."""
    target = -1 - e   # power of zeta needed from P(x + shift) E
    n = len(mono); out = {}
    jmax = max(0, sum(k for _, k in mono) + target) if True else 0
    E = Eseries(s, max(0, sum(k for _, k in mono) + target))
    for r in range(0, n + 1):
        for sel in combinations(range(n), r):
            fac = F(1); D = 0
            for i in sel:
                a, k = mono[i]; fac *= -F(s[a]) * factorial(k - 1); D += k
            j = target + D
            if j < 0 or j >= len(E) or fac == 0: continue
            rest = tuple(mono[i] for i in range(n) if i not in sel)
            padd(out, pmul({rest: fac}, E[j]))
    return out
def kernel(w, screens, sector=(F(0), F(0))):
    """Basis of the weight-w (polynomial weight) fields P e^{sector Phi} killed by every screening residue."""
    B = basis(w); rows = {}
    for si, s in enumerate(screens):
        e = s[0] * sector[0] + s[1] * sector[1]
        assert e.denominator == 1, 'sector not residue-admissible'
        for j, m in enumerate(B):
            for om, c in residue(m, s, int(e)).items(): rows.setdefault((si, om), {})[j] = c
    return nullspace([rows[k] for k in rows], len(B)), B
def nullspace(rows, n):
    R = [dict(r) for r in rows if r]; piv = []
    M = []
    for r in R:
        r = dict(r)
        for (pc, pr) in M:
            if pc in r:
                f = r[pc]
                for k, v in pr.items():
                    r[k] = r.get(k, 0) - f * v
                    if r[k] == 0: del r[k]
        if r:
            pc = min(r); f = r[pc]; r = {k: v / f for k, v in r.items()}
            M2 = []
            for (qc, qr) in M:
                if pc in qr:
                    g = qr[pc]; qr = dict(qr)
                    for k, v in r.items():
                        qr[k] = qr.get(k, 0) - g * v
                        if qr[k] == 0: del qr[k]
                M2.append((qc, qr))
            M = M2 + [(pc, r)]
    pivs = {pc for pc, _ in M}; free = [j for j in range(n) if j not in pivs]; basisv = []
    for fj in free:
        v = [F(0)] * n; v[fj] = F(1)
        for pc, pr in M: v[pc] = -pr.get(fj, 0)
        basisv.append(v)
    return basisv
def deriv(P):
    out = {}
    for m, c in P.items():
        for i, (a, k) in enumerate(m):
            nm = tuple(sorted(m[:i] + ((a, k + 1),) + m[i + 1:])); out[nm] = out.get(nm, 0) + c
    return {m: c for m, c in out.items() if c}
