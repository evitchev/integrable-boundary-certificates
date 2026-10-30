"""VIR1-R -- free-field Wick OPE, copied VERBATIM from stageWINF1/winf1.py (the section between the markers); letters (i, p) = d^p a_i,
a_i(z)a_j(w) ~ delta_ij (z-w)^-2.  Component -1 is used for X."""
from fractions import Fraction as F
from math import factorial
# ======================= free-field Wick OPE (explicit components) =======================
# letter (i, p) = d^p a_i, weight p + 1; a_i(z) a_j(w) ~ delta_ij (z-w)^-2;  <d^p a(z) d^q a(w)> = (-1)^p (p+q+1)! (z-w)^-(p+q+2)
def fadd(A, B, s=1):
    out = dict(A)
    for m, c in B.items():
        v = out.get(m, 0) + s * c
        if v: out[m] = v
        else: out.pop(m, None)
    return out
def fscale(A, s): return {m: s * c for m, c in A.items() if s * c}
def dfield(A):
    out = {}
    for m, c in A.items():
        for idx in range(len(m)):
            i, p = m[idx]; nm = tuple(sorted(m[:idx] + ((i, p + 1),) + m[idx + 1:]))
            out[nm] = out.get(nm, 0) + c
    return {m: c for m, c in out.items() if c}
def _matchings(A, B):
    """partial matchings of letters of A to letters of B with equal component: yields (pairs, remA, remB)"""
    def rec(ia, usedB, pairs):
        if ia == len(A):
            remA = [A[j] for j in range(len(A)) if j not in {pa for pa, _ in pairs}]
            remB = [B[j] for j in range(len(B)) if j not in usedB]
            yield [(A[pa], B[pb]) for pa, pb in pairs], remA, remB; return
        yield from rec(ia + 1, usedB, pairs)                                      # A[ia] unmatched
        for jb in range(len(B)):
            if jb not in usedB and B[jb][0] == A[ia][0]:
                yield from rec(ia + 1, usedB | {jb}, pairs + [(ia, jb)])
    yield from rec(0, frozenset(), [])
def _taylor(rem, T):
    """distributions of total Taylor order T over the letters rem: yields (shifted letters, coefficient prod 1/n!)"""
    if not rem:
        if T == 0: yield [], F(1)
        return
    (i, p), rest = rem[0], rem[1:]
    for n0 in range(T + 1):
        for sh, cf in _taylor(rest, T - n0):
            yield [(i, p + n0)] + sh, cf / factorial(n0)
def ope(A, B, poles=None):
    """{r: field} for pole orders r >= 0 (r = 0 is the normally ordered product); only the orders in `poles` if given"""
    out = {}
    for mA, cA in A.items():
        for mB, cB in B.items():
            for pairs, remA, remB in _matchings(mA, mB):
                P = 0; fac = F(cA * cB)
                for (i, p), (_, q) in pairs:
                    fac *= (-1)**p * factorial(p + q + 1); P += p + q + 2
                for T in range(P + 1):
                    if poles is not None and (P - T) not in poles: continue
                    for sh, cf in _taylor(remA, T):
                        mono = tuple(sorted(sh + list(remB))); r = P - T
                        d = out.setdefault(r, {}); d[mono] = d.get(mono, 0) + fac * cf
    return {r: {m: c for m, c in f.items() if c} for r, f in out.items()}
