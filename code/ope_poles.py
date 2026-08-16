"""General pole-coefficient extraction for bosonic OPEs (the
higher-pole extension; first consumer: the Getmanov stage-2 work).

For canonical bosonic super-engine monomials P, Q (parity-0 letters),
ope_pole_mono(P, Q, j) returns the local field multiplying
(z-w)^{-j} in the OPE P(z)Q(w), including the c-number (empty
monomial) contributions -- the same Wick enumerator as
super_engine.residue_mono / normal_prod_mono with Taylor order D - j.
j = 1 reproduces the residue; j = 0 the point-split normal product.

Coefficients may be exact Fractions or sympy expressions (scalar
arithmetic is generic)."""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import factorial

from ope_engine import _compositions
from super_engine import _koszul, _pairing, canonicalize


def ope_pole_mono(P, Q, j):
    out = {}
    np_, nq = len(P), len(Q)
    for k in range(0, min(np_, nq) + 1):
        for psel in combinations(range(np_), k):
            for qsel in permutations(range(nq), k):
                num, D = 1, 0
                bad = False
                for a, b in zip(psel, qsel):
                    pr = _pairing(P[a], Q[b])
                    if pr is None:
                        bad = True
                        break
                    num *= pr[0]
                    D += pr[1]
                if bad or not num or D < j:
                    continue
                num *= _koszul(P, Q, psel, qsel)
                restP = [P[a] for a in range(np_) if a not in psel]
                restQ = [Q[b] for b in range(nq) if b not in set(qsel)]
                r = len(restP)
                if r == 0:
                    if D == j:
                        cm, sg = canonicalize(tuple(restQ))
                        if cm is not None:
                            c = out.get(cm, 0) + F(num * sg)
                            if c:
                                out[cm] = c
                            else:
                                out.pop(cm, None)
                    continue
                for ts in _compositions(D - j, r):
                    den = 1
                    shifted = []
                    for (p, m, jj), t in zip(restP, ts):
                        den *= factorial(t)
                        shifted.append((p, m + t, jj))
                    cm, sg = canonicalize(tuple(shifted) + tuple(restQ))
                    if cm is None:
                        continue
                    c = out.get(cm, 0) + F(num * sg, den)
                    if c:
                        out[cm] = c
                    else:
                        out.pop(cm, None)
    return out


def ope_pole(A, B, j):
    """Pole-j coefficient of A(z)B(w) for currents (dicts)."""
    out = {}
    for ma, ca in A.items():
        for mb, cb in B.items():
            for m, c in ope_pole_mono(ma, mb, j).items():
                v = out.get(m, 0) + ca * cb * c
                if v:
                    out[m] = v
                else:
                    out.pop(m, None)
    return out
