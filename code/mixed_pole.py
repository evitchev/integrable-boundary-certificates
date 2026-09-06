"""Pole-k coefficient extraction for the mixed Z2 x O(N) engine -- M1 of the
mode-action layer for the cylindrical oper program.

ope_pole_mixed_mono(P, Q, k) -> {mono: {N_power: Fraction}} is the local
field multiplying (z-w)^{-k} in P(z)Q(w) for mixed monomials P, Q,
including c-number (empty-monomial) contributions; k = 1 reproduces
mixed_engine.residue_mixed_mono, k = 0 the point-split normal product.
Same Wick enumerators as the mixed engine (imported, not copied), Taylor
order D - k instead of D - 1 -- the mixed sibling of ope_poles.ope_pole_mono.
mixed_engine.py itself is untouched: it sits in the source closures of
pinned attested records.  Validation: mixed_pole_validate.py."""
from fractions import Fraction as F
from math import factorial

from mixed_engine import _x_matchings, _y_matchings, canon_mixed, canon_pair
from ope_engine import _compositions


def ope_pole_mixed_mono(P, Q, k):
    pxs, pyp = P
    qxs, qyp = Q
    out = {}
    ymatch = list(_y_matchings(pyp, qyp))
    for xnum, xD, xfreeP, xrestQ in _x_matchings(pxs, qxs):
        for ynum, yD, ncyc, templates in ymatch:
            D = xD + yD
            if D < k:
                continue
            num = xnum * ynum
            slots = list(xfreeP)
            tmark = [(ti, si) for ti, pr in enumerate(templates)
                     for si in (0, 1) if pr[si][1]]
            r = len(slots) + len(tmark)
            if r == 0:
                if D == k:
                    yp2 = [canon_pair(pr[0][0], pr[1][0]) for pr in templates]
                    mono = canon_mixed(list(xrestQ), yp2)
                    slot = out.setdefault(mono, {})
                    slot[ncyc] = slot.get(ncyc, 0) + F(num)
                continue
            for ts in _compositions(D - k, r):
                den = 1
                for t in ts:
                    den *= factorial(t)
                xs2 = list(xrestQ) + [o + t for o, t in zip(slots, ts[:len(slots)])]
                shifts = {tm: t for tm, t in zip(tmark, ts[len(slots):])}
                yp2 = [canon_pair(pr[0][0] + shifts.get((ti, 0), 0),
                                  pr[1][0] + shifts.get((ti, 1), 0))
                       for ti, pr in enumerate(templates)]
                mono = canon_mixed(xs2, yp2)
                slot = out.setdefault(mono, {})
                slot[ncyc] = slot.get(ncyc, 0) + F(num, den)
    return {m: {p: c for p, c in d.items() if c}
            for m, d in out.items() if any(d.values())}


def ope_pole_mixed_cur(Pcur, Qcur, k, Nval):
    """Pole-k coefficient of P(z)Q(w) for mixed currents at rational N."""
    out = {}
    for mp, cp in Pcur.items():
        for mq, cq in Qcur.items():
            for mono, powd in ope_pole_mixed_mono(mp, mq, k).items():
                val = cp * cq * sum(c * Nval ** p for p, c in powd.items())
                if val:
                    v = out.get(mono, 0) + val
                    if v:
                        out[mono] = v
                    else:
                        out.pop(mono, None)
    return out
