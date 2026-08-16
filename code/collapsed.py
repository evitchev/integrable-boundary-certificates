"""Multiplicity-class-collapsed Wick enumeration for the mixed engine.

Drop-in fast replacement for mixed_engine.residue_mixed_mono.  Two
collapses, each exact:

1. X-letter contractions are enumerated as class-count profiles
   {c_ab: # contractions between P-letters of order a and Q-letters of
   order b}, with the raw-matching multiplicity restored by
       prod_a n_a!/((n_a-r_a)! prod_b c_ab!) * prod_b m_b!/(m_b-s_b)!
   (r_a = sum_b c_ab, s_b = sum_a c_ab) --- replacing the
   sum-over-position-matchings, which is exponential in multiplicities.
2. Taylor budgets on free P X-letters are distributed as PARTITIONS per
   order class with weight k!/prod(mult_v!) * prod 1/t_i!, replacing
   ordered compositions (whose lru_cache was the OOM source: up to
   ~1e8 tuples for the spin-21 targets).

Y-pair slots keep the original raw enumeration (bounded: few pairs),
including the index-graph path/cycle analysis.  Output format and
values are identical to mixed_engine.residue_mixed_mono; equality is
checked on retained randomized and medium all-ones cases by
collapsed_validate.py (not an exhaustive proof of equivalence).
"""

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from math import factorial

from mixed_engine import _y_matchings, canon_mixed
from invariant_engine import canon_pair


@lru_cache(maxsize=None)
def partitions_leq(T, k):
    """Sorted (descending) tuples of at most k positive ints summing to T."""
    if T == 0:
        return ((),)
    if k == 0:
        return ()
    out = []
    def rec(rem, kmax, vmax, acc):
        if rem == 0:
            out.append(tuple(acc))
            return
        if kmax == 0:
            return
        for v in range(min(vmax, rem), 0, -1):
            acc.append(v)
            rec(rem - v, kmax - 1, v, acc)
            acc.pop()
    rec(T, k, T, [])
    return tuple(out)


def _xprofiles(pcls, qcls):
    """Yield (weight, D, freeP, restQ) over contraction-class profiles.
    pcls/qcls: tuples of (order, count).  freeP/restQ: same format."""
    B = len(qcls)

    def rec_a(i, qrem, w, D, freeP):
        if i == len(pcls):
            wq = 1
            for (b, mb), rem in zip(qcls, qrem):
                wq *= factorial(mb) // factorial(rem)
            yield (w * wq, D, tuple(freeP),
                   tuple((b, rem) for (b, _), rem in zip(qcls, qrem)
                         if rem))
            return
        a, na = pcls[i]

        def rec_b(j, left, qr, wa, Da, cprod):
            if j == B:
                r_a = na - left
                w_a = factorial(na) // (factorial(left) * cprod)
                nf = freeP + ([(a, left)] if left else [])
                yield from rec_a(i + 1, qr, w * wa * w_a, D + Da, nf)
                return
            b, _ = qcls[j]
            for c in range(0, min(left, qr[j]) + 1):
                qr2 = list(qr)
                qr2[j] -= c
                fac = ((-1) ** (a - 1) * factorial(a + b - 1)) ** c
                yield from rec_b(j + 1, left - c, qr2,
                                 wa * fac, Da + c * (a + b),
                                 cprod * factorial(c))
        yield from rec_b(0, na, list(qrem), 1, 0, 1)

    yield from rec_a(0, [m for _, m in qcls], 1, 0, [])


def _class_taylor(freeP, budget):
    """Distribute `budget` Taylor units over X-letter classes.
    Yields (weight_num, weight_den, letters_tuple)."""
    classes = list(freeP)

    def rec(i, rem, wn, wd, letters):
        if i == len(classes):
            if rem == 0:
                yield (wn, wd, tuple(letters))
            return
        a, k = classes[i]
        top = rem if i == len(classes) - 1 else rem
        for Ta in range(0, top + 1):
            if i == len(classes) - 1 and Ta != rem:
                continue
            for part in partitions_leq(Ta, k):
                mult = Counter(part)
                mult[0] = k - len(part)
                wnum = factorial(k)
                wden = 1
                for v, m in mult.items():
                    wden *= factorial(m)
                for t in part:
                    wden *= factorial(t)
                new = [a + t for t in part] + [a] * (k - len(part))
                yield from rec(i + 1, rem - Ta, wn * wnum, wd * wden,
                               letters + new)
    yield from rec(0, budget, 1, 1, [])


def residue_mixed_mono_fast(Pm, Qm):
    pxs, pyp = Pm
    qxs, qyp = Qm
    pcls = tuple(sorted(Counter(pxs).items()))
    qcls = tuple(sorted(Counter(qxs).items()))
    ymatch = list(_y_matchings(pyp, qyp))
    out = {}
    for xw, xD, freeP, restQ in _xprofiles(pcls, qcls):
        rq_letters = tuple(b for b, m in restQ for _ in range(m))
        for ynum, yD, ncyc, templates in ymatch:
            D = xD + yD
            if D == 0:
                continue
            tmark = [(ti, si) for ti, pr in enumerate(templates)
                     for si in (0, 1) if pr[si][1]]
            base = xw * ynum

            def emit(ts_y, xbudget):
                for wn, wd, xletters in _class_taylor(freeP, xbudget):
                    den = wd
                    for t in ts_y:
                        den *= factorial(t)
                    shifts = dict(zip(tmark, ts_y))
                    yp2 = []
                    for ti, pr in enumerate(templates):
                        o0 = pr[0][0] + shifts.get((ti, 0), 0)
                        o1 = pr[1][0] + shifts.get((ti, 1), 0)
                        yp2.append(canon_pair(o0, o1))
                    mono = canon_mixed(rq_letters + xletters, yp2)
                    c = F(base * wn, den)
                    slot = out.setdefault(mono, {})
                    slot[ncyc] = slot.get(ncyc, 0) + c

            T = D - 1
            ns = len(tmark)
            if ns == 0:
                if freeP:
                    emit((), T)
                continue

            def comp(rem, k, acc):
                if k == 0:
                    emit(tuple(acc), rem)
                    return
                for v in range(rem + 1):
                    acc.append(v)
                    comp(rem - v, k - 1, acc)
                    acc.pop()
            if freeP:
                comp(T, ns, [])
            else:
                # all budget must land on the y-slots exactly
                def comp2(rem, k, acc):
                    if k == 1:
                        acc.append(rem)
                        emit(tuple(acc), 0)
                        acc.pop()
                        return
                    for v in range(rem + 1):
                        acc.append(v)
                        comp2(rem - v, k - 1, acc)
                        acc.pop()
                comp2(T, ns, [])
    return {m: {p: c for p, c in d.items() if c}
            for m, d in out.items() if any(d.values())}


# ---------------------------------------------------------------- router
from math import comb, perm
import mixed_engine as _me

# bind at import time so later monkeypatching of the module attribute
# cannot make the fallback recurse into this router
_ORIG_RESIDUE = _me.residue_mixed_mono

RAW_THRESHOLD = 200_000


def _raw_count(nx, mx):
    return sum(comb(nx, k) * perm(mx, k) for k in range(min(nx, mx) + 1))


@lru_cache(maxsize=4096)
def residue_mixed_mono_smart(Pm, Qm):
    """Collapsed enumeration where raw X-matchings explode; original
    engine elsewhere (each validated exactly against the other on the
    overlap)."""
    nx, mx = len(Pm[0]), len(Qm[0])
    if _raw_count(nx, mx) > RAW_THRESHOLD:
        return residue_mixed_mono_fast(Pm, Qm)
    return _ORIG_RESIDUE(Pm, Qm)
