"""VIR1 library (Fable seat) -- own engine, independent of the repo's.

Z2 x O(N) mixed free-field algebra: one boson X, an O(N) vector Y entering
through pairs (ab) = d^a Y . d^b Y.  Monomial = (xs, yp), xs a sorted tuple of
X derivative orders, yp a sorted tuple of pairs (a, b), a <= b.

Propagator (engine convention, phi(z) phi(w) ~ log(z - w)):
    d^m phi(z) d^n phi(w)  ~  (-1)^(m-1) (m+n-1)! / (z-w)^(m+n).

ope(A, B, power): the coefficient of (z-w)^power in A(z) B(w), ALL powers
(singular and regular), as {mono: {ncycles: Fraction}}; a closed Y index loop
carries one factor of the loop variable.
"""
import sys
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations
from math import factorial

sys.dont_write_bytecode = True


# ---------------------------------------------------------------- monomials
def canon(xs, yp):
    return (tuple(sorted(xs)), tuple(sorted(tuple(sorted(p)) for p in yp)))


def weight(m):
    return sum(m[0]) + sum(a + b for a, b in m[1])


def _partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in _partitions(n - k, k):
            yield (k,) + rest


@lru_cache(maxsize=None)
def y_monos(w):
    """All Y-only monomials (multisets of pairs) of weight w."""
    pairs = [(a, b) for a in range(1, w) for b in range(a, w) if a + b <= w]
    out = set()

    def rec(start, remaining, acc):
        if remaining == 0:
            out.add(tuple(acc))
            return
        for i in range(start, len(pairs)):
            p = pairs[i]
            if p[0] + p[1] <= remaining:
                rec(i, remaining - p[0] - p[1], acc + [p])

    rec(0, w, [])
    return sorted(out)


@lru_cache(maxsize=None)
def x_monos(w, parity=0):
    """X-only monomials of weight w with letter count == parity mod 2."""
    return sorted(tuple(sorted(p)) for p in _partitions(w)
                  if len(p) % 2 == parity)


@lru_cache(maxsize=None)
def basis(w):
    """Even mixed basis of weight w."""
    out = []
    for a in range(0, w + 1):
        for xs in x_monos(a, 0):
            for yp in y_monos(w - a):
                out.append((xs, yp))
    return sorted(out)


def d_mono(m):
    xs, yp = m
    out = {}
    for i in range(len(xs)):
        new = canon(xs[:i] + (xs[i] + 1,) + xs[i + 1:], yp)
        out[new] = out.get(new, 0) + 1
    for i, (a, b) in enumerate(yp):
        for np_ in ((a + 1, b), (a, b + 1)):
            new = canon(xs, yp[:i] + (np_,) + yp[i + 1:])
            out[new] = out.get(new, 0) + 1
    return out


def d_vec(v):
    out = {}
    for m, c in v.items():
        for m2, c2 in d_mono(m).items():
            out[m2] = out.get(m2, 0) + c * c2
    return {m: c for m, c in out.items() if c}


# ---------------------------------------------------------------- Wick / OPE
def _compositions(total, parts):
    if parts == 0:
        if total == 0:
            yield ()
        return
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in _compositions(total - first, parts - 1):
            yield (first,) + rest


def _contr(m, n):
    return (-1) ** (m - 1) * factorial(m + n - 1)


def _x_matchings(pxs, qxs):
    lp, lq = len(pxs), len(qxs)
    for k in range(0, min(lp, lq) + 1):
        for psel in combinations(range(lp), k):
            for qsel in permutations(range(lq), k):
                num, D = 1, 0
                for x, y in zip(psel, qsel):
                    num *= _contr(pxs[x], qxs[y])
                    D += pxs[x] + qxs[y]
                ps, qs = set(psel), set(qsel)
                freeP = [pxs[i] for i in range(lp) if i not in ps]
                restQ = [qxs[j] for j in range(lq) if j not in qs]
                yield num, D, freeP, restQ


def _y_matchings(pyp, qyp):
    """Yield (factor, D, ncycles, paths); a path is a pair of free slots
    ((order, isP), (order, isP)) joined by a chain of contractions."""
    kP, kQ = len(pyp), len(qyp)
    PS = [(i, s) for i in range(kP) for s in (0, 1)]
    QS = [(j, s) for j in range(kQ) for s in (0, 1)]
    for k in range(0, min(len(PS), len(QS)) + 1):
        for psel in combinations(range(len(PS)), k):
            for qsel in permutations(range(len(QS)), k):
                num, D = 1, 0
                link = {}
                for x, y in zip(psel, qsel):
                    pi, ps = PS[x]
                    qj, qs = QS[y]
                    num *= _contr(pyp[pi][ps], qyp[qj][qs])
                    D += pyp[pi][ps] + qyp[qj][qs]
                    link[('P', pi, ps)] = ('Q', qj, qs)
                    link[('Q', qj, qs)] = ('P', pi, ps)
                # walk the index graph
                seen = set()
                ncyc = 0
                paths = []
                nodes = [('P', i) for i in range(kP)] + [('Q', j) for j in range(kQ)]

                def order(side, idx, sl):
                    return (pyp if side == 'P' else qyp)[idx][sl]

                # paths: start from every free slot
                for side, idx in nodes:
                    for sl in (0, 1):
                        if (side, idx, sl) in link or (side, idx, sl) in seen:
                            continue
                        # free slot: walk through the partner slot
                        start = (side, idx, sl)
                        seen.add(start)
                        cur = (side, idx, 1 - sl)
                        while cur in link:
                            seen.add(cur)
                            nxt = link[cur]
                            seen.add(nxt)
                            cur = (nxt[0], nxt[1], 1 - nxt[2])
                        seen.add(cur)
                        paths.append(((order(*start), start[0] == 'P'),
                                      (order(*cur), cur[0] == 'P')))
                # cycles: remaining unseen slots are all contracted
                for side, idx in nodes:
                    for sl in (0, 1):
                        if (side, idx, sl) in seen:
                            continue
                        ncyc += 1
                        cur = (side, idx, sl)
                        while cur not in seen:
                            seen.add(cur)
                            nxt = link[cur]
                            seen.add(nxt)
                            cur = (nxt[0], nxt[1], 1 - nxt[2])
                yield num, D, ncyc, paths


@lru_cache(maxsize=None)
def ope(P, Q, power):
    """Coefficient of (z-w)^power in P(z) Q(w): {mono: {ncyc: Fraction}}."""
    pxs, pyp = P
    qxs, qyp = Q
    out = {}
    ymatch = list(_y_matchings(pyp, qyp))
    for xnum, xD, xfreeP, xrestQ in _x_matchings(pxs, qxs):
        for ynum, yD, ncyc, paths in ymatch:
            D = xD + yD
            T = power + D          # total Taylor shift of the free P letters
            if T < 0:
                continue
            tmark = [(ti, si) for ti, pr in enumerate(paths)
                     for si in (0, 1) if pr[si][1]]
            r = len(xfreeP) + len(tmark)
            if r == 0 and T != 0:
                continue
            num = xnum * ynum
            for ts in _compositions(T, r):
                den = 1
                for t in ts:
                    den *= factorial(t)
                xs2 = list(xrestQ)
                for o, t in zip(xfreeP, ts[:len(xfreeP)]):
                    xs2.append(o + t)
                shifts = dict(zip(tmark, ts[len(xfreeP):]))
                yp2 = []
                for ti, pr in enumerate(paths):
                    yp2.append((pr[0][0] + shifts.get((ti, 0), 0),
                                pr[1][0] + shifts.get((ti, 1), 0)))
                mono = canon(xs2, yp2)
                slot = out.setdefault(mono, {})
                slot[ncyc] = slot.get(ncyc, 0) + F(num, den)
    return {m: {p: c for p, c in d.items() if c}
            for m, d in out.items() if any(d.values())}


def ope_vec(Pv, Qv, power, loop):
    """Coefficient of (z-w)^power in Pv(z) Qv(w), loop variable -> `loop`."""
    out = {}
    for mp, cp in Pv.items():
        for mq, cq in Qv.items():
            for mono, powd in ope(mp, mq, power).items():
                val = cp * cq * sum(c * loop ** p for p, c in powd.items())
                if val:
                    out[mono] = out.get(mono, 0) + val
    return {m: c for m, c in out.items() if c}


def residue(Pv, Qv, loop):
    return ope_vec(Pv, Qv, -1, loop)


T_Y = {((), ((1, 1),)): F(1, 2)}


def L(n, v, loop):
    """Virasoro mode L_n of T_Y = (1/2)(11) acting on the state v."""
    return ope_vec(T_Y, v, -n - 2, loop)


# ---------------------------------------------------------------- Virasoro
def vir_partitions(w):
    """Partitions of w into parts >= 2 (L_{-n1}...L_{-nk}|0>, n1>=...>=nk)."""
    return [p for p in _partitions(w) if all(k >= 2 for k in p)]


VAC = {((), ()): F(1)}


@lru_cache(maxsize=None)
def vir_state(part, loop):
    v = dict(VAC)
    for n in reversed(part):          # L_{-n1} ... L_{-nk} |0>
        v = L(-n, v, loop)
    return tuple(sorted(v.items()))


def vir_basis(w, loop):
    return [dict(vir_state(p, loop)) for p in vir_partitions(w)]


def hv_span(w, loop=F(7, 3)):
    """Spanning set of HV_w = (even X-monomials) x Vir_N at total weight w."""
    out = []
    for a in range(0, w + 1):
        for xs in x_monos(a, 0):
            for v in (vir_basis(w - a, loop) if w - a else [dict(VAC)]):
                out.append({canon(xs, m[1]): c for m, c in v.items()})
    return out


# ---------------------------------------------------------------- classical
def classical_y_gens(w, maxspin):
    """Classical differential-polynomial algebra generated by the singlets
    (kk), 2k <= maxspin: Y-only spanning set at weight w."""
    letters = []                          # (weight, vector) = d^j (kk)
    for k in range(1, maxspin // 2 + 1):
        v = {((), ((k, k),)): F(1)}
        wt = 2 * k
        while wt <= w:
            letters.append((wt, v))
            v = d_vec(v)
            wt += 1
    out = []

    def mult(v1, v2):
        r = {}
        for m1, c1 in v1.items():
            for m2, c2 in v2.items():
                m = canon((), m1[1] + m2[1])
                r[m] = r.get(m, 0) + c1 * c2
        return r

    def rec(start, remaining, acc):
        if remaining == 0:
            out.append(acc)
            return
        for i in range(start, len(letters)):
            wt, v = letters[i]
            if wt <= remaining:
                rec(i, remaining - wt, mult(acc, v))

    rec(0, w, dict(VAC))
    return out


def classical_span(w, maxspin):
    out = []
    for a in range(0, w + 1):
        for xs in x_monos(a, 0):
            ys = classical_y_gens(w - a, maxspin) if w - a else [dict(VAC)]
            for v in ys:
                out.append({canon(xs, m[1]): c for m, c in v.items()})
    return out


# ---------------------------------------------------------------- linear algebra
def to_row(v, index):
    row = [F(0)] * len(index)
    for m, c in v.items():
        row[index[m]] += F(c)
    return row


def rref(rows):
    """Reduced row echelon form over Q; returns (rows, pivots)."""
    rows = [list(r) for r in rows]
    piv = []
    r = 0
    ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        p = None
        for i in range(r, len(rows)):
            if rows[i][c]:
                p = i
                break
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        inv = 1 / rows[r][c]
        rows[r] = [x * inv for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                f = rows[i][c]
                rows[i] = [x - f * y for x, y in zip(rows[i], rows[r])]
        piv.append(c)
        r += 1
        if r == len(rows):
            break
    return rows[:r], piv


def rank(rows):
    return len(rref(rows)[1]) if rows else 0


def nullspace(cols, nrows):
    """Right null space of the matrix whose COLUMNS are `cols` (lists of
    length nrows); returns list of coefficient vectors."""
    n = len(cols)
    rows = [[cols[j][i] for j in range(n)] for i in range(nrows)]
    R, piv = rref(rows) if rows else ([], [])
    free = [j for j in range(n) if j not in piv]
    out = []
    for f in free:
        v = [F(0)] * n
        v[f] = F(1)
        for i, p in enumerate(piv):
            v[p] = -R[i][f]
        out.append(v)
    return out


def left_null(span_rows):
    """Functionals lam with lam . r = 0 for every row r of span_rows, in
    reduced echelon form."""
    n = len(span_rows[0])
    R, piv = rref(span_rows)
    # null space of R as a matrix acting on column vectors
    free = [j for j in range(n) if j not in piv]
    out = []
    for f in free:
        v = [F(0)] * n
        v[f] = F(1)
        for i, p in enumerate(piv):
            v[p] = -R[i][f]
        out.append(v)
    return rref(out)[0] if out else []


def quotient_functionals(w, span):
    """Functionals on V_w vanishing on span + d V_{w-1}."""
    B = basis(w)
    idx = {m: i for i, m in enumerate(B)}
    rows = [to_row(v, idx) for v in span]
    rows += [to_row(d_mono(m), idx) for m in basis(w - 1)]
    return left_null(rows), B, idx


# ---------------------------------------------------------------- the curve
def curve_point(t):
    t = F(t)
    return (t * t - 25) / (t * t - 1), F(-24) * t / (t * t - 1)


def P4(sol, N, s):
    """I_3 densities, transcribed from mixed_engine.cyl_P4 (checked against
    the repo function in v0_selftest.py)."""
    assert s * s == (N - 25) * (N - 1)
    M = canon
    base = {M((), [(1, 1), (1, 1)]): F(1), M((), [(2, 2)]): F(-2)}
    if sol == 1:
        extra = {M((2, 2), ()): (-9 - 3 * N - s) / F(6),
                 M((1, 1), [(1, 1)]): F(6),
                 M((1, 1, 1, 1), ()): (2 + N + s) / F(3)}
    elif sol == 2:
        extra = {M((2, 2), ()): -(409 + 46 * N + N * N + (N - 5) * s) / F(192),
                 M((1, 1), [(1, 1)]): (11 + N + s) / F(8),
                 M((1, 1, 1, 1), ()): (117 + 2 * N + N * N + (15 + N) * s) / F(192)}
    elif sol == 3:
        extra = {M((2, 2), ()): (13 - 27 * N + 2 * N * N + (5 - 2 * N) * s) / F(6),
                 M((1, 1), [(1, 1)]): (7 - N + s) * F(1),
                 M((1, 1, 1, 1), ()): F(1)}
    else:
        raise ValueError(sol)
    base.update(extra)
    return {m: c for m, c in base.items() if c}


# ---------------------------------------------------------------- kernels
class Commutant:
    """ad_{I_3} kernel at weight w modulo total derivatives."""

    def __init__(self, w):
        self.w = w
        self.B = basis(w)
        self.Bt = basis(w + 3)
        self.it = {m: i for i, m in enumerate(self.Bt)}
        # projector onto V_{w+3} / d V_{w+2}: left null functionals of Im d
        drows = [to_row(d_mono(m), self.it) for m in basis(w + 2)]
        self.proj = left_null(drows)
        # complement of d V_{w-1} in V_w: non-pivot monomials
        ib = {m: i for i, m in enumerate(self.B)}
        R, piv = rref([to_row(d_mono(m), ib) for m in basis(w - 1)])
        self.dpiv = set(piv)
        self.comp = [i for i in range(len(self.B)) if i not in self.dpiv]

    def kernel(self, P, loop):
        """Kernel classes: list of vectors over self.B supported on the
        complement monomials (unique representatives mod d V_{w-1})."""
        cols = []
        for i in self.comp:
            r = residue(P, {self.B[i]: F(1)}, loop)
            row = to_row(r, self.it)
            cols.append([sum(a * b for a, b in zip(lam, row) if a and b)
                         for lam in self.proj])
        ns = nullspace(cols, len(self.proj))
        out = []
        for v in ns:
            full = [F(0)] * len(self.B)
            for c, i in zip(v, self.comp):
                full[i] = c
            out.append(full)
        return out


def apply(funcs, vec):
    return [sum(a * b for a, b in zip(lam, vec) if a and b) for lam in funcs]
