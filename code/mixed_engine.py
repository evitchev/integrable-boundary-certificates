"""Symbolic-N engine for the Z2 x O(N) ("cylindrical") sector:
one boson X (plain letters, even count) plus N bosons Y entering via
O(N) dot-product pairs (ab)_Y = d^a Y . d^b Y.

Monomial = (xs, yp): xs = sorted tuple of X derivative orders (even
length), yp = sorted tuple of Y pairs (a, b), a <= b.  Wick rules:
X letters contract only with X letters (plain factor), Y slots only
with Y slots (index graph: paths -> new pairs, cycles -> factors of N).
Residue coefficients are {N_power: Fraction}; derivatives are N-free.

Families: Vitchev hep-th/0404195 eqs. (52)-(58), s = sqrt((N-25)(N-1)),
using 4 T_Y^2 = 2 (11)(11)_Y - 4 (22)_Y and (dX)^2 T_Y = (1/2)(11)_X (11)_Y.
"""

from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations
from math import factorial

from invariant_engine import _compositions, _rref, canon_pair


def canon_mixed(xs, yp):
    return (tuple(sorted(xs)), tuple(sorted(yp)))


def spin_mixed(m):
    xs, yp = m
    return sum(xs) + sum(a + b for a, b in yp)


def gen_mixed_basis(s):
    """Z2 x O(N) monomials of spin s: even # of X letters + Y pairs."""
    # X multisets of even length
    def xparts(total):
        res = [[]] if total == 0 else []
        def rec(mx, remaining, acc):
            if remaining == 0:
                if len(acc) % 2 == 0:
                    res.append(list(acc))
                return
            for v in range(min(mx, remaining), 0, -1):
                acc.append(v)
                rec(v, remaining - v, acc)
                acc.pop()
        rec(total, total, [])
        return [xs for xs in res if len(xs) % 2 == 0]

    from invariant_engine import gen_inv_basis
    out = []
    for xspin in range(0, s + 1):
        for xs in xparts(xspin):
            if xspin == s and xs == []:
                continue
            rem = s - xspin
            if rem == 0:
                if xs:
                    out.append(canon_mixed(xs, ()))
                continue
            for yp in gen_inv_basis(rem):
                out.append(canon_mixed(xs, yp))
    return sorted(set(out))


def d_mixed_mono(m):
    xs, yp = m
    out = {}
    for i in range(len(xs)):
        new = canon_mixed(xs[:i] + (xs[i] + 1,) + xs[i + 1:], yp)
        out[new] = out.get(new, 0) + 1
    for i, (a, b) in enumerate(yp):
        for np_ in (canon_pair(a + 1, b), canon_pair(a, b + 1)):
            new = canon_mixed(xs, yp[:i] + (np_,) + yp[i + 1:])
            out[new] = out.get(new, 0) + 1
    return out


def _x_matchings(pxs, qxs):
    """Yield (factor, D, free_P_orders, rest_Q_orders) over X matchings
    (including the empty one)."""
    lp, lq = len(pxs), len(qxs)
    for k in range(0, min(lp, lq) + 1):
        for psel in combinations(range(lp), k):
            for qsel in permutations(range(lq), k):
                num, D = 1, 0
                for x, y in zip(psel, qsel):
                    m, n = pxs[x], qxs[y]
                    num *= (-1) ** (m - 1) * factorial(m + n - 1)
                    D += m + n
                freeP = [pxs[i] for i in range(lp) if i not in psel]
                restQ = [qxs[j] for j in range(lq) if j not in set(qsel)]
                yield num, D, freeP, restQ


def _y_matchings(pyp, qyp):
    """Yield (factor, D, ncycles, templates) over Y-slot matchings
    (including the empty one); template = ((order, isP), (order, isP))."""
    kP, kQ = len(pyp), len(qyp)
    PS = [(i, s) for i in range(kP) for s in (0, 1)]
    QS = [(j, s) for j in range(kQ) for s in (0, 1)]
    for k in range(0, min(2 * kP, 2 * kQ) + 1):
        for psel in combinations(range(2 * kP), k):
            for qsel in permutations(range(2 * kQ), k):
                num, D = 1, 0
                padj = [[None, None] for _ in range(kP)]
                qadj = [[None, None] for _ in range(kQ)]
                for x, y in zip(psel, qsel):
                    pi, ps = PS[x]
                    qj, qs = QS[y]
                    m = pyp[pi][ps]
                    n = qyp[qj][qs]
                    num *= (-1) ** (m - 1) * factorial(m + n - 1)
                    D += m + n
                    padj[pi][ps] = ('Q', qj)
                    qadj[qj][qs] = ('P', pi)
                visited = set()
                ncycles = 0
                templates = []
                for side, adj, prs in (('P', padj, pyp), ('Q', qadj, qyp)):
                    for v in range(len(prs)):
                        if (side, v) in visited:
                            continue
                        comp, stack = [], [(side, v)]
                        while stack:
                            node = stack.pop()
                            if node in visited or node in comp:
                                continue
                            comp.append(node)
                            sd, idx = node
                            adjacency = padj[idx] if sd == 'P' else qadj[idx]
                            for sl in (0, 1):
                                if adjacency[sl] is not None:
                                    stack.append(
                                        (adjacency[sl][0], adjacency[sl][1]))
                        for node in comp:
                            visited.add(node)
                        free = []
                        for sd, idx in comp:
                            adjacency = padj[idx] if sd == 'P' else qadj[idx]
                            prs2 = pyp if sd == 'P' else qyp
                            for sl in (0, 1):
                                if adjacency[sl] is None:
                                    free.append((prs2[idx][sl], sd == 'P'))
                        if not free:
                            ncycles += 1
                        else:
                            templates.append((free[0], free[1]))
                yield num, D, ncycles, templates


@lru_cache(maxsize=None)
def residue_mixed_mono(P, Q):
    """Res of P(z) Q(w) for mixed monomials; {mono: {N_power: Fraction}}."""
    pxs, pyp = P
    qxs, qyp = Q
    out = {}
    ymatch = list(_y_matchings(pyp, qyp))
    for xnum, xD, xfreeP, xrestQ in _x_matchings(pxs, qxs):
        for ynum, yD, ncyc, templates in ymatch:
            D = xD + yD
            if D == 0:
                continue
            num = xnum * ynum
            # free P letters: X ones, then P-side template slots
            slots = []          # (kind, ref); kind 'x' or ('t', ti, si)
            for o in xfreeP:
                slots.append(o)
            tmark = []
            for ti, pr in enumerate(templates):
                for si in (0, 1):
                    if pr[si][1]:
                        tmark.append((ti, si))
            r = len(slots) + len(tmark)
            if r == 0:
                continue
            for ts in _compositions(D - 1, r):
                den = 1
                for t in ts:
                    den *= factorial(t)
                xs2 = list(xrestQ)
                for o, t in zip(slots, ts[:len(slots)]):
                    xs2.append(o + t)
                shifts = {}
                for (ti, si), t in zip(tmark, ts[len(slots):]):
                    shifts[(ti, si)] = t
                yp2 = []
                for ti, pr in enumerate(templates):
                    o0 = pr[0][0] + shifts.get((ti, 0), 0)
                    o1 = pr[1][0] + shifts.get((ti, 1), 0)
                    yp2.append(canon_pair(o0, o1))
                mono = canon_mixed(xs2, yp2)
                slot = out.setdefault(mono, {})
                slot[ncyc] = slot.get(ncyc, 0) + F(num, den)
    return {m: {p: c for p, c in d.items() if c}
            for m, d in out.items() if any(d.values())}


def residue_mixed_cur(Pcur, Qcur, Nval):
    out = {}
    for mp, cp in Pcur.items():
        for mq, cq in Qcur.items():
            for mono, powd in residue_mixed_mono(mp, mq).items():
                val = cp * cq * sum(c * Nval ** p for p, c in powd.items())
                if val:
                    v = out.get(mono, 0) + val
                    if v:
                        out[mono] = v
                    else:
                        out.pop(mono, None)
    return out


class MixedSector:
    _cache = {}

    def __new__(cls, s):
        if s in cls._cache:
            return cls._cache[s]
        self = super().__new__(cls)
        self.s = s
        self.basis = gen_mixed_basis(s)
        self.index = {m: i for i, m in enumerate(self.basis)}
        cols = []
        for b in gen_mixed_basis(s - 1):
            col = [F(0)] * len(self.basis)
            for mono, c in d_mixed_mono(b).items():
                col[self.index[mono]] = F(c)
            if any(col):
                cols.append(col)
        rows = [list(c) for c in cols]
        rows, piv = _rref(rows, len(self.basis))
        self._img_rows = rows[:len(piv)]
        self._img_pivots = piv
        cls._cache[s] = self
        return self

    def vec(self, cur):
        v = [F(0)] * len(self.basis)
        for mono, c in cur.items():
            v[self.index[mono]] = c
        return v

    def reduce_mod_d(self, v):
        v = list(v)
        for row, p in zip(self._img_rows, self._img_pivots):
            if v[p]:
                f = v[p]
                v = [a - f * b for a, b in zip(v, row)]
        return v

    def is_total_derivative(self, cur):
        return not any(self.reduce_mod_d(self.vec(cur)))

    def q_dim(self):
        return len(self.basis) - len(self._img_pivots)


def genuine_kernel_mixed(P_low, cand_basis, Nval):
    from invariant_engine import nullspace_cols
    s_low = spin_mixed(next(iter(P_low)))
    s_cand = spin_mixed(cand_basis[0])
    tgt = MixedSector(s_low + s_cand - 1)
    cs = MixedSector(s_cand)
    cols = [tgt.reduce_mod_d(tgt.vec(
        residue_mixed_cur(P_low, {b: F(1)}, Nval))) for b in cand_basis]
    ker = nullspace_cols(cols, len(tgt.basis))
    reduced, genuine = [], []
    for v in ker:
        cur = {b: c for b, c in zip(cand_basis, v) if c}
        w = cs.reduce_mod_d(cs.vec(cur))
        if not any(w):
            continue
        rows = [list(x) for x in reduced + [w]]
        _, piv = _rref(rows, len(w))
        if len(piv) == len(reduced) + 1:
            reduced.append(w)
            genuine.append(cur)
    return genuine


# ------------------------------------------------- cylindrical families

def _M(xs, yp=()):
    return canon_mixed(xs, tuple(canon_pair(*p) for p in yp))


def cyl_P4(sol, N, s):
    """I_3 densities of Solutions 1-3, eqs. (53), (55), (57).

    Normalization note: taking the printed "4 T_Y^2" at face value with
    eq. (52) (2 T_Y^2 = (11)(11)_Y - 2 (22)_Y) gives densities whose
    N=1 limit FAILS the paper's own paperclip-reduction claim and whose
    ad-kernel at spin 5 is empty; with the Y-part read as
    (11)(11)_Y - 2 (22)_Y (i.e. half), the N=1 limits of Solutions 1
    and 3 match the paperclip P_4(n=-1) exactly.  We use the latter;
    certified a posteriori by nonempty I_5 kernels (cyl_diagnosis.py).
    (dX)^2 T_Y = (1/2)(11)_X (11)_Y throughout."""
    assert s * s == (N - 25) * (N - 1)
    base = {_M((), [(1, 1), (1, 1)]): F(1), _M((), [(2, 2)]): F(-2)}
    if sol == 1:
        extra = {
            _M((2, 2)): (-9 - 3 * N - s) / F(6),
            _M((1, 1), [(1, 1)]): F(6),
            _M((1, 1, 1, 1)): (2 + N + s) / F(3),
        }
    elif sol == 2:
        extra = {
            _M((2, 2)): -(409 + 46 * N + N * N + (N - 5) * s) / F(192),
            _M((1, 1), [(1, 1)]): (11 + N + s) / F(8),
            _M((1, 1, 1, 1)): (117 + 2 * N + N * N + (15 + N) * s) / F(192),
        }
    elif sol == 3:
        extra = {
            _M((2, 2)): (13 - 27 * N + 2 * N * N + (5 - 2 * N) * s) / F(6),
            _M((1, 1), [(1, 1)]): (7 - N + s) * F(1),
            _M((1, 1, 1, 1)): F(1),
        }
    else:
        raise ValueError(sol)
    base.update(extra)
    return {m: c for m, c in base.items() if c}
