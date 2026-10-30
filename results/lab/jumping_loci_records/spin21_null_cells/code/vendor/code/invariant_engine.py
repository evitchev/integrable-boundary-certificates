"""Symbolic-N OPE engine for O(N)-invariant free-boson currents.

An invariant monomial is a multiset of "pairs" (a, b), a <= b, each pair
meaning the dot product d^a X . d^b X (index summed).  By the first
fundamental theorem for O(N), these generate all invariants, and for
formal N (equivalently N larger than the letter count) monomials are
linearly independent -- the algebra is FORMAL in N.

Wick contractions (paper's eq. (17) pairing) connect letters of P(z) to
letters of Q(w).  Each pair is a valence-2 vertex in the index graph, so
contraction components are paths (two free end-letters -> a new pair) or
cycles (closed index loop -> one factor of N).  Hence the residue of an
OPE of invariant monomials is again a combination of invariant
monomials, with coefficients in Q[N].

Coefficient convention: residues are returned as
    {result_monomial: {N_power: Fraction}}.
Total derivatives of invariant monomials are N-free, so sector reduction
matrices are numeric and shared across all N.
"""

from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations
from math import factorial


# ------------------------------------------------------------- monomials

def canon_pair(a, b):
    return (a, b) if a <= b else (b, a)


def canon(pairs):
    return tuple(sorted(pairs))


def spin_inv(mono):
    return sum(a + b for a, b in mono)


def gen_inv_basis(s):
    """All invariant monomials of spin s: multisets of pairs (a,b),
    1 <= a <= b, summing to s.  Generated in non-increasing pair order."""
    pairs = [(a, b) for tot in range(2, s + 1)
             for a in range(1, tot // 2 + 1) for b in [tot - a]]
    pairs.sort(reverse=True)
    out = []

    def rec(idx, remaining, acc):
        if remaining == 0:
            out.append(tuple(sorted(acc)))
            return
        for i in range(idx, len(pairs)):
            a, b = pairs[i]
            if a + b <= remaining:
                acc.append((a, b))
                rec(i, remaining - a - b, acc)
                acc.pop()

    rec(0, s, [])
    return out


def d_inv_mono(mono):
    """Derivative (Leibniz over pairs and slots), integer coefficients."""
    out = {}
    for i, (a, b) in enumerate(mono):
        for np_ in (canon_pair(a + 1, b), canon_pair(a, b + 1)):
            new = canon(mono[:i] + (np_,) + mono[i + 1:])
            out[new] = out.get(new, 0) + 1
    return out


# ------------------------------------------------------------- residue

from math import comb as _comb


@lru_cache(maxsize=None)
def _compositions_small(total, parts):
    if parts == 0:
        return ((),) if total == 0 else ()
    if parts == 1:
        return ((total,),)
    return tuple((h,) + rest for h in range(total + 1)
                 for rest in _compositions_small(total - h, parts - 1))


def _gen_comp(total, parts):
    if parts == 0:
        if total == 0:
            yield ()
        return
    if parts == 1:
        yield (total,)
        return
    for h in range(total + 1):
        for rest in _gen_comp(total - h, parts - 1):
            yield (h,) + rest


def _compositions(total, parts):
    """Weak compositions of `total` into `parts` parts.  Cached tuples
    for small counts; a streaming generator above 100k (the unbounded
    tuple cache was the machine-freezing memory bomb: ~C(30,11) tuples
    per key at the spin-21 targets)."""
    if parts <= 1 or _comb(total + parts - 1, parts - 1) <= 100_000:
        return _compositions_small(total, parts)
    return _gen_comp(total, parts)


@lru_cache(maxsize=None)
def residue_inv_mono(P, Q):
    """Res_{z->w} of invariant monomials P(z) Q(w).

    Returns {mono: {N_power: Fraction}}.
    Sum over nonempty partial matchings of P letter-slots to Q
    letter-slots; matched pair (m, n) contributes
    (-1)^(m-1) (m+n-1)! and pole m+n plus an index edge between the two
    parent pairs; components of the index graph: cycles give N, paths
    give a new pair from their two free end letters; free P letters are
    Taylor-shifted to soak up the pole down to (z-w)^{-1}.
    """
    kP, kQ = len(P), len(Q)
    # slots: (pair_index, slot) with letter order
    PS = [(i, s, P[i][s]) for i in range(kP) for s in (0, 1)]
    QS = [(j, s, Q[j][s]) for j in range(kQ) for s in (0, 1)]
    out = {}
    for k in range(1, min(2 * kP, 2 * kQ) + 1):
        for psel in combinations(range(2 * kP), k):
            for qsel in permutations(range(2 * kQ), k):
                num = 1
                D = 0
                for x, y in zip(psel, qsel):
                    m = PS[x][2]
                    n = QS[y][2]
                    num *= (-1) ** (m - 1) * factorial(m + n - 1)
                    D += m + n
                # index graph: vertices P0..P(kP-1), Q0..Q(kQ-1)
                # adjacency per vertex slot: matched -> partner vertex
                padj = [[None, None] for _ in range(kP)]
                qadj = [[None, None] for _ in range(kQ)]
                for x, y in zip(psel, qsel):
                    pi, ps, _ = PS[x]
                    qj, qs, _ = QS[y]
                    padj[pi][ps] = ('Q', qj, qs)
                    qadj[qj][qs] = ('P', pi, ps)
                # walk components
                visited = set()
                ncycles = 0
                # template pairs: entries are (order, is_P_side) for free ends
                template = []   # list of ((o1, p1?), (o2, p2?))
                for side, adj, pairs_ in (('P', padj, P), ('Q', qadj, Q)):
                    for v in range(len(pairs_)):
                        if (side, v) in visited:
                            continue
                        # explore component containing (side, v)
                        comp = []
                        stack = [(side, v)]
                        cvis = set()
                        while stack:
                            node = stack.pop()
                            if node in cvis:
                                continue
                            cvis.add(node)
                            comp.append(node)
                            sd, idx = node
                            adjacency = padj[idx] if sd == 'P' else qadj[idx]
                            for sl in (0, 1):
                                nb = adjacency[sl]
                                if nb is not None:
                                    stack.append((nb[0], nb[1]))
                        visited |= cvis
                        # free slots in component
                        free = []
                        for sd, idx in comp:
                            adjacency = padj[idx] if sd == 'P' else qadj[idx]
                            prs = P if sd == 'P' else Q
                            for sl in (0, 1):
                                if adjacency[sl] is None:
                                    free.append((prs[idx][sl], sd == 'P'))
                        if not free:
                            ncycles += 1
                        else:
                            # path: exactly two free ends
                            template.append((free[0], free[1]))
                free_p_positions = []
                flat = []
                for pr in template:
                    for o, isp in pr:
                        if isp:
                            free_p_positions.append(len(flat))
                        flat.append(o)
                r = len(free_p_positions)
                if r == 0:
                    continue  # D >= 2 always
                for ts in _compositions(D - 1, r):
                    den = 1
                    orders = list(flat)
                    for pos, t in zip(free_p_positions, ts):
                        den *= factorial(t)
                        orders[pos] += t
                    mono = canon(tuple(
                        canon_pair(orders[2 * i], orders[2 * i + 1])
                        for i in range(len(template))))
                    c = F(num, den)
                    slot = out.setdefault(mono, {})
                    slot[ncycles] = slot.get(ncycles, 0) + c
    return {m: {p: c for p, c in d.items() if c}
            for m, d in out.items()
            if any(c for c in d.values())}


def residue_inv_cur(Pcur, Qcur, Nval=None):
    """Residue of OPE of invariant currents.  Coefficients of the inputs
    and output are Fractions if Nval (a Fraction) is given; otherwise
    output coefficients are {N_power: Fraction} convolved with numeric
    input coefficients (inputs must then be Fractions too)."""
    out = {}
    for mp, cp in Pcur.items():
        for mq, cq in Qcur.items():
            for mono, powd in residue_inv_mono(mp, mq).items():
                if Nval is not None:
                    val = cp * cq * sum(c * Nval ** p
                                        for p, c in powd.items())
                    if val:
                        slot = out.get(mono, 0) + val
                        if slot:
                            out[mono] = slot
                        else:
                            out.pop(mono, None)
                else:
                    slot = out.setdefault(mono, {})
                    for p, c in powd.items():
                        v = slot.get(p, 0) + cp * cq * c
                        if v:
                            slot[p] = v
                        else:
                            slot.pop(p, None)
    if Nval is None:
        return {m: d for m, d in out.items() if d}
    return out


# ------------------------------------------------- sector mod derivatives

def _rref(rows, ncols):
    pivots = []
    rank = 0
    for col in range(ncols):
        piv = next((r for r in range(rank, len(rows)) if rows[r][col]), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        inv = F(1) / rows[rank][col]
        rows[rank] = [x * inv for x in rows[rank]]
        for r in range(len(rows)):
            if r != rank and rows[r][col]:
                f = rows[r][col]
                rows[r] = [a - f * b for a, b in zip(rows[r], rows[rank])]
        pivots.append(col)
        rank += 1
    return rows, pivots


class InvSector:
    """Invariant densities of spin s modulo total derivatives.
    The derivative matrix is integer (N-free): one reduction for all N."""

    _cache = {}

    def __new__(cls, s):
        if s in cls._cache:
            return cls._cache[s]
        self = super().__new__(cls)
        self.s = s
        self.basis = gen_inv_basis(s)
        self.index = {m: i for i, m in enumerate(self.basis)}
        cols = []
        for b in gen_inv_basis(s - 1):
            col = [F(0)] * len(self.basis)
            for mono, c in d_inv_mono(b).items():
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


def nullspace_cols(cols, nrows):
    ncols = len(cols)
    rows = [[cols[c][r] for c in range(ncols)] for r in range(nrows)]
    rows, pivots = _rref(rows, ncols)
    free = [c for c in range(ncols) if c not in pivots]
    basis = []
    for fc in free:
        v = [F(0)] * ncols
        v[fc] = F(1)
        for i, pc in enumerate(pivots):
            v[pc] = -rows[i][fc]
        basis.append(v)
    return basis


def genuine_kernel_inv(P_low, cand_basis, Nval):
    """Genuine charges (mod d) in span(cand_basis) commuting with
    oint P_low, at numeric N = Nval.  Returns list of currents."""
    s_low = spin_inv(next(iter(P_low)))
    s_cand = spin_inv(cand_basis[0])
    tgt = InvSector(s_low + s_cand - 1)
    cs = InvSector(s_cand)
    cols = [tgt.reduce_mod_d(tgt.vec(
        residue_inv_cur(P_low, {b: F(1)}, Nval))) for b in cand_basis]
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


# ------------------------------------------------------------- families

def fam_current(coeff_pairs):
    """Build a current from [(coeff, [pairs...]), ...] with accumulation
    (guards against monomial collisions)."""
    out = {}
    for c, pairs in coeff_pairs:
        m = canon(tuple(canon_pair(*p) for p in pairs))
        out[m] = out.get(m, 0) + c
    return {m: c for m, c in out.items() if c}


def kdv_P4():
    return fam_current([(F(-2), [(2, 2)]), (F(1), [(1, 1), (1, 1)])])


def kdv_P6(N):
    return fam_current([
        (F(56) / 18 + F(N, 18), [(3, 3)]),
        (F(-6), [(2, 2), (1, 1)]),
        (F(-2, 3) * (20 + N), [(1, 2), (1, 2)]),
        (F(1), [(1, 1), (1, 1), (1, 1)]),
    ])


def oN_P4(N):
    return fam_current([(F(-2, 3) * (2 + N), [(2, 2)]),
                        (F(1), [(1, 1), (1, 1)])])


def oN_P6(N):
    return fam_current([
        (F(1, 150) * (4 + N) * (59 + 36 * N), [(3, 3)]),
        (F(-6, 5) * (4 + N), [(2, 2), (1, 1)]),
        (F(-14, 5) * (4 + N), [(1, 2), (1, 2)]),
        (F(1), [(1, 1), (1, 1), (1, 1)]),
    ])


def curve_point(t):
    """Rational point of s^2 = (N-25)(N-1):
    N = (t^2-25)/(t^2-1), s = -24 t/(t^2-1).
    t = 0 -> N = 25 (branch point); t -> inf -> N = 1 (branch point);
    deck transformation s -> -s is t -> -t."""
    t = F(t)
    N = (t * t - 25) / (t * t - 1)
    s = F(-24) * t / (t * t - 1)
    assert s * s == (N - 25) * (N - 1)
    return N, s


def a2_P6(N, s):
    """A_2^(2) family I_5 density, eq. (34); (N, s) on the double cover."""
    assert s * s == (N - 25) * (N - 1)
    return fam_current([
        (F(1, 96) * (107 + 17 * N + 15 * s), [(3, 3)]),
        (F(-6), [(2, 2), (1, 1)]),
        (F(-1, 8) * (-85 + 17 * N + 15 * s), [(1, 2), (1, 2)]),
        (F(1), [(1, 1), (1, 1), (1, 1)]),
    ])


def fmt_inv(cur):
    parts = []
    for mono in sorted(cur):
        pairs = "".join(f"({a}{b})" for a, b in mono)
        parts.append(f"{cur[mono]} * {pairs}")
    return "  +  ".join(parts)
