"""Exact OPE/commutator engine for chiral free bosons.

Conventions follow Vitchev, hep-th/0404195:
  - fields X^j(v), letters are d^m X^j represented as pairs (m, j), m >= 1;
  - a normal-ordered monomial current is a sorted tuple of letters;
  - a current is a dict {monomial: Fraction coefficient};
  - Wick pairing (paper's eq. 17):
        <dX^i(z) dX^j(w)> = delta^{ij} (z-w)^{-2}
    which propagates to
        <d^m X^i(z) d^n X^j(w)> = delta^{ij} (-1)^(m-1) (m+n-1)! (z-w)^{-(m+n)}.

A charge I = oint P dz vanishes iff its density P is a total derivative, so
the commutator [oint P, oint Q] = oint Res_{z->w} P(z)Q(w) vanishes iff
Res P(z)Q(w) lies in the image of d (paper's eqs. 10-16).

All arithmetic is exact (fractions.Fraction).  Rational functions of the
paperclip parameter n are handled by evaluating at sample rational values
of n and, where needed, verifying interpolated results symbolically with
sympy (numeric coefficient matrices keep that cheap).
"""

from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import factorial


# ---------------------------------------------------------------- currents

def cur_add(a, b, coeff=1):
    """a += coeff * b, in place; returns a."""
    for mono, c in b.items():
        v = a.get(mono, 0) + coeff * c
        if v:
            a[mono] = v
        else:
            a.pop(mono, None)
    return a


def cur_scale(a, coeff):
    return {m: coeff * c for m, c in a.items() if coeff * c}


def d_mono(mono):
    """Derivative of a monomial current, as a current."""
    out = {}
    for i in range(len(mono)):
        m, j = mono[i]
        new = tuple(sorted(mono[:i] + ((m + 1, j),) + mono[i + 1:]))
        out[new] = out.get(new, 0) + 1
    return {m: c for m, c in out.items() if c}


def d_cur(cur):
    out = {}
    for mono, c in cur.items():
        cur_add(out, d_mono(mono), c)
    return out


def spin(mono):
    return sum(m for m, _ in mono)


# ------------------------------------------------------------------ bases

def gen_basis(s, nfields, z2=True):
    """All monomials of total spin s in `nfields` fields; if z2, keep only
    those with an even number of letters of each field (Z2^nfields sector).
    Letters (m, j) are generated in non-increasing order to avoid dups."""
    out = []

    def rec(remaining, maxletter, acc):
        if remaining == 0:
            if not z2 or all(
                sum(1 for _, j in acc if j == f) % 2 == 0
                for f in range(nfields)
            ):
                out.append(tuple(sorted(acc)))
            return
        m0, j0 = maxletter
        for m in range(min(remaining, m0), 0, -1):
            jmax = j0 if m == m0 else nfields - 1
            for j in range(jmax, -1, -1):
                rec(remaining - m, (m, j), acc + [(m, j)])

    rec(s, (s, nfields - 1), [])
    return out


# -------------------------------------------------------- residue of OPE

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
def residue_mono(P, Q):
    """Res_{z->w} of P(z) Q(w) for monomials P, Q, as a current (dict).

    Sum over all nonempty partial matchings of letter *positions* of P
    to letter positions of Q (each occurrence contracts individually);
    uncontracted P letters are Taylor-expanded around w to soak up the
    pole (z-w)^{-D} down to (z-w)^{-1}.
    """
    out = {}
    np_, nq = len(P), len(Q)
    for k in range(1, min(np_, nq) + 1):
        for psel in combinations(range(np_), k):
            for qsel in permutations(range(nq), k):
                num, den, D = 1, 1, 0
                for a, b in zip(psel, qsel):
                    m, i = P[a]
                    n, j = Q[b]
                    if i != j:
                        num = 0
                        break
                    num *= (-1) ** (m - 1) * factorial(m + n - 1)
                    D += m + n
                if not num:
                    continue
                rest_p = [P[a] for a in range(np_) if a not in psel]
                rest_q = tuple(Q[b] for b in range(nq) if b not in set(qsel))
                r = len(rest_p)
                if r == 0:
                    continue  # D >= 2 always, no (z-w)^{D-1} available
                for ts in _compositions(D - 1, r):
                    dd = den
                    for t in ts:
                        dd *= factorial(t)
                    mono = tuple(sorted(
                        rest_q
                        + tuple((mu + t, j) for (mu, j), t in zip(rest_p, ts))
                    ))
                    c = out.get(mono, 0) + Fraction(num, dd)
                    if c:
                        out[mono] = c
                    else:
                        out.pop(mono, None)
    return out


def residue_cur(P, Q):
    """Residue of the OPE of two currents (dicts)."""
    out = {}
    for mp, cp in P.items():
        for mq, cq in Q.items():
            r = residue_mono(mp, mq)
            if r:
                cur_add(out, r, cp * cq)
    return out


# ------------------------------------------------- exact linear algebra

def rref(rows, ncols):
    """Row-reduce a list of Fraction rows in place; return (rows, pivots)."""
    pivots = []
    rank = 0
    for col in range(ncols):
        piv = next((r for r in range(rank, len(rows)) if rows[r][col]), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        inv = Fraction(1, 1) / rows[rank][col]
        rows[rank] = [x * inv for x in rows[rank]]
        for r in range(len(rows)):
            if r != rank and rows[r][col]:
                f = rows[r][col]
                rows[r] = [a - f * b for a, b in zip(rows[r], rows[rank])]
        pivots.append(col)
        rank += 1
    return rows, pivots


def nullspace(matrix_cols, nrows):
    """Nullspace basis of the matrix whose COLUMNS are `matrix_cols`
    (each a list of length nrows).  Returns list of coefficient vectors."""
    ncols = len(matrix_cols)
    rows = [[matrix_cols[c][r] for c in range(ncols)] for r in range(nrows)]
    rows, pivots = rref(rows, ncols)
    free = [c for c in range(ncols) if c not in pivots]
    basis = []
    for fc in free:
        v = [Fraction(0)] * ncols
        v[fc] = Fraction(1)
        for i, pc in enumerate(pivots):
            v[pc] = -rows[i][fc]
        basis.append(v)
    return basis


class Sector:
    """Linear-algebra helper for the space B^s of densities of spin s
    (nfields free bosons, optional Z2^nfields symmetry) with its
    total-derivative subspace d(B^{s-1})."""

    def __init__(self, s, nfields, z2=True):
        self.s = s
        self.basis = gen_basis(s, nfields, z2)
        self.index = {m: i for i, m in enumerate(self.basis)}
        lower = gen_basis(s - 1, nfields, z2=False)
        self.dcols = []
        for b in lower:
            col = [Fraction(0)] * len(self.basis)
            for mono, c in d_mono(b).items():
                if mono in self.index:  # derivative preserves parity sector
                    col[self.index[mono]] = Fraction(c)
            if any(col):
                self.dcols.append(col)
        # RREF of the image of d, kept as reduced row vectors for membership
        rows = [list(c) for c in self.dcols]
        self._img_rows, self._img_pivots = rref(rows, len(self.basis))
        self._img_rows = self._img_rows[:len(self._img_pivots)]

    def vec(self, cur):
        v = [Fraction(0)] * len(self.basis)
        for mono, c in cur.items():
            if mono not in self.index:
                raise ValueError(f"monomial {mono} outside sector spin={self.s}")
            v[self.index[mono]] = Fraction(c)
        return v

    def reduce_mod_d(self, v):
        """Reduce a vector modulo the image of d (v untouched)."""
        v = list(v)
        for row, p in zip(self._img_rows, self._img_pivots):
            if v[p]:
                f = v[p]
                v = [a - f * b for a, b in zip(v, row)]
        return v

    def is_total_derivative(self, cur):
        return not any(self.reduce_mod_d(self.vec(cur)))

    def q_dim(self):
        """dim Q^{s-1} for this symmetry sector: densities mod d."""
        return len(self.basis) - len(self._img_pivots)


def commuting_kernel(P_low, cand_basis, nfields, z2=True):
    """Charges oint J with J in span(cand_basis) commuting with oint P_low.

    Returns (kernel_vectors, n_trivial): coefficient vectors c such that
    Res(P_low, J(c)) is a total derivative, and the number of those that
    are themselves total derivatives (trivial charges).  The count of
    genuine commuting charges is len(kernel_vectors) - n_trivial.
    """
    s_low = spin(next(iter(P_low)))
    s_cand = spin(cand_basis[0])
    target = Sector(s_low + s_cand - 1, nfields, z2)
    cols = []
    for b in cand_basis:
        r = residue_cur(P_low, {b: Fraction(1)})
        cols.append(target.reduce_mod_d(target.vec(r)))
    kernel = nullspace(cols, len(target.basis))
    cand_sec = Sector(s_cand, nfields, z2)
    n_trivial = 0
    if kernel:
        triv_cols = []
        for v in kernel:
            cur = {b: c for b, c in zip(cand_basis, v) if c}
            triv_cols.append(cand_sec.reduce_mod_d(cand_sec.vec(cur)))
        # trivial directions = nullspace of the reduced vectors
        n_trivial = len(nullspace(triv_cols, len(cand_sec.basis)))
    return kernel, n_trivial


# ------------------------------------------------------- pretty printing

def fmt_mono(mono, nfields=2):
    """Vitchev's notation: (n1 n2 ... | m1 m2 ...) for two fields,
    (n1 n2 ...) for one field."""
    if nfields == 1:
        return "(" + "".join(str(m) for m, _ in mono) + ")"
    xs = "".join(str(m) for m, j in mono if j == 0)
    ys = "".join(str(m) for m, j in mono if j == 1)
    return f"({xs}|{ys})"


def fmt_cur(cur, nfields=2):
    parts = []
    for mono in sorted(cur):
        parts.append(f"{cur[mono]} * {fmt_mono(mono, nfields)}")
    return "\n".join(parts)
