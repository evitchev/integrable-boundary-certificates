"""The family matrix M(t) of ad I_3 on spin-(W-1) charge classes, exact over Q(t) (used by the certificate
cyl_stratification_w1014; written for an independent open-tool recomputation, 2026-10-29):

  * domain = the non-pivot monomials of MixedSector(W) (a basis of the classes V_W / d V_{W-1});
  * target = the non-pivot coordinates of MixedSector(W+3) (reduce_mod_d zeroes every pivot coordinate);
  * entry  = class coordinate of reduce_mod_d( Res_{z->w} I_3-density(z) b(w) ), summed over the I_3 monomials, on
             N = (t^2-25)/(t^2-1), s = -24t/(t^2-1).
I_3 is the engine's own cyl_P4(sol, N, s) evaluated with N, s exact elements of Q(t) (class QT), so no coefficient is transcribed.
Residues are the engine's residue_mixed_mono ({mono: {N_power: Fraction}}). Since total derivatives commute with oint I_3,
dim ker M(t0) = len(genuine_kernel_mixed(...)) at every t0 != +-1.

Speed devices (the original run compared every entry with the engine's numeric pipeline at 84 points; the certificate pins
the rebuilt matrices by sha256 and re-runs check_reduction):
  * the reduction mod total derivatives uses the engine's RREF rows (_img_rows, _img_pivots) in one pass,
    red(v) = v - sum_p v[p] row_p, restricted to the class coordinates.  This equals MixedSector.reduce_mod_d because the rows
    are in reduced row echelon form (row_p[p] = 1, row_p[q] = 0 at every other pivot q); target_reduction() asserts that, and
    check_reduction() compares the two on the actual residue vectors;
  * entries are assembled as integer polynomials directly: c_m(t) N(t)^q = num_m(t) (t^2-25)^q / (t^2-1)^(e_m+q);
  * the matrix is stored sparse.
Clearing (as before): M_poly = (t^2-1)^k * c0 * M with the least k making every entry a polynomial and c0 the least positive
integer making every coefficient an integer; (t^2-1)^k c0 vanishes only at t = +-1 (excluded).
"""
import hashlib
import json
import os
import sys
from fractions import Fraction as F
from math import lcm

from sympy import QQ
from sympy.polys.fields import field

import mixed_engine as me           # noqa: E402  (code/, the directory of this file)

K, T = field("t", QQ)
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def lift(x):
    """int / Fraction / QT -> element of K = Q(t)."""
    if isinstance(x, QT):
        return x.v
    if isinstance(x, bool):
        raise TypeError("bool")
    if isinstance(x, int):
        return K(x)
    if isinstance(x, F):
        return K(QQ(x.numerator, x.denominator))
    raise TypeError(f"cannot lift {type(x)}")


class QT:
    """Exact element of Q(t), with the arithmetic the engine's cyl_P4 uses (int and Fraction operands)."""
    __slots__ = ("v",)

    def __init__(self, v):
        self.v = v

    def __add__(self, o):
        return QT(self.v + lift(o))

    __radd__ = __add__

    def __sub__(self, o):
        return QT(self.v - lift(o))

    def __rsub__(self, o):
        return QT(lift(o) - self.v)

    def __mul__(self, o):
        return QT(self.v * lift(o))

    __rmul__ = __mul__

    def __truediv__(self, o):
        return QT(self.v / lift(o))

    def __rtruediv__(self, o):
        return QT(lift(o) / self.v)

    def __neg__(self):
        return QT(-self.v)

    def __eq__(self, o):
        d = self.v - lift(o)
        return d.numer == 0

    def __bool__(self):
        return self.v.numer != 0

    __hash__ = None


def curve_N_s():
    N = (T**2 - 25) / (T**2 - 1)
    s = -24 * T / (T**2 - 1)
    return N, s


def p4_family(sol):
    """The engine's I_3 density of sheet `sol` with exact N(t), s(t): {monomial: element of K}."""
    N, s = curve_N_s()
    P4 = me.cyl_P4(sol, QT(N), QT(s))      # the engine's assert s*s == (N-25)(N-1) runs exactly in Q(t)
    return {m: lift(c) for m, c in P4.items()}


# ----------------------------------------------------------------- polynomials (Fraction or int lists, low -> high)

Q1 = [-1, 0, 1]      # t^2 - 1
Q25 = [-25, 0, 1]    # t^2 - 25


def ptrim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def padd(a, b):
    n = max(len(a), len(b))
    return ptrim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])


def pmul(a, b):
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return ptrim(out)


def ppow(a, e):
    out = [1]
    for _ in range(e):
        out = pmul(out, a)
    return out


def pdivmod(a, b):
    """Quotient and remainder of polynomials over Q."""
    a = [F(x) for x in a]
    q = [F(0)] * max(len(a) - len(b) + 1, 0)
    lb = F(b[-1])
    for i in range(len(a) - len(b), -1, -1):
        c = a[i + len(b) - 1] / lb
        q[i] = c
        for j, y in enumerate(b):
            a[i + j] -= c * y
    return ptrim(q), ptrim(a)


def pdivexact(a, b):
    """Exact division a / b of polynomials over Q (asserts zero remainder)."""
    q, r = pdivmod(a, b)
    assert not r, "inexact division"
    return q


def peval(a, t0):
    acc = 0
    for c in reversed(a):
        acc = acc * t0 + c
    return acc


def kpoly_coeffs(p):
    """sympy PolyElement over QQ in t -> Fraction list low -> high."""
    if p == 0:
        return []
    out = [F(0)] * (p.degree() + 1)
    for (d,), c in p.terms():
        out[d] = F(int(c.numerator), int(c.denominator))
    return out


def coeff_split(c):
    """c in K whose reduced denominator divides (t^2-1)^e -> (least such e, numerator list) with c = num / (t^2-1)^e."""
    num = kpoly_coeffs(c.numer)
    den = kpoly_coeffs(c.denom)
    e = 0
    while True:
        qq, rr = pdivmod(ppow(Q1, e), den)
        if not rr:
            break
        e += 1
        assert e <= 10, "denominator is not a divisor of a power of t^2 - 1"
    return e, pmul(num, qq)


# ----------------------------------------------------------------- sectors

def domain_target(W):
    cs = me.MixedSector(W)
    tgt = me.MixedSector(W + 3)
    cpiv, tpiv = set(cs._img_pivots), set(tgt._img_pivots)
    cols = [j for j in range(len(cs.basis)) if j not in cpiv]
    rows = [i for i in range(len(tgt.basis)) if i not in tpiv]
    return cs, tgt, cols, rows


def target_reduction(tgt, rows):
    """{pivot p: [(class row position, row_p[i]) for nonzero i in rows]} from the engine's RREF rows (asserted RREF)."""
    rowpos = {i: ii for ii, i in enumerate(rows)}
    piv = list(tgt._img_pivots)
    red = {}
    for row, p in zip(tgt._img_rows, piv):
        assert row[p] == 1, "pivot entry != 1"
        assert all(row[q] == 0 for q in piv if q != p), "not reduced at another pivot"
        red[p] = [(rowpos[i], row[i]) for i in rows if row[i]]
    return rowpos, red


def reduce_sparse(v, rowpos, red):
    """Class coordinates of reduce_mod_d(v) for a sparse vector v = {target basis index: Fraction}."""
    out = {}
    for i, c in v.items():
        ii = rowpos.get(i)
        if ii is not None:
            out[ii] = out.get(ii, 0) + c
    for p, c in v.items():
        rp = red.get(p)
        if rp is not None:
            for ii, x in rp:
                out[ii] = out.get(ii, 0) - c * x
    return {ii: c for ii, c in out.items() if c}


def residue_columns(P4keys, W, cs, tgt, cols, rowpos, red):
    """R[(m, q)] = {(ii, jj): Fraction}: class coordinates of the N^q part of Res(m(z) b_jj(w)), reduced mod d."""
    R = {}
    for jj, j in enumerate(cols):
        b = cs.basis[j]
        for m in P4keys:
            res = me.residue_mixed_mono(m, b)
            byq = {}
            for mono, d in res.items():
                for q, c in d.items():
                    byq.setdefault(q, {})
                    idx = tgt.index[mono]
                    byq[q][idx] = byq[q].get(idx, 0) + c
            for q, v in byq.items():
                for ii, c in reduce_sparse(v, rowpos, red).items():
                    R.setdefault((m, q), {})[(ii, jj)] = c
    return R


def check_reduction(P4keys, W, cs, tgt, cols, rows, rowpos, red, ncols=None):
    """Compare reduce_sparse with the engine's reduce_mod_d on the residue vectors of the first ncols domain columns."""
    nchk = 0
    for j in cols[: (ncols if ncols is not None else len(cols))]:
        b = cs.basis[j]
        for m in P4keys:
            res = me.residue_mixed_mono(m, b)
            qs = sorted({q for d in res.values() for q in d})
            for q in qs:
                v = [F(0)] * len(tgt.basis)
                vs = {}
                for mono, d in res.items():
                    if q in d:
                        v[tgt.index[mono]] += d[q]
                        vs[tgt.index[mono]] = vs.get(tgt.index[mono], 0) + d[q]
                full = tgt.reduce_mod_d(v)
                assert not any(full[p] for p in tgt._img_pivots)
                fast = reduce_sparse(vs, rowpos, red)
                assert all(full[i] == fast.get(ii, 0) for ii, i in enumerate(rows)), f"reduction mismatch at {b} {m} {q}"
                nchk += 1
    return nchk


def assemble(P4, R):
    """Integer polynomial entries of M_poly = (t^2-1)^k c0 M from the I_3 coefficients and the residue parts R[(m, q)]."""
    split = {m: coeff_split(c) for m, c in P4.items()}
    used = sorted(R)
    Kmax = max(split[m][0] + q for m, q in used)
    fac = {}
    for m, q in used:
        e, num = split[m]
        fac[(m, q)] = pmul(pmul(num, ppow(Q25, q)), ppow(Q1, Kmax - e - q))
    ent = {}
    for mq in used:
        f = fac[mq]
        for key, c in R[mq].items():
            ent[key] = padd(ent.get(key, []), [c * x for x in f])
    ent = {key: p for key, p in ent.items() if p}
    k = Kmax
    while k > 0 and all(peval(p, 1) == 0 and peval(p, -1) == 0 for p in ent.values()):
        ent = {key: pdivexact(p, Q1) for key, p in ent.items()}
        k -= 1
    c0 = 1
    for p in ent.values():
        for x in p:
            c0 = lcm(c0, F(x).denominator)
    out = {}
    for key, p in ent.items():
        ip = [F(x) * c0 for x in p]
        assert all(x.denominator == 1 for x in ip)
        out[key] = [int(x) for x in ip]
    return k, c0, out, split


def mono_str(m):
    xs, yp = m
    sx = "X(" + ",".join(map(str, xs)) + ")" if xs else ""
    sy = "".join(f"Y({a}{b})" for a, b in yp)
    return (sx + sy) or "1"


def save(path, sheet, W, k, c0, m, n, ent, rowlab, collab, note=""):
    obj = {
        "sheet": sheet, "W": W, "k": k, "c0": c0, "scalar": f"M_poly = (t^2-1)^{k} * {c0} * M",
        "m": m, "n": n, "rows": [mono_str(x) for x in rowlab], "cols": [mono_str(x) for x in collab],
        "entries": [[ii, jj, ent[(ii, jj)]] for (ii, jj) in sorted(ent)], "note": note,
    }
    with open(path, "w") as fh:
        json.dump(obj, fh, separators=(",", ":"))
        fh.write("\n")


def load(path, expect_sha=None):
    raw = open(path, "rb").read()
    if expect_sha is not None:
        h = hashlib.sha256(raw).hexdigest()
        assert h == expect_sha, f"sha256 mismatch for {path}: {h}"
    return json.loads(raw)


def dense(d):
    """m x n list of integer coefficient lists (low -> high); [] = 0."""
    M = [[[] for _ in range(d["n"])] for _ in range(d["m"])]
    for ii, jj, p in d["entries"]:
        M[ii][jj] = p
    return M
