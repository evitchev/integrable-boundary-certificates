"""stratification_w6_w8 S1: the family matrix M(t) of ad I_3 on spin-W charge classes, exact over Q(t).

Construction (mine; only the engine's Wick, reduction and I_3 functions are imported):
  * domain  = spin-W classes (densities mod total derivatives): the non-pivot monomials of MixedSector(W), a basis of the
              quotient V_W / d V_{W-1} (the engine's RREF of the derivative image has those pivots);
  * target  = spin-(W+3) classes: the non-pivot coordinates of MixedSector(W+3) (reduce_mod_d zeroes every pivot coordinate);
  * entry   = coordinate i of reduce_mod_d( Res_{z->w} I_3-density(z) b(w) ), summed over the I_3 monomials, with
              N = (t^2-25)/(t^2-1), s = -24t/(t^2-1).
I_3 is the engine's own cyl_P4(sol, N, s), evaluated with N, s as exact elements of Q(t) (class QT below), so no coefficient
is transcribed. Residues are the engine's residue_mixed_mono (exact polynomials in N, by number of index cycles); the
reduction mod total derivatives is linear and N-free, so each N-power is reduced separately and then recombined in Q(t).
Since total derivatives commute with oint I_3 (Res(P(z) dQ(w)) = d_w Res(P(z) Q(w))), the restriction to the non-pivot
monomials is the induced map on classes, and dim ker M(t0) = len(genuine_kernel_mixed(...)) at every t0 != +-1.

Clearing: M_poly = (t^2-1)^k * c0 * M with the least k making every entry a polynomial and c0 a positive integer making
every coefficient an integer; (t^2-1)^k * c0 vanishes only at t = +-1 (excluded).
"""
import json
import os
import sys
from fractions import Fraction as F
from math import lcm

from sympy import QQ
from sympy.polys.fields import field

# the repository's code/ directory (this file sits in <repository>/results/lab/jumping_loci_records/<name>/code/)
REPO_CODE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "..", "code"))
sys.path.append(REPO_CODE)          # appended: stdlib and site-packages keep precedence
import mixed_engine as me           # noqa: E402

K, T = field("t", QQ)
BASE = os.getcwd()                  # output root (data/, results/): the current directory; run.sh runs in its new output directory


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


def domain_target(W):
    cs = me.MixedSector(W)
    tgt = me.MixedSector(W + 3)
    cpiv, tpiv = set(cs._img_pivots), set(tgt._img_pivots)
    cols = [j for j in range(len(cs.basis)) if j not in cpiv]
    rows = [i for i in range(len(tgt.basis)) if i not in tpiv]
    return cs, tgt, cols, rows


def build_M(P4, W):
    """M(t) over K: rows = target class coordinates, columns = domain class monomials."""
    cs, tgt, cols, rows = domain_target(W)
    N, _ = curve_N_s()
    tpiv = tgt._img_pivots
    M = [[K(0) for _ in cols] for _ in rows]
    for jj, j in enumerate(cols):
        b = cs.basis[j]
        for m, c in P4.items():
            res = me.residue_mixed_mono(m, b)
            powers = sorted({p for d in res.values() for p in d})
            for p in powers:
                v = [F(0)] * len(tgt.basis)
                for mono, d in res.items():
                    if p in d:
                        v[tgt.index[mono]] += d[p]
                red = tgt.reduce_mod_d(v)
                assert not any(red[i] for i in tpiv), "reduce_mod_d left a pivot coordinate"
                fac = c * N**p
                for ii, i in enumerate(rows):
                    if red[i]:
                        M[ii][jj] += fac * lift(red[i])
    return M, [tgt.basis[i] for i in rows], [cs.basis[j] for j in cols]


def clear(M):
    """Return (k, c0, Mpoly) with Mpoly = (t^2-1)^k * c0 * M, Mpoly entries as integer coefficient lists (low -> high)."""
    q = T**2 - 1
    k = 0
    while True:
        ok = all((e * q**k).denom.degree() <= 0 for row in M for e in row)
        if ok:
            break
        k += 1
        assert k <= 20, "denominators are not powers of t^2 - 1"
    polys = []
    den_lcm = 1
    for row in M:
        prow = []
        for e in row:
            g = e * q**k
            num, den = g.numer, g.denom          # den is a nonzero constant here
            c_den = den.LC
            coeffs = [F(0)] * (num.degree() + 1) if num != 0 else []
            for (deg,), c in num.terms():
                coeffs[deg] = F(int(c.numerator), int(c.denominator)) / F(int(c_den.numerator), int(c_den.denominator))
            for c in coeffs:
                den_lcm = lcm(den_lcm, c.denominator)
            prow.append(coeffs)
        polys.append(prow)
    c0 = den_lcm
    Mpoly = [[[int(c * c0) for c in e] for e in row] for row in polys]
    for row, prow in zip(Mpoly, polys):
        for e, pe in zip(row, prow):
            assert all(F(a) == b * c0 for a, b in zip(e, pe))
    return k, c0, Mpoly


def ev(poly, t0):
    """Horner evaluation of an integer coefficient list (low -> high) at a Fraction t0."""
    acc = F(0)
    for c in reversed(poly):
        acc = acc * t0 + c
    return acc


def ev_matrix(Mpoly, t0):
    return [[ev(e, t0) for e in row] for row in Mpoly]


def mono_str(m):
    xs, yp = m
    sx = "X(" + ",".join(map(str, xs)) + ")" if xs else ""
    sy = "".join(f"Y({a}{b})" for a, b in yp)
    return (sx + sy) or "1"


def save(path, sheet, W, k, c0, Mpoly, rowlab, collab, note=""):
    obj = {
        "sheet": sheet, "W": W, "k": k, "c0": c0,
        "scalar": f"M_poly = (t^2-1)^{k} * {c0} * M",
        "rows": [mono_str(m) for m in rowlab], "cols": [mono_str(m) for m in collab],
        "Mpoly": Mpoly, "note": note,
    }
    with open(path, "w") as fh:
        json.dump(obj, fh, separators=(",", ":"))
        fh.write("\n")


def load(path):
    with open(path) as fh:
        return json.load(fh)
