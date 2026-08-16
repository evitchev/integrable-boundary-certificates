"""Hunt for screening operators annihilating the commuting families.

A screening charge Q_alpha = oint V_alpha, V_alpha = :e^{alpha.X}:,
commutes with a charge oint P iff Res_{z->w} P(z) V_alpha(w) is a
*twisted* total derivative:  Res = d_alpha(J) V_alpha with
d_alpha J := dJ + (alpha . dX) J   (Feigin-Frenkel screening formalism).

With the paper's pairing <dX^i dX^j> = delta^{ij}/(z-w)^2 we have
<X^i(z) X^j(w)> = delta^{ij} ln(z-w), so contracting d^m X^j(z) with
e^{alpha.X(w)} gives  alpha_j (-1)^{m-1} (m-1)! (z-w)^{-m}.

Residue of P(z)V(w): sum over nonempty subsets S of P's letters
contracted with the exponential (pole D = sum of their m's), remaining
letters Taylor-expanded around w soaking up (z-w)^{D-1}.

Everything is polynomial in the components of alpha, handled with sympy.

Calibration: for one boson the quantum KdV family must be annihilated by
two screenings; then we scan for the paperclip's screening vectors.
"""

from fractions import Fraction as F
from itertools import combinations
from math import factorial

import sympy as sp

from ope_engine import gen_basis, d_mono, spin, _compositions
from families import L, KDV_P4, KDV_P6, paperclip_P4, paperclip_P6


def residue_screen(P, alpha):
    """Res of P(z) V_alpha(w) as {monomial: sympy coeff} (V implicit).
    P: current dict with Fraction coeffs; alpha: tuple of sympy symbols."""
    out = {}
    for mono, cP in P.items():
        cP = sp.Rational(cP.numerator, cP.denominator) \
            if isinstance(cP, F) else cP
        nlet = len(mono)
        for k in range(1, nlet + 1):
            for sel in combinations(range(nlet), k):
                num = sp.Integer(1)
                D = 0
                for a in sel:
                    m, j = mono[a]
                    num *= alpha[j] * (-1) ** (m - 1) * factorial(m - 1)
                    D += m
                rest = [mono[a] for a in range(nlet) if a not in sel]
                r = len(rest)
                if r == 0:
                    if D == 1:  # bare 1/(z-w): contributes V itself
                        out[()] = out.get((), 0) + cP * num
                    continue
                for ts in _compositions(D - 1, r):
                    den = 1
                    for t in ts:
                        den *= factorial(t)
                    new = tuple(sorted(
                        (mu + t, j) for (mu, j), t in zip(rest, ts)))
                    out[new] = out.get(new, 0) + cP * num / den
    return {m: sp.expand(c) for m, c in out.items() if sp.expand(c) != 0}


def twisted_d(mono, alpha, nfields):
    """d_alpha applied to a monomial (times V): dJ + (alpha.dX) J."""
    out = {m: sp.Integer(c) for m, c in d_mono(mono).items()}
    for j in range(nfields):
        new = tuple(sorted(mono + ((1, j),)))
        out[new] = out.get(new, 0) + alpha[j]
    return out


def screening_conditions(P, alpha, nfields):
    """Polynomial conditions on alpha for Res(P V_alpha) to be a twisted
    total derivative.  Returns a list of sympy polynomials."""
    res = residue_screen(P, alpha)
    s = spin(next(iter(P))) - 1  # polynomial spin of the residue
    tgt = []
    for t in range(0, s + 1):
        tgt += gen_basis(t, nfields, z2=False) if t else [()]
    # twisted image: d_alpha of every monomial of poly-spin s-1
    lower = [()] if s - 1 == 0 else gen_basis(s - 1, nfields, z2=False)
    index = {m: i for i, m in enumerate(tgt)}
    cols = []
    for u in lower:
        col = [sp.Integer(0)] * len(tgt)
        for m, c in twisted_d(u, alpha, nfields).items():
            col[index[m]] = c
        cols.append(col)
    vec = [sp.Integer(0)] * len(tgt)
    for m, c in res.items():
        vec[index[m]] = c
    # Gaussian elimination pivoting ONLY on coefficient columns (never on
    # the augmented one), over QQ(alpha); leftover augmented entries in
    # pivot-free rows are the (generic) screening conditions.
    ncols = len(cols)
    rows = [[cols[c][r] for c in range(ncols)] + [vec[r]]
            for r in range(len(tgt))]
    used = set()
    for c in range(ncols):
        piv, best = None, None
        for r in range(len(rows)):
            if r in used:
                continue
            e = sp.cancel(rows[r][c])
            rows[r][c] = e
            if e != 0:
                if not e.free_symbols:
                    piv = r
                    break
                if best is None:
                    best = r
        if piv is None:
            piv = best
        if piv is None:
            continue
        used.add(piv)
        pe = rows[piv][c]
        for r in range(len(rows)):
            if r == piv:
                continue
            f = sp.cancel(rows[r][c])
            if f != 0:
                rows[r] = [sp.cancel(x - f * y / pe)
                           for x, y in zip(rows[r], rows[piv])]
    conds = []
    for r in range(len(rows)):
        if r in used:
            continue
        if all(sp.cancel(rows[r][c]) == 0 for c in range(ncols)):
            e = sp.cancel(rows[r][ncols])
            if e != 0:
                conds.append(sp.factor(sp.numer(sp.together(e))))
    return conds


if __name__ == "__main__":
    # ---------------- calibration: one boson, quantum KdV --------------
    a = sp.Symbol('a')
    conds = screening_conditions(KDV_P4, (a,), 1)
    sols3 = set()
    for c in conds:
        sols3.update(sp.solve(c, a))
    print(f"KdV I_3 screening condition roots: {sorted(sols3, key=str)}")
    conds5 = screening_conditions(KDV_P6, (a,), 1)
    sols5 = set()
    for c in conds5:
        for root in sp.solve(c, a):
            sols5.add(sp.simplify(root))
    print(f"KdV I_5 screening condition roots: {sorted(sols5, key=str)}")
    both = {s for s in sols3
            if any(sp.simplify(s - t) == 0 for t in sols5)}
    print(f"common roots (screenings of the KdV family): "
          f"{sorted(both, key=str)}")

    # ---------------- paperclip: two bosons ---------------------------
    b = sp.Symbol('b')
    for nval in [F(2), F(3)]:
        c4 = screening_conditions(paperclip_P4(nval), (a, b), 2)
        c6 = screening_conditions(paperclip_P6(nval), (a, b), 2)
        print(f"\nn = {nval}:")
        print(f"  I_3 gives {len(c4)} condition(s), "
              f"I_5 gives {len(c6)} condition(s)")
        sols = sp.solve(c4 + c6, [a, b], dict=True)
        for s in sols:
            print(f"  alpha = {tuple(sp.simplify(s.get(v, v)) for v in (a, b))}")
