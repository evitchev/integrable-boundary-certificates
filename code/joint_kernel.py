"""Joint first-order-conservation kernel over Q(sqrt2).

Shared machinery (extracted 2026-08-18 from c1_four_condition_char.py
via kdv_tower_classes.py so certificates can import it without
executing another certificate).  joint(w, moms, cand_z2) returns the
Q-dimension of the joint kernel of the first-order conservation
conditions at density weight w for the momentum vectors moms
(entries in Q(sqrt2)); a K-line over Q(sqrt2) counts as Q-dim 2.
cand_z2=True restricts CANDIDATES to the Z2-even sector."""
from fractions import Fraction as F

import sympy as sp

from ope_engine import d_mono, gen_basis
from screening_fit import residue_screen, twisted_d
from sparse_linalg import SparsePivots, lcm

SQ2_ = sp.sqrt(2)
NB = 2


def split2(expr):
    e = sp.expand(expr)
    c1, c2 = e.coeff(SQ2_, 0), e.coeff(SQ2_, 1)
    if sp.simplify(e - c1 - c2 * SQ2_) != 0:
        raise ValueError(f"split reconstruction failed: {e}")
    r1, r2 = sp.Rational(c1), sp.Rational(c2)
    return (F(int(r1.p), int(r1.q)), F(int(r2.p), int(r2.q)))


def mul2(pair, bi):
    c1, c2 = pair
    return (c1, c2) if bi == 0 else (2 * c2, c1)


def int_vec(fr):
    Lc = 1
    for v in fr.values():
        Lc = lcm(Lc, v.denominator)
    return {i: int(v * Lc) for i, v in fr.items() if v}, Lc


def joint(w, moms, cand_z2=False):
    basisw = gen_basis(w, 1, z2=cand_z2)
    tgt = gen_basis(w - 1, 1, z2=False)
    tidx = {m: i for i, m in enumerate(tgt)}
    lower = [()] if w == 2 else gen_basis(w - 2, 1, z2=False)
    nt = len(tgt)
    R = [[{tidx[m]: split2(c)
           for m, c in residue_screen({b: F(1)}, al).items()}
          for b in basisw] for al in moms]
    T = [[{tidx[m]: split2(c) for m, c in twisted_d(u, al, 1).items()}
          for u in lower] for al in moms]
    cols, scales = [], []
    nx = len(basisw) * NB
    for j in range(len(basisw)):
        for bi in range(NB):
            col = {}
            for smom in range(len(moms)):
                base = smom * nt * NB
                for r, e in R[smom][j].items():
                    for beta, cv in enumerate(mul2(e, bi)):
                        if cv:
                            col[base + r * NB + beta] = \
                                col.get(base + r * NB + beta, F(0)) + cv
            iv, Lc = int_vec(col)
            cols.append(iv)
            scales.append(Lc)
    for smom in range(len(moms)):
        base = smom * nt * NB
        for u in range(len(lower)):
            for bi in range(NB):
                col = {}
                for r, e in T[smom][u].items():
                    for beta, cv in enumerate(mul2(e, bi)):
                        if cv:
                            col[base + r * NB + beta] = \
                                col.get(base + r * NB + beta, F(0)) + cv
                iv, Lc = int_vec(col)
                cols.append(iv)
                scales.append(Lc)
    nc = len(cols)
    piv = SparsePivots()
    kers = []
    for jc, col in enumerate(cols):
        aug = {i + nc: v for i, v in col.items()}
        aug[jc] = 1
        r = piv.reduce(aug)
        if r and max(r) < nc:
            kers.append(r)
        elif r:
            piv.rows[max(r)] = r
    vecs = []
    for kv in kers:
        xv = {c: co * scales[c] for c, co in kv.items()
              if c < nx and co * scales[c]}
        if xv:
            vecs.append(xv)
    fullb = gen_basis(w, 1, z2=False)
    fidx = {m: i for i, m in enumerate(fullb)}
    dp = SparsePivots()
    for bm in gen_basis(w - 1, 1, z2=False):
        col = {}
        for m, c in d_mono(bm).items():
            col[fidx[m]] = col.get(fidx[m], 0) + c
        col = {i: c for i, c in col.items() if c}
        if col:
            for beta in range(NB):
                dp.insert({i * NB + beta: c for i, c in col.items()})
    qd = 0
    for v in vecs:
        mapped = {}
        for c, co in v.items():
            j, beta = divmod(c, NB)
            mapped[fidx[basisw[j]] * NB + beta] = co
        if mapped and dp.insert(mapped):
            qd += 1
    return qd
