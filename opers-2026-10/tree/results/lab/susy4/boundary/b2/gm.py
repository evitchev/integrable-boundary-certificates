"""gm.py -- B2 (Fable seat): free-field super-OPE layer with a DIAGONAL RATIONAL metric.

NEW FILE in the seat's workspace; the registered engines (code/super_engine.py, code/super2_screening.py,
code/super_cyl_vev.py) are imported read-only and never modified (sys.dont_write_bytecode set before import).

Letters as in super_engine: (0, m, j) = d^m X_j (m >= 1), (1, a, j) = d^a psi_j (a >= 0).  Currents = {monomial: coeff}.
Metric eta = tuple of Fractions (or sympy expressions): <dX_j dX_j> = eta_j/(z-w)^2, <psi_j psi_j> = eta_j/(z-w).
Every Wick contraction removes one field-j letter from each side; the emitted monomial's metric factor is
prod_j eta_j^{n_j}, n_j = (#field-j letters in P + in Q - in the output)/2, read off from letter counts
(the same device as super2_screening._metric_factor, generalised from +-1 to any eta).

Also: all-pole OPE coefficients (for V1), T and G with background charge, a joint_kernel with metric-aware
current conditions, and the NS cylinder vev with metric (wrapper over super_cyl_vev.vev_monomial).
"""
import sys
sys.dont_write_bytecode = True
import os as _os_kit; sys.path.insert(0, _os_kit.path.join(_os_kit.path.dirname(_os_kit.path.abspath(__file__)), '..', '..', '..', '..', '..', 'code'))  # KIT PATCH: repo code/ from __file__ (was an absolute path)
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations
from math import factorial

from ope_engine import _compositions
from sparse_linalg import SparsePivots, _to_int_vec
from super_engine import (_koszul, _pairing, canonicalize, cur_add, d_mono, gen_basis_super, is_odd,
                          mono_weight2, normal_prod_mono, residue_mono)
from super2_screening import (Vertex, basis, field_counts, is_twisted_exact, twisted_d, vertex_residue)
import super_cyl_vev as SCV

NF = 2


def X(m, j):
    return (0, m, j)


def PSI(a, j):
    return (1, a, j)


def is_zero(c):
    try:
        return c == 0
    except Exception:
        return False


def clean(cur):
    out = {}
    for m, c in cur.items():
        if hasattr(c, 'expand'):
            c = c.expand()
        if not is_zero(c):
            out[m] = c
    return out


def add(*terms):
    """add((coef, cur), (coef, cur), ...)"""
    out = {}
    for coef, cur in terms:
        for m, c in cur.items():
            out[m] = out.get(m, 0) + coef * c
    return clean(out)


def scale(cur, k):
    return clean({m: k * c for m, c in cur.items()})


def d_cur(cur):
    out = {}
    for m, c in cur.items():
        for mm, cc in d_mono(m).items():
            out[mm] = out.get(mm, 0) + c * cc
    return clean(out)


def counts(mono, nf):
    c = [0] * nf
    for _, _, j in mono:
        c[j] += 1
    return c


def metric_factor(cP, cQ, cO, eta):
    f = 1
    for j, e in enumerate(eta):
        n2 = cP[j] + cQ[j] - cO[j]
        assert n2 % 2 == 0 and n2 >= 0
        n = n2 // 2
        if n:
            f = f * e ** n
    return f


def _apply_mono_op(P, Q, eta, mono_op):
    nf = len(eta)
    out = {}
    for mp, cp in P.items():
        fp = counts(mp, nf)
        for mq, cq in Q.items():
            r = mono_op(mp, mq)
            if not r:
                continue
            fq = counts(mq, nf)
            for mo, c in r.items():
                s = metric_factor(fp, fq, counts(mo, nf), eta)
                out[mo] = out.get(mo, 0) + s * c * cp * cq
    return clean(out)


def residue(P, Q, eta):
    """Res_{z->w} P(z) Q(w)."""
    return _apply_mono_op(P, Q, eta, residue_mono)


def normal_prod(P, Q, eta):
    """Point-split :PQ:(w) (the vertex-algebra (-1)-product)."""
    return _apply_mono_op(P, Q, eta, normal_prod_mono)


@lru_cache(maxsize=200_000)
def ope_coeff_mono(P, Q, n):
    """Coefficient of (z-w)^(-n) in P(z) Q(w) for canonical monomials at metric +1 (n >= 1: poles; n <= 0:
    regular terms, including the contraction-free Taylor term).  Same Wick enumeration as residue_mono with
    Taylor order D - n."""
    out = {}
    np_, nq = len(P), len(Q)
    for k in range(0, min(np_, nq) + 1):
        for psel in combinations(range(np_), k):
            for qsel in permutations(range(nq), k):
                num, D = 1, 0
                for a, b in zip(psel, qsel):
                    pr = _pairing(P[a], Q[b])
                    if pr is None:
                        num = 0
                        break
                    num *= pr[0]
                    D += pr[1]
                if not num or D < n:
                    continue
                if k:
                    num *= _koszul(P, Q, psel, qsel)
                restP = [P[a] for a in range(np_) if a not in psel]
                restQ = [Q[b] for b in range(nq) if b not in set(qsel)]
                r = len(restP)
                if r == 0:
                    if D == n:
                        cm, sg = canonicalize(tuple(restQ))
                        if cm is not None:
                            out[cm] = out.get(cm, 0) + F(num * sg)
                    continue
                for ts in _compositions(D - n, r):
                    den = 1
                    sh = []
                    for (p, m, j), t in zip(restP, ts):
                        den *= factorial(t)
                        sh.append((p, m + t, j))
                    cm, sg = canonicalize(tuple(sh) + tuple(restQ))
                    if cm is None:
                        continue
                    out[cm] = out.get(cm, 0) + F(num * sg, den)
    return {m: c for m, c in out.items() if c}


def ope_coeff(P, Q, n, eta):
    return _apply_mono_op(P, Q, eta, lambda a, b: ope_coeff_mono(a, b, n))


def ope_all(P, Q, eta, maxpole=6):
    return {n: ope_coeff(P, Q, n, eta) for n in range(1, maxpole + 1)}


# ------------------------------------------------------------------ super-Virasoro with background charge

def T_eta(eta, lam=None):
    """T = sum_j (1/eta_j)[1/2 (dX_j)^2 - 1/2 psi_j dpsi_j] + sum_j lam_j d^2 X_j."""
    nf = len(eta)
    lam = lam or [0] * nf
    T = {}
    for j, e in enumerate(eta):
        T[(X(1, j), X(1, j))] = F(1, 2) / e
        T[(PSI(0, j), PSI(1, j))] = F(-1, 2) / e
        if not is_zero(lam[j]):
            T[(X(2, j),)] = lam[j]
    return clean(T)


def G_eta(eta, lam=None):
    """G = sum_j (1/eta_j) psi_j dX_j + 2 sum_j lam_j dpsi_j  (Res(G, G) = 2T asserted in n2frame)."""
    nf = len(eta)
    lam = lam or [0] * nf
    G = {}
    for j, e in enumerate(eta):
        G[(X(1, j), PSI(0, j))] = F(1) / e
        if not is_zero(lam[j]):
            G[(PSI(1, j),)] = 2 * lam[j]
    return clean(G)


# ------------------------------------------------------------------ total derivatives, screening, kernels

def is_total_derivative(cur, nf=NF):
    if not cur:
        return True
    w2 = mono_weight2(next(iter(cur)))
    index = {m: i for i, m in enumerate(basis(w2, nf))}
    piv = SparsePivots()
    for u in basis(w2 - 2, nf):
        iv, _ = _to_int_vec(d_mono(u), index)
        if iv:
            piv.insert(iv)
    iv, _ = _to_int_vec(cur, index)
    return not piv.reduce(iv)


def screens(P, V, eta):
    """[oint V, oint P] = 0 exactly (rational momenta): Res(P, V) twisted-exact."""
    w2P = mono_weight2(next(iter(P)))
    return is_twisted_exact(vertex_residue(P, V, eta), V.b, w2P + V.w2W - 2, len(eta))


# ------------------------------------------------------------------ classical (single-contraction) brackets

@lru_cache(maxsize=200_000)
def residue1_mono(P, Q):
    """Res_{z->w} P(z) Q(w) keeping EXACTLY ONE Wick contraction: the classical Poisson bracket of the local
    functionals (the k-contraction terms carry eta^k; with eta = s eta-hat the k = 1 term is the s -> 0 leading part)."""
    out = {}
    np_, nq = len(P), len(Q)
    for psel in combinations(range(np_), 1):
        for qsel in permutations(range(nq), 1):
            pr = _pairing(P[psel[0]], Q[qsel[0]])
            if pr is None:
                continue
            num, D = pr
            num *= _koszul(P, Q, psel, qsel)
            restP = [P[a] for a in range(np_) if a not in psel]
            restQ = [Q[b] for b in range(nq) if b not in set(qsel)]
            r = len(restP)
            if r == 0:
                if D == 1:
                    cm, sg = canonicalize(tuple(restQ))
                    if cm is not None:
                        out[cm] = out.get(cm, 0) + F(num * sg)
                continue
            for ts in _compositions(D - 1, r):
                den = 1
                sh = []
                for (p, m, j), t in zip(restP, ts):
                    den *= factorial(t)
                    sh.append((p, m + t, j))
                cm, sg = canonicalize(tuple(sh) + tuple(restQ))
                if cm is None:
                    continue
                out[cm] = out.get(cm, 0) + F(num * sg, den)
    return {m: c for m, c in out.items() if c}


def residue1(P, Q, eta):
    return _apply_mono_op(P, Q, eta, residue1_mono)


def _vertex_mono_residue1(Pm, Wm, b, eta):
    """Single-contraction part of Res P(z) [W e^{b.X}](w): one boson letter of P hits the exponential, or one
    letter of P pairs with one letter of W; the rest Taylor-expanded (classical bracket with a vertex charge)."""
    out = {}
    n = len(Pm)
    # (a) one boson letter hits the exponential
    for i in range(n):
        p_, m, j = Pm[i]
        if p_ or not b[j]:
            continue
        num = eta[j] * b[j] * ((-1) ** (m - 1) * factorial(m - 1))
        D = m
        restP = [Pm[k] for k in range(n) if k != i]
        if not restP:
            if D == 1:
                cm, sg = canonicalize(tuple(Wm))
                if cm is not None:
                    out[cm] = out.get(cm, 0) + sg * num
            continue
        for ts in _compositions(D - 1, len(restP)):
            den = 1
            sh = []
            for (p2, m2, j2), t in zip(restP, ts):
                den *= factorial(t)
                sh.append((p2, m2 + t, j2))
            cm, sg = canonicalize(tuple(sh) + tuple(Wm))
            if cm is None:
                continue
            out[cm] = out.get(cm, 0) + F(sg * num, den) if isinstance(num, (int, F)) else out.get(cm, 0) + sg * num / den
    # (b) one letter of P pairs with one letter of W
    nw = len(Wm)
    for i in range(n):
        for c_ in range(nw):
            pr = _pairing(Pm[i], Wm[c_])
            if pr is None:
                continue
            num = pr[0] * eta[Pm[i][2]] * _koszul(Pm, Wm, (i,), (c_,))
            D = pr[1]
            restP = [Pm[k] for k in range(n) if k != i]
            restW = [Wm[k] for k in range(nw) if k != c_]
            if not restP:
                if D == 1:
                    cm, sg = canonicalize(tuple(restW))
                    if cm is not None:
                        out[cm] = out.get(cm, 0) + sg * num
                continue
            for ts in _compositions(D - 1, len(restP)):
                den = 1
                sh = []
                for (p2, m2, j2), t in zip(restP, ts):
                    den *= factorial(t)
                    sh.append((p2, m2 + t, j2))
                cm, sg = canonicalize(tuple(sh) + tuple(restW))
                if cm is None:
                    continue
                out[cm] = out.get(cm, 0) + F(sg * num, den) if isinstance(num, (int, F)) else out.get(cm, 0) + sg * num / den
    return {m: c for m, c in out.items() if c}


def vertex_residue1(P, V, eta):
    out = {}
    for pm, cp in P.items():
        for wm, cw in V.W.items():
            for mo, c in _vertex_mono_residue1(pm, wm, V.b, eta).items():
                out[mo] = out.get(mo, 0) + c * cp * cw
    return clean(out)


def joint_kernel(cands, conditions, eta, stats=None, classical=False):
    """Genuine charges in span(cands) annihilated by every condition (Vertex: twisted-exact residue; current Q:
    Res(Q, cand) in Im d).  Exactly super2_screening.joint_kernel with the metric-aware residue for currents."""
    nf = len(eta)
    if not cands:
        return []
    nc = len(cands)
    w2c = mono_weight2(next(iter(cands[0])))
    blocks, offset = [], nc
    apiv = SparsePivots()
    for cond in conditions:
        if isinstance(cond, Vertex):
            w2t = w2c + cond.w2W - 2
            tgt = basis(w2t, nf)
            idx = {m: i + offset for i, m in enumerate(tgt)}
            for u in basis(w2t - 2, nf):
                col = {}
                for m, c in twisted_d(u, cond.b).items():
                    col[idx[m]] = col.get(idx[m], 0) + c
                iv, _ = _to_int_vec({k: v for k, v in col.items() if v}, {k: k for k in col})
                if iv:
                    apiv.insert(iv)
        else:
            w2t = w2c + mono_weight2(next(iter(cond))) - 2
            tgt = basis(w2t, nf)
            idx = {m: i + offset for i, m in enumerate(tgt)}
            for u in basis(w2t - 2, nf):
                col = {}
                for m, c in d_mono(u).items():
                    col[idx[m]] = col.get(idx[m], 0) + c
                col = {k: v for k, v in col.items() if v}
                if col:
                    apiv.insert(col)
        blocks.append((cond, idx))
        offset += len(tgt)
    kernel_tags, scales = [], []
    for jc, cand in enumerate(cands):
        vec = {}
        for cond, idx in blocks:
            if isinstance(cond, Vertex):
                R = vertex_residue1(cand, cond, eta) if classical else vertex_residue(cand, cond, eta)
            else:
                R = residue1(cond, cand, eta) if classical else residue(cond, cand, eta)
            for m, c in R.items():
                vec[idx[m]] = vec.get(idx[m], 0) + c
        vec = {k: v for k, v in vec.items() if v}
        iv, L = _to_int_vec(vec, {k: k for k in vec})
        scales.append(L)
        iv[jc] = 1
        r = apiv.reduce(iv)
        if r and max(r) < nc:
            kernel_tags.append(r)
        elif r:
            apiv.rows[max(r)] = r
    cs_index = {m: i for i, m in enumerate(basis(w2c, nf))}
    span = SparsePivots()
    for u in basis(w2c - 2, nf):
        col = {}
        for m, c in d_mono(u).items():
            col[cs_index[m]] = col.get(cs_index[m], 0) + c
        col = {k: v for k, v in col.items() if v}
        if col:
            span.insert(col)
    genuine = []
    for t in kernel_tags:
        cur = {}
        for j, c in t.items():
            cur_add(cur, cands[j], F(c) * scales[j])
        cur = {m: c for m, c in cur.items() if c}
        if not cur:
            continue
        iv, _ = _to_int_vec(cur, cs_index)
        if span.insert(iv):
            genuine.append(cur)
    if stats is not None:
        stats.update(n_cand=nc, n_rows=offset - nc, n_pivots=len(apiv.rows))
    return genuine


def mod_d_coords(cur, nf=NF):
    """Reduce a current modulo total derivatives to a canonical representative (integer-scaled); returns
    (reduced sparse vector dict, basis index) -- for span tests of densities mod d."""
    w2 = mono_weight2(next(iter(cur)))
    index = {m: i for i, m in enumerate(basis(w2, nf))}
    piv = SparsePivots()
    for u in basis(w2 - 2, nf):
        iv, _ = _to_int_vec(d_mono(u), index)
        if iv:
            piv.insert(iv)
    iv, L = _to_int_vec(cur, index)
    return piv.reduce(iv), piv, index


def in_span_mod_d(cur, gens, nf=NF):
    """Is cur in span(gens) + Im d ?  (all of one doubled weight; exact)"""
    if not cur:
        return True
    w2 = mono_weight2(next(iter(cur)))
    index = {m: i for i, m in enumerate(basis(w2, nf))}
    piv = SparsePivots()
    for u in basis(w2 - 2, nf):
        iv, _ = _to_int_vec(d_mono(u), index)
        if iv:
            piv.insert(iv)
    for g in gens:
        iv, _ = _to_int_vec(g, index)
        if iv:
            piv.insert(iv)
    iv, _ = _to_int_vec(cur, index)
    return not piv.reduce(iv)


def rank_mod_d(gens, nf=NF):
    if not gens:
        return 0
    w2 = mono_weight2(next(iter(gens[0])))
    index = {m: i for i, m in enumerate(basis(w2, nf))}
    piv = SparsePivots()
    for u in basis(w2 - 2, nf):
        iv, _ = _to_int_vec(d_mono(u), index)
        if iv:
            piv.insert(iv)
    n0 = len(piv.rows)
    for g in gens:
        iv, _ = _to_int_vec(g, index)
        if iv:
            piv.insert(iv)
    return len(piv.rows) - n0


# ------------------------------------------------------------------ cylinder vev, NS sector, metric eta

def vev_ns(dens, eta, P):
    """NS cylinder eigenvalue of a density on the Fock highest-weight state with momenta P (list, per field,
    sympy or Fractions): super_cyl_vev's validated recipe per monomial times eta_j^(#contractions of field j).
    Bosons: contractions = (#boson letters of field j - exponent of P_j)/2; fermions (all paired in NS): #/2."""
    import sympy as sp
    nf = len(eta)
    total = sp.Integer(0)
    for mono, c in dens.items():
        c = sp.Rational(c.numerator, c.denominator) if isinstance(c, F) else sp.sympify(c)
        nb = [0] * nf
        nfm = [0] * nf
        for (p, m, j) in mono:
            if p == 0:
                nb[j] += 1
            else:
                nfm[j] += 1
        r = SCV.vev_monomial(mono, 'NS', nf)
        for zw, dct in r.items():
            assert zw == (), zw
            for ex, v in dct.items():
                term = sp.Rational(v.numerator, v.denominator)
                for j in range(nf):
                    nc = (nb[j] - ex[j])
                    assert nc % 2 == 0
                    term = term * sp.sympify(eta[j]) ** (nc // 2 + nfm[j] // 2) * sp.sympify(P[j]) ** ex[j]
                total += c * term
    return sp.expand(total)
