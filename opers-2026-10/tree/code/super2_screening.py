"""Two-superfield screening layer (SUSY Phase 4.1, Track G, step G1).

Two free N=1 superfields Phi_j = X_j + theta psi_j (j = 0, 1 in code;
"field 1, field 2" in prose), letters as in super_engine.py:
(parity, order, field).  NEW FILE: no existing engine is modified; all
OPE residues between local currents come from super_engine.residue_mono
/ normal_prod_mono, post-multiplied by the metric factor (below).

METRIC.  Field j carries a sign eta_j in its propagators,
    <d X_j d X_j> = eta_j (z-w)^-2,   <psi_j psi_j> = eta_j (z-w)^-1.
eta = (1, 1) is the engine's convention.  eta_0 = -1 is the exact
rational stand-in for an imaginary momentum component: substituting
X_0 = i Y_0, psi_0 = i chi_0 and b_0 = i p maps the Euclidean theory at
complex momentum onto the eta = (-1, 1) theory at real momentum, so the
paperclip-like regime n > 0 (bosonic: alpha_1 = sqrt(-n/2) imaginary)
is reachable in exact Fraction arithmetic.  Every Wick contraction
removes one field-j letter from each side, so the metric factor of an
emitted monomial is prod_j eta_j^{#field-j contractions}, read off from
letter counts (residue_cur_eta).  In the Z2 x Z2-even sector every
monomial has an even number of field-0 letters, so the substitution is
real there.

    T_eta = sum_j eta_j [1/2 (dX_j)^2 - 1/2 psi_j dpsi_j]
    G_eta = sum_j eta_j psi_j dX_j        (Res(G, G) = 2T, validated)

VERTICES.  V = W e^{b.X} with b a momentum vector (Fractions or sympy)
and W a prefactor current (polynomial in letters; W = {(): 1} for a
bosonic vertex, W = u.psi for the fermionic one).  Contractions of a
density P(z) against V(w):
  - any subset of P's BOSON letters hits the exponential,
    d^m X_j -> eta_j b_j (-1)^(m-1) (m-1)! (z-w)^-m;
  - any partial matching of P's remaining letters with W's letters,
    with the engine pairing times eta, and the engine's Koszul sign;
  - P's surviving letters are Taylor-expanded at w; W's never are.
oint V commutes with oint P iff Res_{z->w} P(z) V(w) is a TWISTED total
derivative: Res = D_b(J) e^{b.X},  D_b J = dJ + (b.dX) J  (chain rule,
no metric).  The exponential's conformal weight never enters this
algebra; weight conditions arise only when one demands that V screen
particular composites (derived, not assumed: super2_screening_validate
derives the quadric from the total super-Virasoro composites).

SECTORS.  Z2 x Z2 acts by (X_j, psi_j) -> -(X_j, psi_j) separately; the
even sector = monomials with an even number of field-j letters for each
j.  d preserves letter counts, so the sector is closed under d and
"densities mod d" restricts cleanly.  Odd doubled weight = odd number
of fermions = fermionic charges, scanned alongside the even ones.

G-SINGLET.  oint P is a G-singlet iff Res(G_eta, P) is in Im d.

JOINT KERNEL.  joint_kernel(cands, conditions): one augmented integer
echelon per run, exactly the tag scheme of
super_engine.genuine_kernel_super, with one target block per condition
(screening vertex: twisted-derivative witnesses; current: ordinary d
witnesses).  Genuineness = new pivot beyond the d-span of the
candidate sector.  Exact arithmetic throughout (Fractions -> integers).
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import factorial

from ope_engine import _compositions
from sparse_linalg import SparsePivots, _to_int_vec
from super_engine import (_koszul, _pairing, canonicalize, cur_add,
                          d_mono, gen_basis_super, is_odd, mono_weight2,
                          normal_prod_mono, residue_mono)

NF = 2


def X(m, j):
    return (0, m, j)


def PSI(a, j):
    return (1, a, j)


# ------------------------------------------------------------ basics

def basis(w2, nf=NF):
    """Monomial basis at doubled weight w2, including the empty
    monomial at w2 = 0 (needed for witnesses / bosonic-vertex targets)."""
    if w2 < 0:
        return []
    if w2 == 0:
        return [()]
    return gen_basis_super(w2, nf)


def field_counts(mono, nf=NF):
    c = [0] * nf
    for _, _, j in mono:
        c[j] += 1
    return tuple(c)


def is_even_sector(mono, nf=NF):
    return all(k % 2 == 0 for k in field_counts(mono, nf))


def even_basis(w2, nf=NF):
    return [m for m in basis(w2, nf) if is_even_sector(m, nf)]


def _metric_factor(cP, cQ, cO, eta):
    f = 1
    for j, e in enumerate(eta):
        n2 = cP[j] + cQ[j] - cO[j]
        assert n2 % 2 == 0 and n2 >= 0
        if e == -1 and (n2 // 2) % 2:
            f = -f
    return f


def residue_cur_eta(P, Q, eta=(1, 1)):
    """Res_{z->w} P(z) Q(w) with metric eta (engine residues times the
    metric factor counted from letter numbers)."""
    out = {}
    trivial = all(e == 1 for e in eta)
    for mp, cp in P.items():
        fp = field_counts(mp)
        for mq, cq in Q.items():
            r = residue_mono(mp, mq)
            if not r:
                continue
            fq = field_counts(mq)
            for mo, c in r.items():
                s = 1 if trivial else _metric_factor(fp, fq,
                                                     field_counts(mo), eta)
                v = out.get(mo, 0) + s * c * cp * cq
                if v:
                    out[mo] = v
                else:
                    out.pop(mo, None)
    return out


def normal_prod_eta(A, B, eta=(1, 1)):
    out = {}
    trivial = all(e == 1 for e in eta)
    for ma, ca in A.items():
        fa = field_counts(ma)
        for mb, cb in B.items():
            fb = field_counts(mb)
            for mo, c in normal_prod_mono(ma, mb).items():
                s = 1 if trivial else _metric_factor(fa, fb,
                                                     field_counts(mo), eta)
                v = out.get(mo, 0) + s * c * ca * cb
                if v:
                    out[mo] = v
                else:
                    out.pop(mo, None)
    return out


def d_cur(cur):
    out = {}
    for mono, c in cur.items():
        cur_add(out, d_mono(mono), c)
    return out


def T_eta(eta=(1, 1)):
    T = {}
    for j, e in enumerate(eta):
        T[(X(1, j), X(1, j))] = F(e, 2)
        T[(PSI(0, j), PSI(1, j))] = F(-e, 2)
    return T


def G_eta(eta=(1, 1)):
    return {(X(1, j), PSI(0, j)): F(e) for j, e in enumerate(eta)}


def I3_total(eta=(1, 1)):
    """:TT: - 1/4 :G dG: of the TOTAL (c = 3) super-Virasoro."""
    T, G = T_eta(eta), G_eta(eta)
    out = normal_prod_eta(T, T, eta)
    cur_add(out, normal_prod_eta(G, d_cur(G), eta), F(-1, 4))
    return {m: c for m, c in out.items() if c}


# ------------------------------------------------------------ vertices

class Vertex:
    """V = W e^{b.X}.  b: tuple of numbers; W: prefactor current."""

    def __init__(self, b, W=None, label=""):
        self.b = tuple(b)
        self.W = W if W is not None else {(): F(1)}
        self.label = label or f"V{self.b}"
        w2s = {mono_weight2(m) for m in self.W}
        assert len(w2s) == 1, "prefactor must be weight-homogeneous"
        self.w2W = w2s.pop()

    def __repr__(self):
        return self.label


def susy_vertex(b, eta=(1, 1), label=""):
    """(b.psi) e^{b.X} = sum_j b_j psi_j e^{b.X}: the top component of
    the superfield exponential (dtheta-integral of e^{b.Phi}).  The
    proportionality u ~ b (with NO metric factor, for every eta) is
    DERIVED in the validation from G-covariance of the vertex:
    Res(G_eta, (u.psi) e^{b.X}) = (u.dX) + (b0 u1 - b1 u0) psi0 psi1,
    which lies in Im D_b iff u is parallel to b.  `eta` is accepted for
    call-site symmetry only."""
    W = {(PSI(0, j),): b[j] for j in range(len(b)) if b[j]}
    return Vertex(b, W, label or f"(b.psi)e^(b.X), b={tuple(map(str, b))}")


def fermion_vertex(b, u, label=""):
    W = {(PSI(0, j),): u[j] for j in range(len(u)) if u[j]}
    return Vertex(b, W, label or f"(u.psi)e^(b.X), b={b}, u={u}")


def boson_vertex(b, label=""):
    return Vertex(b, {(): F(1)}, label or f"e^(b.X), b={b}")


def _div(x, den):
    return F(x, den) if isinstance(x, int) else x / den


def _vertex_mono_residue(Pm, Wm, b, eta, shift=0):
    """Coefficient of (z-w)^(-1-shift) in Pm(z) [Wm e^{b.X}](w), as
    {monomial: coeff} (exponential implicit).  Generic in the number
    type of b."""
    out = {}
    n = len(Pm)
    bos = [i for i in range(n) if not is_odd(Pm[i])]
    for kb in range(len(bos) + 1):
        for bsel in combinations(bos, kb):
            num_e, D_e = 1, 0
            zero = False
            for i in bsel:
                _, m, j = Pm[i]
                if not b[j]:
                    zero = True
                    break
                num_e = num_e * eta[j] * b[j] * ((-1) ** (m - 1)
                                                 * factorial(m - 1))
                D_e += m
            if zero:
                continue
            restP0 = tuple(Pm[i] for i in range(n) if i not in bsel)
            n0, nw = len(restP0), len(Wm)
            for k in range(0, min(n0, nw) + 1):
                for psel in combinations(range(n0), k):
                    for qsel in permutations(range(nw), k):
                        num, D = num_e, D_e
                        ok = True
                        for a_, c_ in zip(psel, qsel):
                            pr = _pairing(restP0[a_], Wm[c_])
                            if pr is None:
                                ok = False
                                break
                            num = num * pr[0] * eta[restP0[a_][2]]
                            D += pr[1]
                        if not ok:
                            continue
                        if k:
                            num = num * _koszul(restP0, Wm, psel, qsel)
                        if D - shift <= 0:
                            continue
                        restP = [restP0[i] for i in range(n0)
                                 if i not in psel]
                        qs = set(qsel)
                        restW = [Wm[i] for i in range(nw) if i not in qs]
                        r = len(restP)
                        if r == 0:
                            if D - shift == 1:
                                cm, sg = canonicalize(tuple(restW))
                                if cm is not None:
                                    out[cm] = out.get(cm, 0) + sg * num
                            continue
                        for ts in _compositions(D - 1 - shift, r):
                            den = 1
                            shifted = []
                            for (p, m, j), t in zip(restP, ts):
                                den *= factorial(t)
                                shifted.append((p, m + t, j))
                            cm, sg = canonicalize(tuple(shifted)
                                                  + tuple(restW))
                            if cm is None:
                                continue
                            out[cm] = out.get(cm, 0) + _div(sg * num, den)
    return out


def vertex_residue(P, V, eta=(1, 1), shift=0):
    """Res_{z->w} P(z) V(w) (shift = 0) or the coefficient of
    (z-w)^(-1-shift), polynomial part only.  Coefficients exact
    (Fractions) or sympy, following the input types."""
    out = {}
    for pm, cp in P.items():
        for wm, cw in V.W.items():
            for mo, c in _vertex_mono_residue(pm, wm, V.b, eta,
                                              shift).items():
                out[mo] = out.get(mo, 0) + c * cp * cw
    try:
        import sympy as sp
        return {m: (sp.expand(c) if isinstance(c, sp.Basic) else c)
                for m, c in out.items()
                if (sp.expand(c) if isinstance(c, sp.Basic) else c) != 0}
    except ImportError:
        return {m: c for m, c in out.items() if c}


def twisted_d(mono, b):
    """D_b J = dJ + (b.dX) J on a monomial (J times e^{b.X})."""
    out = {m: c for m, c in d_mono(mono).items()}
    for j, bj in enumerate(b):
        if not bj:
            continue
        cm, sg = canonicalize(mono + (X(1, j),))
        if cm is not None:
            out[cm] = out.get(cm, 0) + sg * bj
    return {m: c for m, c in out.items() if c}


def is_twisted_exact(R, b, w2, nf=NF):
    """Exact membership of R (polynomial part, doubled weight w2) in
    Im D_b at a rational b."""
    if not R:
        return True
    tgt = basis(w2, nf)
    index = {m: i for i, m in enumerate(tgt)}
    piv = SparsePivots()
    for u in basis(w2 - 2, nf):
        iv, _ = _to_int_vec(twisted_d(u, b), index)
        if iv:
            piv.insert(iv)
    iv, _ = _to_int_vec(R, index)
    return not piv.reduce(iv)


def screens(P, V, eta=(1, 1)):
    """True iff [oint V, oint P] = 0 (exact, rational momenta)."""
    w2P = mono_weight2(next(iter(P)))
    return is_twisted_exact(vertex_residue(P, V, eta), V.b,
                            w2P + V.w2W - 2)


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


def g_singlet(P, eta=(1, 1)):
    return is_total_derivative(residue_cur_eta(G_eta(eta), P, eta))


# ------------------------------------------------------- joint kernel

def joint_kernel(cands, conditions, eta=(1, 1), nf=NF, stats=None):
    """Genuine charges in span(cands) (weight-homogeneous currents of
    one doubled weight) annihilated by every condition.  A condition is
    a Vertex (screening: Res(P, V) in Im D_b) or a current Q
    (Res(Q, P) in Im d, e.g. Q = G_eta for the G-singlet condition).
    Returns a list of currents (Fractions), a basis of the genuine
    joint kernel modulo total derivatives."""
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
                iv, _ = _to_int_vec({k: v for k, v in col.items() if v},
                                    {k: k for k in col})
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
                R = vertex_residue(cand, cond, eta)
            else:
                R = residue_cur_eta(cond, cand, eta)
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
        stats.update(n_cand=nc, n_rows=offset - nc,
                     n_pivots=len(apiv.rows))
    return genuine


def sector_candidates(w2, sector="even", nf=NF):
    """Monomial candidates at doubled weight w2: 'even' = Z2 x Z2-even,
    'all' = every monomial."""
    mons = even_basis(w2, nf) if sector == "even" else basis(w2, nf)
    return [{m: F(1)} for m in mons]


def charge_space_dim(w2, sector="even", nf=NF):
    """dim of densities mod d in the sector."""
    mons = even_basis(w2, nf) if sector == "even" else basis(w2, nf)
    index = {m: i for i, m in enumerate(basis(w2, nf))}
    piv = SparsePivots()
    for u in basis(w2 - 2, nf):
        if sector == "even" and not is_even_sector(u, nf):
            continue
        iv, _ = _to_int_vec(d_mono(u), index)
        if iv:
            piv.insert(iv)
    return len(mons) - len(piv.rows)


def four_screenings(b, eta=(1, 1), kind="susy"):
    """The Z2 x Z2 orbit of one vertex: momenta (+-b0, +-b1); for the
    SUSY vertex the prefactor u = b flips with b (the orbit of
    (b.psi) e^{b.X} under (X_j, psi_j) -> -(X_j, psi_j) up to an
    overall sign, irrelevant for a screening condition)."""
    out = []
    for s0 in (1, -1):
        for s1 in (1, -1):
            bb = (s0 * b[0], s1 * b[1])
            out.append(susy_vertex(bb, eta) if kind == "susy"
                       else boson_vertex(bb))
    return out
