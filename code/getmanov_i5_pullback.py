"""Getmanov I5 pullback test -- the decisive identification
computation for the enhanced c=1 tower (successor to
getmanov_stage2.py).

Question: the phi-axis pullback of the Getmanov I3 matches the
enhanced I3' ray at k = -1 and k = -16 (certified).  Does the
Getmanov I5 -- the weight-6 charge of the structure, derived as the
joint kernel of the three Getmanov screenings on the full two-boson
weight-6 space mod total derivatives -- pull back to the enhanced
tower's I5 class (the A2(2) I5 at N=1) at those points?

Arithmetic stack (v3): a first implementation ran the kernel
intersections in sympy exact algebraic arithmetic and spent 9.7 h
inside one nullspace before being killed.  This version does
RESTRICTION OF SCALARS: every screening-system entry lives in an
explicit number field K (Q(sqrt2) at k=-1; Q(sqrt2, i sqrt7) at
k=-16; Q(sqrt2, i) at the generic validation point), split into
exact rational coordinates over the field basis (per-entry
reconstruction asserts), with all linear algebra as sparse
fraction-free INTEGER elimination on the split system.

Correlation discipline (audit round: a v2 separated per-basis
component currents, which ENLARGES the span -- the negative verdict
survives a fortiori, but dims and positive validations were
computed in the wrong space): here solution vectors are kept as
FULL split vectors (all field coordinates jointly), whose Q-span is
exactly the restriction of scalars of the K-solution space
(B-stable by K-linearity); a rational target t lies in the K-span
iff its block-0 embedding lies in the Q-span of the split vectors.
Pullback and mod-d act blockwise and commute with the
multiplication-by-basis action, so every space tested is B-stable
and its K-dimension is Q-dim/|B| (divisibility asserted).

Scope of the verdict: the two RATIONAL phi-axis matching points
k=-1, -16 only.  The theta-axis matching locus (irrational cubic
roots) and tilted embeddings are NOT tested here; special-fiber
kernel dims are reported as such (they do not by themselves prove
regular specialization of the generic I5^G).  The committed verdict
is baked as a fail-closed regression check."""
import sys
from fractions import Fraction as F

import sympy as sp

from families import KDV_P6
from ope_engine import d_mono, gen_basis
from ope_poles import ope_pole
from screening_fit import residue_screen, twisted_d
from sparse_linalg import SparsePivots, lcm
from super_engine import d_cur as _sd_cur
from super_engine import normal_prod

failures = []


def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok:
        failures.append(name)


# ------- W-current assembly (construction certified in stage 2) -----
rk = sp.Symbol('rk')                     # unrestricted branch symbol
k = rk**2
p = sp.I / sp.sqrt(2)
t = sp.I / sp.sqrt(2)


def X(m, j):
    return (0, m, j)


def scur(d):
    out = {}
    for m, c in d.items():
        key = tuple(sorted(m))
        out[key] = sp.expand(out.get(key, 0) + c)
    return {m: c for m, c in out.items() if c != 0}


def add(A, B, co=1):
    out = dict(A)
    for m, c in B.items():
        out[m] = sp.expand(out.get(m, 0) + co * c)
    return {m: c for m, c in out.items() if c != 0}


def dsc(A):
    return {m: sp.expand(c) for m, c in _sd_cur(A).items()
            if sp.expand(c) != 0}


W2 = scur({(X(1, 0), X(1, 0)): -p**2,
           (X(1, 1), X(1, 1)): -t**2,
           (X(2, 0),): -p / rk})
W3 = scur({(X(1, 1), X(1, 1), X(1, 1)): -sp.I * (6*k + 4) / (3*k) * t**3,
           (X(1, 0), X(1, 0), X(1, 1)): -sp.I * 2 * p**2 * t,
           (X(2, 0), X(1, 1)): sp.I * rk * p * t,
           (X(1, 0), X(2, 1)): -sp.I * (k + 2) / rk * p * t,
           (X(3, 1),): -sp.I * (k + 2) / (6 * k) * t})
cval = 2 + 6 / k
C4 = (cval + 7) * (2 * cval - 1) / (3 * (cval - 2))
mu = (2 * cval**2 + 22 * cval - 25) / (30 * (cval - 2))
W2sq = {m: sp.expand(c) for m, c in normal_prod(W2, W2).items()
        if sp.expand(c) != 0}
W4 = add(add(ope_pole(W3, W3, 2), W2sq, -2), dsc(dsc(W2)),
         2 * mu - C4 / 2)
W4 = {m: sp.expand(c / 2) for m, c in W4.items() if sp.expand(c) != 0}
I3G = add({m: (k + 6) * c for m, c in W4.items()}, W2sq, 4 * (2*k + 1))


def bos(cur):
    return {tuple((m, j) for _, m, j in mono): c
            for mono, c in cur.items()}


def at_rk(cur, rkval):
    out = {}
    for m, c in cur.items():
        v = sp.nsimplify(sp.simplify(c.subs(rk, rkval)))
        if v != 0:
            out[m] = v
    return out


def momenta(rkval):
    aphi = sp.nsimplify(sp.simplify(-sp.I * rkval / sp.sqrt(2)))
    ath = sp.nsimplify(sp.simplify(sp.sqrt((rkval**2 + 2) / 2)))
    amup = sp.nsimplify(sp.simplify(sp.I * 4 / (sp.sqrt(2) * rkval)))
    return [(aphi, ath), (aphi, -ath), (amup, sp.Integer(0))]


# ------- number-field splitting -------------------------------------
SQ2, SQ7, SQ14 = sp.sqrt(2), sp.sqrt(7), sp.sqrt(14)


def split(expr, basis):
    """Exact rational coordinates of expr over the field basis, with
    reconstruction assert (fail-closed)."""
    e = sp.expand(sp.radsimp(sp.nsimplify(expr)))
    re_, im_ = [sp.expand(z) for z in e.as_real_imag()]
    out, recon = [], sp.Integer(0)
    for b in basis:
        if b == 1:
            src, r = re_, sp.Integer(1)
        elif b.has(sp.I):
            src, r = im_, sp.expand(b / sp.I)
        else:
            src, r = re_, b
        if r == 1:
            c = src.coeff(SQ2, 0).coeff(SQ7, 0).coeff(SQ14, 0)
        else:
            c = src.coeff(r, 1)
        c = sp.nsimplify(c)
        if not c.is_Rational:
            raise ValueError(f"non-rational coord {c} of {e} on {b}")
        out.append(F(int(sp.numer(c)), int(sp.denom(c))))
        recon += c * b
    if sp.simplify(e - recon) != 0:
        raise ValueError(f"split reconstruction failed: {e} vs {recon}")
    return out


def int_vec(fr):
    """Fraction dict -> (integer dict, scale L) with int = L * frac."""
    L = 1
    for v in fr.values():
        L = lcm(L, v.denominator)
    return {i: int(v * L) for i, v in fr.items() if v}, L


# ------- charge-space linear algebra --------------------------------
def dspan(w, nf, nb=1):
    """d-image pivots in split coordinates (nb blocks; nb=1 plain)."""
    basis = gen_basis(w, nf, z2=False)
    idx = {m: i for i, m in enumerate(basis)}
    piv = SparsePivots()
    for bmono in gen_basis(w - 1, nf, z2=False):
        col = {}
        for m, c in d_mono(bmono).items():
            col[idx[m]] = col.get(idx[m], 0) + c
        col = {i: c for i, c in col.items() if c}
        if col:
            for beta in range(nb):
                piv.insert({i * nb + beta: c for i, c in col.items()})
    return basis, idx, piv


def clone(piv):
    q = SparsePivots()
    q.rows = dict(piv.rows)
    return q


def joint_split(w, moms, basis_field):
    """Kernel of the joint screening system as FULL split vectors
    {j*nb+beta: int} -- the exact restriction of scalars of the
    K-solution space (correlation-preserving)."""
    nb = len(basis_field)
    basisw = gen_basis(w, 2, z2=False)
    tgt = gen_basis(w - 1, 2, z2=False)
    tidx = {m: i for i, m in enumerate(tgt)}
    lower = gen_basis(w - 2, 2, z2=False)
    nt = len(tgt)

    Rcols_s, Tcols_s = [], []
    for al in moms:
        Rcols_s.append([{tidx[m]: c for m, c in
                         residue_screen({b: F(1)}, al).items()}
                        for b in basisw])
        Tcols_s.append([{tidx[m]: c for m, c in
                         twisted_d(u, al, 2).items()} for u in lower])

    bigcols, scales, xmeta = [], [], []
    for j in range(len(basisw)):
        for bi, b in enumerate(basis_field):
            col = {}
            for s in range(len(moms)):
                base = s * nt * nb
                for r, e in Rcols_s[s][j].items():
                    for beta, cv in enumerate(split(sp.expand(e * b),
                                                    basis_field)):
                        if cv:
                            kkey = base + r * nb + beta
                            col[kkey] = col.get(kkey, F(0)) + cv
            iv, L = int_vec(col)
            bigcols.append(iv)
            scales.append(L)
            xmeta.append((j, bi))
    nx = len(bigcols)
    for s in range(len(moms)):
        base = s * nt * nb
        for u in range(len(lower)):
            for b in basis_field:
                col = {}
                for r, e in Tcols_s[s][u].items():
                    for beta, cv in enumerate(split(sp.expand(e * b),
                                                    basis_field)):
                        if cv:
                            kkey = base + r * nb + beta
                            col[kkey] = col.get(kkey, F(0)) + cv
                iv, L = int_vec(col)
                bigcols.append(iv)
                scales.append(L)

    nc = len(bigcols)
    piv = SparsePivots()
    kers = []
    for jcol, col in enumerate(bigcols):
        aug = {i + nc: v for i, v in col.items()}
        aug[jcol] = 1
        r = piv.reduce(aug)
        if r and max(r) < nc:
            kers.append(r)
        elif r:
            piv.rows[max(r)] = r

    vecs = []
    for kv in kers:
        xv = {}
        for cidx, coef in kv.items():
            if cidx < nx:
                v = coef * scales[cidx]
                if v:
                    xv[cidx] = v          # cidx = j*nb + beta
        if xv:
            vecs.append(xv)
    return basisw, vecs


def qdim_mod_d(vecs, dpiv):
    q = clone(dpiv)
    return sum(1 for v in vecs if v and q.insert(dict(v)))


def member_mod_d(target, vecs, dpiv):
    q = clone(dpiv)
    for v in vecs:
        if v:
            q.insert(dict(v))
    return not q.reduce(dict(target))


def pull_phi(splitvec, basisw, idx1, nb):
    out = {}
    for key, c in splitvec.items():
        j, beta = divmod(key, nb)
        m = basisw[j]
        if all(fj == 0 for _, fj in m):
            kk = idx1[m] * nb + beta
            out[kk] = out.get(kk, 0) + c
    return {a: v for a, v in out.items() if v}


def to_int_target(cur, idx, nb=1):
    """Rational target embedded in block 0 of the split space."""
    fr = {}
    for m, c in cur.items():
        if not isinstance(c, F):
            c = F(int(sp.numer(c)), int(sp.denom(c)))
        fr[idx[m] * nb] = fr.get(idx[m] * nb, F(0)) + c
    iv, _ = int_vec(fr)
    return iv


def kdim(vecs, dpiv, nb, label):
    qd = qdim_mod_d(vecs, dpiv)
    if qd % nb:
        check(f"{label}: Q-dim {qd} divisible by |B|={nb} "
              f"(B-stability)", False)
    return qd // nb


def field_split_vec(cur_bos, basis_field, basisw):
    """A field-valued two-boson current as ONE full split vector."""
    nb = len(basis_field)
    idxw = {m: i for i, m in enumerate(basisw)}
    fr = {}
    for m, c in cur_bos.items():
        for beta, cv in enumerate(split(c, basis_field)):
            if cv:
                fr[idxw[m] * nb + beta] = cv
    iv, _ = int_vec(fr)
    return iv


# ------- targets ----------------------------------------------------
A2_N1 = {((3, 0), (3, 0)): sp.Rational(31, 24),
         ((1, 0), (1, 0), (2, 0), (2, 0)): sp.Rational(5, 2),
         ((1, 0),) * 6: sp.Integer(1)}
KDV6 = {m: sp.Rational(c.numerator, c.denominator)
        for m, c in KDV_P6.items()}
I3P_RAY = {((1, 0),) * 4: sp.Integer(2),
           ((1, 0), (3, 0)): sp.Integer(-5)}

B_GEN = [sp.Integer(1), SQ2, sp.I, sp.I * SQ2]
B_M1 = [sp.Integer(1), SQ2]
B_M16 = [sp.Integer(1), SQ2, sp.I * SQ7, sp.I * SQ14]

# plain rational comparator check (nb=1 space)
_, idxc, dpc = dspan(6, 1, 1)
check("A2(N=1) and KdV weight-6 classes independent (comparators ok)",
      qdim_mod_d([to_int_target(A2_N1, idxc, 1),
                  to_int_target(KDV6, idxc, 1)], dpc) == 2)

# ------- generic-point validation (k = 1/4, rk = 1/2) ---------------
nbg = len(B_GEN)
mg = momenta(sp.Rational(1, 2))
bw4, vecs4 = joint_split(4, mg, B_GEN)
_, _, dp24g = dspan(4, 2, nbg)
kd4 = kdim(vecs4, dp24g, nbg, "generic w4")
check("generic k=1/4: weight-4 screened charge K-dim = 1", kd4 == 1)
i3g_vec = field_split_vec(bos(at_rk(I3G, sp.Rational(1, 2))), B_GEN,
                          bw4)
check("generic: assembled I3^G lies in the kernel span (ray equality "
      "given K-dim 1)", member_mod_d(i3g_vec, vecs4, dp24g))

bw6g, vecs6g = joint_split(6, mg, B_GEN)
_, _, dp26g = dspan(6, 2, nbg)
kd6 = kdim(vecs6g, dp26g, nbg, "generic w6")
print(f"generic k=1/4: weight-6 charge K-dim = {kd6}", flush=True)
check("generic k=1/4: weight-6 screened charge K-dim = 1 (I5^G "
      "exists, multiplicity one)", kd6 == 1)

# ------- the special points -----------------------------------------
REVIEW_DIMS = {-1: (2, 1), -16: (1, 1)}   # (joint w6, pullback) K-dims
verdict = {}
for kval, rkval, BF in ((-1, sp.I, B_M1), (-16, 4 * sp.I, B_M16)):
    nb = len(BF)
    mo = momenta(rkval)
    b14, idx14, dp14 = dspan(4, 1, nb)
    bw4p, vecs4p = joint_split(4, mo, BF)
    pb4 = [pull_phi(v, bw4p, idx14, nb) for v in vecs4p]
    check(f"k={kval}: weight-4 kernel pullback contains the I3' ray "
          f"(coherence)",
          member_mod_d(to_int_target(I3P_RAY, idx14, nb), pb4, dp14))

    b16s, idx16s, dp16s = dspan(6, 1, nb)
    bw6p, vecs6p = joint_split(6, mo, BF)
    _, _, dp26p = dspan(6, 2, nb)
    kd6p = kdim(vecs6p, dp26p, nb, f"k={kval} w6")
    pb6 = [pull_phi(v, bw6p, idx16s, nb) for v in vecs6p]
    pbk = kdim(pb6, dp16s, nb, f"k={kval} w6 pullback")
    inA2 = member_mod_d(to_int_target(A2_N1, idx16s, nb), pb6, dp16s)
    inKdV = member_mod_d(to_int_target(KDV6, idx16s, nb), pb6, dp16s)
    verdict[kval] = (kd6p, pbk, inA2, inKdV)
    print(f"k={kval}: weight-6 charge K-dim {kd6p}, pullback K-dim "
          f"{pbk}; enhanced-A2 class in pullback span: {inA2}; KdV "
          f"class in span: {inKdV}", flush=True)
    check(f"k={kval}: special-fiber weight-6 kernel nonzero", kd6p >= 1)
    check(f"k={kval}: K-dims match the review-verified grouped values "
          f"{REVIEW_DIMS[kval]}", (kd6p, pbk) == REVIEW_DIMS[kval])
    check(f"k={kval}: KdV class absent from the pullback span (the "
          f"printed 'nor is KdV' is certified)", not inKdV)

print()
hitA2 = all(verdict[kv][2] for kv in (-1, -16))
missA2 = all(not verdict[kv][2] for kv in (-1, -16))
# committed verdict baked as a regression check (audit finding 4)
check("committed verdict holds: I5 pullback misses at BOTH points",
      missA2)

if failures:
    print(f"FAILURES: {failures}")
    print("GETMANOV I5 PULLBACK TEST FAILED")
    sys.exit(1)

print("VERDICT: I5 PULLBACK FAILS at both points -- the enhanced I5 "
      "class is NOT in the pulled-back Getmanov weight-6 charge "
      "space at k=-1 or k=-16 (nor is KdV).  The I3-ray coincidence "
      "does not extend: the two RATIONAL phi-axis Getmanov "
      "candidates are excluded at I5.  Scope: theta-axis (irrational "
      "cubic roots) and tilted embeddings remain untested; AKNS "
      "remains disfavored, not excluded.")
print("GETMANOV I5 PULLBACK TEST COMPLETE")
