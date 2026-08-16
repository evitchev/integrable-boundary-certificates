"""Getmanov theta-axis test -- closing the last named-candidate
locus for the enhanced c=1 tower.

The phi-axis story is closed (getmanov_i5_pullback.py): the two
rational matching points k = -1, -16 are excluded at I5.  The
remaining axis-aligned locus is the THETA-axis I3 matching cubic
k^3 + 23 k^2 + 55 k + 42 = 0 (one real root ~ -20.41 and a complex
pair).  This certificate tests it exactly, at ALL THREE roots and
both sqrt branches simultaneously, by working over the etale
algebra

    A = Q[k, u, v] / (C(k), u^2 + k/2, v^2 - (k+2)/2),

dim_Q A = 12, where u = -i sqrt(k/2), v = sqrt((k+2)/2) carry the
Getmanov screening momenta: cos-pair (u, +-v), mu' vertex (2/u, 0)
with 2/u polynomial in A (u is invertible since C(0) != 0).  The
relations' leading monomials u^2, v^2, k^3 are pairwise coprime, so
they form a Groebner basis and sympy's reduced() gives a normal
form on the 12-monomial basis k^a u^b v^c.

Method: the weight-w Getmanov charge space is the joint kernel of
the three screening conditions on the full two-boson space mod
total derivatives -- the span-based construction validated at the
rational points -- computed by restriction of scalars over A with
sparse fraction-free integer elimination (correlation-preserving
full split vectors, as audit-corrected in the phi certificate).
Because A is a PRODUCT of number fields (one factor per Galois
orbit of roots/branches), a single membership answer could differ
between factors; the certificate therefore computes the MATCH
IDEAL

    E(t) = { a in A : a * t  in  A-span(pullback classes) mod d },

whose Q-dimension counts matched components with multiplicity:
dim E = 12 means the target class matches at EVERY root and branch,
dim E = 0 means it matches at NONE (a clean global verdict), and
anything between flags root-dependent structure explicitly instead
of hiding it.  Computationally E is the kernel of the map
a_j -> residual of (t embedded in block j) against the span pivots.
For the production cubic, A is in fact a degree-12 FIELD -- z = u+v
satisfies an irreducible degree-12 polynomial (asserted;
review-supplied), and a field of dimension dim_Q A inside A forces
A itself to be that field -- so E is necessarily 0 or all of A,
Q-dimension-12 spaces are exactly A-rank one, and the weight-4
joint and pullback dimensions are ASSERTED at 12 so that E4 = 12
means genuine RAY EQUALITY (membership in a one-dimensional
pullback line), not mere membership.  The E-mechanism and the rank
caveats remain relevant for general relations.

Validation ladder (all against previously certified facts, run
through THIS generalized code with the degree-1 relation k + 1):
phi-axis at k = -1 must give dim E = 4 (full) for the I3' ray at
weight 4 and dim E = 0 for the enhanced I5 at weight 6; theta-axis
at k = -1 must give dim E = 0 at weight 4.  Ring sanity (normal
forms, u-inverse, mu'-vertex identity) is asserted.

Opposite-branch coverage is intrinsic: the algebra contains both
signs of u and v as components, and the cos-pair is v-symmetric.
Scope: axis-aligned theta embeddings at the I3-matching cubic;
tilted embeddings remain outside this certificate."""
import sys
from fractions import Fraction as F

import sympy as sp

from families import KDV_P6
from ope_engine import d_mono, gen_basis
from screening_fit import residue_screen, twisted_d
from sparse_linalg import SparsePivots, lcm

failures = []


def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok:
        failures.append(name)


ks, us, vs = sp.symbols('ks us vs')


class Algebra:
    """Q[k,u,v]/(C(k), u^2 + k/2, v^2 - (k+2)/2) on the monomial
    basis k^a u^b v^c (a < deg C; basis[0] = 1)."""

    def __init__(self, Cpoly):
        self.rels = [sp.expand(Cpoly), us**2 + ks / 2,
                     vs**2 - (ks + 2) / 2]
        dk = sp.degree(Cpoly, ks)
        self.abc = [(a, b, c) for a in range(dk)
                    for b in (0, 1) for c in (0, 1)]
        self.idx = {t: i for i, t in enumerate(self.abc)}
        self.nb = len(self.abc)
        self.monos = [ks**a * us**b * vs**c for a, b, c in self.abc]
        self.dk = dk
        self.mt = [[self.coords(self.monos[i] * self.monos[j])
                    for j in range(self.nb)] for i in range(self.nb)]

    def nf(self, expr):
        return sp.reduced(sp.expand(expr), self.rels, us, vs, ks,
                          order='lex')[1]

    def coords(self, expr):
        out = [F(0)] * self.nb
        e = self.nf(expr)
        if e == 0:
            return out
        for (b, c, a), coef in sp.Poly(e, us, vs, ks).as_dict().items():
            r = sp.Rational(coef)
            out[self.idx[(a, b, c)]] = F(int(r.p), int(r.q))
        return out

    def mul_basis(self, x, j):
        """coords(x * basis_j) from the multiplication table."""
        out = [F(0)] * self.nb
        for i, xi in enumerate(x):
            if xi:
                row = self.mt[i][j]
                for t in range(self.nb):
                    if row[t]:
                        out[t] += xi * row[t]
        return out

    def inv(self, expr):
        """Inverse of expr in A as a sympy polynomial (raises if the
        element is a zero divisor)."""
        x = self.coords(expr)
        M = sp.Matrix([[sp.Rational(self.mul_basis(x, j)[t].numerator,
                                    self.mul_basis(x, j)[t].denominator)
                        for j in range(self.nb)]
                       for t in range(self.nb)])
        rhs = sp.Matrix([1] + [0] * (self.nb - 1))
        y = M.solve(rhs)
        return sp.expand(sum(y[j] * self.monos[j]
                             for j in range(self.nb)))


def int_vec(fr):
    L = 1
    for v in fr.values():
        L = lcm(L, v.denominator)
    return {i: int(v * L) for i, v in fr.items() if v}, L


def clone(piv):
    q = SparsePivots()
    q.rows = dict(piv.rows)
    return q


def dspan(w, nf_, nb):
    basis = gen_basis(w, nf_, z2=False)
    idx = {m: i for i, m in enumerate(basis)}
    piv = SparsePivots()
    for bmono in gen_basis(w - 1, nf_, z2=False):
        col = {}
        for m, c in d_mono(bmono).items():
            col[idx[m]] = col.get(idx[m], 0) + c
        col = {i: c for i, c in col.items() if c}
        if col:
            for beta in range(nb):
                piv.insert({i * nb + beta: c for i, c in col.items()})
    return basis, idx, piv


def joint_split(w, moms, alg):
    """Joint screening kernel as full split vectors {j*nb+beta: int}
    over the algebra (restriction of scalars)."""
    nb = alg.nb
    basisw = gen_basis(w, 2, z2=False)
    tgt = gen_basis(w - 1, 2, z2=False)
    tidx = {m: i for i, m in enumerate(tgt)}
    lower = gen_basis(w - 2, 2, z2=False)
    nt = len(tgt)

    Rcols_s, Tcols_s = [], []
    for al in moms:
        Rcols_s.append([{tidx[m]: alg.coords(c) for m, c in
                         residue_screen({b: F(1)}, al).items()}
                        for b in basisw])
        Tcols_s.append([{tidx[m]: alg.coords(c) for m, c in
                         twisted_d(u_, al, 2).items()} for u_ in lower])

    bigcols, scales = [], []
    nx = len(basisw) * nb
    for j in range(len(basisw)):
        for bi in range(nb):
            col = {}
            for s in range(len(moms)):
                base = s * nt * nb
                for r, e in Rcols_s[s][j].items():
                    for beta, cv in enumerate(alg.mul_basis(e, bi)):
                        if cv:
                            kkey = base + r * nb + beta
                            col[kkey] = col.get(kkey, F(0)) + cv
            iv, L = int_vec(col)
            bigcols.append(iv)
            scales.append(L)
    for s in range(len(moms)):
        base = s * nt * nb
        for u_ in range(len(lower)):
            for bi in range(nb):
                col = {}
                for r, e in Tcols_s[s][u_].items():
                    for beta, cv in enumerate(alg.mul_basis(e, bi)):
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
                val = coef * scales[cidx]
                if val:
                    xv[cidx] = val
        if xv:
            vecs.append(xv)
    return basisw, vecs


def pull_axis(splitvec, basisw, idx1, nb, axis):
    out = {}
    for key, c in splitvec.items():
        j, beta = divmod(key, nb)
        mono = basisw[j]
        if all(jj == axis for _, jj in mono):
            m1 = tuple((m, 0) for m, _ in mono)
            kk = idx1[m1] * nb + beta
            out[kk] = out.get(kk, 0) + c
    return {a: v for a, v in out.items() if v}


def qdim(vecs, dpiv):
    q = clone(dpiv)
    return sum(1 for v in vecs if v and q.insert(dict(v)))


def match_ideal_dim(t_base, span_vecs, dpiv, nb):
    """dim_Q of E(t) = {a in A: a*t in span mod d}: nb minus the rank
    of the per-block residuals of t."""
    q = clone(dpiv)
    for v in span_vecs:
        if v:
            q.insert(dict(v))
    rp = SparsePivots()
    rank = 0
    for beta in range(nb):
        r = q.reduce({i * nb + beta: c for i, c in t_base.items()})
        if r and rp.insert(dict(r)):
            rank += 1
    return nb - rank


def target_base(cur, idx):
    fr = {}
    for m, c in cur.items():
        if not isinstance(c, F):
            c = F(int(sp.numer(c)), int(sp.denom(c)))
        fr[idx[m]] = fr.get(idx[m], F(0)) + c
    iv, _ = int_vec(fr)
    return iv


# ------- targets (rational one-boson classes) -----------------------
I3P_RAY = {((1, 0),) * 4: sp.Integer(2), ((1, 0), (3, 0)): sp.Integer(-5)}
A2_N1 = {((3, 0), (3, 0)): sp.Rational(31, 24),
         ((1, 0), (1, 0), (2, 0), (2, 0)): sp.Rational(5, 2),
         ((1, 0),) * 6: sp.Integer(1)}
KDV6 = {m: sp.Rational(c.numerator, c.denominator)
        for m, c in KDV_P6.items()}


def run_point(alg, axis, label):
    """Weight-4 and weight-6 kernels, axis pullback, match-ideal dims
    for I3' (w4) and A2/KdV (w6).  Returns (E4, E6_A2, E6_KdV,
    w4_joint_qdim, w4_pullback_qdim, w6_joint_qdim)."""
    nb = alg.nb
    mu_al = 2 * alg.inv(us)
    check(f"{label}: mu' vertex identity (2/u) * u = 2 in A",
          alg.coords(sp.expand(mu_al * us)) == alg.coords(2))
    moms = [(us, vs), (us, -vs), (sp.expand(mu_al), sp.Integer(0))]

    b14, idx14, dp14 = dspan(4, 1, nb)
    bw4, vecs4 = joint_split(4, moms, alg)
    pb4 = [pull_axis(v, bw4, idx14, nb, axis) for v in vecs4]
    E4 = match_ideal_dim(target_base(I3P_RAY, idx14), pb4, dp14, nb)
    _, _, dp24 = dspan(4, 2, nb)
    qd4 = qdim(vecs4, dp24)
    pbd4 = qdim(pb4, dp14)

    b16, idx16, dp16 = dspan(6, 1, nb)
    bw6, vecs6 = joint_split(6, moms, alg)
    _, _, dp26 = dspan(6, 2, nb)
    qd6 = qdim(vecs6, dp26)
    pb6 = [pull_axis(v, bw6, idx16, nb, axis) for v in vecs6]
    tA2 = target_base(A2_N1, idx16)
    tKdV = target_base(KDV6, idx16)
    E6A = match_ideal_dim(tA2, pb6, dp16, nb)
    E6K = match_ideal_dim(tKdV, pb6, dp16, nb)
    print(f"{label}: E4(I3')={E4}/{nb}, E6(A2)={E6A}/{nb}, "
          f"E6(KdV)={E6K}/{nb}; w4 joint/pullback Q-dims "
          f"{qd4}/{pbd4}; w6 kernel Q-dim {qd6}", flush=True)
    return E4, E6A, E6K, qd4, pbd4, qd6


# ------- ring sanity + validation controls at k = -1 ----------------
A1 = Algebra(ks + 1)
check("control algebra: dim 4, C(k)=0, u^2=-k/2, v^2=(k+2)/2 in "
      "normal form",
      A1.nb == 4 and A1.coords(ks + 1) == [F(0)] * 4
      and A1.coords(us**2 + ks / 2) == [F(0)] * 4
      and A1.coords(vs**2 - (ks + 2) / 2) == [F(0)] * 4)

E4p, E6Ap, E6Kp, _, _, qd6p = run_point(A1, 0,
                                        "control k=-1 phi-axis")
check("control: k=-1 weight-6 kernel Q-dim 8 = A-rank 2 (the "
      "enrichment independently verified in the phi certificate)",
      qd6p == 8)
check("control: phi-axis I3' ray matches at ALL components (E4 = 4; "
      "certified fact via independent stack)", E4p == 4)
check("control: phi-axis enhanced I5 absent at ALL components "
      "(E6 = 0; certified fact)", E6Ap == 0)
check("control: phi-axis KdV absent (E6_KdV = 0; certified fact)",
      E6Kp == 0)

E4t, _, _, _, _, _ = run_point(A1, 1, "control k=-1 theta-axis")
check("control: theta-axis I3' ray does NOT match at k=-1 (E4 = 0)",
      E4t == 0)

# ------- production: the theta-axis matching cubic ------------------
CUBIC = ks**3 + 23 * ks**2 + 55 * ks + 42
A3 = Algebra(CUBIC)
check("cubic algebra: dim 12; relations reduce to 0",
      A3.nb == 12 and A3.coords(CUBIC) == [F(0)] * 12
      and A3.coords(us**2 + ks / 2) == [F(0)] * 12)
check("cubic is irreducible over Q (no rational root escapes the "
      "algebra treatment)",
      sp.Poly(CUBIC, ks).is_irreducible)
zz = sp.Symbol('zz')
Z12 = (zz**12 - 6 * zz**10 + 388 * zz**8 - 1512 * zz**6
       + 1288 * zz**4 + 432 * zz**2 + 81)
check("A is a degree-12 FIELD: z = u+v satisfies the irreducible "
      "degree-12 polynomial (review-supplied), and dim matching "
      "forces A = Q(z) -- so E is 0 or A, and Q-dim 12 = A-rank 1",
      sp.Poly(Z12, zz).is_irreducible
      and all(c == 0 for c in A3.coords(Z12.subs(zz, us + vs))))

E4, E6A, E6K, qd4, pbd4, qd6 = run_point(A3, 1, "THETA at the cubic")
check("cubic: weight-4 joint kernel Q-dim 12 = A-rank 1 (so E4 = 12 "
      "is RAY EQUALITY, not mere membership)", qd4 == 12)
check("cubic: weight-4 pullback Q-dim 12 = A-rank 1 (the pullback "
      "is a line over A)", pbd4 == 12)

# committed verdict baked fail-closed (regression discipline): ray
# match on the full cubic, I5 miss everywhere, KdV miss everywhere,
# w6 fiber A-rank 1
check("committed verdict holds: (E4, E6_A2, E6_KdV) = (12, 0, 0)",
      (E4, E6A, E6K) == (12, 0, 0))
check("w6 kernel A-rank 1 on the cubic (Q-dim 12)", qd6 == 12)

print()
if failures:
    print(f"FAILURES: {failures}")
    print("GETMANOV THETA-AXIS TEST FAILED")
    sys.exit(1)

if E4 == 0:
    print("VERDICT: the theta-axis I3 candidacy DIES AT I3 -- the "
          "ratio-matching cubic does not deliver ray equality at any "
          "root or branch (E4 = 0).")
elif E4 == 12 and E6A == 0:
    print("VERDICT: theta-axis EXCLUDED AT I5 at every root and "
          "branch -- the I3' ray does match on the full cubic "
          f"(E4 = 12/12), but the enhanced I5 class is absent from "
          f"the weight-6 pullback span at every component "
          f"(E6 = 0/12; KdV also absent, {E6K}/12).  Combined with "
          "the phi-axis exclusion, NO axis-aligned Getmanov "
          "embedding reproduces the enhanced tower at two-charge "
          "depth; only tilted embeddings remain untested.")
elif E4 == 12 and E6A == 12:
    print("VERDICT: theta-axis MATCHES AT I5 at every root and "
          "branch -- investigate immediately.")
else:
    print(f"VERDICT: ROOT-DEPENDENT STRUCTURE -- E4 = {E4}/12, "
          f"E6(A2) = {E6A}/12, E6(KdV) = {E6K}/12; decompose the "
          "match ideal to identify which components carry it.")
print("GETMANOV THETA-AXIS TEST COMPLETE")
