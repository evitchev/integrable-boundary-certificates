"""Item 1 of the post-publication program: the hidden-KdV / sinh-Gordon
scan for the three Z2 x O(N) cylindrical families, ALONG the curve.

The submitted main paper states the caution (the c=1 lesson): an
all-odd multiplicity-one pattern is compatible with several affine
exponent sets, so coincidence with a background-charge KdV /
sinh-Gordon tower on a distinguished boson must be excluded before a
new label is assigned.  This certificate runs that exclusion for the
cylindrical families (mixed_engine.cyl_P4, Solutions 1-3) along the
double cover s^2 = (N-25)(N-1), both sheets -- and finds the
exclusion holds GENERICALLY, with an exactly computed finite
jumping locus where it fails in structured ways.

O(N)-invariance leaves exactly one distinguished boson, X, and one
invariant Feigin--Fuchs family T_Q = 1/2(dX)^2 + Q d^2X + T_Y, so the
only invariant hidden-KdV candidate is the ray 4 oint :T_Q^2: and the
only invariant sinh-Gordon pair is {e^{+-aX}}.

Part A (KdV ray).  The full-parity mixed class of 4:T_Q^2: mod d is
computed in-engine (component normal product at N = 1 and N = 2,
projected onto the invariant basis and reduced; the two agree).  The
background-charge term gives it a NONZERO X-odd class (Q d^2X dY.dY
is not a total derivative), so a Z2-even family can meet the ray only
at Q = 0; the even-sector cross term (11)_X(11)_Y then forces N = 0.
At every sampled curve point with N != 0 no (lam, Q) solves
cyl = lam * KdV(Q); at the formal point N = 0 Solutions 2 and 3 ARE
the c=1 KdV ray (positive control, Q = 0 found by the solver).

Part B (sinh-Gordon on X).  First-order conservation of I_3 under
V_a = e^{aX} in the mixed representation (X-letter contractions, Y
pairs as Taylor-shifted spectators, twisted derivative d + a dX, full
X-parity basis).  Symbolically in (N, s) each solution yields exactly
two condition polynomials; intersected with the curve they have a
FINITE common zero set (baked, regression-checked):
  Sol 1: (1,0), (28,-9)                           at a^2 = 2
  Sol 2: (0,5), (28,9), (-25/2,-45/2)             at a^2 = 2
  Sol 3: (-2,-9) at a in {+-1,+-2} -- the c_{1,2} quartet;
         (-25/2,-45/2), (0,-5), (1,0)              at a^2 = 2
Generic sampled points: no nonzero root.  Controls: the c=1 KdV ray
returns a^2 = 2; the cross term (11)_X(11)_Y alone is ALSO screened at
a^2 = 2 (double contraction; recorded as a fact); a KdV ray at Q^2 = 1
plus the cross term has no common root (negative control).

Part C (tower level).  At each special point the spin-5 kernel of
ad_{I_3} is computed and the X-screened SUBSPACE (screening is linear
in the density) is measured at the special a:
  * a^2 = 2 points with kernel dim 1 -- I_5 is screened too: the tower
    through spin 5 sits in the self-dual X-kernel, i.e. X enters only
    through T_X (the companion's L(1,0) structure on this curve);
    N = 28 is a physical integer.
  * Sol 3 at (-2,-9): kernel dim 2 = {pure-Y I_5} + {X-part screened
    at +-1, +-2}: the family DECOUPLES as c=-2 KdV on X (+) KdV on Y.
    Sol 2 at the SAME point also decouples, differently: cross term and
    (dX)^4 vanish, X-part = -2 (d2X)^2, the free quadratic hierarchy
    (found by the Groebner completeness check: a double root a = 0).
  * (-25/2,-45/2): kernel dim 2 (multiplicity jump) AND screened
    subspace dim 1 (common quotient coordinates; the combination
    P_1 - (48/17) P_0 of the stored kernel basis): one spin-5 class
    continues the self-dual hierarchy, the other is new -- a
    tower-level embedding point that is also a jump point (Codex
    audit correction 2026-08-21: a per-current numerator rescaling
    had reported 0).

Fail-closed: exits nonzero unless every check passes.
"""
import sys
from fractions import Fraction as F
from itertools import combinations
from math import factorial

import sympy as sp

from invariant_engine import _rref, canon_pair, curve_point, gen_inv_basis
from mixed_engine import (MixedSector, _M, canon_mixed, cyl_P4,
                          d_mixed_mono, gen_mixed_basis)
from ope_engine import Sector, d_mono, gen_basis
from super_engine import normal_prod

failures = []


def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok:
        failures.append(name)


def mixed_to_comp(mono, N):
    xs, yp = mono
    cur = {tuple(sorted((o, 0) for o in xs)): F(1)}
    for a, b in yp:
        new = {}
        for m, c in cur.items():
            for mu in range(1, N + 1):
                mm = tuple(sorted(m + ((a, mu), (b, mu))))
                new[mm] = new.get(mm, 0) + c
        cur = new
    return cur


# ====================================================================
# Part A: the invariant KdV ray, computed in-engine
# ====================================================================
Q = sp.Symbol('Q')


def gen_mixed_full(s):
    """All mixed monomials of spin s (any number of X letters)."""
    def xparts(total):
        res = [[]] if total == 0 else []

        def rec(mx, remaining, acc):
            if remaining == 0:
                res.append(list(acc))
                return
            for v in range(min(mx, remaining), 0, -1):
                acc.append(v)
                rec(v, remaining - v, acc)
                acc.pop()
        rec(total, total, [])
        return res
    if s == 0:
        return [canon_mixed((), ())]
    out = []
    for xspin in range(0, s + 1):
        for xs in xparts(xspin):
            rem = s - xspin
            if rem == 0:
                if xs:
                    out.append(canon_mixed(tuple(xs), ()))
                continue
            for yp in gen_inv_basis(rem):
                out.append(canon_mixed(tuple(xs), yp))
    return sorted(set(out))


class FullSector:
    """Mixed monomials of both X-parities, modulo total derivatives."""
    def __init__(self, s):
        self.basis = gen_mixed_full(s)
        self.index = {m: i for i, m in enumerate(self.basis)}
        cols = []
        for b in gen_mixed_full(s - 1):
            col = [F(0)] * len(self.basis)
            for mono, c in d_mixed_mono(b).items():
                col[self.index[mono]] += F(c)
            if any(col):
                cols.append(col)
        rows, piv = _rref([list(c) for c in cols], len(self.basis))
        self.rows = rows[:len(piv)]
        self.piv = piv

    def vec(self, cur):
        v = [sp.Integer(0)] * len(self.basis)
        for m, c in cur.items():
            v[self.index[m]] += sp.sympify(c)
        return v

    def reduce(self, v):
        v = list(v)
        for row, p in zip(self.rows, self.piv):
            if v[p] != 0:
                f = v[p]
                v = [sp.expand(x - f * sp.Rational(y.numerator, y.denominator))
                     for x, y in zip(v, row)]
        return v


full4 = FullSector(4)


def kdv_mixed_class(N):
    """Reduced full-mixed coordinates of 4:T_Q^2: mod d, from the
    component normal product at N Y-bosons (X = field 0)."""
    nf = N + 1
    T = {}
    for j in range(nf):
        T[((0, 1, j), (0, 1, j))] = sp.Rational(1, 2)
    T[((0, 2, 0),)] = Q
    TT = normal_prod(T, T)
    comp = {}
    for mono, c in TT.items():
        key = tuple(sorted((m, j) for _, m, j in mono))
        comp[key] = sp.expand(comp.get(key, 0) + 4 * c)
    basis4 = gen_basis(4, nf, z2=False)
    idx = {m: i for i, m in enumerate(basis4)}
    target = [sp.Integer(0)] * len(basis4)
    for m, c in comp.items():
        target[idx[m]] += c
    cols, names = [], []
    for mm in full4.basis:
        col = [sp.Integer(0)] * len(basis4)
        for m, c in mixed_to_comp(mm, N).items():
            col[idx[tuple(sorted(m))]] += sp.Integer(c)
        cols.append(col)
        names.append(('c', mm))
    for u in gen_basis(3, nf, z2=False):
        col = [sp.Integer(0)] * len(basis4)
        for m, c in d_mono(u).items():
            col[idx[tuple(sorted(m))]] += sp.Integer(c)
        cols.append(col)
        names.append(('d', u))
    A = sp.Matrix([[cols[j][i] for j in range(len(cols))]
                   for i in range(len(basis4))])
    b = sp.Matrix(target)
    syms = sp.symbols(f'x0:{len(cols)}')
    sol = sp.linsolve((A, b), *syms)
    if sol == sp.EmptySet:
        return None
    vec = list(sol)[0]
    free = set().union(*[sp.sympify(x).free_symbols for x in vec]) - {Q}
    cur = {}
    for (kind, key), val in zip(names, vec):
        if kind == 'c':
            v = sp.expand(sp.sympify(val).subs({f: 0 for f in free}))
            if v != 0:
                cur[key] = v
    return full4.reduce(full4.vec(cur))


k1 = kdv_mixed_class(1)
k2 = kdv_mixed_class(2)
check("A0: 4:T_Q^2: projects onto the invariant sector at N=1 and N=2",
      k1 is not None and k2 is not None)
check("A1: KdV reduced class N-independent (N=1 == N=2, symbolic in Q)",
      k1 == k2)
odd_idx = [i for i, m in enumerate(full4.basis) if len(m[0]) % 2 == 1]
even_idx = [i for i, m in enumerate(full4.basis) if len(m[0]) % 2 == 0]
odd_part = [k1[i] for i in odd_idx]
check("A2: the odd-in-X part of the KdV class is exactly Q times a "
      "NONZERO class (background charge breaks X -> -X once Y is coupled)",
      any(odd_part) and all(sp.expand(c / Q).free_symbols - {Q} == set()
                            and sp.expand(c / Q).is_polynomial(Q)
                            and sp.Poly(sp.expand(c / Q), Q).degree() == 0
                            for c in odd_part if c != 0))
print("      odd part of 4:T_Q^2: mod d:",
      {str(full4.basis[i]): str(k1[i]) for i in odd_idx if k1[i] != 0})
print("      even part:",
      {str(full4.basis[i]): str(k1[i]) for i in even_idx if k1[i] != 0})


def kdv_match(cur):
    """Solve reduce(cur) = lam * kdv(Q) exactly in (lam, Q); list of dicts."""
    vc = full4.reduce(full4.vec(cur))
    lam = sp.Symbol('lam')
    eqs = [sp.expand(lam * k - c) for k, c in zip(k1, vc)]
    eqs = [e for e in eqs if e != 0]
    return sp.solve(eqs, [lam, Q], dict=True)


SAMPLE_T = [F(2), F(3), F(5), F(7), F(1, 2), F(1, 3), F(2, 3), F(-2),
            F(-3), F(-5), F(-1, 2), F(-2, 3), F(0)]
points = [curve_point(t) for t in SAMPLE_T] + [(F(1), F(0))]
for sol_id in (1, 2, 3):
    for N, s in points:
        if N == 0:
            continue      # the formal N = 0 points are the A4 controls
        m = kdv_match(cyl_P4(sol_id, N, s))
        check(f"A3: Solution {sol_id} at N={N}, s={s}: NOT the KdV ray "
              f"for any (lam, Q)", m == [])

# positive controls at the formal N = 0 point (no Y bosons)
for sol_id, s0 in ((2, F(5)), (3, F(-5))):
    m = kdv_match(cyl_P4(sol_id, F(0), s0))
    ok = bool(m) and all(d[Q] == 0 for d in m)
    check(f"A4 control: Solution {sol_id} at N=0 IS the c=1 KdV ray "
          f"(solver finds it, Q = 0: {m})", ok)
m = kdv_match(cyl_P4(1, F(0), F(5)))
check("A4 control: Solution 1 at N=0 is NOT KdV (cross term 6 != 2)",
      m == [])

# analytic obstruction, symbolic on the curve
Ns, ss = sp.symbols('N s')
curve = sp.Eq(ss**2, (Ns - 25) * (Ns - 1))
c2 = sp.Eq((11 + Ns + ss) / 8, 2)
c3 = sp.Eq(7 - Ns + ss, 2)
sol2 = sp.solve([curve, c2], [Ns, ss], dict=True)
sol3 = sp.solve([curve, c3], [Ns, ss], dict=True)
check("A5: even-sector cross-term coincidence (coefficient 2) on the "
      "curve forces N = 0 for Solutions 2 and 3",
      bool(sol2) and bool(sol3)
      and all(d[Ns] == 0 for d in sol2) and all(d[Ns] == 0 for d in sol3))

# ====================================================================
# Part B: sinh-Gordon on the distinguished boson X
# ====================================================================
def apply_d(cur):
    out = {}
    for m, c in cur.items():
        for mm, k in d_mixed_mono(m).items():
            out[mm] = out.get(mm, 0) + c * k
    return {m: c for m, c in out.items() if c}


def residue_x(P, a):
    """Res_{z->w} P(z) e^{aX}(w) in the mixed representation (V implicit)."""
    out = {}
    for (xs, yp), cP in P.items():
        cP = (sp.Rational(cP.numerator, cP.denominator)
              if isinstance(cP, F) else sp.sympify(cP))
        n = len(xs)
        for k in range(1, n + 1):
            for sel in combinations(range(n), k):
                num = sp.Integer(1)
                D = 0
                for i in sel:
                    m = xs[i]
                    num *= a * (-1) ** (m - 1) * factorial(m - 1)
                    D += m
                rest_xs = tuple(xs[i] for i in range(n) if i not in sel)
                rest = canon_mixed(rest_xs, yp)
                T = D - 1
                cur = {rest: sp.Integer(1)}
                for _ in range(T):
                    cur = apply_d(cur)
                for m, c in cur.items():
                    out[m] = out.get(m, 0) + cP * num * c / factorial(T)
    return {m: sp.expand(c) for m, c in out.items() if sp.expand(c) != 0}


def twisted_d_mixed(mono, a):
    out = {m: sp.Integer(c) for m, c in d_mixed_mono(mono).items()}
    xs, yp = mono
    new = canon_mixed(xs + (1,), yp)
    out[new] = out.get(new, 0) + a
    return out


def x_screening_conditions(P, a, w):
    """Polynomial conditions in a for Res(P V_a) to be a twisted total
    derivative; P a mixed current of spin w."""
    res = residue_x(P, a)
    s = w - 1
    tgt = []
    for t in range(0, s + 1):
        tgt += gen_mixed_full(t)
    lower = gen_mixed_full(s - 1)
    index = {m: i for i, m in enumerate(tgt)}
    cols = []
    for u in lower:
        col = [sp.Integer(0)] * len(tgt)
        for m, c in twisted_d_mixed(u, a).items():
            col[index[m]] += c
        cols.append(col)
    vec = [sp.Integer(0)] * len(tgt)
    for m, c in res.items():
        vec[index[m]] += c
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
                conds.append(sp.numer(sp.together(e)))
    return conds


def nonzero_roots(conds, a):
    """Common roots of the conditions, a != 0."""
    if not conds:
        return 'ALL'   # no condition at all: every a works
    g = conds[0]
    for c in conds[1:]:
        g = sp.gcd(g, c)
    g = sp.factor(g)
    roots = [r for r in sp.solve(g, a) if r != 0]
    return roots


a = sp.Symbol('a')
# controls
P_kdv = {_M((1, 1, 1, 1)): F(1), _M((2, 2)): F(-2)}
r = nonzero_roots(x_screening_conditions(P_kdv, a, 4), a)
check("B0 control: c=1 KdV ray (dX)^4 - 2(d2X)^2 is screened by "
      "e^{aX} exactly at a^2 = 2", set(r) == {sp.sqrt(2), -sp.sqrt(2)})
P_sum = dict(P_kdv)
P_sum[_M((), [(1, 1), (1, 1)])] = F(1)
P_sum[_M((), [(2, 2)])] = F(-2)
r = nonzero_roots(x_screening_conditions(P_sum, a, 4), a)
check("B1 control: decoupled X-KdV + Y-KdV keeps a^2 = 2 (Y spectators "
      "handled)", set(r) == {sp.sqrt(2), -sp.sqrt(2)})
P_cross_only = {_M((1, 1), [(1, 1)]): F(1)}
r = nonzero_roots(x_screening_conditions(P_cross_only, a, 4), a)
check("B2 fact: the cross term (11)_X(11)_Y alone is screened by e^{aX} "
      "at a^2 = 2 (double contraction)", set(r) == {sp.sqrt(2), -sp.sqrt(2)})
P_neg = {_M((1, 1, 1, 1)): F(1), _M((2, 2)): F(2),      # KdV ray at Q^2 = 1
         _M((1, 1), [(1, 1)]): F(1)}
r = nonzero_roots(x_screening_conditions(P_neg, a, 4), a)
check("B2 negative control: KdV ray at Q^2 = 1 plus the cross term has "
      "no common screening", r == [])

# exact locus, symbolic in (N, s) on the curve
Ns, ss = sp.symbols('N s')


def cyl_sym(sol):
    base = {_M((), [(1, 1), (1, 1)]): sp.Integer(1),
            _M((), [(2, 2)]): sp.Integer(-2)}
    if sol == 1:
        extra = {_M((2, 2)): (-9 - 3 * Ns - ss) / 6,
                 _M((1, 1), [(1, 1)]): sp.Integer(6),
                 _M((1, 1, 1, 1)): (2 + Ns + ss) / 3}
    elif sol == 2:
        extra = {_M((2, 2)): -(409 + 46 * Ns + Ns**2 + (Ns - 5) * ss) / 192,
                 _M((1, 1), [(1, 1)]): (11 + Ns + ss) / 8,
                 _M((1, 1, 1, 1)): (117 + 2 * Ns + Ns**2 + (15 + Ns) * ss) / 192}
    else:
        extra = {_M((2, 2)): (13 - 27 * Ns + 2 * Ns**2 + (5 - 2 * Ns) * ss) / 6,
                 _M((1, 1), [(1, 1)]): (7 - Ns + ss),
                 _M((1, 1, 1, 1)): sp.Integer(1)}
    base.update(extra)
    return base


# transcription cross-check of the symbolic copies against cyl_P4
okc = True
for sol_id in (1, 2, 3):
    for N, s in ((F(28), F(9)), (F(-7), F(-16))):
        ref = cyl_P4(sol_id, N, s)
        sym = {m: sp.nsimplify(c.subs({Ns: sp.Rational(N.numerator, N.denominator),
                                        ss: sp.Rational(s.numerator, s.denominator)}))
               for m, c in cyl_sym(sol_id).items()}
        okc &= all(sym.get(m, 0) == sp.Rational(ref[m].numerator, ref[m].denominator)
                   for m in ref) and all(sym[m] == 0 for m in set(sym) - set(ref))
check("B3: symbolic cyl_P4 copies agree with mixed_engine.cyl_P4", okc)

curve = ss**2 - (Ns - 25) * (Ns - 1)


def linear_roots(poly, var):
    """Roots of a univariate polynomial read off from its linear factors
    over Q; returns (roots, complete) with complete=False if any factor
    is nonlinear (enumeration would then be incomplete -- fail closed)."""
    _, factors = sp.factor_list(sp.Poly(poly, var))
    roots, complete = [], True
    for f, mult in factors:
        f = sp.Poly(f, var)
        if f.degree() == 1:
            c1, c0 = f.all_coeffs()
            roots.append(sp.Rational(-c0, c1) if c1 != 0 else None)
        elif f.degree() >= 2:
            complete = False
    return roots, complete


def enumerate_locus(G):
    """Explicit back-substitution through a lex Groebner basis in
    (a, N, s): no root-finding solver is called anywhere."""
    exprs = list(G.exprs)
    s_polys = [g for g in exprs if g.free_symbols == {ss}]
    n_polys = [g for g in exprs if g.free_symbols <= {Ns, ss} and Ns in g.free_symbols]
    a_polys = [g for g in exprs if a in g.free_symbols]
    if len(s_polys) != 1 or not n_polys or not a_polys:
        return {}, False
    s_roots, complete = linear_roots(s_polys[0], ss)
    if not complete:
        return {}, False
    loc = {}
    for sv in s_roots:
        # N: every N-element must be linear in N at this s, all agreeing
        nvals = set()
        for g in n_polys:
            pg = sp.Poly(sp.expand(g.subs(ss, sv)), Ns)
            if pg.degree() != 1:
                return {}, False
            c1, c0 = pg.all_coeffs()
            nvals.add(sp.Rational(-c0, c1))
        if len(nvals) != 1:
            return {}, False
        nv = nvals.pop()
        if sp.expand(curve.subs({Ns: nv, ss: sv})) != 0:
            return {}, False
        # a: common roots of the a-elements, as polynomials in b = a^2
        b = sp.Symbol('b')
        gcd = None
        for g in a_polys:
            pa = sp.Poly(sp.expand(g.subs({Ns: nv, ss: sv})), a)
            if pa.is_zero:
                continue
            if any(m[0] % 2 for m in pa.monoms()):
                return {}, False        # not even in a
            pb = sp.Poly(sp.expand(pa.as_expr().subs(a**2, b)), b)
            gcd = pb if gcd is None else sp.gcd(gcd, pb)
        if gcd is None:
            return {}, False            # no a-constraint recovered: fail closed
        if gcd.degree() < 1:
            continue                    # nonzero constant: no a-root here
        b_roots, complete = linear_roots(gcd.as_expr(), b)
        if not complete:
            return {}, False
        # every nonzero a^2 counts -- negative b is a valid complex
        # screening momentum a = +-sqrt(b); only a = 0 (the trivial
        # vertex operator) is discarded, so no component can hide
        vals = {bv for bv in b_roots if bv is not None and bv != 0}
        if vals:
            loc[(nv, sv)] = vals
    return loc, True


EXPECTED_LOCUS = {
    1: {(1, 0): {2}, (28, -9): {2}},
    2: {(0, 5): {2}, (28, 9): {2}, (sp.Rational(-25, 2), sp.Rational(-45, 2)): {2}},
    3: {(-2, -9): {1, 4}, (sp.Rational(-25, 2), sp.Rational(-45, 2)): {2},
        (0, -5): {2}, (1, 0): {2}},
}
found_locus = {}
for sol_id in (1, 2, 3):
    conds = [sp.factor(c) for c in
             x_screening_conditions(cyl_sym(sol_id), a, 4)]
    check(f"B4: Solution {sol_id} yields exactly two nontrivial condition "
          f"polynomials in a over Q(N,s)", len(conds) == 2
          and all(sp.Poly(c, a).degree() >= 1 for c in conds))
    sols = sp.solve(conds + [curve], [a, Ns, ss], dict=True)
    loc = {}
    for d in sols:
        if d[a] == 0:
            continue
        key = (sp.nsimplify(d[Ns]), sp.nsimplify(d[ss]))
        loc.setdefault(key, set()).add(sp.nsimplify(d[a]**2))
    found_locus[sol_id] = loc
    check(f"B5: Solution {sol_id} exact X-screening locus on the curve = "
          f"{ {k: sorted(v) for k, v in loc.items()} }",
          loc == EXPECTED_LOCUS[sol_id])
    # branch-complete, solver-independent: the ideal <c1/a, c2/a, curve>
    # (a != 0 divided out) is zero-dimensional, and its lex Groebner
    # basis is triangular -- every point is read off, none can hide.
    red = [sp.cancel(c / a) for c in conds]
    check(f"B5g: Solution {sol_id}: both conditions carry the factor a "
          f"(a = 0 trivial) and are polynomial after division",
          all(sp.denom(sp.together(r)) == 1 for r in red))
    G = sp.groebner(red + [curve], a, Ns, ss, order='lex')
    gloc, shape_ok = enumerate_locus(G)
    check(f"B5g: Solution {sol_id}: ideal zero-dimensional, Groebner basis "
          f"triangular (terminal s-polynomial splits into linear factors "
          f"over Q; N-element linear; a-elements even in a), and explicit "
          f"back-substitution reproduces the baked locus",
          G.is_zero_dimensional and shape_ok
          and gloc == EXPECTED_LOCUS[sol_id])

# generic sampled points (special points excluded) -- no nonzero root
special_pts = {(sp.nsimplify(sp.Rational(N.numerator, N.denominator)),
                sp.nsimplify(sp.Rational(s.numerator, s.denominator)))
               for N, s in points}
for sol_id in (1, 2, 3):
    for N, s in points:
        key = (sp.Rational(N.numerator, N.denominator),
               sp.Rational(s.numerator, s.denominator))
        if key in EXPECTED_LOCUS[sol_id]:
            continue
        r = nonzero_roots(x_screening_conditions(cyl_P4(sol_id, N, s), a, 4),
                          a)
        check(f"B6: Solution {sol_id} at N={N}, s={s} (generic): no "
              f"nonzero e^{{aX}} screening of I_3", r == [])

# ====================================================================
# Part C: tower level at the special points
# ====================================================================
from mixed_engine import genuine_kernel_mixed


def screened_subspace_dim(currents, aval, w):
    """Dimension of the subspace of span(currents) conserved under
    e^{aval X}.  Screening is linear in the density, so build the
    twisted-derivative matrix ONCE at a = aval, take a basis of its left
    nullspace (common quotient coordinates), apply it to every current's
    residue vector, and count the nullity of the stacked coordinates.
    (A per-current numerator rescaling does NOT preserve dependence --
    the v1 bug caught by the Codex audit.)"""
    s = w - 1
    tgt = []
    for tt in range(0, s + 1):
        tgt += gen_mixed_full(tt)
    lower = gen_mixed_full(s - 1)
    index = {m: i for i, m in enumerate(tgt)}
    cols = []
    for u in lower:
        col = [sp.Integer(0)] * len(tgt)
        for m, c in twisted_d_mixed(u, a).items():
            col[index[m]] += sp.sympify(c).subs(a, aval)
        cols.append(col)
    M = sp.Matrix(cols).T
    Y = M.T.nullspace()
    coords = []
    for P in currents:
        r = [sp.Integer(0)] * len(tgt)
        for m, c in residue_x(P, a).items():
            r[index[m]] += c.subs(a, aval)
        r = sp.Matrix(r)
        coords.append([sp.simplify((y.T * r)[0]) for y in Y])
    if not coords:
        return 0
    return len(currents) - sp.Matrix(coords).T.rank()


C_EXPECT = {   # (sol, N, s): (I5 kernel dim, screened dim at a=sqrt2 or at 1 and 2)
    (1, 1, 0): (1, 1), (1, 28, -9): (1, 1),
    (2, 0, 5): (1, 1), (2, 28, 9): (1, 1),
    (3, 0, -5): (1, 1), (3, 1, 0): (1, 1),
}
for (sol_id, N, s), (kdim, sdim) in C_EXPECT.items():
    Nf, sf = F(N), F(s)
    ker = genuine_kernel_mixed(cyl_P4(sol_id, Nf, sf), gen_mixed_basis(6), Nf)
    d = screened_subspace_dim(ker, sp.sqrt(2), 6) if ker else 0
    check(f"C1: Solution {sol_id} at ({N},{s}): I_5 kernel dim {len(ker)}, "
          f"e^{{sqrt2 X}}-screened subspace dim {d} -- tower-level "
          f"self-dual embedding", (len(ker), d) == (kdim, sdim))

# decoupling point: Solution 3 at (-2, -9)
ker = genuine_kernel_mixed(cyl_P4(3, F(-2), F(-9)), gen_mixed_basis(6), F(-2))
pureY = [P for P in ker if all(len(m[0]) == 0 for m in P)]
d1 = screened_subspace_dim(ker, sp.Integer(1), 6)
d2 = screened_subspace_dim(ker, sp.Integer(2), 6)
cross = cyl_P4(3, F(-2), F(-9)).get(_M((1, 1), [(1, 1)]), 0)
check("C2: Solution 3 at (-2,-9) DECOUPLES: cross term 0, I_5 kernel "
      f"dim {len(ker)} containing a pure-Y charge, fully screened at "
      f"a = 1 and a = 2 (dims {d1}, {d2}) -- c=-2 KdV on X (+) KdV on Y",
      cross == 0 and len(ker) == 2 and len(pureY) == 1
      and d1 == 2 and d2 == 2)

# Solution 2 at the same point: a decoupling of a different type
P2 = cyl_P4(2, F(-2), F(-9))
ker2 = genuine_kernel_mixed(P2, gen_mixed_basis(6), F(-2))
pureY2 = [P for P in ker2 if all(len(m[0]) == 0 for m in P)]
pureX2 = [P for P in ker2 if all(len(m[1]) == 0 for m in P)]
check("C2b: Solution 2 at (-2,-9) ALSO decouples: cross term and (dX)^4 "
      "both vanish, X-part = -2 (d2X)^2 (the free quadratic hierarchy), "
      f"I_5 kernel dim {len(ker2)} = pure-Y (+) pure-X",
      P2.get(_M((1, 1), [(1, 1)]), 0) == 0
      and P2.get(_M((1, 1, 1, 1)), 0) == 0
      and P2.get(_M((2, 2))) == -2
      and len(ker2) == 2 and len(pureY2) == 1 and len(pureX2) == 1)

# multiplicity-jump point (-25/2, -45/2): kernel dim 2, screened dim 1
for sol_id in (2, 3):
    ker = genuine_kernel_mixed(cyl_P4(sol_id, F(-25, 2), F(-45, 2)),
                               gen_mixed_basis(6), F(-25, 2))
    d = screened_subspace_dim(ker, sp.sqrt(2), 6) if ker else 0
    check(f"C3: Solution {sol_id} at (-25/2,-45/2): I_5 multiplicity JUMPS "
          f"to 2 and the e^{{sqrt2 X}}-screened subspace has dim 1 "
          f"(tower-level embedding inside an enlarged kernel)",
          (len(ker), d) == (2, 1))
    if len(ker) == 2:
        comb = {}
        for m, c in ker[1].items():
            comb[m] = comb.get(m, F(0)) + c
        for m, c in ker[0].items():
            comb[m] = comb.get(m, F(0)) - F(48, 17) * c
        comb = {m: c for m, c in comb.items() if c}
        check(f"C4: Solution {sol_id}: the combination P_1 - (48/17) P_0 of the "
              f"stored kernel basis is itself screened at a = sqrt2",
              screened_subspace_dim([comb], sp.sqrt(2), 6) == 1)

print()
if failures:
    print(f"CYLINDRICAL HIDDEN-KDV SCAN FAILED ({len(failures)}):")
    for f in failures:
        print(" -", f)
    sys.exit(1)
print("CYLINDRICAL HIDDEN-KDV SCAN CERTIFIED (no invariant KdV ray at any "
      "N != 0; X-screening only on the exact finite locus, with tower-level "
      "self-dual embeddings, the c=-2 decoupling point, and the "
      "multiplicity-jump point; controls positive)")
