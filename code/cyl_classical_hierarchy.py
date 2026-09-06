"""Classical-hierarchy certificate (2026-09-05): the top momentum degree of the certified
cylindrical eigenvalue tables is a two-component hydrodynamic hierarchy determined by its
weight-4 density, and the pole lattices of cyl_pole_lattice are a theorem of that hierarchy
(ansatz note secs. 23-24).

Setting.  Write the top-degree part of the weight-W eigenvalue I(X, Y) (X = P^2, Y = Q^2,
article momenta) as a density F_W(x, r) = sum c_ab x^{2a} r^{2b} in the velocities x = dX,
r = |dY| (rho = r^2 = dY.dY).  O(N)-invariant densities of the velocities see the free-field
bracket as the two-component flat bracket diag(d, d) in (x, r) [diag(d, 4 rho d + 2 rho') in
(x, rho)], with no N dependence.  A density F is conserved by the flow of H iff
    H_xr (F_xx - F_rr) = (H_xx - H_rr) F_xr           (exactness of the flux 1-form).
The weight-4 density is x^4 + b x^2 r^2 + c r^4 after normalization; each family is a curve
(b(t), c(t)) in the (b, c) plane, the paperclip family of hep-th/0404195 eq. (59) another curve.

Checks (exact, sympy over Q(t) / Q(b, c)):
  P1  the top-degree densities of the PINNED tables Poisson-commute pairwise (weights 4-10,
      three solutions): all six pairs {F_a, F_b}, a < b in {4, 6, 8, 10};
  P2  the conservation PDE of the weight-4 density has a ONE-dimensional space of polynomial
      densities at each weight 6-10 for each solution, spanned by the certified density;
  P3  universality, GENERICALLY: over the function field Q(b, c) the conserved space is
      one-dimensional at weights 6-12; every denominator of the universal profiles is a
      polynomial in b alone, with zeros exactly at b = 6j/(3j-s) for j = 1..(s-1)/2
      (s = W-1; the j = s/3 factor is at infinity), i.e. on the Beta lattice s m + 2j = 0 of
      the exponent  m(b) = -2b/(3b - 6);
  P3b the exceptional locus (Codex rounds 25-26), FAIL-CLOSED: the constraint matrix drops rank
      on the line b = 0 (H_xr = 0 there, so the PDE forces F_xr = 0 and the even densities are
      A x^W + B r^W -- nullity exactly 2, checked at weights 6-12) and at finitely many isolated
      points, enumerated EXHAUSTIVELY (the lex Groebner ideal of the maximal minors with the
      factor b removed must be zero-dimensional; its univariate element is factored over Q and
      back-substituted; each point must annihilate every minor with nullity exactly 2):
      weight 6 {(-3, 1)}; weight 8 {(-3/2, 1/8), (-3/2, 1), (-12, 8)}; weight 10 {(-1, 0),
      (-1, 1/4), (-1, 1), (-4, 1), (-4, 4)}.  On the family curves these are the isolated
      NON-DECOUPLING classical jumping points: Solution 3 at t = 5/3 (weight 6), 11/5 (8),
      7/5 and 17/7 (10); Solution 2 at t = 5/3 (6), 2 (8); Solution 1 none.  The decoupling
      component b = 0 (nullity 2 at every weight) is reached by Solution 3 at t = 3, and
      Solution 2 reaches the decoupled point at infinity of the plane (c -> infinity, density
      r^W alone, nullity 2) at t = 3 as well -- the classical face of the all-spin decoupling
      drop of sheets 2 and 3.  The certified non-decoupling quantum jumping points at these
      weights -- t = 5/3 at weight 6 on the merged sheets 2 and 3 and t = 7/5 at weight 10 on
      sheet 3 (cyl_strata_points), none on sheet 1 -- are all among the isolated points;
  P4  the exponent-mixing law: m(u_4) with u_4 the certified weight-4 mixing coefficient of
      each solution (cyl_mo2_profiles' PINS, reconstructed exactly by rational fitting) equals the pinned exponent of cyl_pole_lattice for
      Solutions 1, 3 and its lattice dual -2 - m for Solution 2; for the paperclip family
      (eq. (59), all n) it gives -2 - n;
  P5  decoupling: at b = 0 the universal densities are x^W + c^{W/2-1} r^W (weights 6-12) --
      the t = 3 point of Solution 3 (b = 0, c = 1) is the decoupled system x^W + r^W;
  P6  the (b, c)-plane geometry: the paperclip curve is c = b(b-4)/(4(b-3)) (reproduces
      eq. (59)'s (b(n), c(n))); the family curves are c = 1, c = b/6, c = b^2/(6-b); their
      intersections with the paperclip curve are t = 5 and t -> infinity for Solution 3
      (b = 2, 6: the recorded survivors of the two-term exclusion), t = 1 and infinity for
      Solution 1, t = -5, 1 and infinity for Solution 2;
  P7  negative controls: a wrong bracket coefficient breaks P1; a perturbed exponent law
      breaks P4 (guards live).
Fail-closed: any failure exits 1 without the CERTIFIED marker.  Runtime ~10 min.
"""
import sys, os, time
import sympy as sp
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, HERE)
from cyl_first_order import cyl, X, Y, t, IPC, n, fit_rational
from cyl_mo2_profiles import PINS

def pinned_coefficient(sol, W, i):
    """the certified profile coefficient (a callable of t, Fraction-valued) as an exact rational function of t:
    sampled at eight points, fitted at the lowest degree that reproduces them, verified on four more"""
    fn = PINS[sol](W)[i]; pts = [sp.Rational(v) for v in (2, 7, -4, sp.Rational(13, 3), -sp.Rational(9, 2), sp.Rational(23, 7), 10, -sp.Rational(11, 3))]
    samples = {}
    for tv in pts:
        try:
            samples[tv] = sp.Rational(fn(Fr(int(tv.p), int(tv.q))))
        except ZeroDivisionError:
            pass
    f = next(f_ for f_ in (fit_rational(samples, dg) for dg in range(0, 5)) if f_ is not None)
    for tv in (sp.Rational(17, 5), 12, -sp.Rational(7, 3), sp.Rational(31, 4)):
        try:
            if f.subs(t, tv) != sp.Rational(fn(Fr(int(sp.Rational(tv).p), int(sp.Rational(tv).q)))):
                raise ValueError(f"rational fit of the pinned profile fails a hold-out point: solution {sol}, weight {W}, index {i}")
        except ZeroDivisionError:
            pass
    return f

fails = []
def require(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg, flush=True)
    if not cond:
        fails.append(msg)

x, r, z, b, c = sp.symbols('x r z b c')
M_PINNED = {1: (1 - t) / (t + 3), 2: (2 * t - 14) / (t + 5), 3: (3 - t) / (t - 5)}   # cyl_pole_lattice

def density(sol, W):
    P = sp.Poly(sp.expand(cyl(sol, W)), X, Y); k = W // 2
    F = sum(cf * x**(2 * a) * r**(2 * bb) for (a, bb), cf in P.terms() if a + bb == k)
    return sp.expand(F / P.coeff_monomial(X**k))

def bracket_vanishes(F, G, coef=(1, 1)):
    """{int F, int G} under diag(coef[0] d, coef[1] d) in (x, r): the bracket density is a total derivative
    iff its Euler-Lagrange derivatives vanish."""
    xf, rf = sp.Function('xf')(z), sp.Function('rf')(z); sub = {x: xf, r: rf}
    dens = sp.expand(coef[0] * sp.diff(F, x).subs(sub) * sp.diff(sp.diff(G, x).subs(sub), z) + coef[1] * sp.diff(F, r).subs(sub) * sp.diff(sp.diff(G, r).subs(sub), z))
    from sympy.calculus.euler import euler_equations
    return all(sp.simplify(e.lhs) == 0 for e in euler_equations(dens, [xf, rf], z))

def pde(F, H):
    return sp.expand(sp.diff(H, x, r) * (sp.diff(F, x, 2) - sp.diff(F, r, 2)) - (sp.diff(H, x, 2) - sp.diff(H, r, 2)) * sp.diff(F, x, r))

def conserved(H, W):
    """(dimension of the polynomial conserved space, normalized profile [a_1..a_k] with x^W coefficient 1)"""
    k = W // 2; cs = [sp.Symbol(f'a{i}') for i in range(1, k + 1)]
    F = x**W + sum(cs[i - 1] * x**(W - 2 * i) * r**(2 * i) for i in range(1, k + 1))
    eqs = sp.Poly(pde(F, H), x, r).coeffs()
    gen = sp.linsolve(eqs, cs); (vals,) = gen
    free = {s_ for v in vals for s_ in v.free_symbols if str(s_).startswith('a')}
    sol = sp.solve(eqs, cs, dict=True)
    return len(free) + 1, [sp.factor(sol[0][ci]) for ci in cs] if sol else None

if __name__ == "__main__":
    t0 = time.time()
    dens = {sol: {W: density(sol, W) for W in (4, 6, 8, 10)} for sol in (1, 2, 3)}
    # P1
    pairs = [(a_, b_) for a_ in (4, 6, 8, 10) for b_ in (4, 6, 8, 10) if a_ < b_]
    ok1 = all(bracket_vanishes(dens[sol][a_], dens[sol][b_]) for sol in (1, 2, 3) for a_, b_ in pairs)
    require(ok1, "P1: the top-degree densities of the pinned tables Poisson-commute under diag(d, d) in (x, r): all six pairs {F_a, F_b}, a < b in {4, 6, 8, 10}, three solutions")
    # P2
    ok2 = True
    for sol in (1, 2, 3):
        for W in (6, 8, 10):
            dim, prof = conserved(dens[sol][4], W); k = W // 2
            F = x**W + sum(prof[i - 1] * x**(W - 2 * i) * r**(2 * i) for i in range(1, k + 1))
            ok2 &= (dim == 1) and sp.simplify(sp.expand(F - dens[sol][W])) == 0
    require(ok2, "P2: the weight-4 density's conservation PDE has a one-dimensional polynomial solution space at weights 6, 8, 10, spanned by the certified density (three solutions)")
    # families in the (b, c) plane (used by P3b and P6): the certified weight-4 mixing and r^4 coefficients
    u4 = {sol: pinned_coefficient(sol, 4, 1) for sol in (1, 2, 3)}
    fam = {sol: (u4[sol], pinned_coefficient(sol, 4, 2)) for sol in (3, 1, 2)}
    # P3
    H = x**4 + b * x**2 * r**2 + c * r**4; m_b = -2 * b / (3 * b - 6); ok3 = True; uni = {}
    for W in (6, 8, 10, 12):
        dim, prof = conserved(H, W); uni[W] = prof; s = W - 1
        facs = {f_ for v in prof for f_, _ in sp.factor_list(sp.denom(v))[1] if f_.free_symbols}
        ok3 &= dim == 1 and all(c not in f_.free_symbols for f_ in facs)
        zeros = sorted({sp.nsimplify(bp) for f_ in facs for bp in sp.solve(f_, b)})
        expected = sorted({sp.Rational(6 * j, 3 * j - s) for j in range(1, (s - 1) // 2 + 1) if 3 * j != s})
        ok3 &= zeros == expected
    require(ok3, "P3: universal (b, c) hierarchy, GENERICALLY over Q(b, c) -- one polynomial density per weight 6-12; denominators in b alone with zeros exactly at b = 6j/(3j-s), j = 1..(s-1)/2 (the Beta lattice of m(b) = -2b/(3b-6))")
    # P3b: the exceptional locus (rank drops of the constraint matrix) -- FAIL-CLOSED enumeration (Codex round 26):
    # the lex Groebner ideal of the maximal minors (decoupling factor b removed) must be zero-dimensional; its points are
    # enumerated exhaustively by factoring the univariate element over Q and back-substituting; each point must annihilate
    # every minor and have nullity exactly 2; the decoupling component b = 0 is checked by the general argument at 6-12.
    import itertools
    EXPECT_ISO = {6: {(-3, 1)}, 8: {(sp.Rational(-3, 2), sp.Rational(1, 8)), (sp.Rational(-3, 2), 1), (-12, 8)},
                  10: {(-1, 0), (-1, sp.Rational(1, 4)), (-1, 1), (-4, 1), (-4, 4)}}
    EXPECT_JUMPS = {3: {6: {sp.Rational(5, 3)}, 8: {sp.Rational(11, 5)}, 10: {sp.Rational(7, 5), sp.Rational(17, 7)}},
                    2: {6: {sp.Rational(5, 3)}, 8: {2}, 10: set()}, 1: {6: set(), 8: set(), 10: set()}}
    QUANTUM = {3: {6: {sp.Rational(5, 3)}, 10: {sp.Rational(7, 5)}}, 2: {6: {sp.Rational(5, 3)}}}     # cyl_strata_points (spins 5, 9)
    def cmatrix(W, Hq):
        kk = W // 2; cs0 = [sp.Symbol(f'a{i}') for i in range(0, kk + 1)]
        F0 = sum(cs0[i] * x**(W - 2 * i) * r**(2 * i) for i in range(0, kk + 1))
        Mc, _ = sp.linear_eq_to_matrix(sp.Poly(pde(F0, Hq), x, r).coeffs(), cs0)
        return Mc, kk + 1
    def enumerate_points(minors):
        G = sp.groebner(minors, b, c, order='lex')
        if not G.is_zero_dimensional:
            return None, "not zero-dimensional"
        uni = [g for g in G.exprs if g.free_symbols <= {c}]
        if not uni:
            return None, "no univariate element"
        pc = sp.gcd_list(uni) if len(uni) > 1 else uni[0]; pts = set()
        for fac, _ in sp.factor_list(pc)[1]:
            if sp.Poly(fac, c).degree() != 1:
                return None, f"non-rational factor {fac}"
            cv = sp.solve(fac, c)[0]
            rest = [g_ for g_ in (sp.expand(g.subs(c, cv)) for g in G.exprs if not (g.free_symbols <= {c})) if g_ != 0]
            if not rest:
                return None, f"no equation in b at c = {cv}"
            pb = sp.gcd_list(rest) if len(rest) > 1 else rest[0]
            for fb, _ in sp.factor_list(pb)[1]:
                if sp.Poly(fb, b).degree() != 1:
                    return None, f"non-rational factor {fb}"
                pts.add((sp.solve(fb, b)[0], cv))
        return pts, "ok"
    ok3b = True; why = []
    for W in (6, 8, 10, 12):
        Mc, ncol = cmatrix(W, H)
        # decoupling component: H_xr vanishes identically at b = 0 and H_xx - H_rr does not, so the PDE forces F_xr = 0
        # and the even homogeneous densities are A x^W + B r^W: nullity exactly 2, for symbolic c
        H0 = H.subs(b, 0)
        ok3b &= sp.diff(H0, x, r) == 0 and sp.simplify(sp.diff(H0, x, 2) - sp.diff(H0, r, 2)) != 0 and (ncol - Mc.subs(b, 0).rank()) == 2
        if W == 12:
            continue          # the isolated points at weight 12 are a pinned lab result (rank_drop_locus_w12.txt)
        combos = list(itertools.combinations(range(Mc.cols), Mc.rows))
        dets = [sp.expand(Mc.extract(list(range(Mc.rows)), list(cols)).det()) for cols in combos]
        minors = [sp.expand(sp.cancel(d_ / b)) for d_ in dets]
        ok3b &= all(sp.expand(m_ * b - d_) == 0 for m_, d_ in zip(minors, dets))            # b divides every maximal minor
        pts, msg = enumerate_points(minors)
        if pts is None:
            ok3b = False; why.append(f"w{W}: {msg}"); continue
        ok3b &= pts == {(sp.nsimplify(u_), sp.nsimplify(v_)) for u_, v_ in EXPECT_ISO[W]}
        for (bv, cv) in pts:
            ok3b &= all(m_.subs({b: bv, c: cv}) == 0 for m_ in minors) and (ncol - Mc.subs({b: bv, c: cv}).rank()) == 2
        for sol, (bb, cc) in fam.items():
            hits = {tv for (bv, cv) in pts for tv in sp.solve(sp.Eq(bb, bv), t) if sp.simplify(cc.subs(t, tv) - cv) == 0}
            ok3b &= hits == EXPECT_JUMPS[sol][W]
            ok3b &= QUANTUM.get(sol, {}).get(W, set()) <= hits
    # the families at t = 3: Solution 3 on the decoupling line (b = 0, c = 1); Solution 2 at the decoupled point at infinity:
    # 1/c -> 0 and b/c -> 0 as t -> 3, so H_4/c -> r^4, whose conserved space is again {x^W, r^W} (H_xr = 0)
    ok3b &= sp.simplify(fam[3][0].subs(t, 3)) == 0 and sp.simplify(fam[3][1].subs(t, 3)) == 1
    ok3b &= sp.limit(1 / fam[2][1], t, 3) == 0 and sp.limit(fam[2][0] / fam[2][1], t, 3) == 0
    Hinf = r**4
    ok3b &= sp.diff(Hinf, x, r) == 0 and all((cmatrix(W, Hinf)[1] - cmatrix(W, Hinf)[0].rank()) == 2 for W in (6, 8, 10))
    require(ok3b, "P3b: the exceptional locus, fail-closed -- decoupling line b = 0 with nullity exactly 2 at weights 6-12 by the general argument (H_xr = 0 there); isolated NON-DECOUPLING points enumerated exhaustively (zero-dimensional lex ideal, factored over Q, back-substituted), each annihilating every maximal minor with nullity exactly 2: {(-3,1)} (w6), {(-3/2,1/8),(-3/2,1),(-12,8)} (w8), {(-1,0),(-1,1/4),(-1,1),(-4,1),(-4,4)} (w10); on the families Solution 3 jumps at t = 5/3, 11/5, 7/5, 17/7, Solution 2 at 5/3, 2, Solution 1 never, containing the certified non-decoupling quantum jumps (5/3 at w6 on sheets 2, 3; 7/5 at w10 on sheet 3); Solution 3 reaches b = 0 at t = 3 and Solution 2 the decoupled point at infinity (H_4/c -> r^4, nullity 2) at t = 3"
            + ("" if ok3b else f" -- {why}"))
    # P3c: the double-lattice closed form of the isolated points (sec. 26): at spin s = W - 1 they are EXACTLY
    #   b = 6j/(3j - s),  c = j(3j' - s)/(j'(3j - s)),   j, j' >= 1,  j + j' <= (s-1)/2  (3j != s),
    # verified against the exhaustive enumeration at weights 6-10 (P3b) and the pinned loci at weights 12-22; the family
    # formulas follow: Solution 3 (c = 1) needs j = j' -> t = (3s - 10j)/(s - 2j); Solution 2 (c = b^2/(6-b)) needs
    # 2j + 3j' = s; Solution 1 (c = b/6) needs s = 2j', impossible at odd s -- no classical jump ever.
    import re as _re
    def closed_form(s):
        return {(sp.Rational(6 * j, 3 * j - s), sp.Rational(j * (3 * jp - s), jp * (3 * j - s))) for j in range(1, s) for jp in range(1, s) if j + jp <= (s - 1) // 2 and 3 * j != s}
    def parse_locus(path):
        out = {}
        for line in open(path):
            mm = _re.match(r'weight (\d+): isolated rank-drop points beyond b = 0: \[(.*?)\];', line)
            if mm:
                out[int(mm.group(1))] = {(sp.Rational(a_), sp.Rational(b_)) for a_, b_ in _re.findall(r'\((-?\d+(?:/\d+)?), (-?\d+(?:/\d+)?)\)', mm.group(2))}
        return out
    ok3c = all(closed_form(W - 1) == {(sp.nsimplify(u_), sp.nsimplify(v_)) for u_, v_ in EXPECT_ISO[W]} for W in (6, 8, 10))
    pinned = {}
    for fn in ("rank_drop_locus_w12.txt", "rank_drop_locus_w14_18.txt", "rank_drop_locus_w20_22.txt"):
        pinned.update(parse_locus(os.path.join(ROOT, "results", "lab", "classical", fn)))
    ok3c &= set(pinned) == {12, 14, 16, 18, 20, 22} and all(closed_form(W - 1) == pinned[W] for W in pinned)
    # family formulas against the pinned family lists (parsed the same way) at weights 12-22 and the certified low weights
    def parse_family(path, sol):
        out = {}; W = None
        for line in open(path):
            mm = _re.match(r'weight (\d+):', line)
            if mm: W = int(mm.group(1))
            mm = _re.match(rf'   Solution {sol}: (?:non-decoupling )?classical jumps at t = \[(.*?)\]', line)
            if mm and W: out[W] = {sp.Rational(v_) for v_ in _re.findall(r'-?\d+(?:/\d+)?', mm.group(1))}
        return out
    fam_pinned = {sol: {} for sol in (1, 2, 3)}
    for fn in ("rank_drop_locus_w12.txt", "rank_drop_locus_w14_18.txt", "rank_drop_locus_w20_22.txt"):
        for sol in (1, 2, 3): fam_pinned[sol].update(parse_family(os.path.join(ROOT, "results", "lab", "classical", fn), sol))
    def sol3_cf(s): return {sp.Rational(3 * s - 10 * j, s - 2 * j) for j in range(1, s) if 2 * j <= (s - 1) // 2 and 3 * j != s}
    def sol2_cf(s):
        pts = set()
        for j in range(1, s):
            for jp in range(1, s):
                if 2 * j + 3 * jp == s and j + jp <= (s - 1) // 2 and 3 * j != s:
                    bv = sp.Rational(6 * j, 3 * j - s); pts.add(sp.nsimplify(sp.solve(sp.Eq(fam[2][0], bv), t)[0]))
        # the projective points (3j = s, at infinity of the affine chart, (1/b, c/b) = (0, (3j'-s)/(6j'))): Solution 2's
        # curve reaches (0, -1) at its normalization pole t = 11/5, i.e. s = 9j' with j' odd (Codex round 27)
        if s % 9 == 0 and (s // 9) % 2 == 1:
            pts.add(sp.Rational(11, 5))
        return pts
    for W in (12, 14, 16, 18, 20, 22):
        ok3c &= sol3_cf(W - 1) == fam_pinned[3][W] and sol2_cf(W - 1) == fam_pinned[2][W] and fam_pinned[1][W] == set()
    for W in (6, 8, 10):
        ok3c &= sol3_cf(W - 1) == EXPECT_JUMPS[3][W] and sol2_cf(W - 1) == (EXPECT_JUMPS[2][W] | ({sp.Rational(11, 5)} if W == 10 else set()))
    # Solution 1: c = b/6 with the closed form forces s = 2 j'
    jj, jjp, ss = sp.symbols('j jp s', positive=True)
    ok3c &= sp.simplify(sp.solve(sp.Eq(sp.Rational(1, 6) * 6 * jj / (3 * jj - ss), jj * (3 * jjp - ss) / (jjp * (3 * jj - ss))), ss)[0] - 2 * jjp) == 0
    require(ok3c, "P3c: the double-lattice closed form b = 6j/(3j-s), c = j(3j'-s)/(j'(3j-s)), j + j' <= (s-1)/2, reproduces the exhaustively enumerated points (weights 6-10) and the pinned loci (weights 12-22) exactly; Solution 3's jumps are t = (3s-10j)/(s-2j) (j = j'), Solution 2's are the pairs 2j + 3j' = s, and Solution 1 (c = b/6) would need s = 2j' -- never at odd spin")
    # P3e: the projective completion (Codex round 27).  In the chart (B, C') = (1/b, c/b) the line at infinity B = 0 carries
    # H = x^2 r^2 + C' r^4; its rank drops are the 3j = s lattice points C' = (3j'-s)/(6j'), j + j' <= (s-1)/2.  Solution 2's
    # curve reaches (0, -1) at its normalization pole t = 11/5 ((N, s) = (-21/4, -55/4)): a regular fiber whose top density is
    # -r^2 (x^2 - r^2), with nullity 2 at weights 10 and 28 (s = 9, 27) and 1 at the other weights checked.
    ok3e = True
    for W in (6, 8, 10, 12, 14):
        s_ = W - 1; Mc, ncol = cmatrix(W, x**2 * r**2 + c * r**4)
        dets = [sp.expand(Mc.extract(list(range(Mc.rows)), list(cols)).det()) for cols in itertools.combinations(range(Mc.cols), Mc.rows)]
        g_ = sp.factor(sp.gcd_list(dets)); zeros = {sp.nsimplify(z_) for z_ in sp.solve(g_, c)} if g_.free_symbols else set()
        expect = {sp.Rational(3 * jp - s_, 6 * jp) for j in range(1, s_) for jp in range(1, s_) if 3 * j == s_ and j + jp <= (s_ - 1) // 2}
        ok3e &= zeros == expect
    Ht = -r**2 * (x**2 - r**2)
    from mixed_engine import cyl_P4 as _P4
    Nv, sv = Fr(-21, 4), Fr(-55, 4); P4t = _P4(2, Nv, sv)
    Htop = sum(sp.Rational(cf.numerator, cf.denominator) * x**len(m[0]) * r**(2 * len(m[1])) for m, cf in P4t.items() if all(k_ == 1 for k_ in m[0]) and all(pr == (1, 1) for pr in m[1]))
    ok3e &= sp.expand(Htop - Ht) == 0 or sp.expand(Htop + Ht) == 0
    for W, expect_null in ((6, 1), (8, 1), (10, 2), (12, 1), (28, 2)):
        Mc, ncol = cmatrix(W, Ht); ok3e &= (ncol - Mc.rank()) == expect_null
    require(ok3e, "P3e: projective completion -- on the line at infinity the rank drops are exactly the 3j = s lattice points c/b = (3j'-s)/(6j') (weights 6-14: only c/b = -1 at weight 10); Solution 2 reaches (0, -1) at its normalization pole t = 11/5 (top density -r^2(x^2 - r^2) from the committed I_3), a classical jump with nullity 2 at weights 10 and 28 and 1 at 6, 8, 12 -- included in the family formula (spins 9(2h+1)); the coupling rule leaves it without a quantum jump (beta^2 = 3/8), as certified at spin 9")
    # P3f: the double lattice is a THEOREM (derivation from the Codex round-27 review, verified here): the conservation
    # matrix is bidiagonal, coefficient of x^{W-2i-1} r^{2i+1} = A_i a_i + B_i a_{i+1} with
    # A_i = 4(W-2i)[b(s-3i) + 6ci], B_i = -4(i+1)[2b(2i+1) + (6-b)(W-2i-2)]; all maximal minors vanish iff some A_i and
    # some B_m vanish with i <= m, which is the double lattice with j = k-m-1, j' = i, and the nullity there is exactly 2.
    ok3f = True
    for W in (6, 8, 10, 12, 14):
        k_ = W // 2; s_ = W - 1; a_ = sp.symbols(f'a0:{k_ + 1}'); F0 = sum(a_[i] * x**(W - 2 * i) * r**(2 * i) for i in range(k_ + 1))
        Pp = sp.Poly(pde(F0, H), x, r)
        for i in range(k_):
            A_ = 4 * (W - 2 * i) * (b * (s_ - 3 * i) + 6 * c * i); B_ = -4 * (i + 1) * (2 * b * (2 * i + 1) + (6 - b) * (W - 2 * i - 2))
            ok3f &= sp.expand(Pp.coeff_monomial(x**(W - 2 * i - 1) * r**(2 * i + 1)) - (A_ * a_[i] + B_ * a_[i + 1])) == 0
    require(ok3f, "P3f: the conservation matrix is bidiagonal with A_i = 4(W-2i)[b(s-3i)+6ci], B_i = -4(i+1)[2b(2i+1)+(6-b)(W-2i-2)] (weights 6-14, every coefficient) -- whence the double lattice of P3c with nullity exactly 2 is a theorem in the affine chart (derivation: Codex round 27, verified here)")
    # P3d: the lifting rule (sec. 29): a classical jumping point of the closed-form lattice is a certified quantum jumping
    # point iff beta^2 = (t-1)/(t+1) = 1/(2k); checked against the certified quantum jump table of the strata certificates
    # (cyl_strata_points, _s11, cyl_spin13_sol3_sporadic, cyl_spin15_stratification, cyl_spin17_double_prediction) and the
    # scan nulls through weight 20 (cyl_spin19_scans) on sheets 2 and 3
    QTABLE = {(3, 6): {sp.Rational(5, 3)}, (3, 10): {sp.Rational(7, 5)}, (3, 14): {sp.Rational(9, 7)}, (3, 16): {sp.Rational(5, 3)}, (3, 18): {sp.Rational(11, 9)},
              (2, 6): {sp.Rational(5, 3)}, (2, 12): {sp.Rational(9, 7)}, (2, 16): {sp.Rational(5, 3)}, (2, 18): {sp.Rational(13, 11)}}
    def rule(tv):
        b2 = (tv - 1) / (tv + 1); return b2.p == 1 and b2.q % 2 == 0
    ok3d = True
    for W in range(6, 22, 2):
        s_ = W - 1
        for sheet, pts in ((3, sol3_cf(s_)), (2, sol2_cf(s_))):
            ok3d &= {tv for tv in pts if rule(tv)} == QTABLE.get((sheet, W), set())
    ok3d &= rule(sp.Integer(3)) and not rule(sp.Integer(2)) and not rule(sp.Rational(3, 2))
    require(ok3d, "P3d: the lifting rule -- a closed-form classical jumping point is a certified quantum jumping point iff beta^2 = (t-1)/(t+1) = 1/(2k): reproduces the nine certified quantum jumps and every null at weights 6-20 on sheets 2 and 3 (sixteen cells); the decoupling t = 3 is k = 1")
    # P4
    law = {k_: sp.simplify(m_b.subs(b, v)) for k_, v in u4.items()}
    ok4 = sp.simplify(law[1] - M_PINNED[1]) == 0 and sp.simplify(law[3] - M_PINNED[3]) == 0 and sp.simplify(law[2] - (-2 - M_PINNED[2])) == 0
    P3 = sp.Poly(sp.expand(IPC[3]), X, Y); bP = sp.simplify(P3.coeff_monomial(X * Y) / P3.coeff_monomial(X**2)); cP = sp.simplify(P3.coeff_monomial(Y**2) / P3.coeff_monomial(X**2))
    ok4 &= sp.simplify(m_b.subs(b, bP) - (-2 - n)) == 0
    require(ok4, f"P4: exponent-mixing law m(u_4) = -2u_4/(3u_4 - 6): equals the pinned exponents of Solutions 1, 3 and the lattice dual of Solution 2's ({[sp.factor(v) for v in law.values()]}); on the paperclip family it gives -2 - n for all n")
    # P5
    ok5 = all(all(sp.simplify(v.subs(b, 0)) == 0 for v in uni[W][:-1]) and sp.simplify(uni[W][-1].subs(b, 0) - c**(W // 2 - 1)) == 0 for W in (6, 8, 10, 12))
    require(ok5, "P5: at b = 0 the universal densities are x^W + c^{W/2-1} r^W (weights 6-12): the t = 3 point of Solution 3 is the decoupled system x^W + r^W")
    # P6
    cPb = b * (b - 4) / (4 * (b - 3)); ok6 = sp.simplify(cPb.subs(b, bP) - cP) == 0
    curves = {3: sp.Integer(1), 1: b / 6, 2: b**2 / (6 - b)}
    expect_meet = {3: {5}, 1: {1}, 2: {-5, 1}}
    for sol, (bb, cc) in fam.items():
        tb = sp.solve(sp.Eq(bb, b), t)[0]; ok6 &= sp.simplify(cc.subs(t, tb) - curves[sol]) == 0
        meets = {sp.nsimplify(v) for v in sp.solve(sp.Eq(cc, cPb.subs(b, bb)), t)}
        ok6 &= meets == expect_meet[sol]
        binf, cinf = sp.limit(bb, t, sp.oo), sp.limit(cc, t, sp.oo); ok6 &= sp.simplify(cinf - cPb.subs(b, binf)) == 0
    require(ok6, "P6: (b, c) plane -- paperclip curve c = b(b-4)/(4(b-3)) reproduces eq. (59); family curves c = 1, b/6, b^2/(6-b); intersections at t -> infinity (all three), t = 5 (Sol 3), t = -5 (Sol 2), and at the closure point t = 1 (pole of the curve, not a fiber) for Solutions 1 and 2")
    # P7
    bad1 = not bracket_vanishes(dens[3][4], dens[3][6], coef=(1, 2))
    bad4 = sp.simplify((-2 * b / (3 * b - 5)).subs(b, u4[3]) - M_PINNED[3]) != 0
    require(bad1 and bad4, "P7: negative controls -- the bracket diag(d, 2d) breaks P1 for Solution 3; the law -2b/(3b-5) breaks P4 (guards live)")
    print(f"\n{time.time()-t0:.0f}s")
    if fails:
        print(f"CLASSICAL HIERARCHY: {len(fails)} FAILURE(S)"); sys.exit(1)
    print("CLASSICAL HIERARCHY CERTIFIED (the top-degree eigenvalues form a two-component hydrodynamic hierarchy determined by its weight-4 density "
          "x^4 + b x^2 r^2 + c r^4; the pole lattices are the Beta lattice of m(b) = -2b/(3b-6), a theorem of the universal hierarchy; "
          "the three families are the curves c = 1, b/6, b^2/(6-b) in the (b, c) plane against the paperclip curve c = b(b-4)/(4(b-3)))")
