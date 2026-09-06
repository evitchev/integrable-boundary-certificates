"""First-order response-kernel certificate (2026-09-04): what any scalar operator of the
rational class forces on the first-order data at spins 1-9, and that the cylindrical
Solutions 1 and 3 violate it; the kappa-odd lemma; the density pullback.

Setting (ansatz note sec. 21).  At the paperclip base n = -1,
    Psi'' = [kappa^2/(1+u) + V0] Psi,  V0 = -X/2 u/(1+u) + (Y/2 + 1/4) u/(1+u)^2,
the first-order response of the WKB coefficient S_s (s odd) to a potential deformation
f(u) is  delta S_s = int K_s f dz  with K_s the PINNED adjoint kernel of the linearized
Riccati recursion (results/lab/relations/kernels_final.txt, K_s = (1+u)^{s/2} x rational
in u over Q[X, Y], s = 1..11; Wolfram, lab/relations_kernels_all.wls).  A leading-symbol
deformation delta q responds through K_{s+2}, the exponent direction delta m through
K_{s+2} against log(1+u)/(1+u).

Checks (exact rational arithmetic throughout; Beta integrals at half-integers are rational):
  P1  kernel moments reproduce the PINNED first-order columns of cyl_first_order
      (u/(1+u)^b, b = 1..6; the growth terms u, u^2; the exponent column m) at spins 1-9;
  P2  kappa-odd lemma, recomputed here with code/wkb.py's recursion extended by a
      kappa^{-1} source: the response of S_1, S_3, S_5 to two shapes vanishes;
  P3  the 29 PINNED relation vectors (affine_relations.txt) annihilate every response
      column of the affine rational class u^p/(1+u)^q, 1 <= p <= q+2, q <= 9, with
      coefficients 1, X, Y, and the pinned exponent column; the class has rank 26 in the
      55-dimensional data space and its left null space is exactly 29-dimensional;
  P4  positive control: the paperclip's own n-deformation (delta m = 1, delta A = X/2,
      delta B = Y/2), assembled from the pinned Python columns, satisfies all 29 relations;
  P5  Solutions 1 and 3 (data from cyl_first_order.first_order, spins 3-9): with the ten
      normalization freedoms free, the relations are INCONSISTENT (unknown block rank 7,
      augmented 8), for each solution;
  P6  density pullback: the PINNED first-order density coefficients per mixed-basis
      monomial (density_first_order.json) pushed through the certified VEV recipe at N = 1
      reproduce the pinned tables' base and first order at weights 4, 6, 8, 10, both solutions;
  P7  negative controls: a shifted kernel coefficient breaks P1; a shifted relation vector
      breaks P3 (the guards are live);
  P8  ANY pole order (Codex round 25): with p, q symbolic every entry of the response column
      of u^p/(1+u)^q is Gamma(p) Gamma(q-p-1/2)/Gamma(q-1/2) times a rational function of
      (p, q) (Pochhammer shifts of the Beta functions), and all 29 relations vanish
      IDENTICALLY in (p, q) for the three coefficient types -- the relations hold for
      every pole order, not only the tested ones (for non-integral q as an identity of the
      meromorphic continuation, wherever the moments are defined).
Fail-closed: any failure exits 1 without the CERTIFIED marker.  Runtime ~4-6 min.
"""
import sys, os, json, re, time
from fractions import Fraction as F
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
REL = os.path.join(ROOT, "results", "lab", "relations"); FO = os.path.join(ROOT, "results", "lab", "first_order")
sys.path.insert(0, HERE)
from wkb import riccati, integrate_dz, u, m as msym
from cyl_first_order import first_order, cyl, load_col, X, Y, t
from cyl_vev import vev_coefficients
from mixed_engine import gen_mixed_basis

fails = []
def require(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg, flush=True)
    if not cond:
        fails.append(msg)

def beta_int(a, r):
    """int_0^inf u^{a-1} (1+u)^{-(a+r)} du = Gamma(a)Gamma(r)/Gamma(a+r), a >= 1 integer, r half-integer: rational"""
    den = F(1)
    for i in range(a):
        den *= (r + i)
    v = F(1)
    for i in range(1, a):
        v *= i
    return v / den

def beta_int_dlog(a, r):
    """-(d/dr) of beta_int(a, r): int u^{a-1} (1+u)^{-(a+r)} log(1+u) du = beta_int(a, r) * sum_{i<a} 1/(r+i)"""
    return beta_int(a, r) * sum(F(1) / (r + i) for i in range(a))

def parse_kernels():
    K = {}
    for line in open(os.path.join(REL, "kernels_final.txt")):
        name, expr = line.strip().split(" = ", 1)
        s = int(name.split("_")[1]); K[s] = sp.sympify(expr.replace("^", "**"))
    return K

def components(Ks):
    """{(a, b): (numerator coefficient list, q)} of the rational part: K_s = (1+u)^{s/2} N(u, X, Y)/(1+u)^q"""
    num, den = sp.fraction(sp.together(Ks)); den = sp.Poly(den, u); q = 0
    while den.degree() > 0:
        quo, rem = sp.div(den, sp.Poly(1 + u, u))
        if not rem.is_zero:
            raise ValueError("kernel denominator is not a power of (1+u)")     # explicit (not assert: survives python -O)
        den = quo; q += 1
    c = den.coeffs()[0]; P = sp.Poly(sp.expand(num / c), u, X, Y); out = {}
    for (i, a, b), cf in P.terms():
        out.setdefault((a, b), {})[i] = F(int(cf.p), int(cf.q))
    return {ab: ([d.get(i, F(0)) for i in range(max(d) + 1)], q) for ab, d in out.items()}

def moment(comp, s, p, q2, logpow=0):
    """int (1+u)^{s/2} N(u)/(1+u)^q * u^p/(1+u)^q2 [log(1+u)]^logpow du/u  (N given by coefficient list)"""
    coeffs, q = comp; tot = F(0)
    for i, c in enumerate(coeffs):
        if c == 0: continue
        a = i + p; r = F(q + q2) - F(s, 2) - a          # exponent bookkeeping: u^{a-1} (1+u)^{-(a + r)}
        if a < 1:
            raise ValueError('u^0 term in a moment')
        tot += c * (beta_int(a, r) if logpow == 0 else beta_int_dlog(a, r))
    return tot

def wl(sx): return sp.sympify(sx.replace("{", "[").replace("}", "]"))

if __name__ == "__main__":
    t0 = time.time()
    K = parse_kernels(); comps = {s: components(K[s]) for s in K}
    txt = open(os.path.join(REL, "affine_relations.txt")).read()
    coords = [tuple(int(z) for z in c) for c in wl(re.search(r"coords = (.*)", txt).group(1))]
    rels = [[F(int(sp.Rational(v).p), int(sp.Rational(v).q)) for v in r] for r in (wl(mm) for mm in re.findall(r"^rel = (.*)$", txt, re.M))]
    mcol_pinned = [F(int(sp.Rational(v).p), int(sp.Rational(v).q)) for v in wl(re.search(r"mcol = (.*)", txt).group(1))]
    require(len(coords) == 55 and len(rels) == 29 and len(mcol_pinned) == 55 and set(K) == {1, 3, 5, 7, 9, 11}, "P0: pinned inputs -- 55 coordinates, 29 relations, exponent column, kernels K_1..K_11")
    # response of shape u^p/(1+u)^q2 with coefficient X^i Y^j at coordinate (s, a, b)
    def resp(p, q2, ij, c):
        s, a, b = c; i, j = ij
        if a - i < 0 or b - j < 0 or (a - i, b - j) not in comps[s]: return F(0)
        return moment(comps[s][(a - i, b - j)], s, p, q2)
    def poly_of(vals, s):
        return sp.expand(sum(sp.Rational(v.numerator, v.denominator) * X**a * Y**b for (s_, a, b), v in vals.items() if s_ == s))
    # P1: pinned columns
    cols = {nm: load_col(f"n-1_{nm}") for nm in ("b1", "b2", "b3", "b4", "b5", "b6", "a1", "a2", "m")}
    shapes = {"b1": (1, 1), "b2": (1, 2), "b3": (1, 3), "b4": (1, 4), "b5": (1, 5), "b6": (1, 6), "a1": (1, 0), "a2": (2, 0)}
    ok1 = True
    for nm, (p, q2) in shapes.items():
        for s in (1, 3, 5, 7, 9):
            vals = {c: resp(p, q2, (0, 0), c) for c in coords if c[0] == s}
            ok1 &= sp.expand(poly_of(vals, s) - cols[nm][s]) == 0
    for s in (1, 3, 5, 7, 9):       # exponent: K_{s+2} against log(1+u)/(1+u)
        vals = {c: moment(comps[s + 2][(c[1], c[2])], s + 2, 0, 1, logpow=1) if (c[1], c[2]) in comps[s + 2] else F(0) for c in coords if c[0] == s}
        ok1 &= sp.expand(poly_of(vals, s) - cols["m"][s]) == 0
    require(ok1, "P1: kernel moments reproduce the pinned columns u/(1+u)^b (b = 1..6), the growth terms u and u^2, and the exponent column at spins 1-9")
    # P2: kappa-odd lemma with wkb.py's recursion extended by a kappa^{-1} source (m = -1)
    def odd_response(e_shape, kmax):
        # y^2 + y' = kappa^2 (1+u)^m + V0 + c e/kappa: at order kappa^{-1} the source enters r_2
        c = sp.Symbol('c'); V0 = -X/2 * u/(1+u) + (Y/2 + sp.Rational(1, 4)) * u/(1+u)**2
        r = riccati(V0, kmax, msym)                          # unperturbed, symbolic m
        # perturbed recursion: redo with the source (copy of wkb.riccati with rhs at k = 1)
        def dz(k, rk): return u * (sp.diff(rk, u) - (k * msym / 2) * rk / (1 + u))
        rr = {-1: sp.Integer(1)}; rr[0] = sp.cancel(-dz(-1, rr[-1]) / 2); rr[1] = sp.cancel((V0 - rr[0]**2 - dz(0, rr[0])) / 2)
        for k in range(1, kmax):
            ssum = sum(rr[i] * rr[k - i] for i in range(0, k + 1))
            src = c * e_shape if k == 1 else 0     # kappa^{-1} source e = e_shape (1+u)^{-m/2} (the order-kappa^{-1} equation divided by (1+u)^{-m/2}): rational bookkeeping; e = e_shape (1+u)^{1/2} at m = -1
            rr[k + 1] = sp.cancel(-(dz(k, rr[k]) + ssum - src) / 2)
        out = {}
        for s in range(1, kmax + 1, 2):
            d = sp.diff(rr[s], c).subs(c, 0)
            out[s] = sp.expand(integrate_dz(sp.cancel(d.subs(msym, -1)), s, msym).subs(msym, -1)) if d != 0 else sp.Integer(0)
        return out
    ok2 = True
    for e_shape in (u/(1+u)**2, u**2/(1+u)**3):
        resp_e = odd_response(e_shape, 5)
        ok2 &= all(resp_e[s] == 0 for s in (1, 3, 5))
    require(ok2, "P2: kappa-odd lemma -- a kappa^{-1} source e = (1+u)^{1/2} u/(1+u)^2 and (1+u)^{1/2} u^2/(1+u)^3 leaves S_1, S_3, S_5 unchanged at first order (wkb.py recursion, m = -1)")
    # P3: relations annihilate the affine rational class (q <= 9), rank 26, left null space 29
    shapes_all = [(p, q2) for q2 in range(1, 10) for p in range(1, q2 + 3)]
    columns = []
    for (p, q2) in shapes_all:
        for ij in ((0, 0), (1, 0), (0, 1)):
            columns.append([resp(p, q2, ij, c) for c in coords])
    columns.append(mcol_pinned)
    ok3 = all(sum(r[i] * col[i] for i in range(55)) == 0 for r in rels for col in columns)
    M = sp.Matrix([[sp.Rational(v.numerator, v.denominator) for v in col] for col in columns]).T    # 55 x ncols
    rk = M.rank(); ln = 55 - rk
    require(ok3 and rk == 26 and ln == 29, f"P3: the 29 pinned relations annihilate all {len(columns)} response columns of the affine rational class (q <= 9) and the exponent column; rank {rk} (expected 26), left null space {ln} (expected 29)")
    # P4: positive control -- the paperclip n-deformation from the pinned Python columns
    ok4 = True
    for r in rels:
        tot = F(0)
        for i, (s, a, b) in enumerate(coords):
            v = sp.Poly(sp.expand(cols["m"][s] + X/2 * cols["b1"][s] + Y/2 * cols["b2"][s]), X, Y).coeff_monomial(X**a * Y**b)
            tot += r[i] * F(int(sp.Rational(v).p), int(sp.Rational(v).q))
        ok4 &= (tot == 0)
    require(ok4, "P4: positive control -- the paperclip's n-deformation (delta m = 1, delta A = X/2, delta B = Y/2; pinned Python columns) satisfies all 29 relations")
    # P5: the cylindrical data
    def coeff_map(expr, s): P = sp.Poly(sp.expand(expr), X, Y); return {c: P.coeff_monomial(X**c[1] * Y**c[2]) for c in coords if c[0] == s}
    for sol in (3, 1):
        data = first_order(sol, (3, 5, 7, 9)); unk = []; vec = {c: sp.Integer(0) for c in coords}
        da1, db1 = sp.symbols('dal1 dbe1'); unk += [da1, db1]
        for c, v in coeff_map(da1 * (X + Y) + db1, 1).items(): vec[c] += v
        for s in (3, 5, 7, 9):
            al0, I0, I1 = data[s]; da, db = sp.symbols(f'dal{s} dbe{s}'); unk += [da, db]
            for c, v in coeff_map(al0 * I1 + da * I0 + db, s).items(): vec[c] += v
        eqs = [sp.expand(sum(sp.Rational(r[i].numerator, r[i].denominator) * vec[c] for i, c in enumerate(coords))) for r in rels]
        A, rhs = sp.linear_eq_to_matrix(eqs, unk); rkA, rkAb = A.rank(), sp.Matrix.hstack(A, rhs).rank()
        require(rkA == 7 and rkAb == 8, f"P5: Solution {sol} -- the 29 relations with the ten normalization freedoms are INCONSISTENT (unknown block rank {rkA}, augmented {rkAb}; expected 7, 8)")
    # P6: density pullback
    dj = json.load(open(os.path.join(REL, "density_first_order.json")))
    def poly(v): return sp.expand(sum(sp.Rational(c.numerator, c.denominator) * (-X)**a * (-Y)**b for (a, b), c in v.items()))
    e_ = sp.Symbol('e')
    for sol in (3, 1):
        ok6 = True
        for W in (4, 6, 8, 10):
            BASIS = gen_mixed_basis(W); d = dj[f"sol{sol}_w{W}"]
            dens1 = {b: F(d.get(str(b), {'c1': '0'})['c1']) for b in BASIS}; dens0 = {b: F(d.get(str(b), {'c0': '0'})['c0']) for b in BASIS}
            v1 = poly(vev_coefficients(dens1, 1, W)); v0 = poly(vev_coefficients(dens0, 1, W))
            ser = sp.expand(sp.series(cyl(sol, W).subs(t, 1 / e_), e_, 0, 2).removeO())
            I0 = sp.expand(ser.subs(e_, 0)); I1 = sp.expand(sp.diff(ser, e_).subs(e_, 0))
            ok6 &= sp.expand(v0 - I0) == 0 and sp.expand(v1 - I1) == 0
        require(ok6, f"P6: Solution {sol} -- pinned first-order density coefficients per monomial reproduce the pinned tables' base and first order at weights 4, 6, 8, 10 through the VEV recipe at N = 1")
    # P7: negative controls
    comps_bad = {s: dict(v) for s, v in comps.items()}; cl, q_ = comps_bad[3][(0, 0)]; comps_bad[3][(0, 0)] = ([c + (F(1, 7) if i == 1 else 0) for i, c in enumerate(cl)], q_)
    saved = comps; comps = comps_bad
    vals = {c: resp(1, 1, (0, 0), c) for c in coords if c[0] == 3}; bad1 = sp.expand(poly_of(vals, 3) - cols["b1"][3]) != 0
    comps = saved
    rbad = list(rels[0]); rbad[5] += F(1, 3); bad3 = any(sum(rbad[i] * col[i] for i in range(55)) != 0 for col in columns[:3])
    require(bad1 and bad3, "P7: negative controls -- a kernel coefficient shifted by 1/7 breaks P1 at spin 3; a relation vector shifted by 1/3 breaks P3 (guards live)")
    # P8: the relations as identities in symbolic (p, q)
    ps, qs = sp.symbols('p q'); rho = qs - ps - sp.Rational(1, 2)
    def poch(base, mm):
        mm = int(mm)
        if mm >= 0: return sp.prod([base + i for i in range(mm)]) if mm > 0 else sp.Integer(1)
        return 1 / sp.prod([base + i for i in range(mm, 0)])
    def entry_sym(comp, s):
        coeffs, qden = comp; tot = sp.Integer(0)
        for nn, cn in enumerate(coeffs):
            if cn == 0: continue
            tot += sp.Rational(cn.numerator, cn.denominator) * poch(ps, nn) * poch(rho, qden - (s - 1) // 2 - nn) / poch(qs - sp.Rational(1, 2), qden - (s - 1) // 2)
        return tot
    ok8 = True
    for ij in ((0, 0), (1, 0), (0, 1)):
        colsym = [entry_sym(comps[s][(a - ij[0], b - ij[1])], s) if (a - ij[0] >= 0 and b - ij[1] >= 0 and (a - ij[0], b - ij[1]) in comps[s]) else sp.Integer(0) for (s, a, b) in coords]
        for rr in rels:
            ok8 &= sp.cancel(sp.together(sum(sp.Rational(rr[i].numerator, rr[i].denominator) * colsym[i] for i in range(55)))) == 0
    require(ok8, "P8: with p, q SYMBOLIC (Pochhammer form of the Beta functions) all 29 relations vanish identically for the response column of u^p/(1+u)^q with coefficients 1, X, Y -- any pole order")
    print(f"\n{time.time()-t0:.0f}s")
    if fails:
        print(f"FIRST-ORDER KERNELS: {len(fails)} FAILURE(S)"); sys.exit(1)
    print("FIRST-ORDER KERNELS CERTIFIED (pinned adjoint kernels reproduce the pinned columns; the kappa-odd lemma; the affine rational class "
          "has rank 26 in the 55-dimensional data space of spins 1-9 with 29 exact relations, satisfied by the paperclip's own deformation and "
          "violated by cylindrical Solutions 1 and 3; the first-order density coefficients reproduce the pinned tables at weights 4-10)")
