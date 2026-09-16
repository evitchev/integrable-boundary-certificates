"""THE EXACT POINT: at (N, s) = (28, -9) on Solution 1 and (28, 9) on Solution 2 -- the N = 28 embedding anchor where the two families
merge (cyl_n28_tower_embedding.py) and the Solution 1/2 meeting point of the pillow parameter plane (ansatz note secs. 38(ag)-(ah)) --
the pillow-type cylindrical operator with its exact sub-leading completion (sec. 38(ad); the Lukyanov-Zamolodchikov pillow-brane
vacuum oper at P_1^2 = 1/6, cyl_pillow_identification.py) reproduces the cylindrical vacuum eigenvalues EXACTLY at every spin
3, 5, 7, 9, 11:

    psi_zz = [kappa^2 u^(3/2) (1+u)^(-5/2) - X w/2 - (Y/4 + 5/16) w(1-w) + (1-w)/16] psi,   w = u/(1+u),   X = P^2, Y = Q^2,

i.e. (n, eps, c_X, kappa_s, c0, d2) = (-5/2, -3/2, 1, -1/2, 1/16, -9/16) for BOTH solutions (the same operator), and the operator's
continued-Beta WKB coefficients I_s (Riccati recursion in w, monomials w^i X^a Y^b, integrals (s eps/2)_i/(s n/2)_i) equal, after
normalization on X^k (k = (s+1)/2), the vacuum eigenvalues of the certified kernels -- computed HERE from scratch at the point:
the genuine kernel of ad_{I_3} at each weight W = 4..12 by mixed_engine.genuine_kernel_mixed (dimension 1 at every weight, checked),
its vacuum polynomial by the validated Wick recipe cyl_vev.vev_polynomial -- independent of the lab's rational-in-t tables.
Certified statement: the degree-(k-2) mismatch M_s = e_s/e_s[X^k] - I_s/I_s[X^k] is IDENTICALLY ZERO for s = 3..11 on both sheets.
Controls (run with --control; a separate registry entry, exit 0 only if every control FIRES): (a) the same pipeline at the generic
point Solution 1, (N, s) = (-46/3, -77/3) (t = 11/7) gives M_3 = -25/1404 != 0; (b) the operator with d2 perturbed by 1/100 at the
exact point gives M_3 != 0; (c) the deck partner on the SAME sheet, Solution 1 at (28, +9) (t = 1/3), gives M_3 != 0 -- the
exactness is not a property of N = 28 alone but of the two merging sheets.  Exact rational arithmetic throughout; fail-closed.
Usage: cyl_exact_point.py [--control]   (~15 min; the weight-12 kernels dominate)"""
import sys, os, time
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from mixed_engine import cyl_P4, gen_mixed_basis, genuine_kernel_mixed
import cyl_vev
X, Y = sp.symbols('X Y')
R = lambda v: sp.Rational(v.numerator, v.denominator)
CONTROL = '--control' in sys.argv
failures = []
def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok: failures.append(name)
# ---- the operator (sec. 38(ad)) and its continued-Beta WKB coefficients ----
def params(sol, tv):
    if sol == 1: n, eps, cX, d2 = -(tv + 7) / (tv + 3), F(-4) / (tv + 3), F(1), F(-1) / ((tv + 1) * (tv + 3))
    elif sol == 3: n, eps, cX, d2 = -(tv - 7) / (tv - 5), F(4) / (tv - 5), (tv - 3) / (tv - 5), -(tv - 3) / ((tv - 5) * (tv - 1) * (tv + 1))
    else: n, eps, cX, d2 = 2 * (tv - 7) / (tv + 5), F(-8) / (tv + 5), -2 * (tv - 3) / (tv + 5), F(-4) / ((tv + 1) * (tv + 5))
    return n, eps, -(n + 2) / (n - eps), cX, -eps / 24, d2
def add(p, q, s=F(1)):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, F(0)) + s * v
        if r[k] == 0: del r[k]
    return r
def mul(p, q):
    r = {}
    for (i1, a1, b1), v1 in p.items():
        for (i2, a2, b2), v2 in q.items():
            k = (i1 + i2, a1 + a2, b1 + b2); r[k] = r.get(k, F(0)) + v1 * v2
    return {k: v for k, v in r.items() if v != 0}
def scale(p, s): return {k: s * v for k, v in p.items()}
def dw(p):
    r = {}
    for (i, a, b), v in p.items():
        if i == 0: continue
        for k, s in (((i, a, b), i * v), ((i + 1, a, b), -i * v)): r[k] = r.get(k, F(0)) + s
    return {k: v for k, v in r.items() if v != 0}
def poch(a, i):
    r = F(1)
    for j in range(i): r *= a + j
    return r
def operator_I(prm, smax):
    n, eps, kap, cX, c0, d2 = prm
    onemw = {(0, 0, 0): F(1), (1, 0, 0): F(-1)}
    V = {(1, 1, 0): -cX / 2}; V = add(V, mul({(1, 0, 1): kap * cX / 2, (1, 0, 0): F(1, 4) + d2}, onemw)); V = add(V, scale(onemw, c0))
    ll = {(0, 0, 0): -eps / 2, (1, 0, 0): n / 2}
    def dz(k, r): return add(dw(r), mul(r, ll), F(-k))
    rho = {0: scale(ll, F(-1, 2))}
    rho[1] = scale(add(add(V, mul(rho[0], rho[0]), F(-1)), dz(0, rho[0]), F(-1)), F(1, 2))
    for k in range(1, smax):
        acc = dz(k, rho[k])
        for i in range(0, k + 1): acc = add(acc, mul(rho[i], rho[k - i]))
        rho[k + 1] = scale(acc, F(-1, 2))
    out = {}
    for s in range(3, smax + 1, 2):
        a0, ab = s * eps / 2, s * n / 2; I = {}
        for (i, a, b), v in rho[s].items():
            den = poch(ab, i)
            if den == 0: raise ZeroDivisionError(f"continuation lattice: (s n/2)_{i} = 0 at s = {s}")
            I[(a, b)] = I.get((a, b), F(0)) + v * poch(a0, i) / den
        out[s] = sum(R(v) * X**a * Y**b for (a, b), v in I.items())
    return out
def normalized(expr, k):
    p = sp.Poly(sp.expand(expr), X, Y); top = p.coeff_monomial(X**k)
    if top == 0: raise ValueError("vanishing top coefficient")
    return sp.expand(p.as_expr() / top)
def certified_e(sol, N, s, W):
    """vacuum eigenvalue of the (unique) genuine kernel element of ad_{I_3} at weight W, from scratch at the point"""
    ker = genuine_kernel_mixed(cyl_P4(sol, N, s), gen_mixed_basis(W), N)
    if len(ker) != 1: raise ValueError(f"kernel dimension {len(ker)} != 1 at weight {W}")
    e = cyl_vev.vev_polynomial(ker[0], N, W)
    syms = {str(z): z for z in e.free_symbols}
    return sp.expand(e.subs({syms['P']: sp.sqrt(X), syms['Q']: sp.sqrt(Y)}, simultaneous=True))
def mismatch(sol, tv, W, prm=None):
    N = (tv * tv - 25) / (tv * tv - 1); s = F(-24) * tv / (tv * tv - 1)
    k = W // 2; sW = W - 1
    op = operator_I(prm or params(sol, tv), sW)[sW]
    return sp.expand(normalized(certified_e(sol, N, s, W), k) - normalized(op, k)), N, s
t0 = time.time()
if CONTROL:
    fired = []
    m, N, s = mismatch(1, F(11, 7), 4)
    fired.append(("(a) generic point Solution 1 (N, s) = (-46/3, -77/3): M_3 != 0", m != 0 and m == sp.Rational(-25, 1404)))
    print(f"   control (a): M_3 = {m} at (N, s) = ({N}, {s})", flush=True)
    prm = list(params(1, F(-1, 3))); prm[5] += F(1, 100)
    m, N, s = mismatch(1, F(-1, 3), 4, prm=tuple(prm))
    fired.append(("(b) exact point with d2 perturbed by 1/100: M_3 != 0", m != 0))
    print(f"   control (b): M_3 = {m}", flush=True)
    m, N, s = mismatch(1, F(1, 3), 4)
    fired.append(("(c) deck partner on the same sheet, Solution 1 at (28, +9): M_3 != 0", m != 0 and N == 28 and s == 9))
    print(f"   control (c): M_3 = {m} at (N, s) = ({N}, {s})", flush=True)
    for name, ok in fired: check(f"control fired {name}", ok)
    if failures: print(f"CONTROL NOT FIRED: {len(failures)}"); sys.exit(1)
    print("EXACT POINT CONTROL FIRED: the pipeline reports nonzero mismatches at a generic point, for a perturbed operator, and on the "
          "deck partner of the same sheet; the exactness at (28, -+9) is not an artifact of the pipeline.")
    sys.exit(0)
for sol, tv in ((1, F(-1, 3)), (2, F(1, 3))):
    prm = params(sol, tv)
    check(f"Solution {sol}, t = {tv}: operator parameters (n, eps, kappa_s, c_X, c0, d2) = (-5/2, -3/2, -1/2, 1, 1/16, -9/16)",
          prm == (F(-5, 2), F(-3, 2), F(-1, 2), F(1), F(1, 16), F(-9, 16)))
    for W in (4, 6, 8, 10, 12):
        t1 = time.time()
        try:
            m, N, s = mismatch(sol, tv, W)
        except (ValueError, ZeroDivisionError) as ex:
            check(f"Solution {sol}, weight {W}: pipeline ({ex})", False); continue
        check(f"Solution {sol} at (N, s) = ({N}, {s}), spin {W - 1}: certified vacuum eigenvalue = operator's I_{W - 1} (mismatch {m}; {time.time() - t1:.0f}s)",
              N == 28 and abs(s) == 9 and m == 0)
if failures:
    print(f"FAILS: {len(failures)}"); sys.exit(1)
print(f"EXACT POINT CERTIFIED: at (N, s) = (28, -9) [Solution 1] and (28, 9) [Solution 2] the sec.-38(ad) operator "
      f"psi_zz = [kappa^2 u^(3/2)(1+u)^(-5/2) - X w/2 - (Y/4 + 5/16) w(1-w) + (1-w)/16] psi reproduces the cylindrical vacuum "
      f"eigenvalues exactly at spins 3, 5, 7, 9, 11 (kernels and VEVs recomputed at the point; {time.time() - t0:.0f}s).")
