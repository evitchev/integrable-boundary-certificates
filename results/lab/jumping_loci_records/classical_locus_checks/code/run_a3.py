"""Check A3: the two reviewers' closed forms of the bidiagonal recursion against the full conservation matrix.

(1) W = 4..22 (spins 3..21): the matrix of PDE(F, H) (jl_lib.matrix, built from the PDE; it refuses any monomial outside
    x^(W-2i-1) r^(2i+1)) for H = x^4 + b x^2 r^2 + c r^4 and for the projective H = al x^4 + be x^2 r^2 + ga r^4:
    bidiagonality; comparison with the certificate's P3f forms (code/cyl_classical_hierarchy.py line 293), R2's
    (alpha_k, gamma_k), d = W, and R1's (A_i, B_i), sigma = W - 1, n = W/2 (also in R1's orientation of the PDE).
(2) R2's Hessian form: [Hess F, Hess H] = PDE(F, H) [[0, 1], [-1, 0]] identically.
(3) W = 6..14: every maximal minor (column m deleted) equals prod_{i<m} M[i,i] prod_{i>=m} M[i,i+1].
(4) W = 4..30: the rank-drop criterion min{i: A_i = 0} <= max{m: B_m = 0} against direct exact ranks at every point where
    some A_i (i >= 1) and some B_m vanish together (affine, and at infinity), at points with only one of them zero, on the
    line b = 0, and at generic controls; the maximal nullity seen."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
import sympy as sp
import jl_lib as L
from jl_lib import x, r, b, c, al, be, ga

fails = []
def require(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg, flush=True)
    if not cond:
        fails.append(msg)

H_aff = L.quartic(1, b, c)
H_proj = L.quartic(al, be, ga)
def r1_orient(F, G):
    # R1 sec. 3: (F_xx - F_rr) G_xr + F_xr (G_rr - G_xx) = 0 with F the weight-4 density and G the charge density
    return sp.expand((sp.diff(F, x, 2) - sp.diff(F, r, 2)) * sp.diff(G, x, r) + sp.diff(F, x, r) * (sp.diff(G, r, 2) - sp.diff(G, x, 2)))

print("== (1) the closed forms against the matrix, W = 4..22")
MATS = {}
for W in range(4, 24, 2):
    k = W // 2; s = W - 1; n = W // 2; sig = W - 1; d = W
    M = L.matrix(W, H_aff); MP = L.matrix(W, H_proj); MG = L.matrix(W, H_aff, orient=lambda F, H: r1_orient(H, F))
    MATS[W] = M
    bid = all(M[i, j] == 0 for i in range(k) for j in range(k + 1) if j not in (i, i + 1))
    bidp = all(MP[i, j] == 0 for i in range(k) for j in range(k + 1) if j not in (i, i + 1))
    cert = all(sp.expand(M[i, i] - 4 * (W - 2 * i) * (b * (s - 3 * i) + 6 * c * i)) == 0 and
               sp.expand(M[i, i + 1] - (-4 * (i + 1) * (2 * b * (2 * i + 1) + (6 - b) * (W - 2 * i - 2)))) == 0 for i in range(k))
    r2 = all(sp.expand(M[i, i] - 4 * (d - 2 * i) * (b * (d - 3 * i - 1) + 6 * i * c)) == 0 and
              sp.expand(M[i, i + 1] - (-4 * (i + 1) * (b * (6 * i + 4 - d) + 6 * (d - 2 * i - 2)))) == 0 for i in range(k))
    gA = [(n - i) * (b * (3 * i - sig) - 6 * c * i) for i in range(k)]
    gB = [(i + 1) * (6 * (n - i - 1) + b * (3 * i + 2 - n)) for i in range(k)]
    r1 = all(sp.expand(M[i, i] + 8 * gA[i]) == 0 and sp.expand(M[i, i + 1] + 8 * gB[i]) == 0 for i in range(k))
    r1_or = all(sp.expand(MG[i, i] - 8 * gA[i]) == 0 and sp.expand(MG[i, i + 1] - 8 * gB[i]) == 0 for i in range(k)) and (MG + M).applyfunc(sp.expand).is_zero_matrix
    proj = all(sp.expand(MP[i, i] - 4 * (W - 2 * i) * (be * (s - 3 * i) + 6 * ga * i)) == 0 and
               sp.expand(MP[i, i + 1] - (-4 * (i + 1) * (be * (6 * i + 4 - W) + 6 * al * (W - 2 * i - 2)))) == 0 for i in range(k))
    require(bid and bidp and cert and r2 and r1 and r1_or and proj,
            f"W={W} (spin {s}, {k}x{k + 1}): bidiagonal {bid}/{bidp} (affine/projective); = certificate P3f {cert}; "
            f"= R2 (alpha_k, gamma_k) {r2}; = -8 x R1 (A_i, B_i) {r1}; R1's own orientation gives -M = 8 x (A_i, B_i) "
            f"{r1_or}; projective forms 4(W-2i)[be(s-3i)+6 ga i], -4(i+1)[be(6i+4-W)+6 al(W-2i-2)] {proj}")

print("== (2) the Hessian form")
p = sp.symbols('p0:9'); q = sp.symbols('q0:3')
Fg = sum(p[i] * x**(8 - 2 * i) * r**(2 * i) for i in range(5)) + p[5] * x**7 * r + p[6] * x * r**7 + p[7] * x**3 * r**3   # not even: a general test
Hg = q[0] * x**4 + q[1] * x**2 * r**2 + q[2] * r**4 + p[8] * x**3 * r
HF = sp.hessian(Fg, (x, r)); HH = sp.hessian(Hg, (x, r))
comm = sp.expand(HF * HH - HH * HF); E = L.pde(Fg, Hg)
require(sp.expand(comm[0, 1] - E) == 0 and sp.expand(comm[1, 0] + E) == 0 and comm[0, 0] == 0 and comm[1, 1] == 0,
        "[Hess F, Hess H] = PDE(F, H) [[0, 1], [-1, 0]] identically (a general non-even test pair): R2's Hessian "
        "condition is the repository's conservation condition")

print("== (3) maximal minors = products of the diagonals, W = 6..14")
for W in range(6, 16, 2):
    M = MATS[W]; k = W // 2; ok = True
    for m in range(k + 1):
        cols = [j for j in range(k + 1) if j != m]
        det = sp.expand(M.extract(list(range(k)), cols).det(method='berkowitz'))
        prod = sp.expand(sp.Mul(*[M[i, i] for i in range(m)]) * sp.Mul(*[M[i, i + 1] for i in range(m, k)]))
        ok &= sp.expand(det - prod) == 0
    require(ok, f"W={W}: all {k + 1} maximal minors equal prod_(i<m) A_i prod_(i>=m) B_i")

print("== (4) the criterion against direct exact ranks, W = 4..30")
def nul(Mw, subs):
    rows = [[L.fr(sp.expand(Mw[i, j].subs(subs))) for j in range(Mw.cols)] for i in range(Mw.rows)]
    return Mw.cols - L.rank_fr(rows)
maxnull = 0; tested = 0; agree = 0; summary = []
for W in range(4, 32, 2):
    k = W // 2; s = W - 1
    M = MATS[W] if W in MATS else L.matrix(W, H_aff); MP = L.matrix(W, H_proj)
    nA = 0; nB = 0; nboth_drop = 0; nboth_nodrop = 0; ninf = 0; ok = True
    # affine: B_m = 0 at b_m (m <= k-2 and 6m+4 != W), A_i = 0 at c = b(3i - s)/(6i) (i >= 1)
    bm = {m: Fr(6 * (W - 2 * m - 2), W - 6 * m - 4) for m in range(0, k - 1) if 6 * m + 4 != W}
    for m, bv in bm.items():
        for i in range(1, k):
            cv = bv * Fr(3 * i - s, 6 * i)
            nv = nul(M, {b: sp.Rational(bv.numerator, bv.denominator), c: sp.Rational(cv.numerator, cv.denominator)})
            pred = 2 if i <= m else 1
            ok &= nv == pred; maxnull = max(maxnull, nv); tested += 1; agree += nv == pred
            if i <= m: nboth_drop += 1
            else: nboth_nodrop += 1
        # B_m = 0 alone (c generic)
        nv = nul(M, {b: sp.Rational(bv.numerator, bv.denominator), c: sp.Rational(1, 9973)}); ok &= nv == 1; tested += 1; agree += nv == 1; nB += 1
    # A_i = 0 alone (b generic)
    for i in range(1, k):
        bv = Fr(7, 9967); cv = bv * Fr(3 * i - s, 6 * i)
        nv = nul(M, {b: sp.Rational(bv.numerator, bv.denominator), c: sp.Rational(cv.numerator, cv.denominator)}); ok &= nv == 1; tested += 1; agree += nv == 1; nA += 1
    # at infinity: al = 0, be = 1, B_m = 0 iff 6m + 4 = W
    for m in range(0, k):
        if 6 * m + 4 == W:
            for i in range(1, k):
                cp = Fr(3 * i - s, 6 * i)
                nv = nul(MP, {al: 0, be: 1, ga: sp.Rational(cp.numerator, cp.denominator)})
                pred = 2 if i <= m else 1
                ok &= nv == pred; maxnull = max(maxnull, nv); tested += 1; agree += nv == pred; ninf += pred == 2
    # the line b = 0 (four values of c) and the points x^4, r^4; generic controls
    for cv in (Fr(0), Fr(1), Fr(-7, 3), Fr(5)):
        nv = nul(M, {b: 0, c: sp.Rational(cv.numerator, cv.denominator)}); ok &= nv == 2; tested += 1; agree += nv == 2; maxnull = max(maxnull, nv)
    for pt in ((0, 0, 1), (1, 0, 0)):
        nv = nul(MP, {al: pt[0], be: pt[1], ga: pt[2]}); ok &= nv == 2; tested += 1; agree += nv == 2
    for bv, cv in ((Fr(1, 7), Fr(2, 11)), (Fr(3), Fr(5)), (Fr(-5, 2), Fr(7, 3)), (Fr(13, 17), Fr(-19, 23))):
        nv = nul(M, {b: sp.Rational(bv.numerator, bv.denominator), c: sp.Rational(cv.numerator, cv.denominator)}); ok &= nv == 1; tested += 1; agree += nv == 1
    summary.append((W, nboth_drop, nboth_nodrop, ninf))
    require(ok, f"W={W} (spin {s}): criterion = direct rank at {nboth_drop} affine drop points (i <= m, nullity 2), {nboth_nodrop} "
                f"pairs with i > m (nullity 1), {nB} B-only and {nA} A-only points (nullity 1), {ninf} drop points at infinity, b = 0 "
                f"(nullity 2), x^4 and r^4 (nullity 2), 4 controls (nullity 1)")
require(maxnull == 2, f"the largest nullity at any tested point is {maxnull} (exactly 2 on the drop locus)")
require(all(nd == len(L.lattice_affine(W - 1)) for W, nd, _, _ in summary) and all(ni == len(L.lattice_infinity(W - 1)) for W, _, _, ni in summary),
        "the affine drop points (pairs i <= m) are exactly as many as the eq.-(3) lattice points with 3j != s, and the drop points at "
        "infinity as many as the 3j = s lattice points, at every W = 4..30")
print(f"      {tested} exact rank tests, {agree} as predicted")
print()
if fails:
    print(f"A3: {len(fails)} FAILURE(S)"); sys.exit(1)
print("A3 DONE: all checks pass")
