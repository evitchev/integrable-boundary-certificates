"""Library of the classical locus checks (own code; exact arithmetic only, sympy and fractions.Fraction).

The classical problem, in the repository's convention (code/cyl_classical_hierarchy.py lines 6-13 and 105-106): velocities
x = dX, r = |dY|; the weight-4 symbol is the even quartic H = al x^4 + be x^2 r^2 + ga r^4, a point [al : be : ga] of the
projective plane (affine chart al = 1: (b, c) = (be, ga); chart at infinity be = 1: (1/b, c/b) = (0, C')); an even density
F = sum_{i=0}^{k} a_i x^(W-2i) r^(2i) of weight W = 2k (spin s = W - 1) is conserved iff
    PDE(F, H) = H_xr (F_xx - F_rr) - (H_xx - H_rr) F_xr = 0.
matrix(W, H) is the coefficient matrix of PDE(F, H) in a_0..a_k, built from the PDE itself (no closed form assumed)."""
import sympy as sp
from fractions import Fraction as Fr

x, r = sp.symbols('x r')
al, be, ga = sp.symbols('al be ga')
b, c = sp.symbols('b c')
t = sp.symbols('t')


def pde(F, H):
    return sp.expand(sp.diff(H, x, r) * (sp.diff(F, x, 2) - sp.diff(F, r, 2))
                     - (sp.diff(H, x, 2) - sp.diff(H, r, 2)) * sp.diff(F, x, r))


def density(W):
    k = W // 2
    a = sp.symbols(f'a0:{k + 1}')
    return sum(a[i] * x**(W - 2 * i) * r**(2 * i) for i in range(k + 1)), a


def matrix(W, H, orient=pde):
    """rows: the coefficients of x^(W-2i-1) r^(2i+1), i = 0..k-1, as linear forms in a_0..a_k.  Raises if PDE(F, H) has any
    other monomial or a coefficient that is not linear in the a_j (so the row set is the full coefficient set)."""
    k = W // 2
    F, a = density(W)
    P = sp.Poly(orient(F, H), x, r)
    allowed = {(W - 2 * i - 1, 2 * i + 1) for i in range(k)}
    extra = [m for m in P.monoms() if m not in allowed]
    if extra:
        raise ValueError(f"W={W}: monomials outside x^(W-2i-1) r^(2i+1): {extra}")
    rows = []
    for i in range(k):
        cf = sp.expand(P.coeff_monomial(x**(W - 2 * i - 1) * r**(2 * i + 1)))
        row = [sp.expand(sp.diff(cf, a[j])) for j in range(k + 1)]
        if sp.expand(cf - sum(row[j] * a[j] for j in range(k + 1))) != 0:
            raise ValueError(f"W={W}: row {i} not linear in the a_j")
        rows.append(row)
    return sp.Matrix(rows)


def quartic(A, B, C):
    return A * x**4 + B * x**2 * r**2 + C * r**4


def fr(v):
    """exact conversion of a sympy Rational (or int) to Fraction; refuses anything else"""
    v = sp.sympify(v)
    if not v.is_Rational:
        raise ValueError(f"not rational: {v}")
    return Fr(int(v.p), int(v.q))


def rank_fr(rows):
    """exact rank of a matrix of Fractions (list of lists), Gaussian elimination"""
    M = [list(rw) for rw in rows]
    nr = len(M); nc = len(M[0]) if nr else 0
    rk = 0
    for col in range(nc):
        piv = next((i for i in range(rk, nr) if M[i][col] != 0), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        pv = M[rk][col]
        for i in range(nr):
            if i != rk and M[i][col] != 0:
                f = M[i][col] / pv
                M[i] = [M[i][j] - f * M[rk][j] for j in range(nc)]
        rk += 1
    return rk


def nullity_at(W, A, B, C):
    """nullity of the conservation matrix at the quartic [A : B : C] (rational entries)"""
    M = matrix(W, quartic(sp.Rational(A), sp.Rational(B), sp.Rational(C)))
    rows = [[fr(M[i, j]) for j in range(M.cols)] for i in range(M.rows)]
    return M.cols - rank_fr(rows)


# ----------------------------------------------------------------- the curve and the sheets

def N_of(tv):
    tv = Fr(tv); return (tv * tv - 25) / (tv * tv - 1)


def s_of(tv):
    tv = Fr(tv); return -24 * tv / (tv * tv - 1)


def x_of(tv):
    tv = Fr(tv); return (tv - 1) / (tv + 1)


def t_of_x(xv):
    xv = Fr(xv); return (1 + xv) / (1 - xv)


# the per-sheet maps claimed by the reviewers (review R1, sec. 3 table; check A4); verified against the committed I_3 in run_a4.py
MAPS = {1: (6 * (t - 1) / (t - 9), (t - 1) / (t - 9)),
        2: (12 * (t - 1) / (5 * t - 11), 8 * (t - 1)**2 / ((t - 3) * (5 * t - 11))),
        3: (6 * (t - 3) / (t + 1), sp.Integer(1))}


def homogeneous(sheet):
    """[al(t) : be(t) : ga(t)] with polynomial entries: the claimed map cleared of its common denominator"""
    bb, cc = MAPS[sheet]
    den = sp.lcm(sp.denom(sp.together(bb)), sp.denom(sp.together(cc)))
    A, B, C = (sp.factor(sp.cancel(v * den)) for v in (sp.Integer(1), bb, cc))
    return sp.expand(A), sp.expand(B), sp.expand(C)


def quartic_at(sheet, tv):
    A, B, C = homogeneous(sheet)
    return tuple(fr(v.subs(t, sp.Rational(tv.numerator, tv.denominator))) for v in (A, B, C))


# ----------------------------------------------------------------- the closed forms (checked, not assumed, in the scripts)

def lattice_affine(s):
    """eq. (3) with 3j != s: {(j, j'): (b, c)}"""
    out = {}
    for j in range(1, s):
        for jp in range(1, s):
            if j + jp <= (s - 1) // 2 and 3 * j != s:
                out[(j, jp)] = (Fr(6 * j, 3 * j - s), Fr(j * (3 * jp - s), jp * (3 * j - s)))
    return out


def lattice_infinity(s):
    """the points at infinity (3j = s): {(j, j'): C' = c/b}"""
    out = {}
    if s % 3 == 0:
        j = s // 3
        for jp in range(1, s):
            if j + jp <= (s - 1) // 2:
                out[(j, jp)] = Fr(3 * jp - s, 6 * jp)
    return out


def sheet3_points(s):
    """{j: t}: j = j', t = (3s - 10j)/(s - 2j), 1 <= j <= (s-1)/4"""
    return {j: Fr(3 * s - 10 * j, s - 2 * j) for j in range(1, (s - 1) // 4 + 1)}


def sheet2_points(s):
    """{(j, j'): t}: 2j + 3j' = s; t = (3j + 2j')/(j + 2j'), or 11/5 when j = 3j' (the point at infinity)"""
    out = {}
    for jp in range(1, s):
        if (s - 3 * jp) > 0 and (s - 3 * jp) % 2 == 0:
            j = (s - 3 * jp) // 2
            out[(j, jp)] = Fr(11, 5) if j == 3 * jp else Fr(3 * j + 2 * jp, j + 2 * jp)
    return out


def ladder(k):
    return Fr(2 * k + 3, 2 * k + 1)
