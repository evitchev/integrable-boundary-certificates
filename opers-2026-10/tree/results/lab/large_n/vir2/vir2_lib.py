"""VIR2 library (Fable seat): the two-boson (X, phi) realisation of the cylindrical
charges and the exponential-screening test.

Conventions (sealed): dX dX ~ 1/(z-w)^2, dphi dphi ~ 1/(z-w)^2;
T_phi = (1/2)(dphi)^2 + rho d2phi, c = 1 - 12 rho^2;  V = exp(a X + b phi).
Two-boson monomial = (xs, fs): sorted tuples of derivative orders of X and of phi.
Coefficients are sympy expressions.
"""
import sys
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
from itertools import combinations
import sympy as sp

sys.dont_write_bytecode = True
import vir_lib as V

a, b, rho, u, S = sp.symbols('a b rho u S')     # S stands for sqrt(2)


def R(x):
    return sp.Rational(x.numerator, x.denominator) if isinstance(x, F) else sp.sympify(x)


def canon2(xs, fs):
    return (tuple(sorted(xs)), tuple(sorted(fs)))


def weight2(m):
    return sum(m[0]) + sum(m[1])


def add_to(out, m, c):
    if c == 0:
        return
    out[m] = out.get(m, 0) + c


def clean(v):
    out = {}
    for m, c in v.items():
        c = sp.expand(c)
        if c != 0:
            out[m] = c
    return out


def d2(v):
    out = {}
    for (xs, fs), c in v.items():
        for i in range(len(xs)):
            add_to(out, canon2(xs[:i] + (xs[i] + 1,) + xs[i + 1:], fs), c)
        for i in range(len(fs)):
            add_to(out, canon2(xs, fs[:i] + (fs[i] + 1,) + fs[i + 1:]), c)
    return clean(out)


@lru_cache(maxsize=None)
def ope2_mono(P, Q, power):
    """Coefficient of (z-w)^power in P(z) Q(w), two plain bosons: {mono: Fraction}."""
    pxs, pfs = P
    qxs, qfs = Q
    out = {}
    fm = list(V._x_matchings(pfs, qfs))
    for xnum, xD, xfreeP, xrestQ in V._x_matchings(pxs, qxs):
        for fnum, fD, ffreeP, frestQ in fm:
            D = xD + fD
            T = power + D
            if T < 0:
                continue
            r = len(xfreeP) + len(ffreeP)
            if r == 0 and T != 0:
                continue
            for ts in V._compositions(T, r):
                den = 1
                for t_ in ts:
                    den *= factorial(t_)
                xs2 = list(xrestQ) + [o + t_ for o, t_ in zip(xfreeP, ts[:len(xfreeP)])]
                fs2 = list(frestQ) + [o + t_ for o, t_ in zip(ffreeP, ts[len(xfreeP):])]
                m = canon2(xs2, fs2)
                out[m] = out.get(m, 0) + F(xnum * fnum, den)
    return {m: c for m, c in out.items() if c}


def ope2(Pv, Qv, power):
    out = {}
    for mp, cp in Pv.items():
        for mq, cq in Qv.items():
            for m, c in ope2_mono(mp, mq, power).items():
                add_to(out, m, cp * cq * R(c))
    return clean(out)


T_PHI = {((), (1, 1)): sp.Rational(1, 2), ((), (2,)): rho}
VAC2 = {((), ()): sp.Integer(1)}


def Lphi(n, v):
    return ope2(T_PHI, v, -n - 2)


@lru_cache(maxsize=None)
def _vir_phi(part):
    v = dict(VAC2)
    for n in reversed(part):
        v = Lphi(-n, v)
    return tuple(sorted(v.items(), key=lambda kv: kv[0]))


def vir_phi(part):
    """L_{-n1}...L_{-nk}|0> of T_phi as a phi-polynomial (coefficients in rho)."""
    return dict(_vir_phi(tuple(part)))


# ------------------------------------------------- HV decomposition (constant matrices)
@lru_cache(maxsize=None)
def hv_decomposer(w):
    """Labels and a constant matrix Pi with: for v in HV_w + dV_{w-1},
    v = sum_j Pi[j].v * g_j, g_j the chosen independent generators
    (HV elements first, then total derivatives).  Returns (labels, Pi, B) where
    labels[j] = ('hv', xs, part) or ('d', mono)."""
    B = V.basis(w)
    idx = {m: i for i, m in enumerate(B)}
    gens = []
    # total derivatives FIRST, so that HV elements which are themselves total
    # derivatives (e.g. L_-4|0> = d^2 T / 2) are absorbed there and the HV
    # coefficients are the class modulo d.  (v0 of the gates had HV first and
    # returned 4 L_-2^2 - 4 L_-4 for (11)(11) - 2(22): the same class.)
    for m in V.basis(w - 1):
        gens.append((('d', m), V.to_row(V.d_mono(m), idx)))
    for a_ in range(0, w + 1):
        for xs in V.x_monos(a_, 0):
            parts = V.vir_partitions(w - a_) if w - a_ else [()]
            for p in parts:
                vec = dict(V.vir_state(p, F(7, 3))) if p else dict(V.VAC)
                gens.append((('hv', xs, p), V.to_row({V.canon(xs, m[1]): c for m, c in vec.items()}, idx)))
    # independent subset, total derivatives first
    chosen, rows = [], []
    ech = []                      # (pivot column, reduced row): incremental elimination
    for lab, row in gens:
        red = list(row)
        for pc, pr in ech:
            if red[pc]:
                f = red[pc] / pr[pc]
                red = [x - f * y for x, y in zip(red, pr)]
        pc = next((i for i, x in enumerate(red) if x), None)
        if pc is not None:
            ech.append((pc, red))
            rows.append(row)
            chosen.append(lab)
    rk = len(rows)
    n = len(B)
    # pick rk independent coordinates (columns of `rows`)
    _, piv = V.rref(rows)
    Bm = [[rows[j][c] for j in range(rk)] for c in piv]        # rk x rk, entry [coord][gen]
    aug = [Bm[i] + [F(1) if i == k else F(0) for k in range(rk)] for i in range(rk)]
    Rr, pv = V.rref(aug)
    assert pv[:rk] == list(range(rk))
    inv = [r[rk:] for r in Rr]                                  # inv[gen][coord]
    return chosen, inv, piv, B, rows


def to_two_boson(vec, w, check=True):
    """Image in (X, phi) of a mixed-basis vector (list over V.basis(w), sympy or
    Fraction entries) that lies in HV_w + dV_{w-1}.  Returns (image, hv_coeffs)."""
    labels, inv, piv, B, rows = hv_decomposer(w)
    coords = [R(vec[c]) for c in piv]
    coef = [sp.expand(sum(R(inv[j][k]) * coords[k] for k in range(len(piv)) if inv[j][k]))
            for j in range(len(labels))]
    if check:   # fail-closed: the decomposition reproduces the vector exactly
        for i in range(len(B)):
            rec = sum(R(rows[j][i]) * coef[j] for j in range(len(labels)) if rows[j][i])
            assert sp.cancel(rec - R(vec[i])) == 0, 'vector not in HV + dV'
    out = {}
    hv = {}
    for lab, c in zip(labels, coef):
        if lab[0] != 'hv' or c == 0:
            continue
        hv[(lab[1], lab[2])] = c
        vp = vir_phi(lab[2]) if lab[2] else dict(VAC2)
        for (_, fs), cc in vp.items():
            add_to(out, canon2(lab[1], fs), c * cc)
    return clean(out), hv


# ------------------------------------------------- screening residue and normal form
def screening_residue(Pv, qa=a, qb=b):
    """R with Res_z P(z) V(w) = R(w) V(w), V = exp(qa X + qb phi)."""
    out = {}
    for (xs, fs), c in Pv.items():
        letters = [('x', m) for m in xs] + [('f', m) for m in fs]
        n = len(letters)
        for k in range(1, n + 1):
            for sel in combinations(range(n), k):
                fac, D = 1, 0
                for i in sel:
                    sp_, m = letters[i]
                    fac = fac * (qa if sp_ == 'x' else qb) * (-1) ** (m - 1) * factorial(m - 1)
                    D += m
                rest = [letters[i] for i in range(n) if i not in sel]
                T = D - 1
                if not rest:
                    if T == 0:
                        add_to(out, ((), ()), c * fac)
                    continue
                for ts in V._compositions(T, len(rest)):
                    den = 1
                    for t_ in ts:
                        den *= factorial(t_)
                    xs2 = [m + t_ for (sp_, m), t_ in zip(rest, ts) if sp_ == 'x']
                    fs2 = [m + t_ for (sp_, m), t_ in zip(rest, ts) if sp_ == 'f']
                    add_to(out, canon2(xs2, fs2), c * fac * sp.Rational(1, den))
    return clean(out)


def normal_form(Rv, elim='f', qa=a, qb=b):
    """Unique representative modulo the image of d + qa dX + qb dphi with no bare
    d(phi) letter (elim='f', needs qb != 0) or no bare dX letter (elim='x')."""
    cache = {}

    def nf(m):
        if m in cache:
            return cache[m]
        xs, fs = m
        if elim == 'f':
            if 1 not in fs:
                cache[m] = {m: sp.Integer(1)}
                return cache[m]
            i = fs.index(1)
            mp = (xs, fs[:i] + fs[i + 1:])
            other = canon2(xs + (1,), mp[1])
            div, oq = qb, qa
        else:
            if 1 not in xs:
                cache[m] = {m: sp.Integer(1)}
                return cache[m]
            i = xs.index(1)
            mp = (xs[:i] + xs[i + 1:], fs)
            other = canon2(mp[0], fs + (1,))
            div, oq = qa, qb
        res = {}
        for m2, c2 in d2({mp: sp.Integer(1)}).items():
            for m3, c3 in nf(m2).items():
                add_to(res, m3, -c2 * c3 / div)
        for m3, c3 in nf(other).items():
            add_to(res, m3, -oq * c3 / div)
        cache[m] = res
        return res

    out = {}
    for m, c in Rv.items():
        for m3, c3 in nf(m).items():
            add_to(out, m3, c * c3)
    return {m: sp.cancel(sp.together(c)) for m, c in out.items() if sp.cancel(sp.together(c)) != 0}


def conditions(Pv, elim='f', qa=a, qb=b):
    """Numerators of the normal-form coefficients: the conservation conditions."""
    nf = normal_form(screening_residue(Pv, qa, qb), elim, qa, qb)
    return [sp.numer(sp.together(c)) for c in nf.values()]


# ------------------------------------------------- the curve in u
def N_of_u(uu):
    return 13 - 6 * (uu**2 + uu**-2)


def s_of_u(uu):
    return 6 * (uu**2 - uu**-2)


def t_of_u(uu):
    return (1 + uu**2) / (1 - uu**2)


def rho_of_u(uu, sqrt2=sp.sqrt(2)):
    return (uu - 1 / uu) / sqrt2


def two_boson_basis(w):
    out = []
    for k in range(0, w + 1):
        for xs in V.x_monos(k, 0) + V.x_monos(k, 1):
            for fs in V.x_monos(w - k, 0) + V.x_monos(w - k, 1):
                out.append((tuple(xs), tuple(fs)))
    return sorted(set(out))


def is_total_derivative(v, w):
    """Is the two-boson vector v (weight w) in the image of d?  Exact (sympy)."""
    v = clean(v)
    if not v:
        return True
    B = two_boson_basis(w)
    idx = {m: i for i, m in enumerate(B)}
    Bm = two_boson_basis(w - 1)
    M = sp.zeros(len(B), len(Bm))
    for j, m in enumerate(Bm):
        for m2, c in d2({m: sp.Integer(1)}).items():
            M[idx[m2], j] = c
    rhs = sp.zeros(len(B), 1)
    for m, c in v.items():
        rhs[idx[m], 0] = c
    aug = M.row_join(rhs)
    return M.rank() == aug.rank(iszerofunc=lambda x: sp.cancel(x) == 0)
