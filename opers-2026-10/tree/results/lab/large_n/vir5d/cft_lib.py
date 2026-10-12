"""VIR3, CFT side: vacuum eigenvalues of the Solution-3 charges on the cylinder, from the
two-boson densities of VIR2, exact and symbolic in p.

Rule (cylinder, conformal vector T_X + T_phi): a density monomial prod d^{k_i} chi_i, normal
ordered by subtraction of the singular terms in the cylinder coordinate, has vacuum eigenvalue
   sum over partial pairings of same-species letters  prod C(k_i, k_j)  x  prod_unpaired z(letter),
   C(m, n) = (-1)^(n-1) (m+n-2)! r_{(m+n-2)/2}   (0 if m + n odd),
   1/(4 sinh^2(w/2)) = 1/w^2 + sum_j r_j w^{2j},   r_j = -B_{2j+2} (2j+1)/(2j+2)!,
   z(dX) = PX,  z(dphi) = pi  (= P_phi - rho),  z(higher derivative) = 0.
"""
import sys, pickle, json
from fractions import Fraction as F
import sympy as sp
sys.dont_write_bytecode = True
import vir_lib as V
import vir2_lib as W
from vir2_lib import rho

p, PX, PI, t = sp.symbols('p PX PI t')
tp = (2 + p**2) / (2 - p**2)
rp = p / 2 - 1 / p
Np = sp.cancel(1 - 12 * rp**2)
sP = sp.cancel(3 * p**2 - 12 / p**2)


def r_coef(j):
    return -sp.bernoulli(2 * j + 2) * (2 * j + 1) / sp.factorial(2 * j + 2)


def C(m, n):
    if (m + n) % 2:
        return sp.Integer(0)
    return (-1) ** (n - 1) * sp.factorial(m + n - 2) * r_coef((m + n - 2) // 2)


def _species_poly(orders, z):
    """sum over partial pairings of the letters `orders` (one species)."""
    orders = list(orders)
    if not orders:
        return sp.Integer(1)
    first, rest = orders[0], orders[1:]
    tot = (z if first == 1 else 0) * _species_poly(rest, z)
    for i, o in enumerate(rest):
        c = C(first, o)
        if c != 0:
            tot += c * _species_poly(rest[:i] + rest[i + 1:], z)
    return sp.expand(tot)


def vacuum_eigenvalue(density):
    tot = sp.Integer(0)
    for (xs, fs), c in density.items():
        tot += c * _species_poly(xs, PX) * _species_poly(fs, PI)
    return sp.expand(tot)


B4 = V.basis(4)
i4 = {m: i for i, m in enumerate(B4)}


def sol3_I3_density():
    A2 = sp.cancel((13 - 27 * Np + 2 * Np**2 + (5 - 2 * Np) * sP) / 6)
    Cx = sp.cancel(7 - Np + sP)
    vec = [sp.Integer(0)] * len(B4)
    vec[i4[((), ((1, 1), (1, 1)))]] = 1
    vec[i4[((), ((2, 2),))]] = -2
    vec[i4[((2, 2), ())]] = A2
    vec[i4[((1, 1), ((1, 1),))]] = Cx
    vec[i4[((1, 1, 1, 1), ())]] = 1
    d = W.to_two_boson(vec, 4)[0]
    return {m: sp.cancel(c.subs(rho, rp)) for m, c in d.items()}


def sol3_I5_density():
    I5 = [sp.sympify(x) for x in pickle.load(open('i5_sym.pkl', 'rb'))[3]]
    v5 = [sp.cancel(c.subs(t, tp)) for c in I5]
    d = W.to_two_boson(v5, 6)[0]
    return {m: sp.cancel(c.subs(rho, rp)) for m, c in d.items()}


def sol3_promoted_density(sigma):
    import os as _os_kit; _KITROOT = _os_kit.path.join(_os_kit.path.dirname(_os_kit.path.abspath(__file__)), '..', '..', '..', '..') + '/'  # KIT PATCH: repo root from __file__ (was an absolute path)
    QT = json.load(open(_KITROOT + 'results/qt_promotion.json'))
    data = QT['cyl3_sigma%d' % sigma]
    Ns = sp.Symbol('N')
    Bw = V.basis(sigma + 1)
    vec = [sp.Integer(0)] * len(Bw)
    for i, (A_, B_) in data['coeffs_A_B_of_N'].items():
        vec[int(i)] = sp.cancel(sp.sympify(A_, locals={'N': Ns}).subs(Ns, Np) + sP * sp.sympify(B_, locals={'N': Ns}).subs(Ns, Np))
    d = W.to_two_boson(vec, sigma + 1)[0]
    return {m: sp.cancel(c.subs(rho, rp)) for m, c in d.items()}


def normalise_top(Fpoly, W_):
    """Divide by the coefficient of PX^W (the record's normalisation: plane (dX)^W coefficient 1)."""
    c = sp.Poly(Fpoly, PX, PI).coeff_monomial(PX**W_)
    return sp.Poly(sp.cancel(Fpoly / c), PX, PI)
