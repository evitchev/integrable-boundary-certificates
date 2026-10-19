"""Cylinder eigenvalues of N=1 free-superfield densities on Fock highest-weight
states, NS and R sectors (SUSY Phase 4, Track T1).  Extends the validated
bosonic recipe of cyl_vev.py to the fermion letters of super_engine.py.
Exact arithmetic throughout (Fractions; sympy only for symbolic l and P).

Letters (super_engine conventions): (0, m, j) = d^m X^j (m >= 1),
(1, a, j) = d^a psi^j (a >= 0).  Plane OPEs of the engine:
    dX^i(z) dX^j(w) ~ delta_ij/(z-w)^2,     psi^i(z) psi^j(w) ~ delta_ij/(z-w).

RECIPE (plane density -> cylinder eigenvalue), z = e^u:

  * Read the density as the SAME polynomial in cylinder letters,
        d_z^m X -> d_u^{m-1} J,     d_z^a psi -> d_u^a psi_cyl,
    with the cylinder point-split (vertex-algebra) normal ordering.  This is
    legitimate because the point-split normal-ordered product of free fields
    is the Wick product with respect to the SINGULAR part of the two-point
    function only (1/(z-w)^2, 1/(z-w)), and the cylinder singular parts are
    identical (1/(u-v)^2, 1/(u-v)); hence a composite built from T, G by
    point-split products on the plane has the same free-field coefficients on
    the cylinder.  For T_l, G_l of super_lambda.py the plane->cylinder maps
        T~(u) = z^2 T_l(z) - c/24,   G~(u) = z^{3/2} G_l(z)
    become the same polynomials in J_rec = J_cyl - l and psi_cyl:
    z^2 d_z^2 X = d_u J_cyl - J_cyl and z^{3/2} d_z psi = d_u psi_cyl
    - psi_cyl/2 produce -l J_cyl and -l psi_cyl, which complete
    1/2 J^2 - l J -> 1/2 (J - l)^2 - l^2/2 and J psi - l psi = (J - l) psi;
    the constant l^2/2 + (point-split constants) is exactly -c/24 with
    c(l) = 3/2 - 12 l^2 (checked against the engine in the validation).
    The evaluator's momentum P is the zero mode of J_rec:
        P = p - l,  p = a_0 eigenvalue on the plane Fock state |p>.

  * Wick sum over partial pairings of the letters; a contracted pair (left
    letter at u1, right letter at u2, the left one first in the monomial's
    word order) gives d_{u1}^a d_{u2}^b K(u1 - u2) at u1 = u2, i.e.
    (-1)^b K^{(a+b)}(0), where K = full cylinder two-point function on the
    highest-weight state minus its singular part.  Fermion pairings carry
    the Koszul sign of the permutation bringing each pair adjacent.

  * KERNELS (derived from the mode expansions, not quoted):
      bosons: J(u) = sum_n a_n e^{-nu}, [a_m, a_n] = m delta_{m+n};
              <J(u1) J(u2)>_conn = sum_{n>0} n e^{-n eps} = 1/(4 sinh^2(eps/2)),
              R_B(eps) = 1/(4 sinh^2(eps/2)) - 1/eps^2       (as cyl_vev.py);
              unpaired J -> P; unpaired d^k J (k >= 1) -> 0.
      fermions: psi_cyl(u) = z^{1/2} psi(z) = sum_r psi_r e^{-ru},
              {psi_r, psi_s} = delta_{r+s}, psi_r |hw> = 0 for r > 0.
         NS (r in Z + 1/2, antiperiodic on the cylinder):
              <psi(u1) psi(u2)> = sum_{r = 1/2, 3/2, ...} e^{-r eps}
                                = 1/(2 sinh(eps/2)),
              K_NS(eps) = 1/(2 sinh(eps/2)) - 1/eps   (odd in eps);
              every fermion letter must be paired.
         R (r in Z, periodic): non-zero modes give
              sum_{r >= 1} e^{-r eps} = 1/(e^eps - 1),
              K_R(eps) = 1/(e^eps - 1) - 1/eps = (1/2) coth(eps/2) - 1/eps - 1/2,
              and the zero mode psi_0 is kept as an OPERATOR:
              an unpaired underived psi^j -> psi_0^j, an unpaired d^a psi
              (a >= 1) -> 0.
         R ZERO-MODE CONVENTION (Tanabe 2604.14899 eq. (4.28)): the R ground
              space of one Majorana fermion is two-dimensional, |Lambda, +->,
              psi_0 |+-> = (1/sqrt 2)|-+>, so psi_0^2 = 1/2 (Clifford).  The
              evaluator returns R results as elements of the zero-mode
              Clifford algebra: {(j1 < j2 < ...): value} = value *
              psi_0^{j1} psi_0^{j2} ..., psi_0^j psi_0^k = -psi_0^k psi_0^j
              (j != k), (psi_0^j)^2 = 1/2.  The total R two-point function
              is K_R + 1/eps + 1/2 = (1/2) coth(eps/2), the 1/2 being
              psi_0^2.  For an EVEN density with one fermion species the odd
              element never survives (a canonical monomial holds psi at most
              once), so even eigenvalues are scalars and do not depend on
              the ground-state label +-.
    Bose-Fermi pairings vanish.

Conventions match Tanabe (arXiv:2604.14899) by X = phi, l = -alpha_0/2,
p = Lambda, P = Lambda + alpha_0/2, P^2 = sigma_1; see
results/lab/susy4/tanabe/READING.md.

API
  kernel_table(kind, n)          K^{(k)}(0), k < n, kind in {'B', 'NS', 'R'}
  vev_monomial(mono, sector)     {zm_word: {P-exponent tuple: Fraction}}
  vev_density(dens, sector, P)   {zm_word: sympy expr}; dens coefficients may
                                 be Fractions or sympy expressions (e.g. in l)
  scalar_vev(dens, sector, P)    the even (zm_word = ()) part, asserting that
                                 no odd zero-mode part is present
Options for negative controls: kernels=(dict) overrides, keep_zero_mode.
"""
from fractions import Fraction as F
from functools import lru_cache

import sympy as sp

_eps = sp.Symbol('eps')
_NMAX = 40


def _series_derivs(expr, n):
    s = sp.series(expr, _eps, 0, n + 1).removeO()
    return tuple(F(str(sp.diff(s, _eps, k).subs(_eps, 0))) for k in range(n))


@lru_cache(maxsize=None)
def kernel_table(kind, n=_NMAX):
    if kind == 'B':
        return _series_derivs(1 / (4 * sp.sinh(_eps / 2) ** 2) - 1 / _eps ** 2, n)
    if kind == 'NS':
        return _series_derivs(1 / (2 * sp.sinh(_eps / 2)) - 1 / _eps, n)
    if kind == 'R':
        return _series_derivs(1 / (sp.exp(_eps) - 1) - 1 / _eps, n)
    raise ValueError(kind)


DEFAULT_KERNELS = {'NS': ('B', 'NS'), 'R': ('B', 'R')}


def _matchings(idx, compatible):
    """all partial matchings of the index list idx; pairs (i, k) with i
    before k in idx; compatible(i, k) filters pairs."""
    if not idx:
        yield []
        return
    first, rest = idx[0], idx[1:]
    for m in _matchings(rest, compatible):
        yield m
    for t, other in enumerate(rest):
        if compatible(first, other):
            for m in _matchings(rest[:t] + rest[t + 1:], compatible):
                yield [(first, other)] + m


def _perm_sign(seq):
    inv = 0
    for i in range(len(seq)):
        for k in range(i + 1, len(seq)):
            if seq[i] > seq[k]:
                inv += 1
    return -1 if inv % 2 else 1


def _add(acc, key, val):
    v = acc.get(key, 0) + val
    if v:
        acc[key] = v
    else:
        acc.pop(key, None)


@lru_cache(maxsize=200_000)
def _boson_part(bos, nfields, kb):
    """bos: tuple of (m, j) boson letters.  Returns {P-exponent tuple: Fraction}."""
    KB = kernel_table(kb)
    out = {}
    idx = list(range(len(bos)))
    for match in _matchings(idx, lambda i, k: bos[i][1] == bos[k][1]):
        val = F(1)
        paired = set()
        for i, k in match:
            a, b = bos[i][0] - 1, bos[k][0] - 1
            val *= (-1) ** b * KB[a + b]
            paired.update((i, k))
        if not val:
            continue
        exps = [0] * nfields
        ok = True
        for i in idx:
            if i in paired:
                continue
            if bos[i][0] != 1:
                ok = False
                break
            exps[bos[i][1]] += 1
        if ok:
            _add(out, tuple(exps), val)
    return out


@lru_cache(maxsize=200_000)
def _fermion_part(fer, kf, keep_zero_mode):
    """fer: tuple of (a, j) fermion letters in word order.  Returns
    {zero-mode word (sorted tuple of field indices): Fraction}."""
    KF = kernel_table(kf)
    zero_modes_allowed = (kf == 'R') and keep_zero_mode
    out = {}
    idx = list(range(len(fer)))
    for match in _matchings(idx, lambda i, k: fer[i][1] == fer[k][1]):
        paired = set()
        val = F(1)
        for i, k in match:
            a, b = fer[i][0], fer[k][0]
            val *= (-1) ** b * KF[a + b]
            paired.update((i, k))
        if not val:
            continue
        unp = [i for i in idx if i not in paired]
        if any(fer[i][0] != 0 for i in unp):
            continue
        if unp and not zero_modes_allowed:
            continue
        # Koszul sign: word -> (pair1)(pair2)...(unpaired in order)
        order = [x for pr in match for x in pr] + unp
        sign = _perm_sign(order)
        # zero-mode word: unpaired psi_0^j in word order -> Clifford normal form
        js = [fer[i][1] for i in unp]
        coeff = F(sign) * val
        # bubble-sort with anticommutation; equal neighbours square to 1/2
        js = list(js)
        changed = True
        while changed:
            changed = False
            for t in range(len(js) - 1):
                if js[t] > js[t + 1]:
                    js[t], js[t + 1] = js[t + 1], js[t]
                    coeff = -coeff
                    changed = True
                elif js[t] == js[t + 1]:
                    js = js[:t] + js[t + 2:]
                    coeff = coeff / 2
                    changed = True
                    break
        _add(out, tuple(js), coeff)
    return out


def vev_monomial(mono, sector, nfields=1, kernels=None, keep_zero_mode=True):
    """mono: canonical super_engine monomial (tuple of letters, word order).
    Returns {zm_word: {P-exponent tuple: Fraction}}."""
    kb, kf = (kernels or DEFAULT_KERNELS)[sector]
    # bosons commute with everything: separating them costs no sign
    bos = tuple((m, j) for (p, m, j) in mono if p == 0)
    fer = tuple((m, j) for (p, m, j) in mono if p == 1)
    B = _boson_part(bos, nfields, kb)
    Fp = _fermion_part(fer, kf, keep_zero_mode)
    out = {}
    for zw, fv in Fp.items():
        d = {}
        for ex, bv in B.items():
            _add(d, ex, fv * bv)
        if d:
            out[zw] = d
    return out


def vev_density(dens, sector, P=None, nfields=1, kernels=None,
                keep_zero_mode=True):
    """dens: {monomial: coeff}.  P: list of sympy symbols/values per field
    (default symbols P0, P1, ...).  Returns {zm_word: sympy expr}."""
    if P is None:
        P = sp.symbols(f'P0:{nfields}') if nfields > 1 else [sp.Symbol('P')]
    out = {}
    for mono, c in dens.items():
        c = sp.Rational(c.numerator, c.denominator) if isinstance(c, F) else sp.sympify(c)
        for zw, d in vev_monomial(mono, sector, nfields, kernels, keep_zero_mode).items():
            val = sp.Integer(0)
            for ex, v in d.items():
                term = sp.Rational(v.numerator, v.denominator)
                for Pj, e in zip(P, ex):
                    term *= Pj ** e
                val += term
            out[zw] = out.get(zw, 0) + c * val
    return {zw: sp.expand(v) for zw, v in out.items() if sp.expand(v) != 0}


def scalar_vev(dens, sector, P=None, nfields=1, kernels=None, keep_zero_mode=True):
    r = vev_density(dens, sector, P, nfields, kernels, keep_zero_mode)
    odd = [k for k in r if k != ()]
    if odd:
        raise ValueError(f"density has non-scalar zero-mode part {odd}")
    return r.get((), sp.Integer(0))


def vev_fraction_poly(dens, sector, nfields=1, kernels=None, keep_zero_mode=True):
    """Fast path for Fraction-coefficient densities, one field, scalar result:
    {P exponent: Fraction}."""
    out = {}
    for mono, c in dens.items():
        r = vev_monomial(mono, sector, nfields, kernels, keep_zero_mode)
        for zw, d in r.items():
            if zw != ():
                raise ValueError("non-scalar zero-mode part")
            for ex, v in d.items():
                _add(out, ex[0] if nfields == 1 else ex, F(c) * v)
    return out
