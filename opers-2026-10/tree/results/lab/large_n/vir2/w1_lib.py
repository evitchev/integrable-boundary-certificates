"""VIR2 helpers on top of vir2_lib: densities at a fibre and reduced condition systems."""
import sys, pickle, os
from fractions import Fraction as F
import sympy as sp
sys.dont_write_bytecode = True
import vir_lib as V
import vir2_lib as W
from vir2_lib import a, b, rho, u, S

q, g = sp.symbols('q g')        # q = a^2,  b = sqrt2 * g
C6 = V.Commutant(6)
B4, B6 = V.basis(4), V.basis(6)
i4 = {m: i for i, m in enumerate(B4)}


def fibre(uval):
    """(t, N, s, r) at u; r = u - 1/u so that rho = r/sqrt2."""
    uval = sp.Rational(uval) if not isinstance(uval, sp.Basic) else uval
    return W.t_of_u(uval), W.N_of_u(uval), W.s_of_u(uval), uval - 1 / uval


def frac(x):
    x = sp.Rational(x)
    return F(int(x.p), int(x.q))


def I3_two_boson(sol, uval):
    t, N, s, r = fibre(sp.Rational(uval))
    p4 = V.to_row(V.P4(sol, frac(N), frac(s)), i4)
    return W.to_two_boson(p4, 4)[0]


def I5_mixed(sol, uval, loop_shift=0):
    t, N, s, r = fibre(sp.Rational(uval))
    K = C6.kernel(V.P4(sol, frac(N), frac(s)), frac(N) + loop_shift)
    return K


def I5_two_boson(sol, uval):
    K = I5_mixed(sol, uval)
    assert len(K) == 1, 'kernel dim %d' % len(K)
    return W.to_two_boson(K[0], 6)[0]


def fold_S(expr):
    """Reduce a polynomial in S modulo S^2 = 2; return (part0, part1): expr = part0 + S part1."""
    p = sp.Poly(sp.expand(expr), S)
    p0 = p1 = sp.Integer(0)
    for (k,), c in p.terms():
        if k % 2 == 0:
            p0 += c * 2 ** (k // 2)
        else:
            p1 += c * 2 ** ((k - 1) // 2)
    return sp.expand(p0), sp.expand(p1)


def reduce_conditions(conds, r, rho_sign=1, keep_a=False):
    """conds: polynomials in (a, b, rho).  Substitute rho = rho_sign * r S/2, b = S g,
    fold S, strip the overall S and a parities.  Returns polynomials in (q = a^2, g)."""
    out = []
    for c in conds:
        e = sp.expand(c.subs({rho: rho_sign * r * S / 2, b: S * g}))
        p0, p1 = fold_S(e)
        assert p0 == 0 or p1 == 0, 'mixed sqrt2 parity'
        e = p0 if p1 == 0 else p1
        if e == 0:
            continue
        if keep_a:
            out.append(e)
            continue
        pa = sp.Poly(e, a)
        degs = {k for (k,), _ in pa.terms()}
        par = {k % 2 for k in degs}
        assert len(par) == 1, 'mixed parity in a'
        if par == {1}:
            e = sp.expand(sp.cancel(e / a))
            pa = sp.Poly(e, a)
        e = sum(cf * q ** (k // 2) for (k,), cf in pa.terms())
        out.append(sp.expand(e))
    return out


def system(density, r, **kw):
    return reduce_conditions(W.conditions(density), r, **kw)


def SEALED(sol, uval):
    """Sealed (q, g): Sol 1 a^2 = 1 - u^2/2, beta = -u/sqrt2 (g = -u/2);
    Sol 3 a^2 = 1/(2u^2), beta = 1/(sqrt2 u) (g = 1/(2u))."""
    uval = sp.Rational(uval) if not isinstance(uval, sp.Basic) else uval
    if sol == 1:
        return {q: 1 - uval**2 / 2, g: -uval / 2}
    if sol == 3:
        return {q: 1 / (2 * uval**2), g: 1 / (2 * uval)}
    return None
