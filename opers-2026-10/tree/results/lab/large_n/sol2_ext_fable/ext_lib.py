"""REVIEW (Fable secondary, 2026-10-02): the archived Sol 2 family S2(k) (sol2f.py's run/score, unchanged) at fibres OUTSIDE the
range tested in the record (k in [1, 4]): the exact point t = 1/3 (k = -1/2, n = 2) and others with k < 1 or k < 0.
Monic spins 1..9 against the certified charges.  Tamper: delta + 1/10 must fail.  Exit 0 iff the exact-point fibre matches at
every spin with a non-vanishing WKB charge and the tamper fails there; other fibres are reported."""
import sys, time
from fractions import Fraction as F
import sympy as sp
from common import *
def run(k, M, K=10, delta_shift=0):
    kq = sp.Rational(k.numerator, k.denominator)
    tv = (3 * kq + 2) / (kq + 2); g = classG(2, tv)
    assert g['alpha'] == F(k) and g['gamma'] == 1 and g['D'] == 2 * F(k) + 3, g
    n = g['D']; Cv = 2 * n * (1 + F(M))
    bl = blocks_of(g, M)
    if delta_shift: bl = [(a, sh, d + delta_shift) if a != 1 else (a, sh, d) for (a, sh, d) in bl]
    sv = gwkb.wkb(bl, g['s'], F(M), K); bb = n * F(M); ch = {}
    for kk in range(2, K + 1, 2):
        sk = sv[kk]
        if any(a.denominator == 1 and a <= 0 for (a, c) in sk): pw, poly, nonpole = thwkb.integrate_reg(sk, bb); poly = dict(poly)
        else:
            r = thwkb.integrate(sk, bb); poly = {k_: sp.Rational(v.numerator, v.denominator) for k_, v in r[2].items()} if r else {}
        ch[kk - 1] = sp.expand(sum(v * (sp.Rational(Cv.numerator, Cv.denominator) * PX**2 / g['aY'])**(i // 2) * (sp.Rational(Cv.numerator, Cv.denominator) * PI**2 / g['aX'])**(j // 2) for (i, j), v in poly.items()))
    return tv, g, Cv, ch

def gnorm(expr):
    """reduce every gamma(q), q rational non-integer, to gamma(frac) with 0 < frac < 1 by the functional equation; exact."""
    rep = {}
    for g_ in expr.atoms(sp.gamma):
        q = g_.args[0]
        if not q.is_Rational or q.is_Integer: continue
        fac = sp.Integer(1); x = q
        while x > 1: x -= 1; fac *= x            # gamma(x+1) = x gamma(x)
        while x < 0: fac /= x; x += 1            # gamma(x) = gamma(x+1)/x
        rep[g_] = fac * sp.gamma(x)
    return sp.cancel(sp.together(expr.xreplace(rep)))
def score(ch, cert, spins):
    rep = {}
    for s_ in spins:
        if gnorm(ch[s_]) == 0: rep[s_] = 'WKB charge identically 0 (certified: %s)' % ('ZERO' if cert[s_] == 0 else 'non-zero, %d coeffs' % len(sp.Poly(cert[s_], PX, PI).terms())); continue
        d = sp.expand(gnorm(sp.expand(monic(ch[s_], (s_ + 1) // 2) - cert[s_]))); n_ = len(sp.Poly(cert[s_], PX, PI).terms())   # v2: exact Gamma functional-equation reduction (v0: none; v1: sympy gammasimp, insufficient)
        rep[s_] = True if d == 0 else f'DIFF ({len(sp.Poly(d, PX, PI).terms())} of {n_} coeffs)'
    return rep
