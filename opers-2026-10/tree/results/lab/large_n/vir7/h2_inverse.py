"""VIR7 post-hoc diagnostics (LABELLED post-hoc; after h1_scan.py returned exit 3 = no member of class U passes T1).
v1: the series now starts at K = 1 (beta_0 from the spin-1 charge); v0 mis-indexed beta_m by one (log kept as h2_inverse_v0.log).
(1) Control S repaired: the single-string form of item 293's ODE with the X scale from the crossed-curve rule
    l0^2 : l1^2 = PX^2/a_Y : pi^2/a_X (the first version tried a finite set of scales that did not contain it).
(2) INVERSE first-order problem: for an X-even base symbol of normalisation gamma (n = 2 gamma; powers
    beta = gamma, alpha = gamma a_Y, frozen -2 alpha) and the validated dictionary (C = c2), the UNIQUE even series
    b(z) = sum_m beta_m z^(2m) that reproduces Codex's Sol-1 loss-1 layer at every K, printed in units Xh = l0^2/c2,
    Yh = l1^2/c2.  A sum of strings would give  beta_m = (A1 (2m+1) + A3) Xh^m + (A2 (2m+1) + A4) Yh^m,  no mixed terms."""
import hashlib, json, sys
import sympy as sp
sys.dont_write_bytecode = True
import g1_lib as G
from h1_core import load_law, ward
l0, l1, ZERO, ONE = G.l0, G.l1, G.ZERO, G.ONE
law1 = load_law(1)


def record_to_l(layer_minus_shift, lead, K, CX, CY):
    """polynomial in (l0, l1) whose record image is lead * (layer - shift)"""
    P = ZERO
    for (a, b), v in layer_minus_shift.items():
        P += l0 ** (2 * a) * l1 ** (2 * b) * G.q_(lead * v / ((-CX) ** a * (-CY) ** b))
    return P


def inverse_b(t, gamma, KMAX, single=None):
    t = sp.Rational(t)
    aY = (t - 1) / (t + 3)
    rho2 = 2 / (t * t - 1)
    if single is None:
        n = 2 * sp.Rational(gamma)
        fac = [(l1, gamma), (-l1, gamma), (l0, gamma * aY), (-l0, gamma * aY), (ZERO, -2 * gamma * aY)]
        alpha, beta = gamma * aY, sp.Rational(gamma)
    else:
        a = (t - 1) / 2
        n = 2 + a
        fac = [(l1, 1), (-l1, 1), (l0, a)]
        alpha, beta = aY, sp.Integer(1)          # only the RATIO matters for the top; see below
    M = sp.Integer(1)
    c2 = n * n * (M + 1)
    CX, CY = c2 / alpha, c2 / beta
    if single is not None:
        CX = CY / aY * single                    # crossed-curve ratio times a trial factor
    betas = []
    for K in range(1, KMAX + 1):
        if K == 1:
            # spin 1 is L_0 - c/24 for every solution: PX^2 + pi^2 - 1/6, i.e. X + Y + (1/6 - rho^2) in record variables
            lay, top = {(0, 0): sp.Rational(1, 6) - rho2}, {(1, 0): sp.Integer(1), (0, 1): sp.Integer(1)}
        else:
            lay, top = ward(law1, K, t)
        N = 2 * K
        extra = [ZERO] * (N + 1)
        for m, bm in enumerate(betas):
            extra[2 * m] = bm
        topP, l1P = G.charge(2 * K, n, M, fac, [], extra_b=extra)
        r_shift = G.to_record(topP, ZERO, K, CX, CY, rho2)
        if r_shift is None:
            return None, 'top vanishes at K = %d' % K
        if any(r_shift[1].get(key, 0) != v for key, v in top.items()):
            return None, 'top mismatch at K = %d' % K
        # lead of the top in record variables
        XR, X, Y = sp.polys.rings.ring('X,Y', G.QQ)
        lead = sum(c * (-CX) ** (i // 2) for (i, l), c in topP.terms() if l == 0 and i == 2 * K)
        lead = sp.Rational(int(sp.Rational(str(lead)).p), int(sp.Rational(str(lead)).q)) if not isinstance(lead, sp.Rational) else lead
        r_full = G.to_record(topP, l1P, K, CX, CY, rho2)
        diff = {key: lay[key] - r_full[0].get(key, 0) for key in lay}
        nu = sp.Rational(2 * K - 1) / n
        need = record_to_l(diff, lead, K, CX, CY)       # must equal nu * beta_(K-1)
        betas.append(need * G.q_(1 / nu))
    return (betas, c2), 'ok'


def show(betas, c2):
    for m, bm in enumerate(betas):
        terms = []
        for (i, j), c in sorted(bm.terms(), reverse=True):
            terms.append('%s Xh^%d Yh^%d' % (sp.Rational(int(c.numerator), int(c.denominator)) / c2 ** ((i + j) // 2 + 1), i // 2, j // 2))
        print('      beta_%d / c2 = ' % m + (' ; '.join(terms) if terms else '0'))


if __name__ == '__main__':
    print('== (1) control S, repaired scale')
    for tstr in ('2', '9/4'):
        for fac_trial in (1, -1):
            res, msg = inverse_b(tstr, None, 3, single=fac_trial)
            print('   t = %s, l0^2 scale = %s x (pi-scale/a_Y): %s' % (tstr, fac_trial, msg if res is None else 'top law reproduced for K = 2, 3'))
    print('== (2) inverse problem, even bases')
    for tstr in ('2', '9/4'):
        for gamma in (1, sp.Rational(1, 2), 2):
            res, msg = inverse_b(tstr, gamma, 5)
            print('   t = %s, gamma = %s (n = %s): %s' % (tstr, gamma, 2 * sp.Rational(gamma), msg))
            if res:
                show(*res)
