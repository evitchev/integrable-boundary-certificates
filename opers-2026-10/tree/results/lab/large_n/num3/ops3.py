"""NUM3 operators (sealed): three-term Sol 2 (third order, monic after x(-1)), seventh-order Gamma form at t = 2, Sol 1/Sol 3 for controls; generic exact
E = 0 Q-functions for two-term equations P(theta) y + x^p R(theta) y = 0 (Meijer G^(D,0)_(dR,D))."""
import mpmath as mp
from fractions import Fraction as Fr
from ode import Op
def polymul(a, b):
    out = [mp.mpf(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i + j] += x * y
    return out
def fr2mp(q): return mp.mpf(q.numerator) / q.denominator
def fib2(k):
    k = Fr(k); n = 2 * k + 3; b = k / (k + 1); c = n / b; s = k + 1
    return dict(k=k, n=n, c=c, s=s, b=b)
def sol2_three_term(k, l0, l1, E, wrong_order=False):
    F = fib2(k); s = fr2mp(F['s']); c = F['c']; cf = fr2mp(c)
    P = [-x for x in polymul([s, mp.mpf(-1)], [s * s - l1 * l1, -2 * s, mp.mpf(1)])]      # -(s - th)((s - th)^2 - l1^2): monic
    sh = s if wrong_order else s - cf / 2                                                 # x^c L(T - c/2)  [wrong: x^c L(T)]
    R1 = [sh * sh - l0 * l0, -2 * sh, mp.mpf(1)]                                         # (sh - th)^2 - l0^2
    op = Op(P, [(c, R1), (F['n'], [mp.mpf(E)])])
    thetas = [s, s - l1, s + l1]
    rootsR = [sh - l0, sh + l0]                                                           # roots of R1 in theta
    return op, thetas, F, rootsR
def seventh(l0, l1, EG):
    ex = [mp.mpf(3), mp.mpf(5) / 4 - l0, mp.mpf(5) / 4 + l0, mp.mpf(19) / 4 - l0, mp.mpf(19) / 4 + l0, 3 - l1, 3 + l1]
    P = [mp.mpf(1)]
    for e in ex: P = polymul(P, [-e, mp.mpf(1)])
    return Op(P, [(Fr(21, 2), [mp.mpf(1)]), (Fr(7), [-mp.mpf(EG)])]), ex
def exact_E0(thetas, p, rho, rootsR):
    """two-term P(theta) y + x^p R(theta) y = 0, P monic with roots thetas (degree D), R = rho prod (theta - r_m) (degree dR):
    Q_i(0) for the decaying solution normalised y ~ x^beta exp(...), and beta = p theta_G"""
    D = len(thetas); dR = len(rootsR); pf = fr2mp(Fr(p)); sG = D - dR
    a = [t / pf for t in thetas]; B = [1 + r / pf for r in rootsR]
    C = rho * (-1) ** (dR - D + 1) * pf ** (dR - D)
    assert C > 0
    thG = ((1 - mp.mpf(sG)) / 2 + sum(a) - sum(B)) / sG
    out = []
    for i in range(D):
        v = mp.mpf(1)
        for j in range(D):
            if j != i: v *= mp.gamma(a[j] - a[i])
        for b_ in B: v /= mp.gamma(b_ - a[i])
        out.append(v * C ** a[i] / ((2 * mp.pi) ** ((mp.mpf(sG) - 1) / 2) * mp.mpf(sG) ** (-mp.mpf(1) / 2) * C ** thG))
    return out, pf * thG
