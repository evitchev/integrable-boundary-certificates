"""VIR5c: formula F1 (loss-1 layer of the generic-t operator at every spin index K) and helpers.

Operator variables: x0 = PX^2, x1 = pi^2.  Record variables: X = -PX^2, Y = rho^2 - pi^2.
    Q_K = [w^K] phi^nu + (nu/24) [w^(K-1)] phi^nu (const - 2 D log psi_1 - 4 D^2 log psi_1),
    nu = (2K-1)/n, phi = (1 - x0 w)(1 - x1 w), psi_1 = n phi - 2 D phi, D = w d/dw, const = 4.
All arithmetic exact: the polynomial ring QQ[A, B] (sympy.polys.rings); series in w are lists.
(v0 of this file used sympy expressions; same formulas, slower.  S1 was run with both.)
"""
import sys
import sympy as sp
from sympy.polys.rings import ring
from sympy.polys.domains import QQ

sys.dont_write_bytecode = True
RING, A, B = ring('A,B', QQ)
ONE, ZERO = RING(1), RING(0)


def q_(x):
    x = sp.Rational(x)
    return QQ(int(x.p), int(x.q))


def fibre(t):
    t = sp.Rational(t)
    k = 2 * (t - 1) / (3 - t)
    n = k + 4
    rho2 = 2 / (t * t - 1)
    return k, n, rho2


def ser_mul(u, v, N):
    out = [ZERO] * (N + 1)
    for i, a in enumerate(u):
        if i > N or not a:
            continue
        for j, b in enumerate(v):
            if i + j > N:
                break
            if b:
                out[i + j] += a * b
    return out


def binom_series(a, x, N):
    """(1 + x w)^a to order N; a rational, x a ring element."""
    a = sp.Rational(a)
    out, c, xp = [], sp.Integer(1), ONE
    for j in range(N + 1):
        out.append(xp * q_(c))
        c = c * (a - j) / (j + 1)
        xp = xp * x
    return out


def ser_inv(u, N):
    """1/u for u[0] a nonzero constant ring element."""
    u0 = u[0].coeff(1) if hasattr(u[0], 'coeff') else u[0]
    inv0 = QQ(1) / u0
    out = [ONE * inv0]
    for m in range(1, N + 1):
        s = ZERO
        for j in range(1, min(m, len(u) - 1) + 1):
            s += u[j] * out[m - j]
        out.append(-s * inv0)
    return out


def D(u):
    return [c * QQ(j) for j, c in enumerate(u)]


def loss1_operator(K, n, sgn=-1, const=4):
    """second term of F1 (un-normalised), phi = (1 + sgn A w)(1 + sgn B w); ring element of degree K-1."""
    n = sp.Rational(n)
    nu = sp.Rational(2 * K - 1) / n
    N = K - 1
    phi = [ONE, (A + B) * QQ(sgn), A * B]
    Dphi = D(phi)
    psi = [phi[j] * q_(n) - Dphi[j] * QQ(2) for j in range(3)]
    Dpsi = D(psi)
    inv = ser_inv(psi, N)
    dlog = ser_mul(Dpsi, inv, N)                      # D log psi_1
    d2log = D(dlog)                                   # D^2 log psi_1
    bracket = [dlog[j] * QQ(-2) + d2log[j] * QQ(-4) for j in range(N + 1)]
    bracket[0] = bracket[0] + ONE * q_(const)
    phinu = ser_mul(binom_series(nu, A * QQ(sgn), N), binom_series(nu, B * QQ(sgn), N), N)
    tot = ZERO
    for j in range(N + 1):
        tot += phinu[j] * bracket[N - j]
    return tot * q_(nu / 24)


def top_operator(K, n, sgn=-1):
    """[w^K] ((1 + sgn A w)(1 + sgn B w))^nu"""
    nu = sp.Rational(2 * K - 1) / sp.Rational(n)
    tot = ZERO
    for a in range(K + 1):
        tot += A**a * B**(K - a) * q_(sp.binomial(nu, a) * sp.binomial(nu, K - a) * sgn**K)
    return tot


def as_dict(P, scale=1):
    s = q_(scale)
    return {tuple(m): sp.Rational(int((c * s).numerator), int((c * s).denominator)) for m, c in P.terms()}


def record_loss1_from_F1(K, t, const=4, with_rho=True):
    """record-normalised loss-1 layer as {(a, b): coefficient of X^a Y^b}; None if the X^K coefficient vanishes."""
    k, n, rho2 = fibre(t)
    nu = sp.Rational(2 * K - 1) / n
    lead = sp.binomial(nu, K)                # X^K coefficient of [w^K]((1 + Xw)(1 + Yw))^nu
    if lead == 0:
        return None
    l1 = loss1_operator(K, n, sgn=+1, const=const)
    # pi^2 = rho^2 - Y: phi = (1 + X w)(1 + Y w - rho^2 w); first order in rho^2 of the top
    if with_rho:
        shift = ZERO
        for a in range(K):
            shift += A**a * B**(K - 1 - a) * q_(-rho2 * nu * sp.binomial(nu, a) * sp.binomial(nu - 1, K - 1 - a))
        l1 = l1 + shift
    return as_dict(l1, 1 / lead)


def Kd_coeff(d, p, q, r, a):
    return sp.binomial(d, a) * sp.rf(p, a) * sp.rf(q, d - a) * r**(d - a) / sp.rf(p, d)
