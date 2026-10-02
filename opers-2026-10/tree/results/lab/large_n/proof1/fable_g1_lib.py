"""VIR7: formula G1 -- the loss-1 layer (first quantum order) of ANY Gamma-symbol operator
        H(T) psi = x^n (x^(nM) - E) psi,   H = H_0 (1 + b(1/T)/T^2 + ...),   H_0 = T^n phi(1/T),
    phi(z) = prod_i (1 - e_i z)^(alpha_i)   (sum of alpha_i may be less than n: the rest is a power of T),
    Gamma strings: a string of length a centred at e contributes  c_1(a) sigma^2 / (1 - e z)^2  to b,
    c_1(a) = -a (a^2 - 1)/24, sigma = n M.
With nu = (j-1)/n, D = z d/dz, c2 = n^2 (M+1):
    Q_j = [z^j] phi^nu + (nu c2/24) [z^(j-2)] phi^nu { (n-1) - nu (n L1 - L2 - L1^2) } + nu [z^(j-2)] phi^nu b(z),
    L1 = D log phi, L2 = D L1.                                                                            (G1)
(derived as F1 of VIR5c, in z = 1/T instead of w = 1/T^2, with the integration by parts of VIR5d; no evenness assumed.)
Ring: QQ[l0, l1]; series in z are lists.  Exact."""
import sys
import sympy as sp
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
sys.dont_write_bytecode = True
RING, l0, l1 = ring('l0,l1', QQ)
ONE, ZERO = RING(1), RING(0)


def q_(x):
    x = sp.Rational(x)
    return QQ(int(x.p), int(x.q))


def smul(u, v, N):
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


def sadd(u, v, N, cu=1, cv=1):
    cu, cv = q_(cu), q_(cv)
    return [(u[i] if i < len(u) else ZERO) * cu + (v[i] if i < len(v) else ZERO) * cv for i in range(N + 1)]


def bpow(e, a, N):
    """(1 - e z)^a, e a ring element, a rational"""
    a = sp.Rational(a)
    out, c, xp = [], sp.Integer(1), ONE
    for j in range(N + 1):
        out.append(xp * q_(c))
        c = c * (a - j) / (j + 1)
        xp = xp * (-e)
    return out


def geom(e, N, power=1):
    """1/(1 - e z)^power"""
    return bpow(e, -power, N)


def D(u):
    return [c * QQ(j) for j, c in enumerate(u)]


def charge(j, n, M, factors, strings, extra_b=None, parts=False):
    """un-normalised Q_j truncated at loss 1: returns (top, loss1) as ring elements (degrees j and j-2 in l0, l1).
    factors: list of (e, alpha); strings: list of (e, a); extra_b: optional series (list) added to b(z)."""
    n, M = sp.Rational(n), sp.Rational(M)
    nu = sp.Rational(j - 1) / n
    c2 = n * n * (M + 1)
    sigma2 = (n * M) ** 2
    N = j
    phinu = [ONE] + [ZERO] * N
    L1 = [ZERO] * (N + 1)
    for e, al in factors:
        phinu = smul(phinu, bpow(e, sp.Rational(al) * nu, N), N)
        g = geom(e, N)                       # 1/(1 - e z) = 1 + e z + ...
        # D log (1 - e z)^al = -al e z/(1 - e z) = -al (g - 1)
        L1 = [L1[i] - (g[i] if i else ZERO) * q_(al) for i in range(N + 1)]
    L2 = D(L1)
    L1sq = smul(L1, L1, N)
    br = [L1[i] * q_(-nu * n) + L2[i] * q_(nu) + L1sq[i] * q_(nu) for i in range(N + 1)]
    br[0] = br[0] + ONE * q_(n - 1)
    b = [ZERO] * (N + 1)
    for e, a in strings:
        a = sp.Rational(a)
        c1 = -a * (a * a - 1) / 24
        g2 = geom(e, N, 2)
        b = [b[i] + g2[i] * q_(c1 * sigma2) for i in range(N + 1)]
    if extra_b is not None:
        b = [b[i] + (extra_b[i] if i < len(extra_b) else ZERO) for i in range(N + 1)]
    tot_br = [br[i] * q_(nu * c2 / 24) + b[i] * q_(nu) for i in range(N + 1)]
    top = phinu[j]
    l1 = ZERO
    for i in range(j - 1):
        l1 += phinu[i] * tot_br[j - 2 - i]
    return top, l1


def to_record(top, l1, K, CX, CY, rho2):
    """(PX, pi) -> record variables: l0^2 = CX PX^2, l1^2 = CY pi^2, PX^2 = -X, pi^2 = rho^2 - Y.
    Returns the record-normalised loss-1 layer {(a, b): coeff of X^a Y^b} (a + b = K-1), or None if the X^K coefficient vanishes."""
    XR, X, Y = ring('X,Y', QQ)
    def conv(P):
        out = XR(0)
        for (i, l), c in P.terms():
            assert i % 2 == 0 and l % 2 == 0, 'odd monomial at an even order'
            out += (X * q_(-CX)) ** (i // 2) * ((XR(1) * q_(rho2) - Y) * q_(CY)) ** (l // 2) * c
        return out
    T = conv(top)
    lead = T.coeff(X ** K)
    if lead == 0:
        return None
    full = (T + conv(l1)) * (1 / lead)
    return {(a, b): sp.Rational(int(c.numerator), int(c.denominator)) for (a, b), c in full.terms() if a + b == K - 1}, \
           {(a, b): sp.Rational(int(c.numerator), int(c.denominator)) for (a, b), c in full.terms() if a + b == K}
