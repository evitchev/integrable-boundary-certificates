"""VIR3: large-E WKB of the DDMST pseudo-differential equations (rank 2), exact in (M, g).

Equation (3.20) of hep-th/0612298 at n = 2, K = 1 (transcribed from the rendered page 5):
    ( D_2(g^dagger) D_2(g) + sqrt(P) (d/dx) sqrt(P) ) psi = 0,    P = x^{h M} - E,  h = 3,
    D_2(g) = D(g_1 - 1) D(g_0),  D_2(g^dagger) = D(-g_0) D(1 - g_1),  D(g) = d/dx - g/x.
With E -> -E, x = E^{1/(hM)} y, eps = E^{-(M+1)/(hM)}:
    DD psi + eps^{-3} (p d/dy + p'/2) psi = 0,   p = 1 + y^{hM}.
WKB: psi = exp(int sigma dy), sigma = sum_k eps^{k-1} S_k, decaying branch S_0 = -p^{1/3}.
Elements of the coefficient ring: sums of c(M, g) y^a p^b, a integer, b in (1/3) Z.
J_k = int_0^infty S_k dy = B(alpha, beta0)/(hM) * R_k,  alpha = -(k-1)/(hM),
beta0 = (k-1)(M+1)/(hM),  R_k = sum_j c_j (beta0)_j / ((k-1)/3)_j  (rational in M, g).

Generalised: operator = product of first-order factors D(c_j) (list `factors`, rightmost
first) plus eps^{-r} times a first-order "potential" operator; see `Problem`.
"""
import sys
from fractions import Fraction as F
import sympy as sp

sys.dont_write_bytecode = True
M, g0, g1 = sp.symbols('M g0 g1')
l0, l1, sh, m1, ma, mb = sp.symbols('l0 l1 sh m1 ma mb')        # VIR4: centred exponents and frozen parameters
GENS = (M, g0, g1, l0, l1, sh, m1, ma, mb)


def P_(x):
    return sp.Poly(x, *GENS, domain='QQ')


ZERO = P_(0)
ONE = P_(1)


class El(dict):
    """sum of coeff * y^a p^b ; keys (a, b) with a int, b Fraction; values sympy Poly."""

    def copy(self):
        return El(self)


def el_add(u, v, cv=ONE):
    out = El(u)
    for k, c in v.items():
        n = out.get(k, ZERO) + c * cv
        if n.is_zero:
            out.pop(k, None)
        else:
            out[k] = n
    return out


def el_mul(u, v):
    out = El()
    for (a1, b1), c1 in u.items():
        for (a2, b2), c2 in v.items():
            k = (a1 + a2, b1 + b2)
            n = out.get(k, ZERO) + c1 * c2
            if n.is_zero:
                out.pop(k, None)
            else:
                out[k] = n
    return out


def el_scale(u, c):
    c = P_(c)
    return El({k: v * c for k, v in u.items() if not (v * c).is_zero})


def el_diff(u, h):
    """d/dy [y^a p^b] = y^{a-1} [ (a + hM b) p^b - hM b p^{b-1} ],  p = 1 + y^{hM}."""
    out = El()
    for (a, b), c in u.items():
        bq = sp.Rational(b.numerator, b.denominator)
        t1 = c * P_(a + h * M * bq)
        t2 = c * P_(-h * M * bq)
        for k, t in (((a - 1, b), t1), ((a - 1, b - 1), t2)):
            if t.is_zero:
                continue
            n = out.get(k, ZERO) + t
            if n.is_zero:
                out.pop(k, None)
            else:
                out[k] = n
    return out


def ser_add(A, B):
    out = dict(A)
    for n, e in B.items():
        out[n] = el_add(out.get(n, El()), e)
    return out


def ser_mul(A, B, nmax):
    out = {}
    for n1, e1 in A.items():
        for n2, e2 in B.items():
            if n1 + n2 > nmax:
                continue
            out[n1 + n2] = el_add(out.get(n1 + n2, El()), el_mul(e1, e2))
    return out


class Problem:
    """prod_j D(c_j) psi + sign * eps^{-r} * (potential operator) psi = 0 after rescaling.

    factors: list of c_j, RIGHTMOST factor first.  The WKB series index n counts powers
    eps^{n - order} in  (D...D psi)/psi.
    potential: 'sqrtP d sqrtP'  ->  p sigma + p'/2      (order r = len(factors) - 1)
               'P'              ->  p                   (order r = len(factors))
    h: exponent of the potential p = 1 + y^{hM}.
    """

    def __init__(self, factors, potential, sign, h, root):
        self.factors = [sp.sympify(c) for c in factors]
        self.order = len(factors)
        self.potential = potential
        self.sign = sign
        self.h = h
        self.S = {}
        self.root = root      # (power of p in S_0 as Fraction, coefficient)

    def A_series(self, nmax):
        """(D...D [A0])/psi as a series: dict n -> El, meaning eps^{n - J}, J = number of
        sigma-factors at leading order (order, or order + 1 for 'P dinv P')."""
        sig = {k: e for k, e in self.S.items()}          # sigma: index k <-> eps^{k-1}
        y_inv = lambda c: El({(-1, F(0)): P_(-c)})
        if self.potential == 'P dinv P':
            # unknown chi with psi = chi'/P: start from sigma/p  (index k <-> eps^{k-1})
            A = {k: el_mul(e, El({(0, F(-1)): ONE})) for k, e in sig.items() if k <= nmax}
        else:
            A = {0: El({(0, F(0)): ONE})}                # the function 1: eps^0
        for c in self.factors:
            new = {}
            for n, e in A.items():                       # derivative and -c/y: index n -> n + 1
                if n + 1 <= nmax:
                    new[n + 1] = el_add(new.get(n + 1, El()), el_diff(e, self.h))
                    new[n + 1] = el_add(new[n + 1], el_mul(y_inv(c), e))
            prod = ser_mul(sig, A, nmax)                 # sigma * A: index k + n
            new = ser_add(new, prod)
            A = {n: e for n, e in new.items() if e}
        return A

    def solve(self, K):
        b0, c0 = self.root
        self.S = {0: El({(0, b0): P_(c0)})}
        h = self.h
        pprime = El({(-1, F(1)): P_(h * M), (-1, F(0)): P_(-h * M)})       # p' = hM (p - 1)/y
        for n in range(1, K + 1):
            A = self.A_series(n)
            rest = A.get(n, El())
            if self.potential == 'sqrtP d sqrtP':
                # + sign (p sigma + p'/2); coefficient of S_n: order S_0^{order-1} + sign p
                if n == 1:
                    rest = el_add(rest, el_scale(pprime, sp.Rational(self.sign, 2)))
                coef = self.order * c0 ** (self.order - 1) + self.sign
                denom_b = b0 * (self.order - 1)
                assert denom_b == 1
            elif self.potential == 'P':
                coef = self.order * c0 ** (self.order - 1)
                denom_b = b0 * (self.order - 1)
            elif self.potential == 'P dinv P':
                # D_5[(1/p) chi'] - eps^{-6} p chi = 0: coefficient of S_n: (order+1) S_0^order / p
                coef = (self.order + 1) * c0 ** self.order
                denom_b = b0 * self.order - 1
            else:
                raise ValueError
            inv = El({(0, -denom_b): P_(sp.Rational(-1) / coef)})
            self.S[n] = el_mul(rest, inv)
        return self.S

    def R(self, k):
        """J_k = B(alpha, beta0)/(hM) * R_k  for S_k = sum c y^{-k} p^{b}; returns (R_k, ok)
        with ok = every term has y-power -k."""
        e = self.S[k]
        b0, _ = self.root
        beta0 = sp.Rational(k - 1) * (M + 1) / (self.h * M)
        tot = sp.Rational(k - 1, 1) * sp.Rational(b0.numerator, b0.denominator)   # alpha + beta0
        base_b = F(1 - k) * b0                                        # p-power with j = 0
        R = sp.Integer(0)
        ok = True
        for (a, b), c in e.items():
            if a != -k:
                ok = False
            j = base_b - b
            if j.denominator != 1:
                ok = False
                continue
            j = int(j)
            # B(alpha, beta0 + j)/B(alpha, beta0) = (beta0)_j / (alpha + beta0)_j  (j may be negative)
            if j >= 0:
                fac = sp.rf(beta0, j) / sp.rf(tot, j)
            else:
                fac = sp.rf(tot + j, -j) / sp.rf(beta0 + j, -j)
            R += c.as_expr() * fac
        return sp.factor(sp.cancel(R)), ok


def ddmst_320():
    """(3.20) at n = 2: D(-g0) D(1-g1) D(g1-1) D(g0) + sqrtP d sqrtP,  h = 3,  S_0^3 = -p."""
    return Problem([g0, g1 - 1, 1 - g1, -g0], 'sqrtP d sqrtP', +1, 3, (F(1, 3), -1))


def ddmst_321():
    """(3.21) at n = 2: ( D(-g0) D(1-g1) (d/dx) D(g1-1) D(g0) - P (d/dx)^{-1} P ) psi = 0, h = 3.
    Solved for chi = (d/dx)^{-1} P psi:  D_5[(1/P) chi'] = P chi;  S_0^6 = p^2, S_0 = -p^{1/3}."""
    return Problem([g0, g1 - 1, 0, 1 - g1, -g0], 'P dinv P', -1, 3, (F(1, 3), -1))


def dual_4th(h=4):
    """NOT in DDMST: D(-g0) D(1-g1) D(g1-1) D(g0) psi = P psi, P = x^{hM} - E (the A_3 form (3.18)
    with a symplectic g); S_0^4 = p."""
    return Problem([g0, g1 - 1, 1 - g1, -g0], 'P', -1, h, (F(1, 4), -1))


class Chain(Problem):
    """General chain:  [ops applied right-to-left] psi + tailsign * eps^{-J} * tail(psi) = 0, where
    ops is a list (rightmost first) of ('D', c) [d/dy - c/y] and ('p', b) [multiplication by p^b],
    J = number of D's, and tail in {'1' (psi), 'd' (psi')}.  Series index n <-> eps^{n-J}.
    The total p-power of the chain is pw; leading order: S_0^J p^pw + tailsign * (1 or S_0) = 0."""

    def __init__(self, ops, tail, tailsign, h, root):
        self.ops = [(k, sp.sympify(v) if k == 'D' else F(v)) for k, v in ops]
        self.J = sum(1 for k, _ in self.ops if k == 'D')
        self.pw = sum(v for k, v in self.ops if k == 'p')
        self.tail = tail
        self.tailsign = tailsign
        self.h = h
        self.root = root
        self.S = {}

    def A_series(self, nmax):
        sig = {k: e for k, e in self.S.items()}
        y_inv = lambda c: El({(-1, F(0)): P_(-c)})
        A = {0: El({(0, F(0)): ONE})}
        for kind, v in self.ops:
            if kind == 'p':
                A = {n: el_mul(e, El({(0, v): ONE})) for n, e in A.items()}
                continue
            new = {}
            for n, e in A.items():
                if n + 1 <= nmax:
                    new[n + 1] = el_add(new.get(n + 1, El()), el_diff(e, self.h))
                    new[n + 1] = el_add(new[n + 1], el_mul(y_inv(v), e))
            new = ser_add(new, ser_mul(sig, A, nmax))
            A = {n: e for n, e in new.items() if e}
        return A

    def solve(self, K):
        b0, c0 = self.root
        self.S = {0: El({(0, b0): P_(c0)})}
        J = self.J
        # leading-order check and the coefficient of S_n:  J S_0^{J-1} p^pw (+ tailsign for tail 'd')
        lead_b = b0 * (J - 1) + self.pw
        coef = J * c0 ** (J - 1)
        if self.tail == 'd':
            assert lead_b == 0, 'coefficient of S_n is not a monomial in p'
            coef = coef + self.tailsign
            assert b0 * J + self.pw == b0 and c0 ** J + self.tailsign * c0 == 0, 'S_0 does not solve the leading order'
        else:
            assert b0 * J + self.pw == 0 and c0 ** J + self.tailsign == 0, 'S_0 does not solve the leading order'
        inv = El({(0, -lead_b): P_(sp.Rational(-1) / coef)})
        for n in range(1, K + 1):
            A = self.A_series(n)
            rest = A.get(n, El())
            # tail 'd' contributes tailsign * S_n at index n (already separated into `coef`); tail '1' only at n = 0
            self.S[n] = el_mul(rest, inv)
        return self.S


def cand_a():
    """(a) = DDMST (3.20), n = 2, written as a chain with tail sqrtP d sqrtP (kept as Problem)."""
    return ddmst_320()


def cand_b():
    """(b) D_3^(2)-type connection with the potential on the MIDDLE (long) node:
    D(-g0) D(1-g1) p^{-1} D(g1-1) D(g0) psi + psi' = 0,  p = x^{3M} - E (normalised), S_0^3 = -p."""
    return Chain([('D', g0), ('D', g1 - 1), ('p', -1), ('D', 1 - g1), ('D', -g0)], 'd', +1, 3, (F(1, 3), -1))


def cand_d():
    """(d) C_2^(1)-type connection with the potential on the MIDDLE (short) node, i.e. two coupled
    second-order equations:  D(-g0) p^{-1} D(1-g1) D(g1-1) p^{-1} D(g0) psi - psi = 0,
    p = x^{2M} - E,  S_0^4 = p^2,  S_0 = -p^{1/2}."""
    return Chain([('D', g0), ('p', -1), ('D', g1 - 1), ('D', 1 - g1), ('p', -1), ('D', -g0)], '1', -1, 2, (F(1, 2), -1))


def cand_c():
    """(c) C_2^(1)-type connection with the potential on an END (long) node = the A_3 equation (3.18)
    with a symplectic g:  D(-g0) D(1-g1) D(g1-1) D(g0) psi - p psi = 0, p = x^{4M} - E, S_0 = -p^{1/4}."""
    return Chain([('D', g0), ('D', g1 - 1), ('D', 1 - g1), ('D', -g0), ('p', -1)], '1', -1, 4, (F(1, 4), -1))


# ---------------------------------------------------------------- VIR4 classes (Euler form)
class ChainTail(Chain):
    """prod D(c_j) psi + eps^{-(J-1)} [ p psi' - ((ma (p - 1) + mb)/y) psi ] = 0   (after E -> -E):
    in Euler form  prod (theta - e_i) psi + x^3 [ x^{3M} (theta - ma) + E (theta - mb) ] psi = 0."""

    def __init__(self, ops, h, root, ma_, mb_):
        Chain.__init__(self, ops, 'd', +1, h, root)
        self.ma_, self.mb_ = sp.sympify(ma_), sp.sympify(mb_)

    def solve(self, K):
        b0, c0 = self.root
        self.S = {0: El({(0, b0): P_(c0)})}
        J = self.J
        assert b0 * J == b0 + 1 and c0 ** J + c0 == 0          # S_0^J + p S_0 = 0
        coef = J * c0 ** (J - 1) + 1                             # (J S_0^{J-1} + p) = coef * p
        inv = El({(0, F(-1)): P_(sp.Rational(-1) / coef)})
        extra = El({(-1, F(1)): P_(-self.ma_), (-1, F(0)): P_(self.ma_ - self.mb_)})   # -(ma (p-1) + mb)/y
        for n in range(1, K + 1):
            A = Chain.A_series(self, n)
            rest = A.get(n, El())
            if n == 1:
                rest = el_add(rest, extra)
            self.S[n] = el_mul(rest, inv)
        return self.S

    def A_series(self, nmax):
        return Chain.A_series(self, nmax)


def class_A6():
    """Sixth order, frame (PX, pi):  prod_{i=1}^6 (theta - e_i) psi = x^6 (x^{6M} - E) psi,
    exponents {sh +- l0, sh +- l1, m1, 15 - 4 sh - m1};  S_0^6 = p."""
    ex = [sh + l0, sh - l0, sh + l1, sh - l1, m1, 15 - 4 * sh - m1]
    ops = [('D', ex[i] - i) for i in range(6)] + [('p', -1)]
    return Chain(ops, '1', -1, 6, (F(1, 6), -1))


def class_B4():
    """Fourth order, frame (pi_U, pi_W):  prod (theta - 3/2 -+ l_i) psi + x^3 [x^{3M}(theta - ma) - E (theta - mb)] psi = 0."""
    ex = [sp.Rational(3, 2) + l0, sp.Rational(3, 2) - l0, sp.Rational(3, 2) + l1, sp.Rational(3, 2) - l1]
    ops = [('D', ex[i] - i) for i in range(4)]
    return ChainTail(ops, 3, (F(1, 3), -1), ma, mb)


def class_B4b():
    """Candidate (b) of VIR3 in centred variables: D4 D3 p^{-1} D2 D1 psi + psi' = 0, g = l."""
    return Chain([('D', l0), ('D', l1 - 1), ('p', -1), ('D', 1 - l1), ('D', -l0)], 'd', +1, 3, (F(1, 3), -1))
