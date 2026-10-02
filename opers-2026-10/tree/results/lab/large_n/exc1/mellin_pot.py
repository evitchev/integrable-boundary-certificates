"""EXC1: copy of mellin_gen.py with (a) four generators, (b) POTENTIAL TERMS: Lambda = y^n p (1 + sum_w V_w), V_w a ring\nelement at order w (solve: W_N = (V_N - rest_N)/n), (c) R() generalised to y-powers a = -i + sigma m (integer m).\nVIR7: GENERAL copy of mellin_wkb.py: the symbol is Ht(T) = T^n sum_j h_j T^(-j) with ALL integer j (hlist indexed by j);\nonly stage() differs (h index, Z power j + m, alpha = n - j - m, ff(n - j, m)).\nVIR5: Mellin-WKB for an Euler-type operator with a general symbol.

Equation (theta = x d/dx = d/dv):      H(theta) psi = Lambda(x) psi,
    H(theta) = Ht(T),  T = s - theta,   Ht(T) = T^n sum_{j>=0} h_j T^(-2j)   (formal, n real),
    Lambda = x^n (x^{nM} + E)           (E -> -E already done; decaying branch T ~ +Lambda^(1/n)).
psi = exp(S), theta S = s - T(v).  Conjugation:  e^{-S} H(theta) e^{S} . 1
    = sum_m Ht^(m)(T) C_m,   sum_m mu^m C_m = exp( sum_{r>=2} (-1)^(r+1) mu^r T^(r-1)/r! ),
where T^(r) = theta^r T.  After x = E^(1/(nM)) y:  Lambda = y^n p, p = 1 + y^{nM}, T_0 = y p^(1/n),
T = T_0 (1 + W), W = sum_{i>=1} W_i eps^i.  Ring: sums of c(l0, l1) y^a p^b (b rational).
The local-charge polynomial of spin i - 1 is R_i, from J_i = -int_0^infty T_i dy/y
    = B(alpha, beta0)/(nM) * R_i,  alpha = -(i-1)/(nM),  beta0 = (i-1)(M+1)/(nM).

For the family C_k:  n = k + 4, M = 1/k, s = (n-1)/2,
    Ht(T) = (T^2 - l0^2)(T^2 - l1^2) (n/k)^k Gamma(kT/n + (k+1)/2) / Gamma(kT/n - (k-1)/2),
whose large-T expansion is  T^n (1 - l0^2/T^2)(1 - l1^2/T^2) sum_j c_j(k) (n/(kT))^(2j),
    c_j(k) = binom(k, 2j) (2j)! [t^(2j)] ( t / (2 sinh(t/2)) )^(k+1)      (Tricomi-Erdelyi).
For integer k the Gamma ratio is the finite product over the frozen exponents.
"""
import sys
from fractions import Fraction as F
from math import factorial
import sympy as sp

sys.dont_write_bytecode = True
l0, l1, z1, z2 = sp.symbols('l0 l1 z1 z2')      # EXC1: l0 = alpha, l1 = lambda, z1, z2 = apparent singularities
GENS = (l0, l1, z1, z2)


def P_(x):
    return sp.Poly(x, *GENS, domain='QQ')


ZERO, ONE = P_(0), P_(1)


def R_(x):
    return sp.Rational(x.numerator, x.denominator) if isinstance(x, F) else sp.Rational(x)


def el_add(u, v):
    out = dict(u)
    for k_, c in v.items():
        n_ = out.get(k_, ZERO) + c
        if n_.is_zero:
            out.pop(k_, None)
        else:
            out[k_] = n_
    return out


def el_mul(u, v):
    out = {}
    for (a1, b1), c1 in u.items():
        for (a2, b2), c2 in v.items():
            k_ = (a1 + a2, b1 + b2)
            n_ = out.get(k_, ZERO) + c1 * c2
            if n_.is_zero:
                out.pop(k_, None)
            else:
                out[k_] = n_
    return out


def el_scale(u, c):
    c = P_(c)
    if c.is_zero:
        return {}
    return {k_: v * c for k_, v in u.items()}


def ser_add(A, B):
    out = dict(A)
    for e, x in B.items():
        r = el_add(out.get(e, {}), x)
        if r:
            out[e] = r
        else:
            out.pop(e, None)
    return out


def ser_mul(A, B, emax):
    out = {}
    for e1, x1 in A.items():
        for e2, x2 in B.items():
            if e1 + e2 > emax:
                continue
            r = el_add(out.get(e1 + e2, {}), el_mul(x1, x2))
            if r:
                out[e1 + e2] = r
            else:
                out.pop(e1 + e2, None)
    return out


def ser_scale(A, c):
    return {e: el_scale(x, c) for e, x in A.items() if el_scale(x, c)}


def gamma_ratio_coeffs(k, jmax):
    """c_j(k), j = 0..jmax: Gamma(u + (k+1)/2)/Gamma(u - (k-1)/2) = u^k sum_j c_j u^(-2j)."""
    t = sp.Symbol('t')
    k = sp.Rational(k)
    # log( t / (2 sinh(t/2)) ) = - sum_{m>=1} B_{2m} t^(2m) / (2m (2m)!)
    L = sum(-sp.bernoulli(2 * m) * t**(2 * m) / (2 * m * sp.factorial(2 * m)) for m in range(1, jmax + 1))
    ser = sp.series(sp.exp((k + 1) * L), t, 0, 2 * jmax + 1).removeO()
    return [sp.nsimplify(sp.binomial(k, 2 * j) * sp.factorial(2 * j) * ser.coeff(t, 2 * j), rational=True) if False else
            sp.simplify(sp.ff(k, 2 * j) * ser.coeff(t, 2 * j)) for j in range(jmax + 1)]


def family_symbol(k, jmax):
    """h_j for C_k as Poly in (l0, l1): Ht(T) = T^n sum_j h_j T^(-2j)."""
    k = sp.Rational(k)
    n = k + 4
    c = gamma_ratio_coeffs(k, jmax)
    T2 = sp.Symbol('T2')          # stands for T^(-2)
    expr = sp.expand((1 - l0**2 * T2) * (1 - l1**2 * T2) * sum(c[j] * (n / k) ** (2 * j) * T2**j for j in range(jmax + 1)))
    return [P_(expr.coeff(T2, j)) for j in range(jmax + 1)]


class MellinWKB:
    def __init__(self, n, M, hlist):
        self.n = F(sp.Rational(n).p, sp.Rational(n).q)
        self.M = F(sp.Rational(M).p, sp.Rational(M).q)
        self.h = hlist                       # list of Poly, h[0] = 1
        self.hM = self.n * self.M            # exponent of the potential
        self.W = {}                          # i -> El

    def theta(self, u):
        """theta = y d/dy on y^a p^b:  y^a [ (a + hM b) p^b - hM b p^(b-1) ]."""
        out = {}
        for (a, b), c in u.items():
            t1 = c * P_(R_(F(a) + self.hM * b))
            t2 = c * P_(R_(-self.hM * b))
            for key, tt in (((a, b), t1), ((a, b - 1), t2)):
                if tt.is_zero:
                    continue
                n_ = out.get(key, ZERO) + tt
                if n_.is_zero:
                    out.pop(key, None)
                else:
                    out[key] = n_
        return out

    def stage(self, N):
        """coefficient of eps^N of  sum_{m,j} h_j ff(n-2j, m) Z^(2j+m) (1+W)^(n-2j-m) C_m  with W_N = 0."""
        n = self.n
        T0 = {(1, 1 / n): ONE}
        Wser = {i: w for i, w in self.W.items() if i < N}
        # T = T0 (1 + W): exponents of eps: -1 + i
        Tser = {-1: T0}
        for i, w in Wser.items():
            Tser[i - 1] = el_mul(T0, w)
        # derivatives
        derivs = [Tser]
        mmax = 2 * N
        for r in range(1, mmax):
            derivs.append({e: self.theta(x) for e, x in derivs[-1].items()})
            derivs[-1] = {e: x for e, x in derivs[-1].items() if x}
        # C_m by the exponential recursion  m C_m = sum_{r=2}^m r a_r C_{m-r},  a_r = (-1)^(r+1) T^(r-1)/r!
        C = [{0: {(0, F(0)): ONE}}]
        # C_m has minimal eps-exponent -floor(m/2); it multiplies Z^(2j+m) (eps^(2j+m)): keep exponents <= N - m
        for m in range(1, mmax + 1):
            acc = {}
            for r in range(2, m + 1):
                a_r = ser_scale(derivs[r - 1], sp.Rational((-1) ** (r + 1) * r, factorial(r)))
                acc = ser_add(acc, ser_mul(a_r, C[m - r], N - m))
            C.append(ser_scale(acc, sp.Rational(1, m)))
        # powers of W
        Wpow = [{0: {(0, F(0)): ONE}}]
        for q_ in range(1, N + 1):
            Wpow.append(ser_mul(Wpow[-1], Wser, N))
        total = {}
        for j in range(0, len(self.h)):          # h[j] multiplies T^(n-j)   (GENERAL: odd j allowed)
            if j > N:
                break
            if self.h[j].is_zero:
                continue
            for m in range(0, mmax + 1):
                emin = j + m - (m // 2)
                if emin > N or (m == 1):
                    continue
                alpha = n - j - m
                ffac = sp.ff(R_(n - j), m)
                if ffac == 0:
                    continue
                pw = {}
                for q_ in range(0, N + 1):
                    bq = sp.binomial(R_(alpha), q_)
                    if bq == 0:
                        continue
                    pw = ser_add(pw, ser_scale(Wpow[q_], bq))
                Z = {j + m: {(-(j + m), -F(j + m) / n): ONE}}
                term = ser_mul(ser_mul(Z, C[m], N), pw, N)
                term = {e: el_mul(x, {(0, F(0)): self.h[j] * P_(ffac)}) for e, x in term.items()}
                total = ser_add(total, term)
        return total.get(N, {})

    def solve(self, K, V=None):
        V = V or {}
        self.W = {}
        for N in range(1, K + 1):
            rest = self.stage(N)
            # coefficient of W_N in the eps^N term is n (from (1+W)^n);  rest_N + n W_N = V_N
            tot = el_add(el_scale(rest, -1), V.get(N, {}))
            self.W[N] = el_scale(tot, sp.Rational(1) / R_(self.n))
        return self.W

    def R(self, i):
        """R_i with J_i = -int T_i dy/y = B(alpha, beta0)/(nM) R_i; T_i = T0 W_i.  General y-power a = -i + sigma m (m integer):
        int y^(a+1) p^(b+1/n) dy/y = (1/sigma) B(alpha + m, beta0 + j - m),  j = (1-i)/n - (b + 1/n)."""
        n, M = self.n, self.M
        sigma = n * M
        alpha = R_(F(1 - i) / sigma)
        beta0 = R_(F(i - 1) * (M + 1) / (n * M))
        tot = R_(F(i - 1) / n)
        base_b = F(1 - i) / n
        def poch(x, m):
            return sp.rf(x, m) if m >= 0 else 1 / sp.rf(x + m, -m)
        out = sp.Integer(0)
        half = sp.Integer(0)       # terms with half-integer m: relative to B(alpha + 1/2, beta0 - 1/2); a total derivative at even orders
        ok = True
        for (a, b), c in self.W[i].items():
            mm = (F(a) + i) / sigma
            j = base_b - (b + 1 / n)
            if (2 * mm).denominator != 1 or (j - mm).denominator != 1:
                ok = False
                continue
            if mm.denominator == 1:
                mm, j = int(mm), int(j)
                out += -c.as_expr() * poch(alpha, mm) * poch(beta0, j - mm) / poch(tot, j)
            else:
                # half-integer m and j: B(alpha + m, beta0 + j - m) relative to B(alpha + 1/2, beta0)
                mh, jh = int(mm - F(1, 2)), int(j - F(1, 2))
                half += -c.as_expr() * poch(alpha + sp.Rational(1, 2), mh) * poch(beta0, jh - mh) / poch(tot + sp.Rational(1, 2), jh)
        self.last_half = sp.expand(half)
        return sp.expand(out), ok


def family_monic(k, K, PX, PI):
    """monic (in PX) charge polynomials of C_k from the Mellin-WKB, spins 1, 3, ..., K-1, and the odd orders."""
    k = sp.Rational(k)
    n = k + 4
    eng = MellinWKB(n, 1 / k, family_symbol(k, K // 2 + 1))
    eng.solve(K)
    c2 = n**2 * (k + 1) / k
    out, odd = {}, {}
    for i in range(2, K + 1):
        Ri, ok = eng.R(i)
        assert ok, 'inhomogeneous W_%d' % i
        if i % 2 == 1:
            odd[i] = sp.simplify(Ri) == 0
            continue
        Pl = sp.Poly(Ri, l0, l1)
        assert all(a % 2 == 0 and b % 2 == 0 for (a, b), _ in Pl.terms())
        Rm = sp.Poly(sp.expand(sum(cf * c2 ** ((a + b) // 2) * PX**a * PI**b for (a, b), cf in Pl.terms())), PX, PI)
        lead = Rm.coeff_monomial(PX**i)
        out[i - 1] = (None if lead == 0 else sp.Poly(sp.expand(Rm.as_expr() / lead), PX, PI), lead)
    return out, odd


def string_series(a, e, sigma, jmax):
    """expansion of a centred string  sigma^a R_a((T - e)/sigma) / T^a  in powers of 1/T, as a list of sympy expressions
    (index = power of 1/T):  (1 - e/T)^a sum_i c_i(a) sigma^(2i) (T - e)^(-2i)."""
    u = sp.Symbol('u')           # u = 1/T
    a = sp.Rational(a)
    c = gamma_ratio_coeffs(a, jmax // 2 + 1)
    expr = 0
    for i in range(jmax // 2 + 1):
        expr += c[i] * sigma ** (2 * i) * u ** (2 * i) * sp.series((1 - e * u) ** (a - 2 * i), u, 0, jmax + 1).removeO()
    expr = sp.expand(expr)
    return [expr.coeff(u, j) for j in range(jmax + 1)]


def symbol_from_blocks(blocks, jmax):
    """blocks: list of series (lists indexed by power of 1/T); returns hlist of Poly in (l0, l1)."""
    u = sp.Symbol('u')
    tot = sp.Integer(1)
    for b in blocks:
        tot = sp.expand(tot * sum(cf * u ** j for j, cf in enumerate(b)))
        tot = sum(tot.coeff(u, j) * u ** j for j in range(jmax + 1))
    return [P_(sp.expand(tot.coeff(u, j))) for j in range(jmax + 1)]
