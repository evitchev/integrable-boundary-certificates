"""EXC2: wkb3 EXTENDED by a first-order term  + x^(-c) Qt(T) phi  (Qt = sum_j qcoef[j] T^j, coefficients linear in the
markers B0, B1, B2 times U, where U stands for u^(-1)); everything is truncated at first order in U.
  x^(-c) -> u^(-1) tau_0^(dL-dP);  theta(u^(-1)) = -[n sigma (1-u)/w] u^(-1);
  integrals with u^(-1): M'_l (A0 -> A0 - 1), reduced with  dP M'_l - M'_(l+1) = dL M_l  and  M'_0 = M_0 (A0+B0)/(A0-1).
VIR8: WKB for the THREE-TERM equation (class U, un-dualised form)

      [ P(T) - x^(c/2) L(T) x^(c/2) - E x^n ] phi = 0,     T = s - theta,  theta = x d/dx,  c = n (1 + M),

P, L monic of degrees dP > dL >= 1 (coefficients polynomial in l0, l1), n = dP + dL/M.  (E -> -E already done.)
phi = exp(S), theta S = s - Tt.  Conjugation (same identity as mellin_wkb):
      sum_m [ P^(m)(Tt) - x^c Lc^(m)(Tt) ] C_m = E x^n,     Lc(T) = L(T - c/2),
      sum_m mu^m C_m = exp( sum_{r>=2} (-1)^(r+1) mu^r Tt^(r-1)/r! ).
Leading order (momenta off): tau^dP - y^c tau^dL = y^n after x = E^(1/sigma) y, sigma = nM.  Uniformisation:
      u = y^c tau^(dL-dP):   tau_0 = u^(1/sigma) (1-u)^(-c/(n sigma)),   y^n = tau_0^dP (1-u),   u in (0, 1),
      theta = y d/dy = [n sigma (1-u)/w] u d/du,   w = dP - dL u,   theta tau_0 = tau_0 (n + sigma u)/w.
Everything is tau_0^(integer) x Laurent polynomial in w:  Tt = tau_0 sum_i tau_0^(-i) W_i(w), W_0 = 1.
Charges:  J_i = int Tt_i dy/y = (1/(n sigma)) int u^(A0-1) (1-u)^(B0) w W_i(w) du,
          A0 = (1-i)/sigma,  B0 = (i-1) c/(n sigma) - 1   (the SAME Beta exponents as the Gamma form),
reduced EXACTLY to two masters  M_0 = B(A0, B0+1)  and  M_(-1) = int u^(A0-1)(1-u)^B0 / w du  by the recurrence
      (A0+B0+1-l) M_(1-l) + [ -A0 e - (B0+1) dP + l (dP+e) ] M_(-l) - l dP e M_(-l-1) = 0,   e = dP - dL.
"""
import sys
from fractions import Fraction as F
from math import factorial
import sympy as sp

sys.dont_write_bytecode = True
l0, l1, A2, A3, A4, B0, B1, B2, U = sp.symbols('l0 l1 A2 A3 A4 B0 B1 B2 U')
GENS = (l0, l1, A2, A3, A4, B0, B1, B2, U)
UIDX = 8


def P_(x):
    return sp.Poly(x, *GENS, domain='QQ')


ZERO, ONE = P_(0), P_(1)


def trunc(p_):
    """first order in U"""
    if p_.degree(U) < 2:
        return p_
    return sp.Poly.from_dict({m: c for m, c in p_.terms() if m[UIDX] < 2}, *GENS, domain='QQ')


def upart(p_):
    """the U^1 part"""
    return sp.Poly.from_dict({m: c for m, c in p_.terms() if m[UIDX] == 1}, *GENS, domain='QQ') if p_.degree(U) >= 1 else ZERO


def R_(x):
    return sp.Rational(x.numerator, x.denominator) if isinstance(x, F) else sp.Rational(x)


# ---- Laurent polynomials in w: dict {power: Poly}
def la_add(A, B):
    out = dict(A)
    for k_, c in B.items():
        v = out.get(k_, ZERO) + c
        if v.is_zero:
            out.pop(k_, None)
        else:
            out[k_] = v
    return out


def la_mul(A, B):
    out = {}
    for k1, c1 in A.items():
        for k2, c2 in B.items():
            v = out.get(k1 + k2, ZERO) + trunc(c1 * c2)
            if v.is_zero:
                out.pop(k1 + k2, None)
            else:
                out[k1 + k2] = v
    return out


def la_scale(A, c):
    c = P_(c) if not isinstance(c, sp.Poly) else c
    if c.is_zero:
        return {}
    return {k_: w_ for k_, w_ in ((k_, trunc(v * c)) for k_, v in A.items()) if not w_.is_zero}


def la_dw(A):
    return {k_ - 1: v * P_(k_) for k_, v in A.items() if k_ != 0}


# ---- graded series: dict {order: Laurent}
def gr_add(A, B):
    out = dict(A)
    for o, x in B.items():
        r = la_add(out.get(o, {}), x)
        if r:
            out[o] = r
        else:
            out.pop(o, None)
    return out


def gr_mul(A, B, omax):
    out = {}
    for o1, x1 in A.items():
        for o2, x2 in B.items():
            if o1 + o2 > omax:
                continue
            r = la_add(out.get(o1 + o2, {}), la_mul(x1, x2))
            if r:
                out[o1 + o2] = r
            else:
                out.pop(o1 + o2, None)
    return out


def gr_scale(A, c):
    out = {}
    for o, x in A.items():
        r = la_scale(x, c)
        if r:
            out[o] = r
    return out


class WKB3:
    def __init__(self, dP, dL, M, pcoef, lcoef_c, qcoef=None):
        """pcoef[j]: coefficient of T^(dP-j) in P (Poly, pcoef[0] = 1); lcoef_c[j]: coefficient of T^(dL-j) in Lc(T) = L(T - c/2)."""
        self.dP, self.dL = dP, dL
        self.M = sp.Rational(M)
        self.n = dP + dL / self.M
        self.sigma = self.n * self.M
        self.c = self.n + self.sigma
        self.e = dP - dL
        self.p = [P_(x) if not isinstance(x, sp.Poly) else x for x in pcoef]
        self.l = [P_(x) if not isinstance(x, sp.Poly) else x for x in lcoef_c]
        self.q = [P_(x) if not isinstance(x, sp.Poly) else x for x in (qcoef or [])]      # qcoef[j] multiplies T^j
        assert self.p[0] == ONE and self.l[0] == ONE
        dPq, dLq = sp.Integer(dP), sp.Integer(dL)
        self.u = {0: P_(dPq / dLq), 1: P_(-1 / dLq)}                       # u = (dP - w)/dL
        self.g1 = {0: P_(self.n + self.sigma * dPq / dLq), 1: P_(-self.sigma / dLq)}     # n + sigma u
        self.g2 = {0: P_(-self.n * self.sigma * self.e / dLq), 1: P_(self.n * self.sigma / dLq)}   # n sigma (1 - u)
        self.dPmw = {0: P_(dPq), 1: P_(-1)}                               # dP - w = dL u
        self.W = {}

    def theta(self, G):
        """theta on a graded object sum_i tau_0^(1-i) X_i (key = i): (1/w)[(1-i) g1 X - g2 (dP - w) X']"""
        out = {}
        for i, X in G.items():
            t1 = la_mul(la_scale(self.g1, 1 - i), X)
            t2 = la_mul(la_mul(self.g2, self.dPmw), la_dw(X))
            r = la_add(t1, la_scale(t2, -1))
            XU = {k_: upart(v) for k_, v in X.items()}
            XU = {k_: v for k_, v in XU.items() if not v.is_zero}
            if XU:
                r = la_add(r, la_scale(la_mul(self.g2, XU), -1))                 # theta(u^-1) = -(g2/w) u^-1
            r = {k_ - 1: v for k_, v in r.items()}
            if r:
                out[i] = r
        return out

    def stage(self, N):
        Tser = {0: {0: ONE}}
        for i, w_ in self.W.items():
            if i < N:
                Tser[i] = w_
        derivs = [Tser]
        mmax = 2 * N
        for r in range(1, mmax):
            derivs.append(self.theta(derivs[-1]))
        # a_r = (-1)^(r+1) Tt^(r-1)/r!  with Tt^(r) = tau_0 sum_i tau_0^(-i) (...): graded order i - 1
        def shifted(G):
            return {i - 1: X for i, X in G.items()}
        C = [{0: {0: ONE}}]
        for m in range(1, mmax + 1):
            acc = {}
            for r in range(2, m + 1):
                a_r = gr_scale(shifted(derivs[r - 1]), sp.Rational((-1) ** (r + 1) * r, factorial(r)))
                acc = gr_add(acc, gr_mul(a_r, C[m - r], N - m))
            C.append(gr_scale(acc, sp.Rational(1, m)))
        Wser = {i: w_ for i, w_ in self.W.items() if i < N}
        Wpow = [{0: {0: ONE}}]
        for q_ in range(1, N + 1):
            Wpow.append(gr_mul(Wpow[-1], Wser, N))
        total = {}
        for which, coefs, deg in (('P', self.p, self.dP), ('L', self.l, self.dL)):
            for j in range(0, deg + 1):
                if j > N or coefs[j].is_zero:
                    continue
                for m in range(0, mmax + 1):
                    if m == 1 or j + m - (m // 2) > N:
                        continue
                    ffac = sp.ff(sp.Integer(deg - j), m)
                    if ffac == 0:
                        continue
                    alpha = deg - j - m
                    pw = {}
                    for q_ in range(0, N + 1):
                        bq = sp.binomial(alpha, q_)
                        if bq == 0:
                            continue
                        pw = gr_add(pw, gr_scale(Wpow[q_], bq))
                    Z = {j + m: {0: ONE}}
                    term = gr_mul(gr_mul(Z, C[m], N), pw, N)
                    term = gr_scale(term, coefs[j] * P_(ffac))
                    if which == 'L':
                        term = {o: la_scale(la_mul(self.u, x), -1) for o, x in term.items()}
                    total = gr_add(total, term)
        # first-order term  + x^(-c) Qt(T):  after division by tau_0^dP it is  U tau_0^(dL - 2 dP) sum_j q_j Tt^j
        base = 2 * self.dP - self.dL
        for j, qj in enumerate(self.q):
            if qj.is_zero:
                continue
            for m in range(0, min(j, mmax) + 1):
                if m == 1:
                    continue
                order0 = base - j + m - (m // 2)
                if order0 > N:
                    continue
                ffac = sp.ff(sp.Integer(j), m)
                if ffac == 0:
                    continue
                alpha = j - m
                pw = {}
                for q_ in range(0, N + 1):
                    bq = sp.binomial(alpha, q_)
                    if bq == 0:
                        continue
                    pw = gr_add(pw, gr_scale(Wpow[q_], bq))
                Z = {base - j + m: {0: ONE}}
                term = gr_mul(gr_mul(Z, C[m], N), pw, N)
                term = gr_scale(term, qj * P_(ffac) * P_(U))
                total = gr_add(total, term)
        return total.get(N, {})

    def solve(self, K):
        self.W = {}
        for N in range(1, K + 1):
            rest = self.stage(N)
            self.W[N] = {k_ - 1: v * P_(-1) for k_, v in rest.items()}        # W_N = -rest/w
        return self.W

    def masters(self, i):
        """returns function red(l) -> (alpha_l, beta_l) with M_l = alpha_l M_0 + beta_l M_(-1), exact rationals."""
        A0 = sp.Rational(1 - i) / self.sigma
        B0 = sp.Rational(i - 1) * self.c / (self.n * self.sigma) - 1
        dP, e = sp.Integer(self.dP), sp.Integer(self.e)
        cache = {0: (sp.Integer(1), sp.Integer(0)), -1: (sp.Integer(0), sp.Integer(1))}
        def rec_coeffs(l):       # (A0+B0+1-l) M_(1-l) + mid(l) M_(-l) - l dP e M_(-l-1) = 0
            return A0 + B0 + 1 - l, -A0 * e - (B0 + 1) * dP + l * (dP + e), -l * dP * e
        def red(l):
            if l in cache:
                return cache[l]
            if l > 0:
                # use the relation with index l' = 1 - l  (so that M_(1-l') = M_l):  hi M_l + mid M_(l-1) + lo M_(l-2) = 0
                hi, mid, lo = rec_coeffs(1 - l)
                a1, b1 = red(l - 1)
                a2, b2 = red(l - 2) if lo != 0 else (0, 0)
                if hi == 0:
                    raise ZeroDivisionError('master recurrence singular (upward) at l = %d' % l)
                val = (-(mid * a1 + lo * a2) / hi, -(mid * b1 + lo * b2) / hi)
            else:
                # l <= -2: relation with index l' = -l - 1:  hi M_(l+2) + mid M_(l+1) + lo M_l = 0
                hi, mid, lo = rec_coeffs(-l - 1)
                a1, b1 = red(l + 1)
                a2, b2 = red(l + 2)
                if lo == 0:
                    raise ZeroDivisionError('master recurrence singular (downward) at l = %d' % l)
                val = (-(mid * a1 + hi * a2) / lo, -(mid * b1 + hi * b2) / lo)
            cache[l] = val
            return val
        return red, A0, B0

    def J(self, i):
        """J_i = (1/(n sigma)) [ cB M_0 + cG M_(-1) ]: returns (cB, cG) as Polys (U set to 1 after the reduction)."""
        red, A0, B0 = self.masters(i)
        dP, dL = sp.Integer(self.dP), sp.Integer(self.dL)
        # primed family: exponent A0 - 1
        save = None
        redp, _, _ = self.masters_shifted(i)
        ratio = (A0 + B0) / (A0 - 1)                 # M'_0 = B(A0-1, B0+1) = ratio * M_0
        cB, cG = ZERO, ZERO
        for lw, cf in self.W[i].items():
            c0 = sp.Poly.from_dict({m: c for m, c in cf.terms() if m[UIDX] == 0}, *GENS, domain='QQ') if not cf.is_zero else ZERO
            c1 = upart(cf)
            if not c0.is_zero:
                al, be = red(lw + 1)
                cB = cB + c0 * P_(al); cG = cG + c0 * P_(be)
            if not c1.is_zero:
                alp, bep = redp(lw + 1)            # M'_(lw+1) = alp M'_0 + bep M'_(-1);  M'_(-1) = (M'_0 + dL M_(-1))/dP
                cB = cB + c1 * P_((alp + bep / dP) * ratio)
                cG = cG + c1 * P_(bep * dL / dP)
        sub = lambda p_: sp.Poly(p_.as_expr().subs(U, 1), *GENS, domain='QQ')
        return sub(cB), sub(cG)

    def masters_shifted(self, i):
        """the same reduction with A0 -> A0 - 1"""
        A0 = sp.Rational(1 - i) / self.sigma - 1
        B0 = sp.Rational(i - 1) * self.c / (self.n * self.sigma) - 1
        dP, e = sp.Integer(self.dP), sp.Integer(self.e)
        cache = {0: (sp.Integer(1), sp.Integer(0)), -1: (sp.Integer(0), sp.Integer(1))}
        def rec_coeffs(l):
            return A0 + B0 + 1 - l, -A0 * e - (B0 + 1) * dP + l * (dP + e), -l * dP * e
        def red(l):
            if l in cache:
                return cache[l]
            if l > 0:
                hi, mid, lo = rec_coeffs(1 - l)
                a1, b1 = red(l - 1)
                a2, b2 = red(l - 2) if lo != 0 else (0, 0)
                if hi == 0:
                    raise ZeroDivisionError('shifted master recurrence singular (upward) at l = %d' % l)
                val = (-(mid * a1 + lo * a2) / hi, -(mid * b1 + lo * b2) / hi)
            else:
                hi, mid, lo = rec_coeffs(-l - 1)
                a1, b1 = red(l + 1)
                a2, b2 = red(l + 2)
                if lo == 0:
                    raise ZeroDivisionError('shifted master recurrence singular (downward) at l = %d' % l)
                val = (-(mid * a1 + hi * a2) / lo, -(mid * b1 + hi * b2) / lo)
            cache[l] = val
            return val
        return red, A0, B0
