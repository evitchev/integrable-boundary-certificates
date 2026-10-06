"""EXC3: wkb3x EXTENDED to SECOND order in U = u^(-1) (x^(-c)):  + x^(-c) Q1(T) + x^(-2c) Q2(T),  Q1 = sum qcoef[j] T^j,
Q2 = sum q2coef[j] T^j.  Everything truncated at U^3.  theta(u^(-m)) = -m [n sigma (1-u)/w] u^(-m).
Integrals with u^(-m): M^(m)_l (A0 -> A0 - m), reduced with  dP M^(m)_l - M^(m)_(l+1) = dL M^(m-1)_l  and
M^(m)_0 = M^(m-1)_0 (A0 + B0 - m + 1)/(A0 - m).   With q2coef = None this must reproduce wkb3x exactly (regression in t0d).
Everything else as in wkb3x.py / VIR8 wkb3.py (see those docstrings)."""
import sys
from fractions import Fraction as F
from math import factorial
import sympy as sp

sys.dont_write_bytecode = True
l0, l1, A2, A3, A4, B0, B1, B2, C0, C1, U = sp.symbols('l0 l1 A2 A3 A4 B0 B1 B2 C0 C1 U')
GENS = (l0, l1, A2, A3, A4, B0, B1, B2, C0, C1, U)      # C0, C1: second-order markers (added after t0d's first run; t0d re-run)
UIDX = 10
MAXU = 2


def P_(x):
    return sp.Poly(x, *GENS, domain='QQ')


ZERO, ONE = P_(0), P_(1)


def trunc(p_):
    if p_.degree(U) <= MAXU:
        return p_
    return sp.Poly.from_dict({m: c for m, c in p_.terms() if m[UIDX] <= MAXU}, *GENS, domain='QQ')


def upart(p_, m):
    """the U^m part"""
    if p_.is_zero or p_.degree(U) < m:
        return ZERO
    d = {mm: c for mm, c in p_.terms() if mm[UIDX] == m}
    return sp.Poly.from_dict(d, *GENS, domain='QQ') if d else ZERO


def R_(x):
    return sp.Rational(x.numerator, x.denominator) if isinstance(x, F) else sp.Rational(x)


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
    def __init__(self, dP, dL, M, pcoef, lcoef_c, qcoef=None, q2coef=None):
        """pcoef[j]: coefficient of T^(dP-j) in P (pcoef[0] = 1); lcoef_c[j]: coefficient of T^(dL-j) in Lc(T) = L(T - c/2);
        qcoef[j] / q2coef[j]: coefficient of T^j in Q1 / Q2."""
        self.dP, self.dL = dP, dL
        self.M = sp.Rational(M)
        self.n = dP + dL / self.M
        self.sigma = self.n * self.M
        self.c = self.n + self.sigma
        self.e = dP - dL
        self.p = [P_(x) if not isinstance(x, sp.Poly) else x for x in pcoef]
        self.l = [P_(x) if not isinstance(x, sp.Poly) else x for x in lcoef_c]
        self.q = [[P_(x) if not isinstance(x, sp.Poly) else x for x in (qcoef or [])],
                  [P_(x) if not isinstance(x, sp.Poly) else x for x in (q2coef or [])]]
        assert self.p[0] == ONE and self.l[0] == ONE
        dPq, dLq = sp.Integer(dP), sp.Integer(dL)
        self.u = {0: P_(dPq / dLq), 1: P_(-1 / dLq)}
        self.g1 = {0: P_(self.n + self.sigma * dPq / dLq), 1: P_(-self.sigma / dLq)}
        self.g2 = {0: P_(-self.n * self.sigma * self.e / dLq), 1: P_(self.n * self.sigma / dLq)}
        self.dPmw = {0: P_(dPq), 1: P_(-1)}
        self.W = {}

    def theta(self, G):
        out = {}
        for i, X in G.items():
            t1 = la_mul(la_scale(self.g1, 1 - i), X)
            t2 = la_mul(la_mul(self.g2, self.dPmw), la_dw(X))
            r = la_add(t1, la_scale(t2, -1))
            for m in range(1, MAXU + 1):
                XU = {k_: upart(v, m) for k_, v in X.items()}
                XU = {k_: v for k_, v in XU.items() if not v.is_zero}
                if XU:
                    r = la_add(r, la_scale(la_mul(self.g2, XU), -m))          # theta(u^-m) = -m (g2/w) u^-m
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
        # deformation terms + x^(-k c) Q_k(T): after division by tau_0^dP:  U^k tau_0^(k dL - (k+1) dP) sum_j q_j Tt^j
        for kq, qlist in enumerate(self.q, start=1):
            base = (kq + 1) * self.dP - kq * self.dL
            for j, qj in enumerate(qlist):
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
                    term = gr_scale(term, qj * P_(ffac) * P_(U**kq))
                    total = gr_add(total, term)
        return total.get(N, {})

    def solve(self, K):
        self.W = {}
        for N in range(1, K + 1):
            rest = self.stage(N)
            self.W[N] = {k_ - 1: v * P_(-1) for k_, v in rest.items()}
        return self.W

    def masters_shifted(self, i, shift):
        """reduction M^(shift)_l = alpha_l M^(shift)_0 + beta_l M^(shift)_(-1) for exponent A0 - shift."""
        A0 = sp.Rational(1 - i) / self.sigma - shift
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
                    raise ZeroDivisionError('master recurrence singular (upward) at l = %d, shift %d' % (l, shift))
                val = (-(mid * a1 + lo * a2) / hi, -(mid * b1 + lo * b2) / hi)
            else:
                hi, mid, lo = rec_coeffs(-l - 1)
                a1, b1 = red(l + 1)
                a2, b2 = red(l + 2)
                if lo == 0:
                    raise ZeroDivisionError('master recurrence singular (downward) at l = %d, shift %d' % (l, shift))
                val = (-(mid * a1 + hi * a2) / lo, -(mid * b1 + hi * b2) / lo)
            cache[l] = val
            return val
        return red, A0, B0

    def J(self, i):
        """J_i = (1/(n sigma)) [cB M_0 + cG M_(-1)]: (cB, cG) as Polys (U set to 1 after the reduction)."""
        dP, dL = sp.Integer(self.dP), sp.Integer(self.dL)
        reds = [self.masters_shifted(i, m) for m in range(0, MAXU + 1)]
        A0, B0 = reds[0][1], reds[0][2]
        # M^(m)_0 and M^(m)_(-1) in terms of (M_0, M_(-1)):  pairs (coefficient of M_0, coefficient of M_(-1))
        m0 = [(sp.Integer(1), sp.Integer(0))]
        mm1 = [(sp.Integer(0), sp.Integer(1))]
        for m in range(1, MAXU + 1):
            ratio = (A0 + B0 - m + 1) / (A0 - m)                      # M^(m)_0 = M^(m-1)_0 * ratio
            m0.append((m0[-1][0] * ratio, m0[-1][1] * ratio))
            # dP M^(m)_(-1) - M^(m)_0 = dL M^(m-1)_(-1)
            mm1.append(((m0[-1][0] + dL * mm1[-1][0]) / dP, (m0[-1][1] + dL * mm1[-1][1]) / dP))
        cB, cG = ZERO, ZERO
        for lw, cf in self.W[i].items():
            for m in range(0, MAXU + 1):
                cm = upart(cf, m)
                if cm.is_zero:
                    continue
                al, be = reds[m][0](lw + 1)                             # M^(m)_(lw+1) = al M^(m)_0 + be M^(m)_(-1)
                cB = cB + cm * P_(al * m0[m][0] + be * mm1[m][0])
                cG = cG + cm * P_(al * m0[m][1] + be * mm1[m][1])
        sub = lambda p_: sp.Poly(p_.as_expr().subs(U, 1), *GENS, domain='QQ')
        return sub(cB), sub(cG)
