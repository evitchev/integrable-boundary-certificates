"""NUM2 generic order-D solver (own code; extends NUM1's spec.py).  ODE in theta = x d/dx (u = log x):
    sum_r x^(p_r) R_r(theta) y = 0,   r = 0: p_0 = 0, R_0 = P monic of degree D;   r >= 1: extra terms (polynomials R_r of degree < D).
Polynomials are lists of mpf coefficients [c_0, c_1, ...] (c_m multiplies theta^m).  Exponents p_r are Fractions.
- asym(): the large-x series of S = theta log y for the decaying solution, S = sum_i a_i x^(gamma - i delta), computed INCREMENTALLY through
  W_m = (theta + S)^m 1 (W_m[i] = coefficient of x^(m gamma - i delta)); log y = sum a_i x^e/e (+ a_i0 u for e = 0), theta^m y = y W_m.
- integrate(): the linear ODE in u by Taylor steps with step doubling and log-rescaling (inward: the decaying solution grows fastest inward).
- frob(): Frobenius series chi_j = x^(theta_j) sum b x^(offset) from P(nu) b(nu) = -sum_r R_r(nu - p_r) b(nu - p_r).
- Q(): projection y = sum_j Q_j chi_j at x0 on (theta^m .)_(m < D)."""
import mpmath as mp
from fractions import Fraction as Fr
from math import gcd
def pev(c, v):
    r = mp.mpf(0)
    for co in reversed(c): r = r * v + co
    return r
def fgcd(a, b):
    a, b = Fr(a), Fr(b)
    return Fr(gcd(a.numerator * b.denominator, b.numerator * a.denominator), a.denominator * b.denominator)
class Op:
    def __init__(s, P, terms):
        """P: list (monic, degree D); terms: list of (p (Fraction > 0), R (list))"""
        s.P = [mp.mpf(x) for x in P]; s.D = len(P) - 1; assert s.P[-1] == 1
        s.terms = [(Fr(p), [mp.mpf(x) for x in R]) for p, R in terms]
        # dominant balance: S^D + rho x^pmax S^dR = 0
        pmax, Rmax = max(s.terms, key=lambda t: t[0])
        dR = max(i for i, x in enumerate(Rmax) if x != 0); rho = Rmax[dR]
        s.gamma = pmax / (s.D - dR); s.dR = dR; s.rho = rho
        q = s.D - dR                                   # a0^q = -rho, decaying root: a0 < 0 real
        cand = -rho
        if q % 2 == 1: a0 = mp.sign(cand) * abs(cand) ** (mp.mpf(1) / q)
        else:
            assert cand > 0, 'no real root'; a0 = -cand ** (mp.mpf(1) / q)
        assert a0 < 0, 'leading coefficient not decaying'
        s.a0 = a0
        d = s.gamma
        for p, R in s.terms: d = fgcd(d, s.D * s.gamma - p)
        s.delta = d; s.g = int(s.gamma / d)
        s.J0 = s.D * a0 ** (s.D - 1) + rho * dR * a0 ** (dR - 1)
        s._a = None
    def asym_coeffs(s, J):
        """a_0..a_J and W_m[0..J] for m = 0..D"""
        D, g = s.D, s.g; gam, dl = s.gamma, s.delta
        a = [s.a0] + [mp.mpf(0)] * J
        W = [[mp.mpf(0)] * (J + 1) for _ in range(D + 1)]
        W[0][0] = mp.mpf(1)
        # offsets: term x^p R_m W_m contributes at exponent m gam - i' dl + p == D gam - i dl  ->  i' = i - (D gam - m gam - p)/dl
        offs = []
        for (p, R) in [(Fr(0), s.P)] + s.terms:
            for m, co in enumerate(R):
                if co == 0: continue
                o = (D * gam - m * gam - p) / dl
                assert o.denominator == 1; offs.append((m, int(o), co))
        for i in range(0, J + 1):
            if i > 0: a[i] = mp.mpf(0)
            for _pass in range(2):
                for m in range(1, D + 1):
                    v = mp.mpf(0)
                    if i - g >= 0: v += (Fr(m - 1) * gam - (i - g) * dl) * W[m - 1][i - g]
                    v += mp.fsum(a[i1] * W[m - 1][i - i1] for i1 in range(0, i + 1))
                    W[m][i] = v
                if i == 0: break
                if _pass == 0:
                    F = mp.fsum(co * W[m][i - o] for (m, o, co) in offs if 0 <= i - o <= J)
                    a[i] = -F / s.J0
        s._a = a; s._W = W; s._J = J
        return a, W
    def asym_eval(s, x, J=None):
        """log y and (theta^m y / y)_(m<D) at x; block-max convergence on the W_m and log-y terms"""
        if s._a is None or (J and J > s._J): s.asym_coeffs(J or 400)
        a, W, J = s._a, s._W, s._J
        L = mp.mpf(0); Wv = [mp.mpf(0)] * s.D; thr = mp.mpf(10) ** (-mp.mp.dps - 5)
        B = max(4 * s.g, 8); bm = mp.mpf(0); prev = mp.inf
        for i in range(J + 1):
            e = s.gamma - i * s.delta
            tl = (a[i] * x ** e / e) if e != 0 else a[i] * mp.log(x)
            L += tl
            tw = [W[m][i] * x ** (m * s.gamma - i * s.delta) for m in range(s.D)]
            for m in range(s.D): Wv[m] += tw[m]
            if i > 2 * s.g:
                bm = max([bm, abs(tl)] + [abs(t) / (abs(Wv[m]) + 1) for m, t in enumerate(tw) if m > 0])
                if (i + 1) % B == 0:
                    if bm > prev and bm > mp.mpf(10) ** (-mp.mp.dps): raise ValueError('asymptotic series not converged')
                    if bm < thr and prev < thr: return L, Wv
                    prev = bm; bm = mp.mpf(0)
        raise ValueError('asymptotic series: J too small')
    def xmax(s, x_start, J=600):
        x = mp.mpf(x_start)
        while True:
            try: s.asym_eval(x, J); return x
            except ValueError: x *= mp.mpf('1.1')
    # ---- inward Taylor integration in u
    def _coef_taylor(s, u0, N):
        """C_m(u0 + h) Taylor coefficients, m < D: C_m = P_m + sum_r R_(r,m) e^(p_r u)"""
        C = [[mp.mpf(0)] * (N + 1) for _ in range(s.D)]
        for m in range(s.D): C[m][0] += s.P[m]
        for p, R in s.terms:
            pf = mp.mpf(p.numerator) / p.denominator; base = mp.exp(pf * u0)
            tj = [base]
            for j in range(1, N + 1): tj.append(tj[-1] * pf / j)
            for m, co in enumerate(R):
                if co == 0: continue
                for j in range(N + 1): C[m][j] += co * tj[j]
        return C
    def _tables(s, N):
        if getattr(s, '_tabN', None) != N:
            D = s.D
            s._ff = [[mp.mpf(mp.rf(n + 1, m)) if m > 0 else mp.mpf(1) for m in range(D + 1)] for n in range(N + 1)]
            s._bin = [[mp.mpf(mp.binomial(n, m)) for m in range(D)] for n in range(N + 1)]
            s._tabN = N
    def _step(s, u0, Y, h, N):
        """Y = Taylor head (Y_0..Y_(D-1)) at u0; returns the head at u0 + h (v2: integer tables precomputed; v1 kept in ode_v1_slow.py)"""
        D = s.D; C = s._coef_taylor(u0, N); s._tables(N); ff = s._ff; bn = s._bin
        c = list(Y) + [mp.mpf(0)] * (N + 1 - D)
        for n in range(0, N + 1 - D):
            acc = mp.mpf(0)
            for m in range(D):
                Cm = C[m]
                for j in range(0, n + 1):
                    if Cm[j] != 0: acc += Cm[j] * ff[n - j][m] * c[n - j + m]
            c[n + D] = -acc / ff[n][D]
        hp = [mp.mpf(1)]
        for i in range(N): hp.append(hp[-1] * h)
        return [mp.fsum(c[n] * bn[n][m] * hp[n - m] for n in range(m, N + 1)) for m in range(D)]
    def integrate(s, u1, Y1, u2, N=40):
        tol = mp.mpf(10) ** (-(mp.mp.dps - 15)); u = mp.mpf(u1); Y = list(Y1); ls = mp.mpf(0)
        sgn = 1 if u2 > u1 else -1; h = sgn * mp.mpf('0.02')
        while (u2 - u) * sgn > 0:
            if (u + h - u2) * sgn > 0: h = u2 - u
            a = s._step(u, Y, h, N)
            b1 = s._step(u, Y, h / 2, N); b = s._step(u + h / 2, b1, h / 2, N)
            sc = sum(abs(x) for x in b) + mp.mpf(10) ** -40
            err = sum(abs(x - y) for x, y in zip(a, b)) / sc
            if err < tol:
                u = u + h; r = sum(abs(x) for x in b); Y = [x / r for x in b]; ls += mp.log(r)
                if err < tol / 1000: h *= mp.mpf('1.5')
            else: h /= 2
        return Y, ls
    # ---- Frobenius
    def frob(s, th, x0, maxterms=200000):
        """chi = x^th sum_off b x^off; returns [theta^m chi](x0), m < D"""
        b = {Fr(0): mp.mpf(1)}; ps = sorted({p for p, _ in s.terms})
        thf = th
        order = [Fr(0)]; seen = {Fr(0)}
        vals = [mp.mpf(0)] * s.D; idx = 0; small = 0
        import heapq
        heap = [Fr(0)]
        while heap:
            off = heapq.heappop(heap)
            if off != 0:
                nu = thf + mp.mpf(off.numerator) / off.denominator
                rhs = mp.mpf(0)
                for p, R in s.terms:
                    if off - p in b: rhs -= pev(R, nu - mp.mpf(p.numerator) / p.denominator) * b[off - p]
                den = pev(s.P, nu); assert abs(den) > mp.mpf(10) ** (-mp.mp.dps // 2), 'resonance'
                b[off] = rhs / den
            nu = thf + mp.mpf(off.numerator) / off.denominator
            term = b[off] * x0 ** nu
            for m in range(s.D): vals[m] += term * nu ** m
            if abs(term) < (abs(vals[0]) + 1e-300) * mp.mpf(10) ** (-mp.mp.dps - 8) and off > 0: small += 1
            else: small = 0
            if small > 40: break
            for p in ps:
                no = off + p
                if no not in seen: seen.add(no); heapq.heappush(heap, no)
            idx += 1
            if idx > maxterms: raise ValueError('Frobenius not converged')
        return vals
    def Q(s, thetas, x0, xm=None, J=600):
        """Q_j for the Frobenius exponents thetas: returns list of log Q_j (complex-safe via mp.log of real values) and Q_j"""
        xm = xm or s.xmax(max(mp.mpf(2), mp.mpf(1)), J)
        L, Wv = s.asym_eval(xm, J)
        # head: Taylor coefficients in u: Y_m = theta^m y / m!
        Y1 = [Wv[m] / mp.factorial(m) for m in range(s.D)]
        Y0, ls = s.integrate(mp.log(xm), Y1, mp.log(x0))
        V = [Y0[m] * mp.factorial(m) for m in range(s.D)]   # theta^m y at x0, times exp(-(L + ls))
        Fm = mp.matrix(s.D, s.D)
        for j, th in enumerate(thetas):
            col = s.frob(th, x0)
            for m in range(s.D): Fm[m, j] = col[m]
        q = mp.lu_solve(Fm, mp.matrix(V))
        return [mp.log(q[j]) + L + ls if q[j] > 0 else mp.log(-q[j]) + L + ls + 1j * mp.pi for j in range(s.D)], [q[j] * mp.exp(L + ls) for j in range(s.D)]
