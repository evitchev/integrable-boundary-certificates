"""NUM1 numerics (mpmath).  -psi'' + (x^(2M) + alpha x^(M-1) + lam/x^2) psi = E psi, lam = l(l+1).
D(E) = W[y, psi_+], y ~ x^(-a_(M+1)) exp(-x^(M+1)/(M+1)) (all decaying corrections), psi_+ = x^(l+1)(1 + ...).
logD_neg(e): E = -e < 0 via the Riccati q = -y'/y (stable inward); D_lin(E): linear (y, y') inward (any E)."""
import mpmath as mp
def asym_coeffs(M, alpha, lam, e, J):
    """q = sum_j a_j x^(M-j) with q^2 - q' = x^(2M) + alpha x^(M-1) + e + lam x^-2"""
    V = {0: mp.mpf(1)}
    V[M + 1] = V.get(M + 1, 0) + alpha; V[2 * M] = V.get(2 * M, 0) + e; V[2 * M + 2] = V.get(2 * M + 2, 0) + lam
    a = [mp.mpf(1)]
    for Jj in range(1, J + 1):
        s = sum(a[i] * a[Jj - i] for i in range(1, Jj))
        der = (2 * M + 1 - Jj) * a[Jj - M - 1] if Jj - M - 1 >= 0 else 0
        a.append((V.get(Jj, 0) - s + der) / 2)
    return a
def asym_y(M, alpha, lam, e, x, J=None):
    """log y(x), q(x) at large x from the asymptotic series (terms kept while decreasing)"""
    J = J or 60 * (M + 1)
    a = asym_coeffs(M, alpha, lam, e, J)
    L = -a[M + 1] * mp.log(x); q = mp.mpf(0); last = mp.inf
    done = False
    for j, aj in enumerate(a):
        tq = aj * x ** (M - j)
        q += tq
        if j != M + 1: L -= aj * x ** (M + 1 - j) / (M + 1 - j)
        if j > 2 * M + 2 and aj != 0:          # zero coefficients (parity) are skipped in the convergence test (v1 bug: a zero term ended the sum)
            if abs(tq) > last and abs(tq) > mp.mpf(10) ** (-mp.mp.dps): raise ValueError('asymptotic series not converged at x_max')
            last = abs(tq)
            if last < mp.mpf(10) ** (-mp.mp.dps - 5): done = True; break
    if not done: raise ValueError('asymptotic series: J too small')
    return L, q
def frob(M, alpha, l, e, x, nmax=100000):
    """psi_+ and psi_+' at x: b_k k(k+2l+1) = e b_(k-2) + alpha b_(k-M-1) + b_(k-2M-2) (E = -e)"""
    b = [mp.mpf(1)]; s = mp.mpf(0); sd = mp.mpf(0); k = 0; small = 0
    while True:
        if k > 0:
            v = 0
            if k - 2 >= 0: v += e * b[k - 2]
            if k - M - 1 >= 0: v += alpha * b[k - M - 1]
            if k - 2 * M - 2 >= 0: v += b[k - 2 * M - 2]
            b.append(v / (k * (k + 2 * l + 1)))
        term = b[k] * x ** (k + l + 1); s += term; sd += b[k] * (k + l + 1) * x ** (k + l)
        if k > 4 * M and abs(term) < abs(s) * mp.mpf(10) ** (-mp.mp.dps - 5):
            small += 1
            if small > 2 * M + 3: break
        else: small = 0
        k += 1
        if k > nmax: raise ValueError('Frobenius not converged')
    return s, sd
def xmax_for(M, alpha, lam, e):
    """smallest x on a 5% grid (starting where x^(2M) = 4|e|) at which the asymptotic series converges to 10^(-dps-5)"""
    x = max(mp.mpf(4 * abs(e) + 4) ** (mp.mpf(1) / (2 * M)), mp.mpf(2))
    while True:
        try: asym_y(M, alpha, lam, e, x); return x
        except ValueError: x *= mp.mpf('1.05')
def taylor_riccati(M, alpha, lam, e, x1, q1, L1, x0, order=None):
    """integrate q' = q^2 - V(x), L' = -q from x1 down to x0 by Taylor series (own; adaptive step from the coefficient decay)"""
    N = order or int(mp.mp.dps * 1.3) + 10
    x = mp.mpf(x1); q = mp.mpf(q1); L = mp.mpf(L1); steps = 0
    while x > x0:
        # Taylor coefficients of V about x in h (x + h): polynomial part exact, lam/x^2 via binomial series
        Vc = [mp.mpf(0)] * (N + 1)
        for p_, c in ((2 * M, mp.mpf(1)), (M - 1, alpha)):
            if c == 0: continue
            for n in range(0, min(p_, N) + 1): Vc[n] += c * mp.binomial(p_, n) * x ** (p_ - n)
        Vc[0] += e
        for n in range(N + 1): Vc[n] += lam * (n + 1) * (-1) ** n * x ** (-2 - n)
        c = [q]
        for n in range(N):
            c.append((sum(c[i] * c[n - i] for i in range(n + 1)) - Vc[n]) / (n + 1))
        # step: ratio-test radius, safety 0.5
        r = min(abs(c[n]) ** (-mp.mpf(1) / n) for n in range(N - 5, N + 1) if c[n] != 0)
        h = -min(r * mp.mpf('0.25'), x - x0)
        Lnew = L - sum(c[n] * h ** (n + 1) / (n + 1) for n in range(N + 1))
        q = sum(c[n] * h ** n for n in range(N + 1)); L = Lnew; x = x + h; steps += 1
        if steps > 200000: raise RuntimeError('too many steps')
    return q, L, steps
def logD_neg(M, alpha, l, e, x0=None, xmax=None):
    e = mp.mpf(e); lam = l * (l + 1)
    xmax = xmax or xmax_for(M, alpha, lam, e); x0 = x0 or 2 / mp.sqrt(e)
    L1, q1 = asym_y(M, alpha, lam, e, xmax)
    q0, L0, steps = taylor_riccati(M, alpha, lam, e, xmax, q1, L1, x0)
    ps, psd = frob(M, alpha, l, e, x0)
    return L0 + mp.log(psd + q0 * ps)
def logD_lin(M, alpha, l, E, x0=mp.mpf('0.25'), xmax=None):
    """linear integration of (log-scaled) y inward: valid for any E (used for N2 and N1 at E > 0)"""
    E = mp.mpf(E); e = -E; lam = l * (l + 1)
    xmax = xmax or xmax_for(M, alpha, lam, e)
    L1, q1 = asym_y(M, alpha, lam, e, xmax)
    V = lambda x: x ** (2 * M) + alpha * x ** (M - 1) + lam / x ** 2 - E
    # u = y / y(xmax): u'' = V u ; d/ds with s = xmax - x
    f = mp.odefun(lambda s, Y: [-Y[1], -V(xmax - s) * Y[0]], 0, [mp.mpf(1), -q1])
    u0, up0 = f(xmax - x0)
    ps, psd = frob(M, alpha, l, e, x0)
    D = (u0 * psd - up0 * ps) * mp.exp(L1)
    return D
