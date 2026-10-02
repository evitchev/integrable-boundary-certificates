"""NUM1 numerics (mpmath).  -psi'' + (x^(2M) + alpha x^(M-1) + lam/x^2) psi = E psi, lam = l(l+1).
D(E) = W[y, psi_+], y ~ x^(-a_(M+1)) exp(-x^(M+1)/(M+1)) (all decaying corrections), psi_+ = x^(l+1)(1 + ...).
logD_neg(e): E = -e < 0 via the Riccati q = -y'/y (stable inward); D_lin(E): linear (y, y') inward (any E)."""
import mpmath as mp
_ACACHE = {}
def _coeff_gen(M, alpha, lam, e):
    """cached incremental a_j of q = sum_j a_j x^(M-j), q^2 - q' = x^(2M) + alpha x^(M-1) + e + lam x^-2"""
    key = (M, alpha, lam, e, mp.mp.dps)
    if key not in _ACACHE:
        V = {0: mp.mpf(1)}
        V[M + 1] = V.get(M + 1, 0) + alpha; V[2 * M] = V.get(2 * M, 0) + e; V[2 * M + 2] = V.get(2 * M + 2, 0) + lam
        _ACACHE[key] = (V, [mp.mpf(1)])
    return _ACACHE[key]
def _coeff(M, alpha, lam, e, j):
    V, a = _coeff_gen(M, alpha, lam, e)
    while len(a) <= j:
        Jj = len(a)
        s_ = mp.fsum(a[i] * a[Jj - i] for i in range(1, Jj))
        der = (2 * M + 1 - Jj) * a[Jj - M - 1] if Jj - M - 1 >= 0 else 0
        a.append((V.get(Jj, 0) - s_ + der) / 2)
    return a[j]
def asym_y(M, alpha, lam, e, x, J=6000):
    """log y(x), q(x) at large x from the asymptotic series.  Convergence judged on BLOCK maxima of 2M+2 consecutive terms (v2: the term-wise test
    raised on the non-monotone pattern and forced x_max ~ 9); raises if a block maximum grows while above 10^(-dps)."""
    L = -_coeff(M, alpha, lam, e, M + 1) * mp.log(x); q = mp.mpf(0)
    B = 2 * M + 2; blockmax = mp.mpf(0); prev = mp.inf; thr = mp.mpf(10) ** (-mp.mp.dps - 5)
    for j in range(J + 1):
        aj = _coeff(M, alpha, lam, e, j)
        tq = aj * x ** (M - j); tl = aj * x ** (M + 1 - j) / (M + 1 - j) if j != M + 1 else 0
        q += tq; L -= tl
        if j > 2 * M + 2:
            blockmax = max(blockmax, abs(tq), abs(tl))
            if (j + 1) % B == 0:
                if blockmax > prev and blockmax > mp.mpf(10) ** (-mp.mp.dps): raise ValueError('asymptotic series not converged at x_max')
                if blockmax < thr and prev < thr: return L, q
                prev = blockmax; blockmax = mp.mpf(0)
    raise ValueError('asymptotic series: J too small')
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
        except ValueError: x *= mp.mpf('1.1')
def _taylor_coeffs(M, alpha, lam, e, x, q, N):
    Vc = [mp.mpf(0)] * (N + 1)
    for p_, c in ((2 * M, mp.mpf(1)), (M - 1, alpha)):
        if c == 0: continue
        for n in range(0, min(p_, N) + 1): Vc[n] += c * mp.binomial(p_, n) * x ** (p_ - n)
    Vc[0] += e
    if lam != 0:
        for n in range(N + 1): Vc[n] += lam * (n + 1) * (-1) ** n * x ** (-2 - n)
    c = [q]
    for n in range(N):
        c.append((mp.fsum(c[i] * c[n - i] for i in range(n + 1)) - Vc[n]) / (n + 1))
    return c
def _step(M, alpha, lam, e, x, q, h, N):
    c = _taylor_coeffs(M, alpha, lam, e, x, q, N)
    qn = mp.fsum(c[n] * h ** n for n in range(N + 1)); dL = -mp.fsum(c[n] * h ** (n + 1) / (n + 1) for n in range(N + 1))
    return qn, dL
def taylor_riccati(M, alpha, lam, e, x1, q1, L1, x0, order=None):
    """q' = q^2 - V(x), L' = -q from x1 DOWN to x0; Taylor steps with step-doubling error control (v2; v1's ratio-test radius was fooled by the
    stiff noise ~ (2q)^n/n! and took ~10^4 steps).  Working precision = caller's dps + 20 (set by the caller)."""
    N = order or 40
    tol = mp.mpf(10) ** (-(mp.mp.dps - 15))
    x = mp.mpf(x1); q = mp.mpf(q1); L = mp.mpf(L1); steps = 0
    h = -min(mp.mpf(N) / (8 * abs(q)), x / 8)
    while x > x0:
        if x + h < x0: h = x0 - x
        q1_, d1 = _step(M, alpha, lam, e, x, q, h, N)
        qa, da = _step(M, alpha, lam, e, x, q, h / 2, N); q2_, db = _step(M, alpha, lam, e, x + h / 2, qa, h / 2, N)
        err = abs(q1_ - q2_) / (abs(q2_) + 1) + abs(d1 - (da + db))
        if err < tol:
            x = x + h; q = q2_; L = L + da + db; steps += 1
            if err < tol / 1000: h = h * mp.mpf('1.5')
        else:
            h = h / 2
        if steps > 100000: raise RuntimeError('too many steps')
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
# ---------------- linear Taylor integrator (any E): u'' = (V - E) u
def _lin_coeffs(M, alpha, lam, E, x, u, up, N):
    Vc = [mp.mpf(0)] * (N + 1)
    for p_, c in ((2 * M, mp.mpf(1)), (M - 1, alpha)):
        if c == 0: continue
        for n in range(0, min(p_, N) + 1): Vc[n] += c * mp.binomial(p_, n) * x ** (p_ - n)
    Vc[0] -= E
    if lam != 0:
        for n in range(N + 1): Vc[n] += lam * (n + 1) * (-1) ** n * x ** (-2 - n)
    c = [u, up]
    for n in range(N - 1):
        c.append(mp.fsum(Vc[i] * c[n - i] for i in range(n + 1)) / ((n + 2) * (n + 1)))
    return c
def _lin_step(M, alpha, lam, E, x, u, up, h, N):
    c = _lin_coeffs(M, alpha, lam, E, x, u, up, N)
    return mp.fsum(c[n] * h ** n for n in range(N + 1)), mp.fsum(n * c[n] * h ** (n - 1) for n in range(1, N + 1))
def lin_integrate(M, alpha, lam, E, x1, u1, up1, x2, N=40):
    """(u, u') from x1 to x2 (either direction), step doubling; values rescaled by a running log factor to avoid overflow; returns (u, up, logscale)"""
    tol = mp.mpf(10) ** (-(mp.mp.dps - 15)); x = mp.mpf(x1); u = mp.mpf(u1); up = mp.mpf(up1); ls = mp.mpf(0)
    sgn = 1 if x2 > x1 else -1; h = sgn * mp.mpf('0.05')
    while (x2 - x) * sgn > 0:
        if (x + h - x2) * sgn > 0: h = x2 - x
        a = _lin_step(M, alpha, lam, E, x, u, up, h, N)
        b1 = _lin_step(M, alpha, lam, E, x, u, up, h / 2, N); b = _lin_step(M, alpha, lam, E, x + h / 2, b1[0], b1[1], h / 2, N)
        sc = abs(b[0]) + abs(b[1]) + mp.mpf(10) ** -30
        err = (abs(a[0] - b[0]) + abs(a[1] - b[1])) / sc
        if err < tol:
            x = x + h; u, up = b
            r = abs(u) + abs(up)
            if r > 0: u, up = u / r, up / r; ls += mp.log(r)
            if err < tol / 1000: h *= mp.mpf('1.5')
        else: h /= 2
    return u, up, ls
def D_lin(M, alpha, l, E, x0=mp.mpf('0.5')):
    """D(E) = W[y, psi_+] for any real E (linear inward integration of y)"""
    E = mp.mpf(E); lam = l * (l + 1); e = -E
    xm = xmax_for(M, alpha, lam, e); L1, q1 = asym_y(M, alpha, lam, e, xm)
    u, up, ls = lin_integrate(M, alpha, lam, E, xm, mp.mpf(1), -q1, x0)
    ps, psd = frob(M, alpha, l, e, x0)
    return (u * psd - up * ps) * mp.exp(ls + L1)
def shoot(M, alpha, l, E, x0=mp.mpf('0.5'), xfar=mp.mpf(3)):
    """psi_+ integrated OUTWARD to xfar (sign only matters: eigenvalues where it changes sign)"""
    E = mp.mpf(E); lam = l * (l + 1); ps, psd = frob(M, alpha, l, -E, x0)
    u, up, ls = lin_integrate(M, alpha, lam, E, x0, ps, psd, xfar)
    return u
