"""WKB for a Heun-type extension at m = -1:
   Psi'' = [kappa^2/(1+u) + V] Psi,  V = A u/(1+u) + B u/(1+u)^2 + D u/(u+v) + E u/(u+v)^2,
regular singular points u = 0, -1, -v, infinity (v > 0 keeps the new point off the path).
Riccati: y_k = (1+u)^{k/2} r_k(u).  S_k = int_0^inf (1+u)^{k/2} r_k du/u, defined -- as for the
paperclip at n = -1 -- by analytic continuation in the (1+u)-exponent: each term
u^{a-1} (1+u)^{-alpha} (u+v)^{-d} integrates to v^{-d} B(a, alpha+d-a) 2F1(d, a; alpha+d; 1-1/v),
which for half-integer alpha is rational in v plus ONE transcendental
T = asin(sqrt(z))/sqrt(z(1-z)), z = 1 - 1/v.  Rationality of S_k in the conic parameter
requires the coefficient of T to vanish -- UNDER a T-free normalization alpha_s(t); without that assumption
the condition is proportionality of the T-part and the rational part (ansatz note sec. 15)."""
import sympy as sp, mpmath as mp
u, v, z, T = sp.symbols('u v z T', positive=True)
A, B, D, E = sp.symbols('A B D E')

def riccati(V, kmax):
    r = {-1: sp.Integer(1)}
    dz = lambda k, rk: u * (sp.diff(rk, u) + sp.Rational(k, 2) * rk / (1 + u))     # m = -1
    r[0] = sp.cancel(-dz(-1, r[-1]) / 2)
    r[1] = sp.cancel((V - r[0] ** 2 - dz(0, r[0])) / 2)
    for k in range(1, kmax):
        r[k + 1] = sp.cancel(-(dz(k, r[k]) + sum(r[i] * r[k - i] for i in range(0, k + 1))) / 2)
    return r

_cache = {}
def euler(a, alpha, d):
    """I(a, alpha, d) = int_0^inf u^{a-1} (1+u)^{-alpha} (u+v)^{-d} du (continued in alpha), as
    rational(v) + T * rational(v).  Reduced to a = 1 by u^{a-1} = u^{a-2}[(u+v) - v]:
    I(a, alpha, d) = I(a-1, alpha, d-1) - v I(a-1, alpha, d); d = 0 is a pure Beta function.
    The pole position v stays SYMBOLIC here (the formal rational + T split needs the symbolic
    argument: at v = 2, 4/3, 1/2 a numeric argument evaluates the inverse trigonometric term to
    pi or a logarithm before it can be identified with T -- Codex round 27); a numeric position is
    substituted afterwards by S_k, where the speed-up actually lives (the denominator division)."""
    key = (a, alpha, d)
    if key in _cache:
        return _cache[key]
    if d == 0:
        val = sp.gamma(a) * sp.gamma(alpha - a) / sp.gamma(alpha)
    elif a == 1:
        pref = v ** (-d) * sp.gamma(alpha + d - 1) / sp.gamma(alpha + d)
        f = sp.hyperexpand(sp.hyper([d, 1], [alpha + d], z))
        f = f.subs(sp.asin(sp.sqrt(z)), T * sp.sqrt(z) * sp.sqrt(1 - z))
        val = pref * f.subs(z, 1 - 1 / v)
    else:
        val = euler(a - 1, alpha, d - 1) - v * euler(a - 1, alpha, d)
    val = sp.expand(val)
    _cache[key] = val
    return val

def S_k(rk, k, vval=None):
    """S_k = int y_k dz as a polynomial in T.  vval: the pole position -- the symbol v (default: the
    certified symbolic path) or an exact rational, in which case the denominator is divided by
    (u + vval) with the number and the symbolic Euler integrals are specialized AFTER their T split
    (identical results, validated against the pinned general-S files at v = 3 and v = 2; the
    specialization makes spin 5 a five-minute computation instead of hours)."""
    if vval is None:
        vv = v
    else:
        vv = sp.nsimplify(vval)
        if not (vv.is_Rational and vv != 0):
            raise ValueError("vval must be an exact nonzero rational")
    rk = sp.together(rk)
    num, den = sp.fraction(rk)
    b = d = 0
    dpoly = sp.Poly(den, u)
    while True:
        q_, rem = sp.div(dpoly, sp.Poly(1 + u, u))
        if rem.is_zero and dpoly.degree() > 0: dpoly = q_; b += 1
        else: break
    while True:
        q_, rem = sp.div(dpoly, sp.Poly(u + vv, u))
        if rem.is_zero and dpoly.degree() > 0: dpoly = q_; d += 1
        else: break
    if dpoly.degree() != 0:
        raise ValueError(("denominator", den))
    const = dpoly.coeffs()[0]
    npoly = sp.Poly(sp.expand(num), u)
    total = sp.Integer(0)
    for (a,), cf in zip(npoly.monoms(), npoly.coeffs()):
        if a < 1:
            raise ValueError("u^0 term (divergent at z -> -inf)")
        term = euler(a, b - sp.Rational(k, 2), d)                # (1+u)^{k/2 - b} = (1+u)^{-alpha}
        total += cf * (term if vval is None else term.subs(v, vv))
    total = sp.expand(total / const)
    return sp.Poly(total, T)

def selfcheck():
    """the Euler formula against numerical integration where it converges"""
    a, alpha, d = 2, sp.Rational(9, 2), 1
    vv = sp.Rational(7, 3)
    sym = euler(a, alpha, d).subs(v, vv)
    Tval = mp.asin(mp.sqrt(1 - 3/mp.mpf(7))) / mp.sqrt((1 - 3/mp.mpf(7)) * (3/mp.mpf(7)))
    numeric = mp.quad(lambda x: x**(a-1) * (1+x)**(-float(alpha)) * (x + 7/mp.mpf(3))**(-d), [0, mp.inf])
    return abs(float(sym.subs(T, sp.Float(str(Tval), 30))) - float(numeric)) < 1e-10

if __name__ == "__main__":
    import time
    t0 = time.time(); print("Euler/2F1 formula self-check vs numerical integration:", selfcheck(), f"({time.time()-t0:.1f}s)", flush=True)
    V = A * u / (1 + u) + B * u / (1 + u) ** 2 + D * u / (u + v) + E * u / (u + v) ** 2
    t0 = time.time(); r = riccati(V, 3); print(f"riccati k<=3: {time.time()-t0:.1f}s", flush=True)
    for k in (1, 3):
        t0 = time.time(); Sk = S_k(r[k], k); print(f"S_{k} integrated: {time.time()-t0:.1f}s", flush=True)
        for (i,), cf in zip(Sk.monoms(), Sk.coeffs()):
            print(f"S_{k}  T^{i}:  {sp.factor(cf)}", flush=True)
