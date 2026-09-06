"""WKB (large-kappa) expansion of the Wronskian for LVZ-shaped ODEs

    Psi'' = [ kappa^2 (1+u)^m + V(u) ] Psi,   u = e^z,   V a finite sum of
    terms  c * u^p / (1+u)^q  (p >= 1, q >= p),

following hep-th/0404195 sec. 5: Psi = exp(kappa S_{-1} + S_0 + S_1/kappa + ...),
Riccati recursion for y_k = S_k', and S_k(-inf) := int_{-inf}^{inf} y_k dz for
odd k, evaluated termwise by Beta functions (exact rational functions of m).
Validated (2026-09-02, certificate cyl_wkb_validation) against eqs. (65)-(68) of the paper for the paperclip
ODE (62) in its CORRECTED form (the (n+2)Q^2/2 + 1/4 term multiplies
u/(1+u)^2, as in LVZ (90)).
"""
import sympy as sp

u, m, kap = sp.symbols('u m kappa')

def riccati(V, kmax, m_sym=m):
    """y_k = (1+u)^{-k m/2} r_k(u), k = -1..kmax; returns {k: r_k} with r_k rational in u over Q(m).
    d/dz [(1+u)^{-km/2} r] = (1+u)^{-km/2} u [r' - (k m/2) r/(1+u)]; all terms of the order-kappa^{-k}
    equation carry (1+u)^{-km/2}, so dividing by 2 y_{-1} = 2 (1+u)^{m/2} keeps the bookkeeping exact."""
    r = {-1: sp.Integer(1)}
    def dz(k, rk):
        return u * (sp.diff(rk, u) - (k * m_sym / 2) * rk / (1 + u))
    r[0] = sp.cancel(-dz(-1, r[-1]) / 2)                        # 2 y_{-1} y_0 + y_{-1}' = 0
    r[1] = sp.cancel((V - r[0] ** 2 - dz(0, r[0])) / 2)          # 2 y_{-1} y_1 + y_0^2 + y_0' = V
    for k in range(1, kmax):
        s = sum(r[i] * r[k - i] for i in range(0, k + 1))
        r[k + 1] = sp.cancel(-(dz(k, r[k]) + s) / 2)
    return r

def beta_int(p, r):
    """int_0^inf u^{p-1} (1+u)^{-(p+r)} du = Gamma(p)Gamma(r)/Gamma(p+r), p a positive integer -> rational in r"""
    p = int(p)
    den = sp.Integer(1)
    for i in range(p):
        den *= (r + i)
    return sp.factorial(p - 1) / den

def integrate_dz(rk, k, m_sym=m):
    """int_{-inf}^{inf} (1+u)^{-k m/2} r_k(u) dz = int_0^inf (1+u)^{-k m/2} r_k(u) du/u,
    r_k = N(u)/(1+u)^b with N a polynomial: termwise Beta functions, exact in m."""
    num, den = sp.fraction(sp.cancel(sp.together(rk)))
    den = sp.Poly(den, u)
    # den must be c (1+u)^b
    b = 0
    while True:
        q, rem = sp.div(den, sp.Poly(1 + u, u))
        if rem.is_zero and den.degree() > 0:
            den = q; b += 1
        else:
            break
    if den.degree() != 0:
        raise ValueError(("denominator not a power of (1+u)", den))
    cden = den.coeffs()[0]
    N = sp.Poly(sp.expand(num), u)
    total = sp.Integer(0)
    for (a,), c in zip(N.monoms(), N.coeffs()):
        # term c u^a (1+u)^{-b - km/2} du/u = c u^{a-1} (1+u)^{-(a + rr)} with rr = b + k m/2 - a
        if a < 1:
            raise ValueError(("u^0 term would diverge at z -> -inf", rk))
        total += c * beta_int(a, b + k * m_sym / 2 - a)
    return sp.factor(sp.cancel(total / cden))

def wronskian_S(V, kmax, m_sym=m):
    """S_k(-inf) for odd k <= kmax (article eq. (64) normalization: S_1 = int y_1 dz)."""
    r = riccati(V, kmax, m_sym)
    return {k: integrate_dz(r[k], k, m_sym) for k in range(1, kmax + 1, 2)}

if __name__ == "__main__":
    n, P, Q = sp.symbols('n P Q')
    # corrected paperclip ODE (62): V = n P^2/2 * u/(1+u) + ((n+2) Q^2/2 + 1/4) * u/(1+u)^2
    V = n * P ** 2 / 2 * u / (1 + u) + ((n + 2) * Q ** 2 / 2 + sp.Rational(1, 4)) * u / (1 + u) ** 2
    S = wronskian_S(V, 7, n)
    art = {
        1: P**2/2 + Q**2/2 + sp.Rational(1, 8),
        3: (-n/(24*(2+3*n))*P**4 - (2+n)/(24*(4+3*n))*Q**4 - n*(2+n)/(4*(2+3*n))*P**2*Q**2
            - n*(3+2*n)/(24*(2+3*n)*(4+3*n))*P**2 - (2+n)*(1+2*n)/(24*(2+3*n)*(4+3*n))*Q**2
            - (2+6*n+3*n**2)/(192*(2+3*n)*(4+3*n))),
        5: (n**2/(40*(2+5*n)*(4+5*n))*P**6 + (2+n)**2/(40*(6+5*n)*(8+5*n))*Q**6
            + 3*n**2*(2+n)/(8*(2+5*n)*(4+5*n)*(6+5*n))*P**4*Q**2 + 3*n*(2+n)**2/(8*(4+5*n)*(6+5*n)*(8+5*n))*P**2*Q**4
            + n**2*(4+3*n)/(16*(2+5*n)*(4+5*n)*(6+5*n))*P**4 + (2+n)**2*(2+3*n)/(16*(4+5*n)*(6+5*n)*(8+5*n))*Q**4
            + 3*n*(2+n)*(2+10*n+5*n**2)/(8*(2+5*n)*(4+5*n)*(6+5*n)*(8+5*n))*P**2*Q**2
            + n*(30+225*n+250*n**2+76*n**3)/(160*(2+5*n)*(4+5*n)*(6+5*n)*(8+5*n))*P**2
            + (2+n)*(28+137*n+206*n**2+76*n**3)/(160*(2+5*n)*(4+5*n)*(6+5*n)*(8+5*n))*Q**2
            + (24+140*n+290*n**2+220*n**3+55*n**4)/(640*(2+5*n)*(4+5*n)*(6+5*n)*(8+5*n))),
        7: (-15*n**3/(448*(2+7*n)*(4+7*n)*(6+7*n))*P**8 - 15*(2+n)**3/(448*(8+7*n)*(10+7*n)*(12+7*n))*Q**8
            - 15*n**3*(2+n)/(16*(2+7*n)*(4+7*n)*(6+7*n)*(8+7*n))*P**6*Q**2
            - 15*n*(2+n)**2/(16*(6+7*n)*(8+7*n)*(10+7*n)*(12+7*n))*P**2*Q**6
            - 45*n**2*(2+n)**2/(32*(4+7*n)*(6+7*n)*(8+7*n)*(10+7*n))*P**4*Q**4
            - 5*n**3*(5+4*n)/(32*(2+7*n)*(4+7*n)*(6+7*n)*(8+7*n))*P**6
            - 5*(2+n)**3*(3+4*n)/(32*(6+7*n)*(8+7*n)*(10+7*n)*(12+7*n))*Q**6
            - 15*n**2*(2+n)*(6+49*n+28*n**2)/(32*(2+7*n)*(4+7*n)*(6+7*n)*(8+7*n)*(10+7*n))*P**4*Q**2
            - 15*n*(2+n)**2*(20+63*n+28*n**2)/(32*(4+7*n)*(6+7*n)*(8+7*n)*(10+7*n)*(12+7*n))*P**2*Q**4
            - n**2*(122+1407*n+1792*n**2+612*n**3)/(128*(2+7*n)*(4+7*n)*(6+7*n)*(8+7*n)*(10+7*n))*P**4
            - (2+n)**2*(420+1583*n+1880*n**2+612*n**3)/(128*(4+7*n)*(6+7*n)*(8+7*n)*(10+7*n)*(12+7*n))*Q**4
            - 21*n*(2+n)*(40+338*n+985*n**2+816*n**3+204*n**4)/(64*(2+7*n)*(4+7*n)*(6+7*n)*(8+7*n)*(10+7*n)*(12+7*n))*P**2*Q**2
            - n*(2520+26166*n+102459*n**2+136612*n**3+74676*n**4+14552*n**5)/(896*(2+7*n)*(4+7*n)*(6+7*n)*(8+7*n)*(10+7*n)*(12+7*n))*P**2
            - (2+n)*(3720+30202*n+89149*n**2+121284*n**3+70844*n**4+14552*n**5)/(896*(2+7*n)*(4+7*n)*(6+7*n)*(8+7*n)*(10+7*n)*(12+7*n))*Q**2
            - (12240+110628*n+386834*n**2+680568*n**3+606452*n**4+261786*n**5+43631*n**6)/(14336*(2+7*n)*(4+7*n)*(6+7*n)*(8+7*n)*(10+7*n)*(12+7*n))),
    }
    for k in (1, 3, 5, 7):
        d = sp.simplify(sp.expand(S[k] - art[k]))
        print(f"S_{k}: {'MATCHES eq. (' + str(64 + k) + ')' if d == 0 else 'DIFFERS: ' + str(sp.factor(d))}")
