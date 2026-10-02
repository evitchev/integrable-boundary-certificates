"""PROOF1 (1): machine check of every algebraic identity in the Mellin proof of
   [P(T) - x^(c/2) L(T) x^(c/2) - E x^n] phi = 0   <=>   P(T) prod_r G_a(T - r) psi = x^n (x^(n/a) + E) psi,
   T = s - theta, sigma = n/a (= n M), c = n + sigma, G_a(T) = sigma^a Gamma(T/sigma + (1+a)/2) / Gamma(T/sigma + (1-a)/2), L(T) = prod_r (T - r).
Multiplier: phihat(nu) = m(nu) psihat(nu),  m(nu) = prod_r sigma^(nu/sigma) / Gamma(gamma_r - nu/sigma),  gamma_r = (s - r)/sigma + (1-a)/2   (ENTIRE in nu).
Identities (each symbolic in a, sigma, s, r, nu, with Gamma kept symbolic):
 (I0) T x^(c/2) = x^(c/2) (T - c/2), hence x^(c/2) L(T) x^(c/2) x^nu = L(s - nu - c/2) x^(nu + c).
 (I1) m(mu)/m(mu - sigma) = prod_r (sigma gamma_r - mu)     [one root factor per r]
 (I2) m(nu)/m(nu - n) = prod_r G_a(s - nu - r)              [n = a sigma; for REAL a, via Gamma]
 (I3) L(s - nu + c/2) m(nu - c)/m(nu - n) = 1
 (I4) the ODE recurrence times -1/m(nu - n) equals the Gamma-form recurrence:  P G psihat(nu) - psihat(nu - c) - E psihat(nu - n) = 0."""
import sys; sys.dont_write_bytecode = True
import sympy as sp
a, sig, s, nu, mu, E, x = sp.symbols('a sigma s nu mu E x', positive=True)
R = sp.symbols('r1:4')                     # up to three roots of L (deg L = 1, 2, 3 cover Sols 1, 2, 3 and beyond)
ok = True
def report(name, val):
    global ok
    print(f'{name}: {val}'); ok &= bool(val)
# (I0)
f = sp.Function('f'); th = lambda g: x * sp.diff(g, x)
cc = sp.Symbol('c', positive=True)
lhs = (s * x**(cc / 2) * f(x) - th(x**(cc / 2) * f(x)))                       # T (x^(c/2) f)
rhs = x**(cc / 2) * ((s - cc / 2) * f(x) - th(f(x)))                            # x^(c/2) (T - c/2) f
report('(I0) T x^(c/2) = x^(c/2)(T - c/2)', sp.simplify(lhs - rhs) == 0)
for deg in (1, 2, 3):
    roots = R[:deg]
    gam = {r: (s - r) / sig + (1 - a) / 2 for r in roots}
    m = lambda v: sp.Mul(*[sig**(v / sig) / sp.gamma(gam[r] - v / sig) for r in roots])
    G = lambda T: sp.Mul(*[sig**a * sp.gamma((T - r) / sig + (1 + a) / 2) / sp.gamma((T - r) / sig + (1 - a) / 2) for r in roots])
    L = lambda T: sp.Mul(*[(T - r) for r in roots])
    n = a * sig; c = n + sig
    i1 = sp.simplify(sp.expand_func(sp.powsimp(sp.gammasimp(m(mu) / m(mu - sig) - sp.Mul(*[(sig * gam[r] - mu) for r in roots])), force=True)))
    report(f'(I1) deg L = {deg}', i1 == 0)
    i2 = sp.simplify(sp.expand_func(sp.powsimp(sp.gammasimp(m(nu) / m(nu - n) / G(s - nu)), force=True)))
    report(f'(I2) deg L = {deg} (ratio == 1)', sp.simplify(sp.powsimp(i2, force=True)) == 1)
    i3 = sp.simplify(sp.expand_func(sp.powsimp(sp.gammasimp(L(s - nu + c / 2) * m(nu - c) / m(nu - n)), force=True)))
    report(f'(I3) deg L = {deg}', sp.simplify(sp.powsimp(i3, force=True)) == 1)
    # (I4): coefficient structure with symbolic P and psihat
    P = sp.Function('P'); ph = sp.Function('psihat')
    ode = P(s - nu) * m(nu) * ph(nu) - L(s - nu + c / 2) * m(nu - c) * ph(nu - c) - E * m(nu - n) * ph(nu - n)
    gam_form = P(s - nu) * G(s - nu) * ph(nu) - ph(nu - c) - E * ph(nu - n)
    i4 = sp.simplify(sp.expand_func(sp.powsimp(sp.gammasimp(sp.expand(ode / m(nu - n)) - gam_form), force=True)))
    report(f'(I4) deg L = {deg}: ODE recurrence / m(nu - n) == Gamma recurrence', i4 == 0)
# NEGATIVE CONTROLS: the non-symmetric ordering x^c L(T) (VIR8's tamper) breaks (I3); a wrong gamma_r shift breaks (I2)
roots = R[:1]; n = a * sig; c = n + sig
gam = {r: (s - r) / sig + (1 - a) / 2 for r in roots}
m = lambda v: sp.Mul(*[sig**(v / sig) / sp.gamma(gam[r] - v / sig) for r in roots])
bad3 = sp.simplify(sp.expand_func(sp.powsimp(sp.gammasimp((s - nu + c - roots[0]) * m(nu - c) / m(nu - n)), force=True)))     # x^c L(T): L evaluated at s - nu + c instead of s - nu + c/2
print('NEGATIVE CONTROL x^c L(T) ordering: (I3) ratio =', sp.simplify(sp.powsimp(bad3, force=True)), '(must NOT be 1)')
mb = lambda v: sig**(v / sig) / sp.gamma(gam[roots[0]] + sp.Rational(1, 3) - v / sig)
Gb = sig**a * sp.gamma((s - nu - roots[0]) / sig + (1 + a) / 2) / sp.gamma((s - nu - roots[0]) / sig + (1 - a) / 2)
bad2 = sp.simplify(sp.powsimp(sp.gammasimp(mb(nu) / mb(nu - n) / Gb), force=True))
print('NEGATIVE CONTROL gamma_r + 1/3: (I2) ratio =', bad2, '(must NOT be 1)')
print('ALL IDENTITIES', 'PASS' if ok else 'FAIL'); sys.exit(0 if ok else 1)
