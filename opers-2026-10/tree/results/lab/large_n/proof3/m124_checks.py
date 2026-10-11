"""PROOF3 M1, M2, M4 (t = 2, k = a = 2, sigma = 3, gamma = 1/3, s = 5/2, n = 6, c = 9).
M1a Hankel realisation: (1/2 pi i) Int_H e^t t^(-gamma) (sigma t)^(nu/sigma) dt == m(nu) = sigma^(nu/sigma)/Gamma(gamma - nu/sigma) (numerical quadrature, complex nu).
M1b integration-by-parts relation K_(g-1) psi = K_g[(g - 1 - theta/sigma) psi] on monomials (symbolic Gamma recurrence) and on a NON-monomial psi = x^nu0 exp(-x^2)
    by Hankel quadrature at one x (numerical).
M1c intertwining L3 K = x^n K x^-n L6 on monomials, symbolic in (l0, l1, E, mu) (from M1b's relations).
M2  frozen exponents of G_a(T - r) and their shift lattice (n, c multiples of sigma) are zeros of m: symbolic, a = 1..4.
M4  Laplace prediction kappa = 3^(1/3)/sqrt(3 pi) = 3^(-1/6) pi^(-1/2) vs M3's 1/kappa, and beta3 = 2 beta6 + 1/2 with NUM2's own asymptotic solver
    (beta = x^0 coefficient of S) at E = 0 and E = 37 (E-independence)."""
import mpmath as mp, sympy as sp
from fractions import Fraction as Fr
mp.mp.dps = 40
sig, gam = mp.mpf(3), mp.mpf(1) / 3
ok = True
# M1a
def hankel(fun, R=mp.mpf(1)):
    # fun(logt): integrand as a function of log t (explicit branch arg t in [-pi, pi]); dt = t dlogt
    # contour: from -inf below the cut (arg -pi) to |t| = R, circle counterclockwise, back to -inf above (arg +pi)
    g = lambda lt: fun(lt) * mp.exp(lt)
    lo = -mp.quad(lambda u: g(mp.log(u) - 1j * mp.pi) / u, [R, mp.inf])   # t = u e^{-i pi}, u: inf -> R
    ci = mp.quad(lambda ph: g(mp.log(R) + 1j * ph) * 1j, [-mp.pi, mp.pi])
    hi = mp.quad(lambda u: g(mp.log(u) + 1j * mp.pi) / u, [R, mp.inf])
    return (lo + ci + hi) / (2j * mp.pi)
for nu in (mp.mpf('0.7'), mp.mpc('2.3', '0.4'), mp.mpf('-1.2')):
    val = hankel(lambda lt: mp.exp(mp.exp(lt) - gam * lt + (nu / sig) * (mp.log(sig) + lt)))
    ex = sig ** (nu / sig) / mp.gamma(gam - nu / sig)
    d = abs(val - ex) / abs(ex); print(f'M1a nu = {nu}: Hankel vs m(nu): rel {mp.nstr(d, 3)}'); ok &= d < mp.mpf(10) ** -25
# M1b symbolic: on x^nu, K_g x^nu = sigma^(nu/sigma)/Gamma(g - nu/sigma) x^nu
nu_, g_, S_ = sp.symbols('nu g sigma', positive=True)
Kg = lambda g, v: S_**(v / S_) / sp.gamma(g - v / S_)
lhs = Kg(g_ - 1, nu_); rhs = Kg(g_, nu_) * (g_ - 1 - nu_ / S_)
print('M1b symbolic (monomials):', sp.simplify(sp.gammasimp(lhs / rhs)) == 1); ok &= sp.simplify(sp.gammasimp(lhs / rhs)) == 1
# M1b numeric on psi = x^nu0 exp(-x^2) (theta psi = (nu0 - 2 x^2) psi), at x = 0.7, g = gamma
x0, nu0 = mp.mpf('0.7'), mp.mpf('0.4')
psi = lambda lX: mp.exp(nu0 * lX - mp.exp(2 * lX))   # psi as a function of log X
thpsi = lambda lX: (nu0 - 2 * mp.exp(2 * lX)) * psi(lX)
def K(g, f):
    return hankel(lambda lt: mp.exp(mp.exp(lt) - g * lt) * f(mp.log(x0) + (mp.log(sig) + lt) / sig))
L = K(gam - 1, psi); Rr = K(gam, lambda X: (gam - 1) * psi(X) - thpsi(X) / sig)
d = abs(L - Rr) / abs(L); print(f'M1b numeric, psi = x^0.4 exp(-x^2): K_(g-1) psi vs K_g[(g-1-theta/sigma) psi]: rel {mp.nstr(d, 3)}'); ok &= d < mp.mpf(10) ** -20
# M1c symbolic intertwining on monomials, k = 2: L3 = P(s-th) - x^c (s - c/2 - th) - E x^n ;  L6 = P(s-th)(th-1)(th-4) - x^c - E x^n (same E)
l0, l1, E, mu = sp.symbols('l0 l1 E mu')
s, n, c, sg, gm = sp.Rational(5, 2), 6, 9, 3, sp.Rational(1, 3)
m = lambda v: sp.Integer(sg)**(v / sg) / sp.gamma(gm - v / sg)
P = lambda T: (T**2 - l0**2) * (T**2 - l1**2)
# L3 K x^mu = m(mu) [P(s-mu) x^mu - (s - c/2 - mu) x^(mu+c) - E x^(mu+n)]
L3K = {0: m(mu) * P(s - mu), c: -m(mu) * (s - sp.Rational(c, 2) - mu), n: -E * m(mu)}
# x^n K x^-n L6 x^mu: L6 x^mu = P(s-mu)(mu-1)(mu-4) x^mu - x^(mu+c) - E x^(mu+n); then x^n K x^-n x^nu = m(nu - n) x^nu
XKX = {0: m(mu - n) * P(s - mu) * (mu - 1) * (mu - 4), c: -m(mu + c - n), n: -E * m(mu)}
okc = all(sp.simplify(sp.gammasimp(sp.expand_func(L3K[h] - XKX[h]))) == 0 for h in (0, c, n))
print('M1c intertwining L3 K = x^n K x^-n L6 on monomials (symbolic):', okc); ok &= okc
# M1c tamper (must FAIL): gamma -> 1/3 + 1/7 in the multiplier
mt = lambda v: sp.Integer(sg)**(v / sg) / sp.gamma(gm + sp.Rational(1, 7) - v / sg)
L3Kt = {0: mt(mu) * P(s - mu), c: -mt(mu) * (s - sp.Rational(c, 2) - mu), n: -E * mt(mu)}
XKXt = {0: mt(mu - n) * P(s - mu) * (mu - 1) * (mu - 4), c: -mt(mu + c - n), n: -E * mt(mu)}
okt = all(sp.simplify(sp.gammasimp(sp.expand_func(L3Kt[h] - XKXt[h]))) == 0 for h in (0, c, n))
print('M1c tamper (gamma + 1/7): intertwining holds?', okt, '(expected False)'); ok &= not okt
# M2 frozen zeros, symbolic a = 1..4, one root r
a_, r_, s_s = sp.symbols('a r s')
for A in range(1, 5):
    sgm = sp.Symbol('sigma', positive=True); gmr = (s_s - r_) / sgm + sp.Rational(1 - A, 2)
    frozen = [s_s - r_ + sgm * (sp.Rational(A - 1, 2) - i) for i in range(A)]       # roots of G_A(T - r) in theta, T = s - theta
    shifts = [A * sgm * p + (A + 1) * sgm * q for p in range(3) for q in range(3)]  # n = A sigma, c = (A+1) sigma
    good = all(((gmr - (f + h) / sgm).simplify()).is_integer and (gmr - (f + h) / sgm).simplify() <= 0 for f in frozen for h in shifts)
    print(f'M2 a = {A}: frozen exponents + shift lattice are zeros of m (gamma_r - nu/sigma in Z<=0): {good}'); ok &= good
# M2 tamper (must FAIL): frozen exponents shifted by sigma/2 are not zeros of m
for A in range(1, 5):
    sgm = sp.Symbol('sigma', positive=True); gmr = (s_s - r_) / sgm + sp.Rational(1 - A, 2)
    bad = all(((gmr - (s_s - r_ + sgm * (sp.Rational(A - 1, 2) - i) + sgm / 2) / sgm).simplify()).is_integer for i in range(A))
    print(f'M2 tamper a = {A} (frozen + sigma/2): zeros of m? {bad} (expected False)'); ok &= not bad
# M4 Laplace
kap_L = mp.mpf(3) ** (mp.mpf(1) / 3) / mp.sqrt(3 * mp.pi)
kap_M3 = 1 / mp.mpf('2.12860533089645112119295775073')
print(f'M4 Laplace kappa = 3^(1/3)/sqrt(3 pi) = {mp.nstr(kap_L, 25)}; from M3: {mp.nstr(kap_M3, 25)}; closed form 3^(-1/6) pi^(-1/2) = {mp.nstr(3 ** (-mp.mpf(1)/6) / mp.sqrt(mp.pi), 25)}')
ok &= abs(kap_L - kap_M3) / kap_M3 < mp.mpf(10) ** -28
import ops
for Ev in (0, 37):
    op3, th, F = ops.three_term(2, mp.mpf(1) / 3, mp.mpf(3) / 5, Ev); a3, _ = op3.asym_coeffs(40); b3 = a3[int(op3.gamma / op3.delta)]
    op6, ex = ops.sixth(mp.mpf(1) / 3, mp.mpf(3) / 5, -Ev); a6, _ = op6.asym_coeffs(40); b6 = a6[int(op6.gamma / op6.delta)]
    print(f'M4 E = {Ev}: beta3 = {mp.nstr(b3, 20)}, beta6 = {mp.nstr(b6, 20)}, 2 beta6 + 1/2 = {mp.nstr(2 * b6 + mp.mpf(1)/2, 20)}; leading a0: {mp.nstr(op3.a0, 10)}, {mp.nstr(op6.a0, 10)}')
    ok &= abs(b3 - (2 * b6 + mp.mpf(1) / 2)) < mp.mpf(10) ** -30 and abs(b3 + 2) < mp.mpf(10) ** -30
print('M1/M2/M4', 'PASS' if ok else 'FAIL')
import sys; sys.exit(0 if ok else 1)
