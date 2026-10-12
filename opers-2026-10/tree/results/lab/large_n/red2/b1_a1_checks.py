"""RED2 (B1) and (A1) machine checks.
B1: (i) Erdelyi-Kober term identity: Gamma(b+mu)/(Gamma(b)Gamma(mu)) B(k+b, mu) == (b)_k/(b+mu)_k, symbolic in b, mu, k = 0..6;
    (ii) spectral determinant of -psi'' + (y^2 + alpha + l(l+1)/y^2) psi = E psi depends on a = (2l + 3 - E + alpha)/4 only: invariance under
         (alpha, l) -> (alpha - 2 delta, l + delta)  [exact: a unchanged];
    (iii) Stirling coefficients B_(n+1)(a) with a = 1/2 + i (pi + sigma P_X)/sqrt2 (SUPER3b dictionary at M = 1: alpha^2 = -8X, (l+1/2)^2 = -2Y,
         X = -P_X^2, Y = -pi^2 [P_art = i P]): the single member is NOT even in P_X; its P_X-even part == (1/2) sum_sigma, n+1 even, n+1 = 2..10.
A1: (i) local orders of Lambda dz^2 at the three special points, Suzuki(M = -2) vs pillow (t = -1/3, Sol 1: eps = -3/2, n = -5/2);
    (ii) covariance: under x = f(z) the kappa^1 part of the Liouville normal form is f'^2 * alpha U(f) (never removed)."""
import sys; sys.dont_write_bytecode = True
import sympy as sp
b, mu, kk = sp.symbols('b mu k', positive=True)
ok = True
for k in range(7):
    lhs = sp.gamma(b + mu) / (sp.gamma(b) * sp.gamma(mu)) * sp.gamma(k + b) * sp.gamma(mu) / sp.gamma(k + b + mu)
    rhs = sp.rf(b, k) / sp.rf(b + mu, k)
    ok &= sp.simplify(sp.gammasimp(lhs - rhs)) == 0
print('B1(i) Erdelyi-Kober term identity k = 0..6:', ok)
al, l, E, dl = sp.symbols('alpha l E delta')
a = (2 * l + 3 - E + al) / 4
print('B1(ii) a invariant under (alpha, l) -> (alpha - 2 delta, l + delta):', sp.simplify(a.subs({al: al - 2 * dl, l: l + dl}, simultaneous=True) - a) == 0)
PX, pi_ = sp.symbols('P_X pi')
okc = True
for m in range(2, 11, 2):
    single = sp.expand(sp.bernoulli(m, sp.Rational(1, 2) + sp.I * (pi_ + PX) / sp.sqrt(2)))
    other = sp.expand(sp.bernoulli(m, sp.Rational(1, 2) + sp.I * (pi_ - PX) / sp.sqrt(2)))
    even = sp.expand((single + single.subs(PX, -PX)) / 2)
    odd = sp.expand((single - single.subs(PX, -PX)) / 2)
    okc &= sp.expand(even - (single + other) / 2) == 0
    print(f'B1(iii) B_{m}: single member P_X-odd part nonzero: {odd != 0};  P_X-even part == (1/2) sum over the doublet: {sp.expand(even - (single + other) / 2) == 0}')
x, z = sp.symbols('x z', positive=True)
# A1(i): orders of Lambda dw^2 in w with w = x^(2M) (Suzuki) and in w = u/(1+u) (pillow)
M = sp.Integer(-2)
def orders(q, w):
    o0 = sp.limit(sp.log(sp.Abs(q.subs(w, sp.Rational(1, 10**6)))) / sp.log(sp.Rational(1, 10**6)), w, 0) if False else None
    return None
w = sp.Symbol('w', positive=True)
qS = (w - 1) * w**(1 / M - 2)                 # Suzuki: Lambda dx^2 in w = x^(2M), up to a constant
eps, n = sp.Rational(-3, 2), sp.Rational(-5, 2)
qP = w**(-eps - 2) * (1 - w)**(eps - n - 2) * 0 + w**sp.Rational(-1, 2) * (1 - w)**(-1)   # pillow Lambda dz^2 in w (derived: w^(3/2)(1-w) * dw^2/(w(1-w))^2)
def ords(q):
    a0 = sp.Poly(sp.expand(sp.powsimp(q * w**10)), w) if False else None
    at0 = sp.simplify(sp.log(q.subs(w, sp.Symbol('e'))).as_leading_term(sp.Symbol('e'))) if False else None
    o0 = sp.limit(w * sp.diff(sp.log(q), w), w, 0)
    o1 = sp.limit((w - 1) * sp.diff(sp.log(q), w), w, 1)
    oinf = -sp.limit(w * sp.diff(sp.log(q), w), w, sp.oo) - 4
    return (o0, o1, oinf)
print('A1(i) orders of Lambda dw^2 at (w = 0, 1, infinity): Suzuki(M = -2):', ords(qS), '  pillow (t = -1/3):', ords(qP))
# A1(ii) covariance of the kappa-odd part
kap = sp.Symbol('kappa'); f = sp.Function('f')(z)
Lam = sp.Function('Lam'); Uf = sp.Function('U'); Vf = sp.Function('V')
Qx = lambda X: kap**2 * Lam(X) + kap * al * Uf(X) + Vf(X)
Qz = sp.diff(f, z)**2 * Qx(f) - sp.Rational(1, 2) * (sp.diff(f, z, 3) / sp.diff(f, z) - sp.Rational(3, 2) * (sp.diff(f, z, 2) / sp.diff(f, z))**2)
odd = sp.expand((sp.expand(Qz) - sp.expand(Qz).subs(kap, -kap)) / 2)
print('A1(ii) kappa-odd part of the transformed normal form:', sp.simplify(odd), ' == kappa alpha f\'^2 U(f):', sp.simplify(odd - kap * al * sp.diff(f, z)**2 * Uf(f)) == 0)
