"""RED2 (B2, second form): tensor product of two COULOMB problems (BLZ at M = -1/2, beta^2 = 2, c = -2: the Sol 3 family's own variable)
L_pm: -psi'' + (1/x + lam_pm/x^2) psi = E psi, i.e. (theta^2 - theta - lam_pm) psi = (x - E x^2) psi.  Tensor product vs Sol 3's k = -2 member S."""
import sys; sys.dont_write_bytecode = True
import sympy as sp
x, D = sp.symbols('x D')
lp, lm, E, l0, l1, C = sp.symbols('lam_p lam_m E l0 l1 C')
def th(f): return sp.expand(x * sp.diff(f, x))
a1 = 1
a0p = lp + x - E * x**2; a0m = lm + x - E * x**2
def theta_vec(v):
    c1, c2, c3, c4 = v
    out = [th(c1), th(c2), th(c3), th(c4)]
    out[1] += c1; out[2] += c1
    out[0] += c2 * a0p; out[1] += c2 * a1; out[3] += c2
    out[0] += c3 * a0m; out[2] += c3 * a1; out[3] += c3
    out[2] += c4 * a0p; out[1] += c4 * a0m; out[3] += c4 * 2 * a1
    return [sp.expand(o) for o in out]
vecs = [[sp.Integer(1), 0, 0, 0]]
for _ in range(4): vecs.append(theta_vec(vecs[-1]))
ns = sp.Matrix(vecs).T.nullspace(); assert len(ns) == 1
rel = [sp.factor(r / ns[0][4]) for r in ns[0]]
den = sp.lcm([sp.denom(sp.together(r)) for r in rel])
T = sp.expand(sum(sp.cancel(r * den) * D**k for k, r in enumerate(rel)))
print('tensor product of two Coulomb problems, denominator', den); Tp = sp.Poly(T, x, D)
for (i, k), c in sorted(Tp.terms()): print(f'   x^{i} D^{k}: {sp.factor(c)}')
S = sp.expand((((D - sp.Rational(1, 2))**2 - l0**2) * ((D - sp.Rational(1, 2))**2 - l1**2)) - C * x * D * (D - 1) + C * E * x**2 * D * (D - 1))
# try: gauge g, scale mu, lam_pm free: x^(-g) T x^g = T(D -> D + g); x -> mu x scales x^i by mu^i
g, mu = sp.symbols('g mu')
Tg = sp.expand(sum(c * mu**i * x**i * (D + g)**k for (i, k), c in Tp.terms()))
diff = sp.Poly(sp.expand(Tg - S), x, D)
eqs = [c for c in diff.coeffs()]
sol = sp.solve(eqs, [g, mu, lp, lm, C], dict=True)
print('exact solve T(gauged, scaled) == S over (g, mu, lam_p, lam_m, C):', sol)
for s_ in sol:
    print('  check residual zero:', sp.expand(Tg.subs(s_) - S) == 0)
