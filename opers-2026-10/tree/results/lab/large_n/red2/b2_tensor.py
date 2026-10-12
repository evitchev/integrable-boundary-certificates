"""RED2 (B2): tensor (Hadamard) product of two radial oscillators vs Sol 3's k = -2 member (item 299 / item 298 at k = -2).
Oscillators in x = y^2:  L_pm = 4 theta^2 - 2 theta - lam_pm + Ep x - x^2   (i.e. -psi_yy + (y^2 + lam/y^2) psi = Ep psi),  theta = x d/dx.
Sol 3 member: S = P4(theta) - C x theta(theta - 1) + C E x^2 theta(theta - 1),  P4 = ((theta-1/2)^2 - l0^2)((theta-1/2)^2 - l1^2).
Exact sympy; operators as polynomials in a commuting symbol D standing for theta with coefficients written to the LEFT (x^j D^k)."""
import sys; sys.dont_write_bytecode = True
import sympy as sp
x, D = sp.symbols('x D')
lp, lm, Ep, l0, l1, C, E, g, mu = sp.symbols('lam_p lam_m Ep l0 l1 C E g mu')
def th(f): return sp.expand(x * sp.diff(f, x))
# theta^2 u = a1 theta u + a0 u
a1 = sp.Rational(1, 2)
a0p = (lp - Ep * x + x**2) / 4; a0m = (lm - Ep * x + x**2) / 4
def theta_vec(v):
    """v = coefficients on (uv, u'v, uv', u'v') -> coefficients of theta(phi-part)"""
    c1, c2, c3, c4 = v
    out = [th(c1), th(c2), th(c3), th(c4)]
    # theta e1 = e2 + e3
    out[1] += c1; out[2] += c1
    # theta e2 = a0p e1 + a1 e2 + e4
    out[0] += c2 * a0p; out[1] += c2 * a1; out[3] += c2
    # theta e3 = a0m e1 + a1 e3 + e4
    out[0] += c3 * a0m; out[2] += c3 * a1; out[3] += c3
    # theta e4 = a0p e3 + a0m e2 + 2 a1 e4
    out[2] += c4 * a0p; out[1] += c4 * a0m; out[3] += c4 * 2 * a1
    return [sp.expand(o) for o in out]
vecs = [[sp.Integer(1), 0, 0, 0]]
for _ in range(4): vecs.append(theta_vec(vecs[-1]))
Mx = sp.Matrix(vecs).T            # 4 x 5
ns = Mx.nullspace()
assert len(ns) == 1
rel = ns[0] / ns[0][4]
rel = [sp.factor(r) for r in rel]
den = sp.lcm([sp.denom(sp.together(r)) for r in rel])
T = sp.expand(sum(sp.cancel(r * den) * D**k for k, r in enumerate(rel)))
print('tensor product T (coeffs left of D = theta), cleared denominator', den, ':')
Tp = sp.Poly(T, x, D)
for (i, k), c in sorted(Tp.terms()): print(f'   x^{i} D^{k}: {sp.factor(c)}')
S = sp.expand((((D - sp.Rational(1, 2))**2 - l0**2) * ((D - sp.Rational(1, 2))**2 - l1**2)) - C * x * D * (D - 1) + C * E * x**2 * D * (D - 1))
print('Sol 3 member S:'); Sp_ = sp.Poly(S, x, D)
for (i, k), c in sorted(Sp_.terms()): print(f'   x^{i} D^{k}: {sp.factor(c)}')
