"""RED2 (A3): is Sol 2's third-order ODE at t = 1/3 a symmetric square?  k = -1/2, n = 2, c = -2, s = 1/2 (item 303/307):
L3 = T(T^2 - l1^2) - x^(-1)(T^2 - l0^2)x^(-1) - E x^2,  T = s - theta,  x^(-1) f(theta) x^(-1) = x^(-2) f(theta - 1).
In z = log x (theta = d/dz) bring to monic, gauge away the d^2 term -> d^3 + p d + q;  Wilczynski invariant W = q - p'/2 is a weight-3 projective
invariant: W == 0 iff L3 is the symmetric square of a second-order operator, in ANY coordinate (point transformations + gauge).
Controls: (i) sym^2 of a generic Schrodinger operator gives W = 0; (ii) Suzuki(M = -2)'s own sym^2 gives W = 0."""
import sys; sys.dont_write_bytecode = True
import sympy as sp
z, E, l0, l1, al, lam = sp.symbols('z E l0 l1 alpha lam')
y = sp.Function('y')
def apply_theta_poly(coeffs, f):   # sum c_k(z) d^k f
    return sum(c * sp.diff(f, z, k) for k, c in enumerate(coeffs))
s = sp.Rational(1, 2)
th = sp.Symbol('th')
P0 = sp.expand((s - th) * ((s - th)**2 - l1**2))           # x^0 part
P2 = sp.expand(-((s - (th - 1))**2 - l0**2))               # x^-2 part: x^(-2) f(theta - 1) with f(T) = T^2 - l0^2, T = s - theta
def poly_coeffs(P, deg): return [sp.Poly(P, th).coeff_monomial(th**k) for k in range(deg + 1)]
c0 = poly_coeffs(P0, 3); c2 = poly_coeffs(P2, 3)
coeffs = [c0[k] + sp.exp(-2 * z) * c2[k] for k in range(4)]
coeffs[0] += -E * sp.exp(2 * z)
lead = coeffs[3]; coeffs = [sp.simplify(c / lead) for c in coeffs]
def wilczynski(coeffs):
    a0, a1, a2, _ = coeffs
    h = sp.exp(-sp.integrate(a2, z) / 3)
    expr = apply_theta_poly(coeffs, h * y(z)) / h
    expr = sp.expand(sp.simplify(expr))
    q = sp.simplify(expr.coeff(y(z)).subs({sp.Derivative(y(z), z, k): 0 for k in (1, 2, 3)}))
    p = sp.simplify(expr.coeff(sp.Derivative(y(z), z)))
    d2 = sp.simplify(expr.coeff(sp.Derivative(y(z), (z, 2))))
    assert d2 == 0
    return sp.simplify(q - sp.diff(p, z) / 2), p, q
W, p, q = wilczynski(coeffs)
print('Sol 2 at t = 1/3: W = q - p\'/2 =', sp.factor(W))
print('   W contains E:', W.has(E), '; W == 0:', W == 0)
# control (i): sym^2 of psi'' = Q psi is y''' - 4 Q y' - 2 Q' y
Q = sp.exp(2 * z) + al * sp.exp(z) + lam + sp.exp(-2 * z)
Wc, _, _ = wilczynski([-2 * sp.diff(Q, z), -4 * Q, sp.Integer(0), sp.Integer(1)])
print('control sym^2 of a generic Schrodinger operator: W =', sp.simplify(Wc))
# control (ii): Suzuki(M = -2) in z = log x: x^2 psi_xx = (theta^2 - theta) psi; normal form psi = x^(1/2) chi: chi_zz = (x^2 Q_x + 1/4) chi
Qs = sp.exp(-2 * z) + al * sp.exp(-z) + lam - E * sp.exp(2 * z) + sp.Rational(1, 4)
Ws, _, _ = wilczynski([-2 * sp.diff(Qs, z), -4 * Qs, sp.Integer(0), sp.Integer(1)])
print('control sym^2 of Suzuki(M = -2): W =', sp.simplify(Ws))
sys.exit(0)
