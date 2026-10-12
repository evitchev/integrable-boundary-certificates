"""RED2 (A2): Laplace transform of Suzuki(M = -2) at the exact point and the pillow.
Suzuki: x^4 psi'' = (1 + alpha x + lam x^2 - E x^4) psi, alpha^2 = -8X, lam = 2 + Y (SUPER3b at M = -2).
Laplace psi = int e^{p x} chi dp:  L4 chi = ((p^2 + E) chi)'''' - lam chi'' + alpha chi' - chi = 0.
(1) indicial exponents of L4 at p = +-i sqrt(E) (its only finite singular points);
(2) pillow at t = -1/3 (Sol 1; 38(ad) parameters): psi_zz = [kappa^2 w^(3/2)(1 - w) + V(w)] psi, d/dz = w(1 - w) d/dw,
    V = X(-w/2) + Y(-w(1-w)/4) + (-5/16) w(1-w) + (1/16)(1-w); exponents at the regular singular points w = 0, w = 1 (in w), and in s = sqrt(w).
A right factor (gauge s^a (1-s^2)^b, Moebius in p) or a symmetric cube must carry the pillow's exponent DIFFERENCES (mod Z, times ramification) at the matched point."""
import sys; sys.dont_write_bytecode = True
import sympy as sp
p, r, E, al, lam, X, Y, w, kap = sp.symbols('p r E alpha lam X Y w kappa')
chi = sp.Function('chi')
q = sp.Symbol('q'); p0 = sp.I * sp.sqrt(E)
# p = p0 + q, p^2 + E = q^2 + 2 p0 q (p0^2 = -E); chi = q^r; lowest power q^(r-3) comes only from the leading term
f = q**r
expr = sp.diff((q**2 + 2 * p0 * q) * f, q, 4) - lam * sp.diff(f, q, 2) + al * sp.diff(f, q) - f
lead = sp.factor(sp.expand(sp.powsimp(sp.expand(expr * q**(3 - r)), force=True)).subs(q, 0))
print('L4 indicial polynomial at p = i sqrt(E):', lead, '-> exponents', sp.solve(lead, r))
# pillow: in z, near w = 0 (z -> -inf, w ~ e^z) and w = 1 (z -> +inf, 1 - w ~ e^-z); Lambda -> 0 at both, so exponents from V(0), V(1)
V = X * (-w / 2) + Y * (-w * (1 - w) / 4) - sp.Rational(5, 16) * w * (1 - w) + sp.Rational(1, 16) * (1 - w)
V0, V1 = sp.simplify(V.subs(w, 0)), sp.simplify(V.subs(w, 1))
print('pillow V(0) =', V0, ' -> psi ~ w^(+-sqrt(V0)) = w^(+-', sp.sqrt(V0), '): exponent difference', 2 * sp.sqrt(V0), '(in s = sqrt(w):', 4 * sp.sqrt(V0), ')')
print('pillow V(1) =', V1, ' -> psi ~ (1-w)^(+-sqrt(V1)): exponent difference', 2 * sp.sqrt(V1), '(momentum-dependent; generic irrational)')
print('CONCLUSION (exponent obstruction): every exponent difference of L4 at its finite singular points is an integer; the pillow X-point difference')
print('  2 sqrt(-X/2) is not an integer for generic X, and no ramified map makes it one. So the pillow is not a right factor of L4, nor is its symmetric')
print('  cube equal to L4, for ANY identification sending the pillow X-point to p = +-i sqrt(E).  Identifications sending it to p = infinity (irregular) NOT tested.')
