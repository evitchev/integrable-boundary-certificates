"""PILLOW1 route (1): the pillow's Mellin-recurrence shift structure.
The branch points of w^(1/2) (w = 0 and w = infinity) must sit at x = 0 / infinity for a power map to give polynomial coefficients, so the admissible
Moebius placements are w = x and w = 1/x; then x = s^2.  Pillow in s (theta = s d/ds), times 4:
    L = (1 - s^2)^2 theta^2 - 2 s^2 (1 - s^2) theta - 4 kappa^2 s^3 (1 - s^2) - 4 V(s^2).
Algebraic gauge (1 - s^2)^h s^g (h may depend on the momenta): conjugate, then divide by the largest power of (1 - s^2) dividing ALL coefficients.
Report the minimal set of s-monomials (the Mellin shifts) and which carry kappa^2 (= E)."""
import sys; sys.dont_write_bytecode = True
import sympy as sp
s, X, Y, kap, h, g = sp.symbols('s X Y kappa h g')
y = sp.Function('y')
def V(w): return -X * w / 2 - Y * w * (1 - w) / 4 - sp.Rational(5, 16) * w * (1 - w) + sp.Rational(1, 16) * (1 - w)
def theta(f): return s * sp.diff(f, s)
def apply_L(f, inv=False):
    if not inv:
        c2, c1, c0 = (1 - s**2)**2, -2 * s**2 * (1 - s**2), -4 * kap**2 * s**3 * (1 - s**2) - 4 * V(s**2)
    else:   # w = 1/x = 1/s^2: redo the transformation from w to s with w = s^-2 (theta_w = -theta_s/2)
        w = s**-2
        c2, c1, c0 = sp.expand((1 - w)**2), sp.expand(-2 * w * (1 - w) * (-1)), sp.expand(-4 * kap**2 * w**sp.Rational(3, 2) * (1 - w) - 4 * V(w))
    return c2 * theta(theta(f)) + c1 * theta(f) + c0 * f
def coeffs_after_gauge(hval, gval, inv=False):
    G = (1 - s**2)**hval * s**gval
    expr = sp.expand(sp.simplify(apply_L(G * y(s), inv) / G))
    d2 = sp.Derivative(y(s), (s, 2)); d1 = sp.Derivative(y(s), s)
    a2 = sp.simplify(expr.coeff(d2)); rest = sp.expand(expr - a2 * d2)
    a1 = sp.simplify(rest.coeff(d1)); a0 = sp.simplify(sp.expand(rest - a1 * d1).coeff(y(s)))
    # back to theta form: a2 y'' + a1 y' + a0 y with s^2 y'' = theta^2 y - theta y, s y' = theta y
    t2 = sp.simplify(a2 / s**2); t1 = sp.simplify(a1 / s - a2 / s**2)
    return [sp.factor(t2), sp.factor(t1), sp.factor(a0)]
print('== placement w = s^2, no gauge:')
c = coeffs_after_gauge(0, 0)
for name, cc in zip(('theta^2', 'theta^1', 'theta^0'), c): print(f'   {name}: {cc}')
# exponents at s = 1 of the ungauged operator: indicial from the (1-s)^r ansatz
r = sp.Symbol('r')
L0 = [sp.expand(cc) for cc in c]
# local exponent at s = 1: theta^2 coeff ~ 4(1-s)^2, theta^1 ~ -4(1-s), theta^0 -> -4 V(1)
ind = sp.factor(sp.limit(((1 - s)**2 * 4) * r * (r - 1) / (1 - s)**2 + (-2 * 2 * (1 - s)) * (-r) / (1 - s) * 0, s, 1))
print('   V(1) =', V(1), '; exponents at s = 1 satisfy 4 r^2 = 4 V(1) -> r = +-sqrt(V(1)) (in the variable 1 - s^2)')
best = None
for hv in (0, sp.sqrt(V(1)), -sp.sqrt(V(1))):
    for gv in (0, sp.Rational(1, 2), -sp.Rational(1, 2), sp.Rational(1, 4), -sp.Rational(1, 4)):
        cc = coeffs_after_gauge(hv, gv)
        # divide by the common (1 - s^2) power
        k = 0
        while all(sp.simplify(sp.cancel(ci / (1 - s**2)**(k + 1))).is_polynomial(s) for ci in cc if ci != 0): k += 1
        red = [sp.expand(sp.cancel(ci / (1 - s**2)**k)) for ci in cc]
        if not all(ri.is_polynomial(s) for ri in red):
            continue
        mons = sorted({m[0] for ri in red for m in sp.Poly(ri, s).monoms()})
        emons = sorted({m[0] for ri in red for (m, cf) in zip(sp.Poly(ri, s).monoms(), sp.Poly(ri, s).coeffs()) if cf.has(kap)})
        print(f'   gauge h = {hv}, g = {gv}: common (1-s^2)^{k} removed; s-shifts {mons} ({len(mons)} terms); kappa^2 carried by shifts {emons}')
        if best is None or len(mons) < best[0]: best = (len(mons), hv, gv, mons, emons, red)
print('MINIMAL:', best[:5])
print('   reduced theta-polynomials per shift:')
red = best[5]
for m in best[3]:
    poly = sum(sp.Poly(ri, s).coeff_monomial(s**m) * sp.Symbol('T')**(2 - j) for j, ri in enumerate(red))
    print(f'     s^{m}: {sp.factor(poly)}')
