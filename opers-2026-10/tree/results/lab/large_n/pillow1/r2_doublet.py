"""PILLOW1 route (2): Suzuki(M = -2) +-alpha DOUBLET tensor product (4th order; charges add -> charge-isospectral with the pillow after monic
normalisation) vs the pillow's 4th-order relatives, by the Mellin shift-pattern invariant (multipliers / power maps / overall shifts / algebraic gauges).
Doublet: (theta^2 - theta - lam) chi = (x^-2 +- alpha x^-1 - E x^2) chi (same lam for both).  Pillow relatives: Sym^3 of the reduced 3-term pillow
(p0 + s^2 p2 + s^3 p3), and pillow (x) (first-order factor theta - f, f = const or shift monomial)."""
import sys; sys.dont_write_bytecode = True
import sympy as sp
D, E, al, lam, X, Y, kap = sp.symbols('D E alpha lam X Y kappa')
def tensor_relation(A1a, A0a, A1b, A0b, var, thetaf):
    """theta^2 u = A1a theta u + A0a u, theta^2 v = A1b theta v + A0b v (coefficients in var); relation for phi = u v."""
    def tv(v):
        c1, c2, c3, c4 = v
        out = [thetaf(c1), thetaf(c2), thetaf(c3), thetaf(c4)]
        out[1] += c1; out[2] += c1
        out[0] += c2 * A0a; out[1] += c2 * A1a; out[3] += c2
        out[0] += c3 * A0b; out[2] += c3 * A1b; out[3] += c3
        out[2] += c4 * A0a; out[1] += c4 * A0b; out[3] += c4 * (A1a + A1b)
        return [sp.expand(o) for o in out]
    vecs = [[sp.Integer(1), 0, 0, 0]]
    for _ in range(4): vecs.append(tv(vecs[-1]))
    ns = sp.Matrix(vecs).T.nullspace(); rel = [sp.cancel(r / ns[0][4]) for r in ns[0]]
    den = sp.lcm([sp.denom(sp.together(r)) for r in rel])
    return [sp.expand(sp.cancel(r * den)) for r in rel], sp.factor(den)
def sym3_relation(A1, A0, var, thetaf):
    """phi = y^3; basis y^3, y^2 y', y y'^2, y'^3."""
    def tv(v):
        c0, c1, c2, c3 = v   # coefficients on y^3, y^2 y', y y'^2, y'^3
        out = [thetaf(c0), thetaf(c1), thetaf(c2), thetaf(c3)]
        # theta(y^3) = 3 y^2 y'
        out[1] += 3 * c0
        # theta(y^2 y') = 2 y y'^2 + y^2 (A1 y' + A0 y)
        out[2] += 2 * c1; out[1] += c1 * A1; out[0] += c1 * A0
        # theta(y y'^2) = y'^3 + 2 y y' (A1 y' + A0 y)
        out[3] += c2; out[2] += 2 * c2 * A1; out[1] += 2 * c2 * A0
        # theta(y'^3) = 3 y'^2 (A1 y' + A0 y)
        out[3] += 3 * c3 * A1; out[2] += 3 * c3 * A0
        return [sp.expand(o) for o in out]
    vecs = [[sp.Integer(1), 0, 0, 0]]
    for _ in range(4): vecs.append(tv(vecs[-1]))
    ns = sp.Matrix(vecs).T.nullspace(); rel = [sp.cancel(r / ns[0][4]) for r in ns[0]]
    den = sp.lcm([sp.denom(sp.together(r)) for r in rel])
    return [sp.expand(sp.cancel(r * den)) for r in rel], sp.factor(den)
def pattern(rel, var):
    P = [sp.Poly(sp.expand(r), var) if sp.expand(r) != 0 else None for r in rel]
    mons = sorted({m[0] for p in P if p for m in p.monoms()})
    emon = sorted({m[0] for p in P if p for m, cf in zip(p.monoms(), p.coeffs()) if cf.has(E) or cf.has(kap)})
    return mons, emon
def normal(mons, emon):
    """affine-invariant normal forms of the (shift set, E-subset): all affine maps sending two chosen shifts to 0, 1."""
    out = set()
    for p in mons:
        for q in mons:
            if p == q: continue
            f = lambda z: sp.Rational(z - p, q - p)
            out.add((tuple(sorted(f(z) for z in mons)), tuple(sorted(f(z) for z in emon))))
    return out
# doublet in x: use x-power symbol; work with z = 1/x so all shifts are nonnegative: x^-2 = z^2, x^-1 = z, x^2 = z^-2 -> multiply by z^2 later
z = sp.Symbol('z')
thz = lambda f: sp.expand(-z * sp.diff(f, z))           # theta_x = -theta_z
# theta_x^2 chi = theta_x chi + (lam + x^-2 +- alpha x^-1 - E x^2) chi
A1 = sp.Integer(1)
A0p = lam + z**2 + al * z - E * z**-2; A0m = lam + z**2 - al * z - E * z**-2
rel, den = tensor_relation(A1, A0p, A1, A0m, z, thz)
rel = [sp.expand(r * z**8) for r in rel]
print('doublet tensor: denominator', den)
dm, de = pattern(rel, z); print('doublet tensor shifts (in z = 1/x):', dm, ' E-carrying:', de)
# pillow reduced 3-term in s: [p0(theta) + s^2 p2(theta) + s^3 p3] y = 0
s = sp.Symbol('s'); ths = lambda f: sp.expand(s * sp.diff(f, s))
sq = sp.sqrt(2) * sp.sqrt(-X)
p0 = lambda T: (2 * T - 1) * (2 * T + 1) / 4
p2 = lambda T: -((2 * T + 2 * sq + 2)**2 - (4 * Y + 9)) / 4
Tsym = sp.Symbol('T')
full = sp.expand(p0(Tsym) + s**2 * p2(Tsym) - 4 * kap**2 * s**3)
c2 = sp.expand(full.coeff(Tsym, 2)); c1 = sp.expand(full.coeff(Tsym, 1)); c0 = sp.expand(full.coeff(Tsym, 0))
A1p, A0pil = sp.cancel(-c1 / c2), sp.cancel(-c0 / c2)
relS, denS = sym3_relation(A1p, A0pil, s, ths)
sm, se = pattern(relS, s); print('Sym^3(pillow) shifts:', sm, ' E-carrying:', se, ' (denominator', denS, ')')
# pillow (x) first-order factor (theta - f): the product y w with theta w = f w is a gauge (w = s^f or exp), so the pattern = pillow's 3-term pattern
pm, pe = [0, 2, 3], [3]
nd, ns_, np_ = normal(dm, de), normal(sm, se), normal(pm, pe)
print('doublet vs Sym^3(pillow) share an affine normal form:', bool(nd & ns_))
print('doublet vs pillow (x) first order (= pillow 3-term pattern) share a normal form:', bool(nd & np_))
print('doublet term count', len(dm), '; Sym^3(pillow) term count', len(sm), '; pillow 3-term', len(pm))
