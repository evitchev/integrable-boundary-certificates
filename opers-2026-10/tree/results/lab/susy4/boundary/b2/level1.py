"""level1.py -- B2 Level I: the J-orthogonal screening data of Z and K versus Solution 1 (PREREG sec. 3, Level I).

Bosonisation of the fermions in the frame eta = (s, -s): psi^+- = psi_0 +- psi_1 are null (<psi^+ psi^+> = s - s = 0)
with <psi^+ psi^-> = 2s; psi^+- = sqrt(2s) e^{+- X_3}, X_3 a third boson with <X_3 X_3> = +log (eta_3 = 1).
A SUSY vertex (m.psi) e^{m.X} with m = m_1 = (1/2, 1/2): m.psi = psi^+/2 -> e^{m.X + X_3}; m_2: psi^-/2 -> e^{m.X - X_3}.
J = (1/s) psi_0 psi_1 - (2/s) dX_1 (n2frame V1); psi_0 psi_1 = -psi^+ psi^-/2 -> -s dX_3; so J -> a.dY with
a = (0, -2/s, -1) in Y = (X_0, X_1, X_3), metric (s, -s, 1); a.a = (s-4)/s = c/3 (J J double pole, V1).
H~ = a.Y; v = b.Y, b = (0, -1, -2) (b.a = 0), v.v = 4 - s = 4 kappa^2; u = X_0.  Every momentum vector m is written in
(u, v, H~) coordinates and its J-charge a.m is printed.  Solution 1's data (opers paper Sec. 4): V_+- = exp(+-kappa X + beta phi),
beta = -p/2, kappa^2 = 1 - p^2/4, Vir screenings e^{p phi}, e^{-(2/p) phi}, metric (1, 1).

All of this is exact linear algebra in sympy with s = p^2 symbolic; printed for PRIMARY.md.
"""
import sympy as sp

s = sp.Symbol('s', positive=True)
p = sp.sqrt(s)
eta3 = (s, -s, sp.Integer(1))


def dot(x, y, eta=eta3):
    return sp.expand(sum(xi * yi * ei for xi, yi, ei in zip(x, y, eta)))


a = (0, -2 / s, sp.Integer(-1))              # J direction (bosonised)
b = (0, sp.Integer(-1), sp.Integer(-2))      # v direction
u = (sp.Integer(1), 0, 0)
assert dot(a, b) == 0 and dot(a, u) == 0 and dot(b, u) == 0
gH, gv, gu = dot(a, a), dot(b, b), dot(u, u)
print(f"J direction a = {a}: a.a = {sp.factor(gH)} (= c/3 = {sp.factor((3 - 12/s)/3)});  v = b.Y, b = {b}: v.v = {gv};  u.u = {gu}")

# coordinates: m.Y = m_u u + m_v v + m_H H~  <=>  m = m_u u_vec + m_v b + m_H a  (as vectors in Y-space)
M = sp.Matrix([[u[i], b[i], a[i]] for i in range(3)])


def coords(m):
    sol = M.solve(sp.Matrix(m))
    return tuple(sp.simplify(x) for x in sol)


half = sp.Rational(1, 2)
vecs = {
    'S1 = (m_1.psi)e^{m_1.X}': (half, half, 1),
    'S2 = (m_2.psi)e^{m_2.X}': (half, -half, -1),
    'S3 = (-m_1.psi)e^{-m_1.X}': (-half, -half, 1),
    'S4 = (-m_2.psi)e^{-m_2.X}': (-half, half, -1),
    'W0 (exponential e^{-X_0})': (sp.Integer(-1), 0, 0),
}
print("\nmomenta in Y = (X_0, X_1, X_3), metric (s, -s, 1); (u, v, H~) coordinates; J-charge a.m; norm:")
for k, m in vecs.items():
    cu, cv, cH = coords(m)
    print(f"  {k:30s}: Y-momentum {m}; (u, v, H~) = ({cu}, {cv}, {cH}); J-charge {sp.simplify(dot(a, m))}; m.m = {sp.simplify(dot(m, m))}")

Z = ['S1 = (m_1.psi)e^{m_1.X}', 'S2 = (m_2.psi)e^{m_2.X}', 'W0 (exponential e^{-X_0})']
print("\nGram matrix of Tower Z's three (J-neutral) vectors (S1, S2, W0):")
GZ = sp.Matrix([[dot(vecs[i], vecs[j]) for j in Z] for i in Z])
sp.pprint(GZ)
kappa2 = 1 - s / 4
beta = -p / 2
kappa = sp.sqrt(kappa2)
sol1 = {'V+': (kappa, beta), 'V-': (-kappa, beta), 'e^{p phi}': (0, p), 'e^{-2phi/p}': (0, -2 / p)}


def dot2(x, y):
    return sp.expand(x[0] * y[0] + x[1] * y[1])


names = ['V+', 'V-', 'e^{p phi}']
G1 = sp.Matrix([[sp.simplify(dot2(sol1[i], sol1[j])) for j in names] for i in names])
print("\nGram matrix of Solution 1's (V+, V-, e^{p phi}) (metric (1,1), opers paper Sec. 4):")
sp.pprint(G1)
print("EQUAL:", sp.simplify(GZ - G1) == sp.zeros(3, 3))
print("Gram with the dual Virasoro screening e^{-2phi/p} instead:", [sp.simplify(dot2(sol1[i], sol1['e^{-2phi/p}'])) for i in names + ['e^{-2phi/p}']])

# the isometry (u, v) -> (X, phi): u = -p phi, v = (kappa/xi) X with xi = v-coordinate of S1 = -1/2:
# check: momentum components m_phi = -p m_u, m_X = (xi/kappa)^{-1}... as vectors: m.Y = m_u u + m_v v = m_u(-p phi) + m_v (kappa/xi) X
xi = coords(vecs['S1 = (m_1.psi)e^{m_1.X}'])[1]
print(f"\nisometry: u = -p phi, v = (kappa/xi) X with xi = {xi}: (kappa/xi)^2 = {sp.simplify(kappa2/xi**2)} = v.v = {gv}: {sp.simplify(kappa2/xi**2 - gv) == 0}")
for k in Z:
    cu, cv, cH = coords(vecs[k])
    mphi, mX = -p * cu, cv * kappa / xi
    print(f"  {k:30s} -> (m_X, m_phi) = ({sp.simplify(mX)}, {sp.simplify(mphi)})")
print("Solution 1:", {k: (sp.simplify(v[0]), sp.simplify(v[1])) for k, v in sol1.items()})

# background charges: coset T_uv inherits lambda along X_0 = u with T ∋ -(1/s) d^2 X_0; weight of e^{m u}: s m^2/2 + m
# (h = 1/2 eta_u m^2 - eta_u lambda_u m = s m^2/2 - s(-1/s) m = s m^2/2 + m); under u = -p phi: e^{m u} = e^{-p m phi}:
# in phi (metric 1) with T_Y' = 1/2 (d phi)^2 + rho' d^2 phi: h = (p m)^2/2 + rho' p m  => rho' p m = m => rho' = 1/p.
rho_opers = p / 2 - 1 / p
print(f"\nbackground charge along phi: coset inherits rho' = 1/p (c_uv = 1 + (1 - 12/p^2) = 2 - 12/p^2 = {sp.simplify(2 - 12/s)});"
      f" opers paper rho = p/2 - 1/p (c = 1 + N = {sp.simplify(1 + 1 - 12*rho_opers**2)}); difference (rho - rho') = {sp.simplify(rho_opers - 1/p)} = p/2 - 2/p")
for k, v in sol1.items():
    h_op = sp.simplify((v[0]**2 + v[1]**2) / 2 - rho_opers * v[1])
    h_cs = sp.simplify((v[0]**2 + v[1]**2) / 2 - (1 / p) * v[1])
    print(f"  weight of {k:12s}: under rho (opers) {h_op}; under rho' = 1/p (coset) {h_cs}")
