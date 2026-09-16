"""Certificate: the pillow-type cylindrical operator with its exact sub-leading completion (ansatz note sec. 38(ad)) IS the
pillow-brane vacuum oper of Lukyanov-Zamolodchikov (arXiv:1208.5259, eqs. (6.1)-(6.2)) with the third momentum P_1^2 = 1/6, and the
cylindrical spin-3 charge differs from that oper's by the closed-form constant Delta(t) (sec. 38(af)).  Exact rational arithmetic.

Operator (z = log u, w = u/(1+u)):  psi_zz = [kappa^2 u^-eps (1+u)^n + V] psi,
   V = -c_X X w/2 + (kappa_s c_X Y/2 + 1/4 + d2) w(1-w) + c0 (1-w),  c0 = -eps/24,  kappa_s = -(n+2)/(n-eps),
   Solution 1: n = -(t+7)/(t+3), eps = -4/(t+3), c_X = 1,             d2 = -1/((t+1)(t+3));
   Solution 3: n = -(t-7)/(t-5), eps =  4/(t-5), c_X = (t-3)/(t-5),   d2 = -(t-3)/((t-5)(t-1)(t+1));
   Solution 2: n = 2(t-7)/(t+5), eps = -8/(t+5), c_X = -2(t-3)/(t+5), d2 = -4/((t+1)(t+5)).
Its large-kappa coefficients I_s = int rho_s Lambda^{-s/2} dz (Riccati polynomials, continued Beta integrals as Pochhammer
ratios (s eps/2)_i / (s n/2)_i) are compared, normalized on X (spin 1) and X^2 (spin 3), with
LZ (5.13):  I_1 = sum P_j^2/4 - 1/8,
            I_3 = sum_j E_j (P_j^4/16 - P_j^2/16 + 1/192) + sum_{m != j} E_mj (P_m^2/4 - 1/24)(P_j^2/4 - 1/24) + sum_j H_j/240,
   (A.3): E_j = -a_j (3 a_m + 2)(3 a_k + 2), E_mj = -3 a_m a_j (3 a_k + 2), H_j = 8 - a_j^2 - 9(a_1 a_2 + a_2 a_3 + a_3 a_1) - 15 a_1 a_2 a_3,
   (4.31): (a_1, a_2, a_3) = (alpha_1^2, alpha_2^2, alpha_3^2) = (n nu, n(1 - nu), -(n + 2)), nu = eps/n,
at the identification of (6.2) [V = kappa^2 e^{-n nu x}(1+e^x)^n - (alpha^2 + beta^2 e^x)/(1+e^x) - (gamma^2 - 1/4) e^x/(1+e^x)^2,
alpha = sqrt(n nu) P_1/2, beta = sqrt(n(1-nu)) P_2/2, gamma = sqrt(n+2) P_3/2] with our V:
   alpha^2 = -c0  =>  P_1^2 = 1/6;   beta^2 = c_X X/2  =>  P_2^2 = 2 c_X X/(n(1-nu));   gamma^2 = -kappa_s c_X Y/2 - d2  =>  P_3^2 = 4 gamma^2/(n+2).
Checks at every sample t (poles of the families excluded):
  C1  (5.13) I_1 == the operator's I_1 (normalized on X) == X + Y + (N+1)/12, N = (t^2-25)/(t^2-1)   [the level-1 law];
  C2  (5.13) I_3 == the operator's I_3 exactly, all coefficients (normalized on X^2)                  [the identification];
  C3  e_3(certified, cyl_first_order) scaled to I_3's X^2 coefficient, minus I_3, is a CONSTANT equal to the closed form
      Delta(t): Sol 1 (t-23)(t-3)(3t+1)/(720(t+1)(t+3)^2); Sol 2 -(t-7)(t-2)(3t-1)(5t-11)/(45(t+1)(t+5)^3);
      Sol 3 (3t^4-22t^3+380t^2-1082t+817)/(720(t-5)^3(t-1))   [the third-sector gap; reconstructed from 25 t, hold-outs here].
Fail-closed: any inequality prints FAIL and the certificate exits 1.  ~20 s."""
import sys, os
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cyl_first_order import X, Y, cyl, t
R = lambda x: sp.Rational(x.numerator, x.denominator)
def params(sol, tv):
    if sol == 1: n, eps, cX, d2 = -(tv + 7) / (tv + 3), F(-4) / (tv + 3), F(1), F(-1) / ((tv + 1) * (tv + 3))
    elif sol == 3: n, eps, cX, d2 = -(tv - 7) / (tv - 5), F(4) / (tv - 5), (tv - 3) / (tv - 5), -(tv - 3) / ((tv - 5) * (tv - 1) * (tv + 1))
    else: n, eps, cX, d2 = 2 * (tv - 7) / (tv + 5), F(-8) / (tv + 5), -2 * (tv - 3) / (tv + 5), F(-4) / ((tv + 1) * (tv + 5))
    return n, eps, -(n + 2) / (n - eps), cX, -eps / 24, d2
# ---- the operator's WKB coefficients: Riccati polynomials in w with monomials (i, a, b) = w^i X^a Y^b
def add(p, q, s=F(1)):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, F(0)) + s * v
        if r[k] == 0: del r[k]
    return r
def mul(p, q):
    r = {}
    for (i1, a1, b1), v1 in p.items():
        for (i2, a2, b2), v2 in q.items():
            k = (i1 + i2, a1 + a2, b1 + b2); r[k] = r.get(k, F(0)) + v1 * v2
    return {k: v for k, v in r.items() if v != 0}
def scale(p, s): return {k: s * v for k, v in p.items()}
def dw(p):
    r = {}
    for (i, a, b), v in p.items():
        if i == 0: continue
        for k, s in (((i, a, b), i * v), ((i + 1, a, b), -i * v)): r[k] = r.get(k, F(0)) + s
    return {k: v for k, v in r.items() if v != 0}
def poch(a, i):
    r = F(1)
    for j in range(i): r *= a + j
    return r
def operator_I(sol, tv):
    n, eps, kap, cX, c0, d2 = params(sol, tv)
    onemw = {(0, 0, 0): F(1), (1, 0, 0): F(-1)}
    V = {(1, 1, 0): -cX / 2}; V = add(V, mul({(1, 0, 1): kap * cX / 2, (1, 0, 0): F(1, 4) + d2}, onemw)); V = add(V, scale(onemw, c0))
    ll = {(0, 0, 0): -eps / 2, (1, 0, 0): n / 2}
    def dz(k, r): return add(dw(r), mul(r, ll), F(-k))
    rho = {0: scale(ll, F(-1, 2))}
    rho[1] = scale(add(add(V, mul(rho[0], rho[0]), F(-1)), dz(0, rho[0]), F(-1)), F(1, 2))
    for k in range(1, 3):
        acc = dz(k, rho[k])
        for i in range(0, k + 1): acc = add(acc, mul(rho[i], rho[k - i]))
        rho[k + 1] = scale(acc, F(-1, 2))
    out = {}
    for s in (1, 3):
        a0, ab = s * eps / 2, s * n / 2; I = {}
        for (i, a, b), v in rho[s].items(): I[(a, b)] = I.get((a, b), F(0)) + v * poch(a0, i) / poch(ab, i)
        out[s] = sum(R(v) * X**a * Y**b for (a, b), v in I.items())
    return out
def pillow_I(sol, tv):
    n, eps, kap, cX, c0, d2 = [R(v) for v in params(sol, tv)]; nu = eps / n
    a = {1: n * nu, 2: n * (1 - nu), 3: -(n + 2)}
    P2 = {1: sp.Rational(1, 6), 2: 2 * cX * X / (n * (1 - nu)), 3: 4 * (-kap * cX * Y / 2 - d2) / (n + 2)}
    oth = lambda j: [m for m in (1, 2, 3) if m != j]
    E = {j: -a[j] * (3 * a[oth(j)[0]] + 2) * (3 * a[oth(j)[1]] + 2) for j in (1, 2, 3)}
    Emj = {(m, j): -3 * a[m] * a[j] * (3 * a[[k for k in (1, 2, 3) if k not in (m, j)][0]] + 2) for m in (1, 2, 3) for j in (1, 2, 3) if m != j}
    S2 = a[1] * a[2] + a[2] * a[3] + a[3] * a[1]; Pr = a[1] * a[2] * a[3]
    H = {j: 8 - a[j]**2 - 9 * S2 - 15 * Pr for j in (1, 2, 3)}
    I1 = sum(P2[j] / 4 for j in (1, 2, 3)) - sp.Rational(1, 8)
    I3 = sum(E[j] * (P2[j]**2 / 16 - P2[j] / 16 + sp.Rational(1, 192)) for j in (1, 2, 3)) \
       + sum(Emj[(m, j)] * (P2[m] / 4 - sp.Rational(1, 24)) * (P2[j] / 4 - sp.Rational(1, 24)) for (m, j) in Emj) + sum(H[j] for j in (1, 2, 3)) / 240
    return {1: sp.expand(I1), 3: sp.expand(I3)}
def Delta_closed(sol, tq):
    return {1: (tq - 23) * (tq - 3) * (3 * tq + 1) / (720 * (tq + 1) * (tq + 3)**2),
            2: -(tq - 7) * (tq - 2) * (3 * tq - 1) * (5 * tq - 11) / (45 * (tq + 1) * (tq + 5)**3),
            3: (3 * tq**4 - 22 * tq**3 + 380 * tq**2 - 1082 * tq + 817) / (720 * (tq - 5)**3 * (tq - 1))}[sol]
def normalized(expr, mono):
    p = sp.Poly(sp.expand(expr), X, Y); return sp.expand(p.as_expr() / p.coeff_monomial(mono))
SAMPLES = {1: [F(11, 7), F(13), F(5), F(3, 4), F(17, 3), F(-9, 2)], 2: [F(7, 5), F(1, 2), F(6, 5), F(9, 5), F(-40, 3)], 3: [F(13, 6), F(9, 5), F(12), F(-3), F(21), F(2, 7)]}
failures = []
def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok: failures.append(name)
for sol in (1, 2, 3):
    for tv in SAMPLES[sol]:
        tq = R(tv); Nq = (tq**2 - 25) / (tq**2 - 1)
        op = operator_I(sol, tv); pl = pillow_I(sol, tv)
        law1 = sp.expand(X + Y + (Nq + 1) / 12)
        check(f"C1 Solution {sol}, t = {tv}: (5.13) I_1 = operator I_1 = X + Y + (N+1)/12",
              sp.expand(normalized(pl[1], X) - normalized(op[1], X)) == 0 and sp.expand(normalized(pl[1], X) - law1) == 0)
        check(f"C2 Solution {sol}, t = {tv}: (5.13) I_3 = operator I_3 exactly (P_1^2 = 1/6)", sp.expand(normalized(pl[3], X**2) - normalized(op[3], X**2)) == 0)
        c3 = sp.Poly(sp.expand(cyl(sol, 4).subs(t, tq)), X, Y); I3p = sp.Poly(pl[3], X, Y)
        e3 = sp.expand(c3.as_expr() * I3p.coeff_monomial(X**2) / c3.coeff_monomial(X**2))
        gap = sp.expand(e3 - pl[3])
        check(f"C3 Solution {sol}, t = {tv}: e_3 - I_3 is the constant Delta(t) = {Delta_closed(sol, tq)}", gap.is_number and sp.simplify(gap - Delta_closed(sol, tq)) == 0)
if failures:
    print(f"FAILS: {len(failures)}"); sys.exit(1)
print("PILLOW IDENTIFICATION CERTIFIED: the sec.-38(ad) operator is the Lukyanov-Zamolodchikov pillow-brane vacuum oper with P_1^2 = 1/6 "
      "(LZ (5.13) equals the operator's I_1 and I_3 exactly at 17 sample t over the three solutions), the spin-1 law is X + Y + (N+1)/12, "
      "and the cylindrical spin-3 gap is the closed-form constant Delta(t) of sec. 38(af).")
