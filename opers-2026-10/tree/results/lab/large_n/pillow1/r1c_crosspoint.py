"""PILLOW1 route (1c, POST-SEAL, labelled): the pillow's 3-term invariant {2/3, -2} equals Sol 1's Gamma-form invariant at t' = +1/3 and t' = -5.
Test: pillow (t = -1/3; 3-term recurrence after the gauge (1-s^2)^(sqrt(-X/2))) = Mellin-multiplier image of Sol 1's Gamma operator at t',
identically in nu, with an affine index map nu = lam mu' + delta and an unknown momentum dictionary (l0, l1).  Cocycle F == 1 (as in SECOND1)."""
import sys; sys.dont_write_bytecode = True
import sympy as sp
tp = sp.Rational(sys.argv[1]) if len(sys.argv) > 1 else sp.Rational(1, 3)
hs = int(sys.argv[2]) if len(sys.argv) > 2 else 1     # gauge sign: h = hs sqrt(-X/2)
nu, X, Y, delta, L0, L1 = sp.symbols('nu X Y delta l0 l1')
sq = sp.sqrt(2) * sp.sqrt(-X)
c0 = sq + 2                                    # from r1_pillow_mellin: s^2 coefficient -[(2T + 2 c0')^2 ...]; here 2T + 2 sqrt2 sqrt(-X) + 2
p0 = lambda T: (2 * T - 1) * (2 * T + 1) / 4
p2 = lambda T: -((2 * T + hs * 2 * sp.sqrt(2) * sp.sqrt(-X) + 2)**2 - (4 * Y + 9)) / 4
p3 = lambda T: sp.Integer(-4)                  # times kappa^2 (absorbed into E by a constant x-scale)
# Sol 1 Gamma form at t': H(T) = (T^2 - l1^2) G_a(T - l0), a = (t'-1)/2, n = (t'+3)/2, sigma = n/a, c = n + sigma, s' = (n - 1)/2 (Gamma-form centre, item 306 T = s - theta)
a = (tp - 1) / 2; n = (tp + 3) / 2; sig = n / a; c = n + sig; sS = (n - 1) / 2
P = lambda T: T**2 - L1**2
Grat = lambda T: ((T - L0) / sig + (a + 1) / 2) / ((T - L0) / sig + (1 - a) / 2)    # G_a(T - l0 + sigma)/G_a(T - l0)
# affine index map nu = lam (mu' - base) with pillow roles: one non-E pillow shift -> Sol-1 base (0), the other -> c, pillow E shift 3 -> n
roles = []
for base, other in ((2, 0), (0, 2)):
    lam_ = c / (other - base)
    if sp.simplify(lam_ * (3 - base) - n) == 0: roles.append((base, other, lam_))
print(f"t' = {tp}: consistent role maps (pillow base, pillow c-shift, lam): {roles}", flush=True)
if not roles: print('RESULT: NONE (no affine shift map; invariant mismatch)'); sys.exit(0)
base, other, lam = roles[0]
mu = lambda v: (v - delta) / lam               # pillow index mu' as a function of nu
# m(nu)/m(nu - c) = f_c = -p0(mu+2) H(sS - nu)/p2(mu);  m(nu)/m(nu - n) = f_n = p3(mu - 1) H(sS - nu)/p2(mu)   (E constant absorbed)
def Hratio(v1, v2):   # H(sS - v1)/H(sS - v2) with v1 - v2 = -sigma or +sigma  (sS - v1 = sS - v2 + sigma when v1 = v2 - sigma)
    T2 = sS - v2; d = sp.simplify((sS - v1) - T2)
    base = P(sS - v1) / P(sS - v2)
    if sp.simplify(d - sig) == 0: return base * Grat(T2)
    if sp.simplify(d + sig) == 0: return base / Grat(sS - v1)
    raise ValueError('shift is not +-sigma')
pc = {0: p0, 2: p2, 3: p3}
fc = lambda v: -pc[other](mu(v) + base - other) / pc[base](mu(v))
fn = lambda v: pc[3](mu(v) + base - 3) / pc[base](mu(v))
# cocycle f_c(nu) f_n(nu - c) = f_n(nu) f_c(nu - n), with the H factors: H(sS-nu)H(sS-nu+c) vs H(sS-nu)H(sS-nu+n) -> Hratio(nu - c, nu - n)
F = fc(nu) * fn(nu - c) / (fn(nu) * fc(nu - n)) * Hratio(nu - c, nu - n)
num = sp.Poly(sp.expand(sp.numer(sp.together(F - 1))), nu)
eqs = [sp.expand(e) for e in num.coeffs()]
print(f"t' = {tp}: a = {a}, n = {n}, sigma = {sig}, c = {c}, lam = {lam}; cocycle numerator degree {num.degree()}, {len(eqs)} equations", flush=True)
sols = sp.solve(eqs, [delta, L0, L1], dict=True)
print('solutions (delta, l0, l1):', sols, flush=True)
for so in sols:
    ok = sp.simplify(F.subs(so) - 1) == 0
    print('  identity check F == 1:', ok, flush=True)
print('RESULT:', 'CLOSES' if any(sp.simplify(F.subs(so) - 1) == 0 for so in sols) else 'NONE')
