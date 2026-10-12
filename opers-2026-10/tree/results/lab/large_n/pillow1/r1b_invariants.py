"""PILLOW1 route (1b): asymmetric algebraic gauges (1-s)^h1 (1+s)^h2 s^g at the pillow's two images of the X-point, and the shift-ratio invariant.
Invariant of a 3-term recurrence with shifts {p, q} (no E) and e (carries E), under Mellin multipliers, power maps (incl. s -> 1/s) and overall shifts:
    the unordered pair { (q - p)/(e - p), (p - q)/(e - q) }."""
import sys; sys.dont_write_bytecode = True
import sympy as sp
from r1_pillow_mellin import apply_L, theta
s, X, Y, kap = sp.symbols('s X Y kappa')
y = sp.Function('y')
hX = sp.sqrt(2) * sp.sqrt(-X) / 2
def reduced(h1, h2, g):
    G = (1 - s)**h1 * (1 + s)**h2 * s**g
    expr = sp.expand(sp.simplify(apply_L(G * y(s)) / G))
    d2 = sp.Derivative(y(s), (s, 2)); d1 = sp.Derivative(y(s), s)
    a2 = sp.simplify(expr.coeff(d2)); rest = sp.expand(expr - a2 * d2)
    a1 = sp.simplify(rest.coeff(d1)); a0 = sp.simplify(sp.expand(rest - a1 * d1).coeff(y(s)))
    cc = [sp.factor(sp.simplify(a2 / s**2)), sp.factor(sp.simplify(a1 / s - a2 / s**2)), sp.factor(a0)]
    for fac in (1 - s, 1 + s):
        while all(sp.cancel(ci / fac).is_polynomial(s) for ci in cc if ci != 0): cc = [sp.cancel(ci / fac) for ci in cc]
    if not all(sp.expand(ci).is_polynomial(s) for ci in cc): return None
    cc = [sp.expand(ci) for ci in cc]
    mons = sorted({m[0] for ci in cc for m in sp.Poly(ci, s).monoms()})
    em = sorted({m[0] for ci in cc for m, cf in zip(sp.Poly(ci, s).monoms(), sp.Poly(ci, s).coeffs()) if cf.has(kap)})
    return mons, em
def inv(shifts, e):
    p, q = [x for x in shifts if x != e]
    return sorted([sp.Rational(q - p, e - p), sp.Rational(p - q, e - q)])
res = set()
for h1 in (0, hX, -hX):
    for h2 in (0, hX, -hX):
        r = reduced(h1, h2, 0)
        if r is None: continue
        mons, em = r
        tag = f'h1 = {h1}, h2 = {h2}: shifts {mons} ({len(mons)} terms), E at {em}'
        if len(mons) == 3 and len(em) == 1: tag += f'  -> invariant {inv(mons, em[0])}'; res.add(tuple(inv(mons, em[0])))
        print(tag, flush=True)
print('PILLOW 3-term invariants found:', res)
# Sol 1 Gamma / three-term form at t = -1/3: shifts 0 (P G), c (L), n (E);  n = (t+3)/2, a = (t-1)/2, sigma = n/a, c = n + sigma
t = sp.Rational(-1, 3); n1 = (t + 3) / 2; a1 = (t - 1) / 2; c1 = n1 + n1 / a1
print('Sol 1 (t = -1/3): shifts {0, c, n} =', (0, c1, n1), ' invariant', inv([0, c1, n1], n1))
# Sol 2 at t = +1/3: k = -1/2, n = 2k + 3, M = 1/k, c = n (1 + M)
k2 = sp.Rational(-1, 2); n2 = 2 * k2 + 3; c2 = n2 * (1 + 1 / k2)
print('Sol 2 (t = +1/3): shifts {0, c, n} =', (0, c2, n2), ' invariant', inv([0, c2, n2], n2))
print('Sol 1 Schrodinger form (Suzuki M = -2): shifts {0, -2, -1, 2 (E)}: 4 terms')
