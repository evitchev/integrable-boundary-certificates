"""SOL12 engine controls: (1) BLZ calibration; (2) gwkb vs thwkb (VIR4-R/VIR5b engine) on Sol 3 k = 2 and k = -6, monic charges spins 1..9."""
import sys, os, time
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
import gwkb, thwkb
LOG = open('run_t_engine.log', 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
fails = 0
l = sp.Symbol('l'); L0, L1 = sp.symbols('l0 l1')
def topoly(expr_poly):
    return sp.expand(sum(sp.Rational(v.numerator, v.denominator) * L0**i * L1**j for (i, j), v in expr_poly.items()))
def charge(sk, bb):
    if any(a.denominator == 1 and a <= 0 for (a, c) in sk):
        pw, poly, nonpole = thwkb.integrate_reg(sk, bb); return sp.expand(sum(v * L0**i * L1**j for (i, j), v in poly.items()))
    r = thwkb.integrate(sk, bb)
    return topoly(r[2]) if r else 0
def monic(e, j):
    c = sp.Poly(e, L0, L1).coeff_monomial(L0**(2 * j)); return sp.expand(e / c)
# (1) BLZ: sigma^2 - (l + 1/2)^2, s = 1/2, n = 2, M = 1
T0 = time.time()
sh = {(1, 0): F(1), (0, 0): F(1, 2)}
blocks = [(1, sh, F(1)), (1, thwkb.pscale(sh, -1), F(1))]
sv = gwkb.wkb(blocks, F(1, 2), F(1), 4)
a_ = (2 * L0 + 3) / 4
rat = []
for k, ex in ((2, -2 * sp.bernoulli(2, a_)), (4, -sp.Rational(16, 3) * sp.bernoulli(4, a_))):
    pw, pref, poly = thwkb.integrate(sv[k], 2)
    rat.append(sp.simplify(pref * topoly(poly) / sp.expand(ex)))
say(f'(1) BLZ calibration ratios {rat} ({time.time()-T0:.1f}s)')
ok1 = all(sp.simplify(r**2 - 1) == 0 for r in rat); fails += not ok1
say(f'    (|ratio| = 1; the overall sign is the branch S = -Sigma): {ok1}')
# (2) Sol 3 family members vs thwkb
for k in (2, -6):
    n = k + 4; s = F(n - 1, 2); d = F(n, k); M = F(1, k)
    blocks = [(1, {(1, 0): F(1)}, d), (1, {(1, 0): F(-1)}, d), (1, {(0, 1): F(1)}, d), (1, {(0, 1): F(-1)}, d), (k, {}, d)]
    T0 = time.time(); sv = gwkb.wkb(blocks, s, M, 10); tg = time.time() - T0
    pairs = [{(0, 0): s, (1, 0): F(1)}, {(0, 0): s, (1, 0): F(-1)}, {(0, 0): s, (0, 1): F(1)}, {(0, 0): s, (0, 1): F(-1)}]
    if k > 0:
        frozen = [s + d * F(2 * jj - (k - 1), 2) for jj in range(k)]
        ref = thwkb.wkb(thwkb.charpoly_coeffs([{(0, 0): fz} for fz in frozen] + pairs, n), n, F(k), d, 10)
    else:
        c = s - d * F(k - 1, 2); right = [{(0, 0): c - j * d} for j in range(1, -k + 1)]
        ref = thwkb.wkb2(thwkb.charpoly_coeffs(pairs, 4), 4, thwkb.charpoly_coeffs(right, len(right)), len(right), F(k), d, 10)
    ok = True; rep = {}
    for kk in range(2, 11, 2):
        a1 = monic(charge(sv[kk], d), kk // 2); a2 = monic(charge(ref[kk], d), kk // 2)
        rep[kk - 1] = sp.expand(a1 - a2) == 0; ok &= rep[kk - 1]
    say(f'(2) Sol 3 k = {k}: gwkb monic charges == thwkb, spins 1..9: {rep} (gwkb {tg:.1f}s)')
    fails += not ok
say('ENGINE CONTROLS', 'FAIL' if fails else 'PASS'); sys.exit(1 if fails else 0)
