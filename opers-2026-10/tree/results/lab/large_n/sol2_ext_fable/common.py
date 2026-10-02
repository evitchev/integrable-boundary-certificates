"""SOL12 shared: certified side (VIR3-R dictionary; spin-1 law; cyl(sol, 4); vev_sol*_w6/8/10), class G blocks, charges in (P_X, pi) with momentum scale C = c^2."""
import os, sys, json, hashlib
from fractions import Fraction as F
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import gwkb, thwkb
assert hashlib.sha256(open(os.path.join(HERE, 'SEAL_SOL2F.md'), 'rb').read()).hexdigest().startswith('80274804')
ORIG = {l.split()[1]: l.split()[0] for l in open(os.path.join(HERE, 'inputs', 'ORIGIN_SHA256'))}
assert all(hashlib.sha256(open(os.path.join(HERE, 'inputs', f_), 'rb').read()).hexdigest() == h for f_, h in ORIG.items())
PX, PI, C = sp.symbols('P_X pi C'); Pa, Qa = sp.symbols('P Q'); tS = sp.Symbol('t')
ROOT = os.path.expanduser('~/opus5-w8-clone/code'); sys.path.insert(0, ROOT); _cwd = os.getcwd(); os.chdir(ROOT)
from cyl_first_order import cyl, X, Y, t as tc
os.chdir(_cwd)
TABLES = {}
for sol in (1, 2, 3):
    for W in (6, 8, 10):
        d = json.load(open(os.path.join(HERE, 'inputs', f'vev_sol{sol}_w{W}.json')))
        TABLES[(sol, W)] = {(int(k_.split()[0][2:]), int(k_.split()[1][2:])): sp.sympify(v, locals={'t': tS}) for k_, v in d['coefficients'].items()}
def aXY(sol, t):
    t = sp.nsimplify(t) if not isinstance(t, sp.Basic) else t
    if sol == 1: return sp.Integer(1), (t - 1) / (t + 3)
    if sol == 2: return -2 * (t - 3) / (t + 5), 4 * (t - 1) / (t + 5)
    return (t - 3) / (t - 5), (t - 3) / (t - 5)
def regular_at(sol, tv):
    for W in (6, 8, 10):
        for c in TABLES[(sol, W)].values():
            if sp.denom(sp.cancel(c)).subs(tS, tv) == 0: return False
    return True
def certified(sol, tv):
    p2 = 2 * (tv - 1) / (tv + 1); Nv = (tv**2 - 25) / (tv**2 - 1); rho2 = (p2 - 2)**2 / (4 * p2)
    def a2i(e): return sp.expand(e.subs({Pa: sp.I * PX}).subs({Qa: sp.sqrt(rho2 - PI**2)}))
    def norm(e, j): return sp.expand(e / sp.Poly(e, PX, PI).coeff_monomial(PX**(2 * j)))
    out = {1: norm(a2i(Pa**2 + Qa**2 + (Nv + 1) / 12), 1),
           3: norm(a2i(sp.expand(cyl(sol, 4).subs(tc, tv).subs({X: Pa**2, Y: Qa**2}))), 2)}
    for W, j in ((6, 3), (8, 4), (10, 5)):
        out[2 * j - 1] = norm(a2i(sp.expand(sum(sp.cancel(c.subs(tS, tv)) * Pa**i * Qa**jj for (i, jj), c in TABLES[(sol, W)].items()))), j)
    return out
def classG(sol, tv, alpha_tamper=0):
    """beta = 1 (sealed): D = n = 2/a_X, alpha = a_Y D/2, gamma = D(1 - a_X - a_Y)"""
    aX, aY = aXY(sol, tv); D = 2 / aX; al = aY * D / 2 + alpha_tamper; ga = D - 2 * al - 2
    tof = lambda x: F(int(sp.fraction(x)[0]), int(sp.fraction(x)[1]))
    return dict(aX=aX, aY=aY, D=tof(D), alpha=tof(al), gamma=tof(ga), s=tof((D - 1) / 2))
def blocks_of(g, M):
    delta = g['D'] * F(M)
    return [(g['alpha'], {(1, 0): F(1)}, delta), (g['alpha'], {(1, 0): F(-1)}, delta),
            (F(1), {(0, 1): F(1)}, delta), (F(1), {(0, 1): F(-1)}, delta), (g['gamma'], {}, delta)]
def charges(g, M, K):
    """monic-ready charges as sympy polys in (P_X, pi, C): l0^2 = C P_X^2/a_Y, l1^2 = C pi^2/a_X"""
    bb = g['D'] * F(M); sv = gwkb.wkb(blocks_of(g, M), g['s'], F(M), K)
    out = {}
    for kk in range(2, K + 1, 2):
        sk = sv[kk]
        if any(a.denominator == 1 and a <= 0 for (a, c) in sk):
            pw, poly, nonpole = thwkb.integrate_reg(sk, bb); poly = dict(poly)
        else:
            r = thwkb.integrate(sk, bb); poly = {k_: sp.Rational(v.numerator, v.denominator) for k_, v in r[2].items()} if r else {}
        assert not [m for m in poly if m[0] % 2 or m[1] % 2]
        out[kk - 1] = sp.expand(sum(v * (C * PX**2 / g['aY'])**(i // 2) * (C * PI**2 / g['aX'])**(j // 2) for (i, j), v in poly.items()))
    return out
def monic(e, j):
    return sp.expand(e / sp.Poly(e, PX, PI).coeff_monomial(PX**(2 * j)))
