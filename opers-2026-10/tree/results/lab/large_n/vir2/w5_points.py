"""VIR2: the screening data at the known points (lead's request b), the Gram matrices
(sealed and computed), and the comparison with the record's exponents (A, B, C).
Post-main reporting script; every conservation statement is recomputed here.
Exit 0 iff all assertions hold."""
import sys, pickle
import os as _os_kit; _KITROOT = _os_kit.path.join(_os_kit.path.dirname(_os_kit.path.abspath(__file__)), '..', '..', '..', '..') + '/'  # KIT PATCH: repo root from __file__ (was an absolute path)
from fractions import Fraction as F
import sympy as sp
sys.dont_write_bytecode = True
import vir_lib as V
import vir2_lib as W
from vir2_lib import a, rho
sys.path.insert(0, _KITROOT + 'code')
import hashlib, json, time
t = sp.Symbol('t'); q = sp.Symbol('q'); p = sp.Symbol('p')
# the shared definitions (Fibre, conserved, PAIRS) are taken verbatim from w4_main.py
src = open('w4_main.py').read()
exec(src[src.index('I5SYM = '):src.index('# ---------------------------------------------------------------- fibres')])
fails = 0

print('== exponents: the record\'s (A, B, C) = (eps, -a_X, -a_Y) against the pairs, symbolic in p ==')
tp = (2 + p**2) / (2 - p**2)
ABC = {1: (-4/(t+3), sp.Integer(-1), -(t-1)/(t+3)),
       2: (-8/(t+5), 2*(t-3)/(t+5), -4*(t-1)/(t+5)),
       3: (4/(t-5), -(t-3)/(t-5), -(t-3)/(t-5))}
for sol in (1, 2, 3):
    A_, B_, C_ = [sp.cancel(x.subs(t, tp)) for x in ABC[sol]]
    assert sp.cancel(A_ + B_ + C_ + 2) == 0
    for name, bv, qv in PAIRS(sol, p):
        ratio_ok = sp.cancel(B_ * bv**2 - C_ * qv) == 0          # B : C = a^2 : beta^2
        k = sp.factor(sp.cancel(-B_ / qv))
        dot = sp.cancel(-qv + bv**2)                              # v_+ . v_-
        norm = sp.cancel(qv + bv**2)
        print('Sol %d %-3s: a^2 = %s, beta = %s;  B : C = a^2 : beta^2  %s;  (B, C) = -k (a^2, beta^2) with k = %s;'
              '  |v|^2 = %s, v_+.v_- = %s;  A = %s' % (sol, name, sp.factor(qv), bv, ratio_ok, k, sp.factor(norm), sp.factor(dot), sp.factor(A_)))
        fails += not ratio_ok
A1 = sp.cancel(ABC[1][0].subs(t, tp)); q1 = 1 - p**2 / 4; b1 = -p / 2
s1 = sp.cancel(A1 - (-q1 + b1**2) / q1) == 0
print('Sealed S-exponents (Sol 1): A = (v_+ . v_-)/a^2: %s ; B = -1 ; C = -beta^2/a^2: %s'
      % (s1, sp.cancel(ABC[1][2].subs(t, tp) + b1**2 / q1) == 0))
fails += not s1

print('\n== Gram matrices of the full screening sets (rows/cols: Vir beta_+, Vir beta_-, then the pairs +,-) ==')
def gram(vs):
    """vs: list of (sign, a^2, beta): vector (sign * a, beta).  Entries sign_i sign_j a_i a_j + beta_i beta_j;
    within one pair a_i a_j = a^2; across the two Sol-3 pairs a_i a_j = (1/p) * p = 1."""
    def xx(v, w):
        if v[0] == 0 or w[0] == 0:
            return 0
        return v[0] * w[0] * sp.sqrt(sp.cancel(v[1] * w[1])) if sp.cancel(v[1] - w[1]) != 0 else v[0] * w[0] * v[1]
    return sp.Matrix(len(vs), len(vs), lambda i, j: sp.factor(sp.cancel(sp.simplify(xx(vs[i], vs[j])) + vs[i][2] * vs[j][2])))
for sol in (1, 2, 3):
    vs = [(0, 0, p), (0, 0, -2 / p)]
    labels = ['v_phi+ = (0, p)', 'v_phi- = (0, -2/p)']
    for name, bv, qv in PAIRS(sol, p):
        vs += [(1, qv, bv), (-1, qv, bv)]
        labels += ['%s+ = (a, %s)' % (name, bv), '%s- = (-a, %s)' % (name, bv)]
    G = gram(vs)
    if sol == 3:      # a_S3 a_S3' = (1/p) p = 1 exactly (no branch ambiguity: sqrt(1/p^2 * p^2) = 1)
        assert G[2, 4] == sp.cancel(1 + (1 / p) * (-p))
    print('Sol %d  [%s]' % (sol, '; '.join(labels)))
    for r in range(G.rows):
        print('     ', [str(x) for x in G.row(r)])
    rk = G.rank(iszerofunc=lambda x: sp.cancel(x) == 0)
    print('      rank over Q(p): %d' % rk)
    fails += rk != 2
u = sp.Symbol('u')
print('Sealed Gram, Sol 1 {v_phi+, v_+, v_-} in u (p^2 = 2u^2): [[2u^2,-u^2,-u^2],[-u^2,1,u^2-1],[-u^2,u^2-1,1]]')
G1 = gram([(0, 0, p), (1, 1 - p**2 / 4, -p / 2), (-1, 1 - p**2 / 4, -p / 2)]).subs(p, sp.sqrt(2) * u).applyfunc(sp.expand)
ok1 = G1 == sp.Matrix([[2*u**2, -u**2, -u**2], [-u**2, 1, u**2 - 1], [-u**2, u**2 - 1, 1]])
print('   computed:', G1.tolist(), ' == sealed:', ok1)
G3 = gram([(0, 0, -2 / p), (1, 1 / p**2, 1 / p), (-1, 1 / p**2, 1 / p)]).subs(p, sp.sqrt(2) * u).applyfunc(sp.expand)
ok3 = (G3 * u**2).applyfunc(sp.expand) == sp.Matrix([[2, -1, -1], [-1, 1, 0], [-1, 0, 1]])
print('Sealed Gram, Sol 3 {v_phi-, v_+, v_-}: (1/u^2)[[2,-1,-1],[-1,1,0],[-1,0,1]];  computed:', G3.tolist(), ' == sealed:', ok3)
fails += (not ok1) + (not ok3)
G2 = gram([(0, 0, p), (1, 1 - p**2, -p), (-1, 1 - p**2, -p)])
print('Sol 2 {v_phi+, S2+, S2-} (not sealed):', G2.tolist())
G3p = gram([(0, 0, p), (1, p**2, -p), (-1, p**2, -p)])
print("Sol 3 {v_phi+, S3'+, S3'-} (not sealed):", G3p.tolist())
print('Linear relations: Sol 1: v_phi+ + S1+ + S1- = 0;  Sol 2: 2 v_phi+ + S2+ + S2- = 0;'
      '  Sol 3: v_phi- + S3+ + S3- = 0 and 2 v_phi+ + S3\'+ + S3\'- = 0;  S3\'+- = -p^2 S3-+ (same lines, dual exponents: product -2)')

print('\n== the known points (values of the formulas; conservation recomputed at each point) ==')
POINTS = [('t = 3 (N = -2)', sp.Integer(1), False), ('t -> infinity (N = 1)', None, True),
          ('t = 5 (N = 0, s = -5)', 2 / sp.sqrt(3), False), ('t = -5 (N = 0, s = +5)', sp.sqrt(3), False),
          ('t = -1/3 (N = 28, s = -9)', 2 * sp.I, False), ('t = +1/3 (N = 28, s = +9)', sp.I, False),
          ('t = 5/3 (N = -25/2)', 1 / sp.sqrt(2), False), ('t = -3 (N = -2, other sheet)', sp.Integer(2), False)]
for name, pv, n1 in POINTS:
    fib = Fibre(N1=True) if n1 else Fibre(pv)
    pv = fib.p
    print('%s: p = beta_+ = %s, beta_- = %s, rho = %s' % (name, pv, sp.simplify(-2 / pv), fib.rho))
    for sol in (1, 2, 3):
        try:
            d3, d5 = fib.I3(sol), fib.I5(sol)
        except Exception as e:
            print('   Sol %d: densities not available here (%s)' % (sol, type(e).__name__)); continue
        for pname, bv, qv in PAIRS(sol, pv):
            qv, bv = sp.simplify(qv), sp.simplify(bv)
            if qv == 0:
                r3, r5 = conserved(d3, bv, 0), conserved(d5, bv, 0)
            else:
                r3, r5 = conserved(d3, bv, qv), conserved(d5, bv, qv)
            h = sp.simplify(bv**2 / 2 - fib.rho * bv)
            print('   Sol %d %-3s: a^2 = %-5s beta = %-12s h_phi = %-6s |v|^2 = %-5s v_+.v_- = %-5s  conserved I_3 %s, I_5 %s'
                  % (sol, pname, qv, bv, h, sp.simplify(qv + bv**2), sp.simplify(-qv + bv**2), r3[0], r5[0]))
            fails += (not r3[0]) or (not r5[0])
print('\n== N = 0: I_3 = 4 (T_X + T_phi)^2 there; is the whole curve a^2/2 + h_phi = 1 conserved? ==')
for name, pv in (('t = 5, Sol 3', 2 / sp.sqrt(3)), ('t = -5, Sol 2', sp.sqrt(3))):
    fib = Fibre(pv)
    sol = 3 if 'Sol 3' in name else 2
    d3, d5 = fib.I3(sol), fib.I5(sol)
    for bv in (sp.Rational(1, 3), sp.Rational(-5, 7)):
        h = sp.simplify(bv**2 / 2 - fib.rho * bv)
        qv = sp.simplify(2 * (1 - h))
        r3, r5 = conserved(d3, bv, qv), conserved(d5, bv, qv)
        print('   %s: generic beta = %s, h = %s, a^2 = 2(1 - h) = %s: conserved I_3 %s, I_5 %s' % (name, bv, h, qv, r3[0], r5[0]))
    r3 = conserved(d3, sp.Rational(1, 3), sp.Rational(5, 7))
    print('   %s: off the curve (beta, a^2) = (1/3, 5/7): I_3 %s' % (name, r3[0]))
print('\nPOINTS', 'FAILED (%d)' % fails if fails else 'OK')
sys.exit(1 if fails else 0)
