"""VIR7 pre-seal validation of G1 on KNOWN operators (no Sol-1 input):
 (1) Sol 3 generic-t operator: G1 == Codex's Sol-3 loss-1 amplitudes;
 (2) Sol 2 family S2(k) (cc, item 301/SOL2F): G1 == Codex's Sol-2 loss-1 amplitudes;
 dictionary l0^2 = C PX^2/a_Y, l1^2 = C pi^2/a_X, C = 2n(1+M).  Exit 0 iff both hold at every tested (K, t)."""
import json, sys
import sympy as sp
sys.dont_write_bytecode = True
import g1_lib as G
AUD = '<repo>/results/lab/audits/codex_wardsol12b_27403fc_2026-09-28/'
ks, ts = sp.symbols('k t')


def load_law(sol):
    law = json.load(open(AUD + 'AMPLITUDES_sol%d_loss1.json' % sol))
    f = lambda d: sp.sympify(d['numerator'], locals={'k': ks, 't': ts}, rational=True) / sp.sympify(d['denominator'], locals={'k': ks, 't': ts}, rational=True)
    return {int(j): f(v) for j, v in law['shift_amplitudes'].items()}, f(law['p']), f(law['q']), f(law['r'])


def Kd(d, p, q, r, a):
    return sp.binomial(d, a) * sp.rf(p, a) * sp.rf(q, d - a) * r ** (d - a) / sp.rf(p, d)


def ward_layer(law, K, t):
    amps, p, q, r = law
    sub = {ks: K, ts: t}
    p, q, r = [sp.Rational(x.subs(sub)) for x in (p, q, r)]
    d = K - 1
    out = {}
    for a in range(d + 1):
        out[(a, d - a)] = sum(sp.Rational(v.subs(sub)) * Kd(d, p + j, q, r, a) for j, v in amps.items())
    return out, (p, q, r)


bad = tot = 0
for sol in (3, 2):
    law = load_law(sol)
    for tstr in ('9/4', '3/2', '2', '4', '-2'):
        t = sp.Rational(tstr)
        k = 2 * (t - 1) / (3 - t)
        rho2 = 2 / (t * t - 1)
        M = 1 / k
        if sol == 3:
            n = k + 4
            aX = aY = (t - 3) / (t - 5)
            C = 2 * n * (1 + M)
            fac = [(G.l0, 1), (-G.l0, 1), (G.l1, 1), (-G.l1, 1)]
            strs = [(G.ZERO, k)]
        else:
            n = 2 * k + 3
            aX = -2 * (t - 3) / (t + 5); aY = 4 * (t - 1) / (t + 5)
            C = 2 * n * (1 + M)
            fac = [(G.l0, k), (-G.l0, k), (G.l1, 1), (-G.l1, 1)]
            strs = [(G.l0, k), (-G.l0, k)]
        fb = ft = 0
        for K in range(2, 21):
            try:
                expect, _ = ward_layer(law, K, t)
            except ZeroDivisionError:
                continue
            if any(v.has(sp.zoo, sp.nan) for v in expect.values()):
                continue
            top, l1 = G.charge(2 * K, n, M, fac, strs)
            got = G.to_record(top, l1, K, C / aY, C / aX, rho2)
            if got is None:
                continue
            got = got[0]
            for key, v in expect.items():
                ft += 1
                fb += got.get(key, 0) != v
        bad += fb; tot += ft
        print('Sol %d t = %-4s (k = %-5s n = %-5s): K = 2..20, %d coefficients, mismatches %d' % (sol, tstr, k, n, ft, fb), flush=True)
print('G1 validation:', 'PASS' if bad == 0 and tot > 0 else 'FAIL', bad, '/', tot)
sys.exit(0 if bad == 0 else 2)
