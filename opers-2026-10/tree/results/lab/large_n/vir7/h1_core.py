"""helpers shared by the VIR7 scripts (copied from h1_scan.py so that importing them does not run the scan)"""
import hashlib, json
import sympy as sp
AUD = '<repo>/results/lab/audits/codex_wardsol12b_27403fc_2026-09-28/'
ks, ts = sp.symbols('k t')


def load_law(sol):
    raw = open(AUD + 'AMPLITUDES_sol%d_loss1.json' % sol, 'rb').read()
    law = json.loads(raw)
    f = lambda d: sp.sympify(d['numerator'], locals={'k': ks, 't': ts}, rational=True) / sp.sympify(d['denominator'], locals={'k': ks, 't': ts}, rational=True)
    return {int(j): f(v) for j, v in law['shift_amplitudes'].items()}, f(law['p']), f(law['q']), f(law['r']), hashlib.sha256(raw).hexdigest()


def Kd(d, p, q, r, a):
    return sp.binomial(d, a) * sp.rf(p, a) * sp.rf(q, d - a) * r ** (d - a) / sp.rf(p, d)


def ward(law, K, t, tamper=False):
    amps, p, q, r, _ = law
    sub = {ks: K, ts: t}
    p, q, r = [sp.Rational(x.subs(sub)) for x in (p, q, r)]
    vals = {j: sp.Rational(v.subs(sub)) for j, v in amps.items()}
    if tamper:
        vals[1] = vals[1] + sp.Rational(1, 1000)
    lay = {(a, K - 1 - a): sum(v * Kd(K - 1, p + j, q, r, a) for j, v in vals.items()) for a in range(K)}
    top = {(a, K - a): Kd(K, p, q, r, a) for a in range(K + 1)}
    return lay, top
