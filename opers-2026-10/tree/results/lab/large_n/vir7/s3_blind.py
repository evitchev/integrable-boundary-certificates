"""VIR7 T2b: BLIND comparison of the registered spin-11 and spin-13 predictions of S1(t) with the certified Solution-1
profiles results/lab/anchor11/vev_profile_sol1_w12.json and anchor13/vev_profile_sol1_w14.json.
To be run only after the lead has registered PREDICTION_VIR7.json (33af7eb7...).
Modes: real | tamper (one table coefficient + 1/1000).  Exit real: 0 = all agree, 2 = mismatch; tamper: 0 = fires."""
import hashlib, json, sys
import sympy as sp
sys.dont_write_bytecode = True
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
PRED = open('PREDICTION_VIR7.sha256').read().split()[0]
assert hashlib.sha256(open('PREDICTION_VIR7.json', 'rb').read()).hexdigest() == PRED and PRED.startswith('33af7eb7')
mode = sys.argv[1] if len(sys.argv) > 1 else 'real'
ROOT = '<repo>/results/lab/'
ts = sp.Symbol('t')
XR, X, Y = ring('X,Y', QQ)


def q_(x):
    x = sp.Rational(x)
    return QQ(int(x.p), int(x.q))


TAB = {}
for s, fn in ((11, 'anchor11/vev_profile_sol1_w12.json'), (13, 'anchor13/vev_profile_sol1_w14.json')):
    raw = open(ROOT + fn, 'rb').read()
    print('spin %d table sha256 %s' % (s, hashlib.sha256(raw).hexdigest()))
    TAB[s] = {}
    for key, val in json.loads(raw)['vev'].items():
        a, b = [int(x) for x in key.split(',')]
        num, den = sp.fraction(sp.together(sp.sympify(val, locals={'t': ts}, rational=True)))
        TAB[s][a, b] = (sp.Poly(num, ts), sp.Poly(den, ts))
pred = json.load(open('PREDICTION_VIR7.json'))['cells']
allbad = alltot = 0
for key in sorted(pred):
    cell = pred[key]
    t = sp.Rational(cell['t'])
    rho2 = 2 / (t * t - 1)
    for s in (11, 13):
        K = (s + 1) // 2
        co = cell['odd_spins'][str(s)]['coefficients']
        if co is None:
            print('t = %-4s spin %d: operator gives no charge' % (cell['t'], s)); continue
        if any(de.eval(t) == 0 for _, de in TAB[s].values()):
            print('t = %-4s spin %d: table pole, UNSCREENED' % (cell['t'], s)); continue
        vals = {m: sp.Rational(nu.eval(t)) / sp.Rational(de.eval(t)) for m, (nu, de) in TAB[s].items()}
        top = vals.get((K, 0), 0)
        if top == 0:
            print('t = %-4s spin %d: table top vanishes, UNSCREENED' % (cell['t'], s)); continue
        vals = {m: v / top for m, v in vals.items()}
        if mode == 'tamper':
            vals[2, 1] = vals.get((2, 1), 0) + sp.Rational(1, 1000)
        P = XR(0)
        for kk, v in co.items():
            a, b = [int(x) for x in kk.split(',')]
            P += (X * QQ(-1)) ** (a // 2) * (XR(1) * q_(rho2) - Y) ** (b // 2) * q_(sp.Rational(v))
        P = P * (1 / P.coeff(X ** K))
        got = {m: sp.Rational(int(c.numerator), int(c.denominator)) for m, c in P.terms()}
        mons = sorted(set(vals) | set(got))
        bad = [m for m in mons if got.get(m, 0) != vals.get(m, 0)]
        allbad += len(bad); alltot += len(mons)
        print('t = %-4s spin %d: %d / %d mismatches' % (cell['t'], s, len(bad), len(mons))
              + ('' if not bad else '  first: X^%d Y^%d predicted %s certified %s' % (bad[0][0], bad[0][1], got.get(bad[0], 0), vals.get(bad[0], 0))), flush=True)
print('T2b mode %s: %d mismatches / %d coefficients' % (mode, allbad, alltot))
if mode == 'real':
    sys.exit(0 if allbad == 0 and alltot > 0 else 2)
sys.exit(0 if allbad > 0 else 1)
