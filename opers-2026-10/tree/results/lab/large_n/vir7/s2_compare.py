"""VIR7 T2a (not blind): the Sol-1 operator S1(t) against (A) Codex's Sol-1 Ward amplitudes, losses 1..4, and
(B) the certified tables vev_sol1_w6/8/10 (spins 5, 7, 9), at every fibre with a prediction file.
Modes: real | tamper_scale (l1 and l0 scales x 11/10 at conversion) | file:<prediction json> (e.g. the tamper_a run).
Exit real: 0 = all agree; 2 = mismatch.  Tamper modes: 0 = fires."""
import glob, hashlib, json, sys
import sympy as sp
sys.dont_write_bytecode = True
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
assert hashlib.sha256(open('SEAL_VIR7_addendum.md', 'rb').read()).hexdigest() == open('SEAL_VIR7_addendum.sha256').read().split()[0]
mode = sys.argv[1] if len(sys.argv) > 1 else 'real'
ROOT = '<repo>/results/lab/'
AUD = ROOT + 'audits/'
SRC = {1: AUD + 'codex_wardsol12b_27403fc_2026-09-28/AMPLITUDES_sol1_loss1.json',
       2: AUD + 'codex_wardsol12_582806a_2026-09-28/AMPLITUDES_sol1_loss2.json',
       3: AUD + 'codex_wardsol12_582806a_2026-09-28/AMPLITUDES_sol1_loss3.json',
       4: AUD + 'codex_wardsol12_582806a_2026-09-28/AMPLITUDES_sol1_loss4.json'}
ks, ts = sp.symbols('k t')
LAW = {}
PQR = None
for m, fn in SRC.items():
    law = json.load(open(fn))
    assert law['solution'] == 1 and law['loss'] == m
    f = lambda d: (sp.Poly(sp.sympify(d['numerator'], locals={'k': ks, 't': ts}, rational=True), ks, ts), sp.Poly(sp.sympify(d['denominator'], locals={'k': ks, 't': ts}, rational=True), ks, ts))
    LAW[m] = {int(j): f(v) for j, v in law['shift_amplitudes'].items()}
    PQR = [f(law[x]) for x in 'pqr']
XR, X, Y = ring('X,Y', QQ)


def q_(x):
    x = sp.Rational(x)
    return QQ(int(x.p), int(x.q))


def ev(pair, K, t):
    return sp.Rational(pair[0].eval({ks: K, ts: t})) / sp.Rational(pair[1].eval({ks: K, ts: t}))


def Kd(d, p, q, r, a):
    return sp.binomial(d, a) * sp.rf(p, a) * sp.rf(q, d - a) * r ** (d - a) / sp.rf(p, d)


def record(co, K, rho2, scale=1):
    P = XR(0)
    for kk, v in co.items():
        a, b = [int(x) for x in kk.split(',')]
        P += (X * q_(-scale)) ** (a // 2) * ((XR(1) * q_(rho2) - Y) * q_(scale)) ** (b // 2) * q_(sp.Rational(v)) * q_(sp.Rational(1) / scale ** ((a + b) // 2))
    lead = P.coeff(X ** K)
    P = P * (1 / lead)
    return {m: sp.Rational(int(c.numerator), int(c.denominator)) for m, c in P.terms()}


def table(s, tv):
    tab = json.load(open(ROOT + 'vev/vev_sol1_w%d.json' % (s + 1)))
    out = {}
    for key, val in tab['coefficients'].items():
        a_, b_ = [int(x.split('^')[1]) for x in key.split()]
        num, den = sp.fraction(sp.together(sp.sympify(val, locals={'t': ts}, rational=True)))
        if den.subs(ts, tv) == 0:
            return None
        out[(a_ // 2, b_ // 2)] = sp.Rational(num.subs(ts, tv)) / sp.Rational(den.subs(ts, tv))
    top = out[((s + 1) // 2, 0)]
    if top == 0:
        return None
    return {m: v / top for m, v in out.items()}


files = sorted(glob.glob('s1_pred_t*.json'))
if mode.startswith('file:'):
    files = [mode[5:]]
else:
    files = [f for f in files if 'tamper' not in f]
tot = bad = 0
for fn in files:
    cell = json.load(open(fn))
    t = sp.Rational(cell['t'])
    rho2 = 2 / (t * t - 1)
    n = sp.Rational(cell['n'])
    print('== t = %s (a = %s, n = %s, M = %s) %s' % (cell['t'], cell['a'], cell['n'], cell['M'], cell['mode']))
    for s in sorted(int(x) for x in cell['odd_spins']):
        if s > 9:
            continue                                   # spins 11, 13 are the blind hold-out (T2b)
        K = (s + 1) // 2
        co = cell['odd_spins'][str(s)]['coefficients']
        if co is None:
            print('   spin %d: operator gives no charge (nu = %s)' % (s, sp.Rational(s) / n)); continue
        if mode == 'tamper_scale':
            # l^2 -> (11/10)^2 l^2: PX^2, pi^2 monomials of degree d get (121/100)^d before normalisation
            co = {kk: str(sp.Rational(v) * sp.Rational(121, 100) ** ((int(kk.split(',')[0]) + int(kk.split(',')[1])) // 2)) for kk, v in co.items()}
        got = record(co, K, rho2)
        rows = []
        if K >= 2:
            p, q, r = [ev(x, K, t) for x in PQR]
            for m in range(1, 5):
                if m > K:
                    continue
                d = K - m
                try:
                    vals = {j: ev(v, K, t) for j, v in LAW[m].items()}
                    fb = 0
                    for a in range(d + 1):
                        expect = sum(v * Kd(d, p + j, q, r, a) for j, v in vals.items())
                        fb += got.get((a, d - a), 0) != expect
                    rows.append('loss %d: %d/%d' % (m, fb, d + 1)); tot += d + 1; bad += fb
                except ZeroDivisionError:
                    rows.append('loss %d: pole' % m)
        line = '   spin %d  vs Ward amplitudes: %s' % (s, '; '.join(rows) if rows else '-')
        if s == 1:
            ok1 = got == {(1, 0): 1, (0, 1): 1, (0, 0): sp.Rational(1, 6) - rho2}
            line += '  spin 1 == X + Y + 1/6 - rho^2: %s' % ok1; tot += 1; bad += not ok1
        if s in (5, 7, 9):
            tb = table(s, t)
            if tb is None:
                line += ' | table: pole/zero top, unscreened'
            else:
                mons = set(tb) | set(got)
                fbad = [mm for mm in mons if got.get(mm, 0) != tb.get(mm, 0)]
                line += ' | full table: %d/%d mismatches' % (len(fbad), len(mons)); tot += len(mons); bad += len(fbad)
                if fbad:
                    mm = sorted(fbad, key=lambda z: -(z[0] + z[1]))[0]
                    line += ' (first, highest degree: X^%d Y^%d predicted %s certified %s)' % (mm[0], mm[1], got.get(mm, 0), tb.get(mm, 0))
        print(line, flush=True)
    ev_sp = {s_: (v['identically_zero'], None if v['coefficients'] is None else len(v['coefficients'])) for s_, v in sorted(cell['even_spins'].items(), key=lambda z: int(z[0]))}
    print('   even spins (identically zero?, #monomials of R/PX): %s' % ev_sp)
print('T2a mode %s: %d mismatches / %d coefficients' % (mode, bad, tot))
if mode == 'real':
    sys.exit(0 if bad == 0 and tot > 0 else 2)
sys.exit(0 if bad > 0 else 1)
