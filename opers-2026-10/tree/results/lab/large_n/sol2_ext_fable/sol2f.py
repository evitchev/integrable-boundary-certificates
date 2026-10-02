"""SOL2F -- the sealed zero-parameter Sol 2 family S2(k) (= SOL12 class G, beta = 1, with M = 1/k, C = 2n(1+M); dual M' = -1/(k+1)) at
k = 2 (regression), 4, 3, 10/3, 1: monic spins 1..9 vs certified (vev_sol2_w6/8/10, cyl(2,4), spin-1 law); missing spins observed; tamper delta + 1/10."""
import sys, time, json
from fractions import Fraction as F
import sympy as sp
from common import *
LOG = open(os.path.join(HERE, 'run_sol2f.log'), 'w'); T00 = time.time()
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
FAIL = []; RES = {}; SEALED_FAIL = []
def tk(k): return sp.Rational(3 * k + 2, 1) / (k + 2)
def run(k, M, label, K=10, delta_shift=0):
    tv = (3 * sp.Rational(k) + 2) / (sp.Rational(k) + 2); g = classG(2, tv)
    assert g['alpha'] == F(k) and g['gamma'] == 1 and g['D'] == 2 * F(k) + 3, g
    n = g['D']; Cv = 2 * n * (1 + F(M))
    bl = blocks_of(g, M)
    if delta_shift: bl = [(a, sh, d + delta_shift) if a != 1 else (a, sh, d) for (a, sh, d) in bl]
    sv = gwkb.wkb(bl, g['s'], F(M), K); bb = n * F(M); ch = {}
    for kk in range(2, K + 1, 2):
        sk = sv[kk]
        if any(a.denominator == 1 and a <= 0 for (a, c) in sk): pw, poly, nonpole = thwkb.integrate_reg(sk, bb); poly = dict(poly)
        else:
            r = thwkb.integrate(sk, bb); poly = {k_: sp.Rational(v.numerator, v.denominator) for k_, v in r[2].items()} if r else {}
        ch[kk - 1] = sp.expand(sum(v * (sp.Rational(Cv.numerator, Cv.denominator) * PX**2 / g['aY'])**(i // 2) * (sp.Rational(Cv.numerator, Cv.denominator) * PI**2 / g['aX'])**(j // 2) for (i, j), v in poly.items()))
    return tv, g, Cv, ch
def score(ch, cert, spins):
    rep = {}
    for s_ in spins:
        if ch[s_] == 0: rep[s_] = 'WKB charge identically 0'; continue
        d = sp.expand(monic(ch[s_], (s_ + 1) // 2) - cert[s_]); rep[s_] = True if d == 0 else f'DIFF ({len(sp.Poly(d, PX, PI).terms())} coeffs)'
    return rep
for k, extra in ((2, 'regression = SOL12 t = 2'), (4, ''), (3, ''), (F(10, 3), 'Gamma X block'), (1, 'coincides with Sol 3 k = 1')):
    k = F(k); tv = (3 * sp.Rational(k.numerator, k.denominator) + 2) / (sp.Rational(k.numerator, k.denominator) + 2)
    if not regular_at(2, tv): say(f'k = {k} t = {tv}: Sol 2 table pole -> fibre DROPPED (pre-registered)'); continue
    cert = certified(2, tv); RES[str(k)] = {}
    for M, nm in ((1 / k, 'M = 1/k'), (-1 / (k + 1), "dual M' = -1/(k+1)")):
        T0 = time.time(); tv_, g, Cv, ch = run(k, M, nm)
        rep = score(ch, cert, (1, 3, 5, 7, 9))
        say(f'k = {k} (t = {tv}, n = {g["D"]}; {extra}) {nm} = {M}, C = {Cv}: {rep}  ({time.time()-T0:.0f}s)')
        RES[str(k)][nm] = {str(s_): (v if v is True else str(v)) for s_, v in rep.items()}
        if nm == 'M = 1/k': RES[str(k)]['charges_M'] = {str(s_): str(v) for s_, v in ch.items()}
    if k == 1:
        c3 = certified(3, tv)
        say(f'   k = 1 (t = 5/3): certified Sol 2 == certified Sol 3 per spin: { {s_: sp.expand(cert[s_] - c3[s_]) == 0 for s_ in (1, 3, 5, 7, 9)} }')
    if k == 2: 
        ok = all(v is True for s_, v in RES['2']['M = 1/k'].items() if s_ != '7')
        if not ok: FAIL.append('regression k = 2')
    if k.denominator == 1 and k > 1:
        _, _, _, cht = run(k, 1 / k, 'tamper', delta_shift=F(1, 10)); rt = score(cht, cert, (3, 5, 7, 9))
        say(f'   TAMPER delta + 1/10 at k = {k}: {rt}')
        if any(v is True for v in rt.values()): FAIL.append(f'tamper k = {k} did not fail')
json.dump(RES, open(os.path.join(HERE, 'sol2f.json'), 'w'), indent=1)
say(f'done ({time.time()-T00:.0f}s)')
if FAIL: say(f'CONTROLS FAILED {FAIL} -> exit 1'); sys.exit(1)
sys.exit(0)
