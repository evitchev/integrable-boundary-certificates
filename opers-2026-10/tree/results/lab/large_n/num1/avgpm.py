"""NUM1 (lead's added test, post-hoc w.r.t. the seal): at A, B compute log D for -alpha (same l), blind-fit both signs, and compare the +-alpha AVERAGE
of each grade with the alpha-even WKB (average of WKB(+alpha), WKB(-alpha) at M = 5 exactly, no limit, no Mellin term)."""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from spec import logD_neg
from rwkb import wkb, integral_coeff
mp.mp.dps = 60
LOG = open('run_avgpm.log', 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
J = json.load(open('num1_dps60.json')); grid = [mp.mpf(g) for g in J['grid']]; NP = len(grid); mu = mp.mpf(3) / 5
S = wkb(5, 24)
def fit(vals, K2):
    with mp.workdps(180):
        rows = [[e ** mu, mp.log(e), 1] + [e ** (mu * (1 - k)) for k in range(2, K2 + 1)] for e in grid]
        x, _ = mp.qr_solve(mp.matrix(rows), mp.matrix(vals))
    out = {0: x[0], 'log': x[1], 'const': x[2]}
    for k in range(2, K2 + 1): out[k] = x[k + 1]
    return out
OUT = {}
for nm, (al, l) in {'A': (mp.mpf(1), mp.mpf(1) / 3), 'B': (mp.mpf(-3) / 2, mp.mpf(7) / 10)}.items():
    lam = l * (l + 1)
    T0 = time.time(); dm = [logD_neg(5, -al, l, e) for e in grid]; say(f'{nm}: -alpha data ({time.time()-T0:.0f}s)')
    dp = [mp.mpf(v) for v in J['data'][nm]]
    fp, fm = fit(dp, 30), fit(dm, 30); fp2, fm2 = fit(dp, 26), fit(dm, 26)
    OUT[nm] = {'minus_alpha_data': [str(v) for v in dm]}
    for k in [0] + list(range(2, 9)):
        wp = -integral_coeff(S[k], 5, al, lam)[0]; wm = -integral_coeff(S[k], 5, -al, lam)[0]
        avg_ext = (fp[k] + fm[k]) / 2; avg_w = (wp + wm) / 2; st = abs((fp[k] + fm[k]) / 2 - (fp2[k] + fm2[k]) / 2)
        if avg_w != 0: dg = -mp.log10(abs(avg_ext - avg_w) / abs(avg_w)); dgs = f'{mp.nstr(dg, 4)} digits'
        else: dgs = f'WKB avg 0, extracted avg {mp.nstr(avg_ext, 3)}'
        say(f'   grade {k}: +-alpha average extracted {mp.nstr(avg_ext, 20)}  alpha-even WKB (M = 5 exactly) {mp.nstr(avg_w, 20)}  {dgs} (fit stability {mp.nstr(-mp.log10(st / (abs(avg_ext) + mp.mpf(10)**-60)), 4)})')
    say(f'   grade-1 log e coefficient: +alpha {mp.nstr(fp["log"], 15)}, -alpha {mp.nstr(fm["log"], 15)}')
json.dump(OUT, open('avgpm.json', 'w'), indent=1)
