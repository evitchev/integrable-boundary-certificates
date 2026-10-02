"""NUM1 post-hoc diagnostic (labelled): at NON-integer M = 51/10 (point A: alpha = 1, l = 1/3), does log D(-e) contain a term e^(-(M+1)/2) OUTSIDE the
WKB family e^(mu(1-k))?  Fit F1: WKB family only (+ log e, const);  F2: F1 + e^(-(M+1)/2).  Compare fitted grades with WKB(M = 51/10) and report the extra term."""
import sys, os, time
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from spec import logD_neg
from rwkb import wkb, integral_coeff
mp.mp.dps = 60
LOG = open('run_diag_generic_M.log', 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
Mq = Fr(51, 10); M = mp.mpf(51) / 10; mu = (M + 1) / (2 * M); al, l = mp.mpf(1), mp.mpf(1) / 3; lam = l * (l + 1)
S = wkb(Mq, 24); pred = {k: -integral_coeff(S[k], Mq, al, lam)[0] for k in range(25) if k != 1}
NP = 40; grid = [mp.mpf(10) ** (3 + mp.mpf(4) * i / (NP - 1)) for i in range(NP)]
T0 = time.time(); data = [logD_neg(Mq, al, l, e) for e in grid]; say(f'data M = 51/10 point A: {NP} points ({time.time()-T0:.0f}s)')
def fit(extra, K2):
    with mp.workdps(180):
        rows = [[e ** mu, mp.log(e), 1] + [e ** (mu * (1 - k)) for k in range(2, K2 + 1)] + [f(e) for f in extra] for e in grid]
        x, res = mp.qr_solve(mp.matrix(rows), mp.matrix(data))
    return x, res
for nm, extra in (('F1 WKB family only', []), ('F2 + e^(-(M+1)/2)', [lambda e: e ** (-(M + 1) / 2)])):
    x1, r1 = fit(extra, 30); x2, r2 = fit(extra, 26)
    say(f'{nm}: lsq residual {mp.nstr(r1, 3)}')
    for k in range(2, 9):
        dg = -mp.log10(abs(x1[k + 1] - pred[k]) / abs(pred[k])); st = -mp.log10(abs(x1[k + 1] - x2[k + 1]) / abs(x1[k + 1]))
        say(f'   grade {k}: fit {mp.nstr(x1[k + 1], 15)}  WKB {mp.nstr(pred[k], 15)}  agree {mp.nstr(dg, 4)} (stability {mp.nstr(st, 4)})')
    if extra: say(f'   extra term e^(-(M+1)/2) coefficient: {mp.nstr(x1[len(x1) - 1], 25)} (stability {mp.nstr(-mp.log10(abs(x1[len(x1)-1] - x2[len(x2)-1]) / abs(x1[len(x1)-1])), 4)})')
