"""NUM3 post-hoc (labelled): t = 2 on the extended grid (E >= 1e4).  (i) the sealed basis WITHOUT the log E terms (registered refinement says R_8 = 0):
fit quality and grades vs WKB; (ii) the E^-3 cell: fitted coefficient (= lim_(k->2) C_8 + non-WKB m = 2) vs lim C_8 (finite_sym, PREDICTION_NUM3_logE.json)."""
import sys, os, json, glob
HERE = os.path.dirname(os.path.abspath(__file__))
import mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 60
LOG = open(os.path.join(HERE, 'run_diag_k8.log'), 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
k = Fr(2); n = 2 * k + 3; pK = lambda K: (k + 1) * (1 - K) / n; fp = lambda q: mp.mpf(q.numerator) / q.denominator
LP = json.load(open(os.path.join(HERE, 'PREDICTION_NUM3_logE.json')))
for pt in ('A', 'B'):
    D = json.load(open(os.path.join(HERE, f'data_2_{pt}_60.json')))
    ext = sorted([json.load(open(f)) for f in glob.glob(os.path.join(HERE, 'points_ext', f'pt_2_{pt}_60_*.json')) if not f.endswith('_wrong.json')], key=lambda p: mp.mpf(p['E']))
    keep = [i for i, e in enumerate(D['grid']) if mp.mpf(e) >= mp.mpf(10) ** 4 - 1]
    grid = [mp.mpf(D['grid'][i]) for i in keep] + [mp.mpf(p['E']) for p in ext]
    Ys = [[mp.mpf(D['logabsQ'][i][j]) for i in keep] + [mp.mpf(p['logabsQ'][j]) for p in ext] for j in range(3)]
    P = json.load(open(os.path.join(HERE, 'wkb_pred3.json')))['2'][pt]['pred']
    limC8 = mp.mpf(LP[pt]['8']['finite_sym'])
    for vname, extra in (('NO log terms', [(-mp.mpf(3) / 2, False), (-mp.mpf(9) / 2, False), (-mp.mpf(15) / 2, False)]),
                         ('with log terms (sealed)', [(-mp.mpf(3) / 2, False), (-mp.mpf(9) / 2, False), (-mp.mpf(15) / 2, False), (mp.mpf(-3), True), (mp.mpf(-6), True)])):
        for j in range(3):
            out = []
            for K2 in (24, 20):
                with mp.workdps(180):
                    cols = [(fp(pK(0)), False), None, (mp.mpf(0), False)] + [(fp(pK(K)), False) for K in range(2, K2 + 1)] + extra
                    A = mp.matrix([[(mp.log(e) if c is None else e ** c[0] * (mp.log(e) if c[1] else 1)) for c in cols] for e in grid])
                    x, res = mp.qr_solve(A, mp.matrix(Ys[j]))
                out.append((x, res))
            (x, res), (x2, _) = out
            dg = [f'K{K} {mp.nstr(-mp.log10(abs(x[K + 1] - mp.mpf(P[str(K)][1])) / abs(mp.mpf(P[str(K)][1]))), 3)}' for K in (2, 4, 6)]
            c8 = x[8 + 1]; st8 = -mp.log10(abs(x[9] - x2[9]) / abs(x[9]))
            say(f'{pt} {vname} Q{j}: res {mp.nstr(res, 2)}; {", ".join(dg)}; E^-3 coefficient {mp.nstr(c8, 16)} [stab {mp.nstr(st8, 3)}] vs lim C_8 {mp.nstr(limC8, 16)}: difference (non-WKB m = 2) {mp.nstr(c8 - limC8, 12)}')
