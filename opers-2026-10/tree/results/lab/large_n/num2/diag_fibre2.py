"""NUM2 post-hoc (labelled): t = 9/4 extended fits.  F3: F1 + E^(-13m/10), m = 1..3;  F4: F1 + E^(-13m/10 - 13j/22), m = 1..2, j = 0..J.
Reports grades K = 0, 2..8 vs WKB and the extra coefficients, per Q_j and for the Weyl pairs."""
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__))
import mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 60
LOG = open(os.path.join(HERE, 'run_diag_fibre2.log'), 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
k = Fr(10, 3); n = k + 4; pK = lambda K: (k + 1) * (1 - K) / n; fp = lambda q: mp.mpf(q.numerator) / q.denominator
for pt in ('A', 'B'):
    D = json.load(open(os.path.join(HERE, f'data_10_3_{pt}_60.json'))); grid = [mp.mpf(e) for e in D['grid']]
    P = json.load(open(os.path.join(HERE, 'wkb_pred.json')))['10/3'][pt]['pred']
    Y = [[mp.mpf(v[j]) for v in D['logabsQ']] for j in range(4)]
    def fit(vals, K2, extra):
        with mp.workdps(180):
            cols = [('K0', fp(pK(0))), ('log', None), ('const', 0)] + [(f'K{K}', fp(pK(K))) for K in range(2, K2 + 1)] + extra
            A = mp.matrix([[(mp.log(e) if p is None else e ** p) for _, p in cols] for e in grid]); x, res = mp.qr_solve(A, mp.matrix(vals))
        return {nm: x[i] for i, (nm, _) in enumerate(cols)}, res
    for vname, extra, K2a, K2b in (('F3: + E^(-13m/10), m=1..3', [(f'N{m}', -mp.mpf(13) * m / 10) for m in (1, 2, 3)], 24, 20),
                                   ('F4: + E^(-13m/10 - 13j/22), m=1..2, j=0..6', [(f'N{m}_{j}', -mp.mpf(13) * m / 10 - mp.mpf(13) * j / 22) for m in (1, 2) for j in range(7)], 18, 15)):
        say(f'\n=== t = 9/4 point {pt}: {vname}')
        for j in range(4):
            f1, r1 = fit(Y[j], K2a, extra); f2, r2 = fit(Y[j], K2b, extra)
            row = []
            for K in [0] + list(range(2, 9)):
                w = mp.mpf(P[str(K)][1]); e1 = f1[f'K{K}']
                if w == 0: row.append(f'K{K}: |fit| {mp.nstr(abs(e1), 2)}')
                else: row.append(f'K{K}: {mp.nstr(-mp.log10(abs(e1 - w) / abs(w)), 3)} d [stab {mp.nstr(-mp.log10(abs(e1 - f2[f"K{K}"]) / abs(w)), 3)}]')
            say(f'  Q{j} (theta = {mp.nstr(mp.mpf(D["thetas"][j]), 6)}): res {mp.nstr(r1, 2)}; ' + '; '.join(row))
            say('     extras: ' + '; '.join(f'{nm}: {mp.nstr(f1[nm], 14)} [stab {mp.nstr(-mp.log10(abs(f1[nm] - f2[nm]) / (abs(f1[nm]) + mp.mpf(10)**-40)), 3)}]' for nm, _ in extra[:4]))
