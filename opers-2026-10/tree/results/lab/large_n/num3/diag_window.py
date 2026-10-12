"""NUM3 post-hoc diagnostic (labelled): does the SEALED-basis fit improve when low energies are dropped (truncation of a slowly converging lattice)?
Fit on grid points with E >= Emin for Emin = 1e3, 1e4, 1e5; report K0, K2, K4, K6 digits vs WKB and N-3/2 stability."""
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__))
import mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 60
tag = sys.argv[1]
D = json.load(open(os.path.join(HERE, f'data_{tag}.json'))); k = Fr(D['k']); n = 2 * k + 3; grid = [mp.mpf(e) for e in D['grid']]
P = json.load(open(os.path.join(HERE, 'wkb_pred3.json')))[str(k)][D['point']]['pred']
pK = lambda K: (k + 1) * (1 - K) / n; fp = lambda q: mp.mpf(q.numerator) / q.denominator
if k == 2: extra = [(-mp.mpf(3) / 2, False), (-mp.mpf(9) / 2, False), (mp.mpf(-3), True)]
else: extra = [(-mp.mpf(13) * m / 10, False) for m in (1, 2)]
y = [mp.mpf(v[0]) for v in D['logabsQ']]
for Emin in (3, 4, 5):
    idx = [i for i, e in enumerate(grid) if e >= mp.mpf(10) ** Emin - 1]
    for K2 in (16, 20):
        with mp.workdps(180):
            cols = [fp(pK(0)), None, mp.mpf(0)] + [fp(pK(K)) for K in range(2, K2 + 1)]
            rows = [[(mp.log(grid[i]) if c is None else grid[i] ** c) for c in cols] + [grid[i] ** p * (mp.log(grid[i]) if lg else 1) for p, lg in extra] for i in idx]
            x, res = mp.qr_solve(mp.matrix(rows), mp.matrix([y[i] for i in idx]))
        dg = []
        for K in (0, 2, 4, 6):
            w = mp.mpf(P[str(K)][1]); e = x[0] if K == 0 else x[K + 1]
            dg.append(f'K{K} {mp.nstr(-mp.log10(abs(e - w) / abs(w)), 3)}')
        print(f'{tag} Emin 1e{Emin} ({len(idx)} pts) K2 {K2}: res {mp.nstr(res, 2)}; ' + ', '.join(dg) + f'; N1 {mp.nstr(x[len(cols)], 12)}', flush=True)
