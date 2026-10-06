"""NUM2 data: log|Q_j(E)| (and signs) on the sealed energy grid E in [1e3, 1e7] (36 points).  Usage: num2_data.py K POINT DPS [wrong]
K = 2 or 10/3 (fibre); POINT = A, B, C.  Writes data_<K>_<POINT>_<DPS>[_wrong].json and a full log."""
import sys, os, json, time, hashlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
assert hashlib.sha256(open(os.path.join(HERE, 'SEAL_NUM2.md'), 'rb').read()).hexdigest().startswith('25ab6266')
import mpmath as mp
from fractions import Fraction as Fr
from ops import three_term
kq = Fr(sys.argv[1]); pt = sys.argv[2]; dps = int(sys.argv[3]); wrong = len(sys.argv) > 4 and sys.argv[4] == 'wrong'
mp.mp.dps = dps
PTS = {'A': (mp.mpf(1) / 3, mp.mpf(3) / 5), 'B': (mp.mpf(9) / 10, mp.mpf(1) / 5), 'C': (mp.mpf(1) / 2, mp.mpf(1) / 4)}
l0, l1 = PTS[pt]
tag = f'{str(kq).replace("/", "_")}_{pt}_{dps}' + ('_wrong' if wrong else '')
LOG = open(os.path.join(HERE, f'run_data_{tag}.log'), 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
NP = 36; grid = [mp.mpf(10) ** (3 + mp.mpf(4) * i / (NP - 1)) for i in range(NP)]
out = {'k': str(kq), 'point': pt, 'l0': str(l0), 'l1': str(l1), 'dps': dps, 'wrong_order': wrong, 'grid': [str(e) for e in grid], 'logabsQ': [], 'sign': []}
T00 = time.time()
for i, E in enumerate(grid):
    op, th, F = three_term(kq, l0, l1, E, wrong_order=wrong)
    n = mp.mpf(F['n'].numerator) / F['n'].denominator
    T0 = time.time(); lq, Q = op.Q(th, 2 * E ** (-1 / n))
    out['logabsQ'].append([str(mp.re(x)) for x in lq]); out['sign'].append([int(mp.sign(q)) for q in Q])
    say(f'{i:2d} E = {mp.nstr(E, 8)}: log|Q| = {[mp.nstr(mp.re(x), 22) for x in lq]} signs {[int(mp.sign(q)) for q in Q]} ({time.time()-T0:.1f}s)')
out['thetas'] = [str(x) for x in th]
json.dump(out, open(os.path.join(HERE, f'data_{tag}.json'), 'w'), indent=1)
say(f'done ({time.time()-T00:.0f}s)')
