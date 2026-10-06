"""NUM2 data, ONE energy (v2 after the v1 runs were killed by their own 'timeout 20000'; v1 logs in v1_timeout/): checkpoint per energy.
Usage: num2_point.py K POINT DPS INDEX  -> points/pt_<K>_<POINT>_<DPS>_<INDEX>.json (same grid, same method as num2_data.py)."""
import sys, os, json, time, hashlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
assert hashlib.sha256(open(os.path.join(HERE, 'SEAL_NUM2.md'), 'rb').read()).hexdigest().startswith('25ab6266')
import mpmath as mp
from fractions import Fraction as Fr
from ops import three_term
kq = Fr(sys.argv[1]); pt = sys.argv[2]; dps = int(sys.argv[3]); i = int(sys.argv[4])
mp.mp.dps = dps
PTS = {'A': (mp.mpf(1) / 3, mp.mpf(3) / 5), 'B': (mp.mpf(9) / 10, mp.mpf(1) / 5), 'C': (mp.mpf(1) / 2, mp.mpf(1) / 4)}
l0, l1 = PTS[pt]; NP = 36; E = mp.mpf(10) ** (3 + mp.mpf(4) * i / (NP - 1))
tag = f'{str(kq).replace("/", "_")}_{pt}_{dps}_{i}'
T0 = time.time(); op, th, F = three_term(kq, l0, l1, E); n = mp.mpf(F['n'].numerator) / F['n'].denominator
lq, Q = op.Q(th, 2 * E ** (-1 / n))
os.makedirs(os.path.join(HERE, 'points'), exist_ok=True)
json.dump({'k': str(kq), 'point': pt, 'dps': dps, 'index': i, 'E': str(E), 'logabsQ': [str(mp.re(x)) for x in lq], 'sign': [int(mp.sign(q)) for q in Q],
           'thetas': [str(x) for x in th], 'seconds': time.time() - T0}, open(os.path.join(HERE, 'points', f'pt_{tag}.json'), 'w'))
print(f'{tag}: E = {mp.nstr(E, 8)} log|Q| = {[mp.nstr(mp.re(x), 22) for x in lq]} ({time.time()-T0:.0f}s)', flush=True)
