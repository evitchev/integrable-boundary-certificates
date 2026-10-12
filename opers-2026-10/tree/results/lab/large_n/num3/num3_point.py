"""NUM3 data, ONE energy (per-energy checkpoint).  Usage: num3_point.py K POINT DPS INDEX [wrong]  ->  points/pt_<K>_<POINT>_<DPS>_<INDEX>[_wrong].json
Grid E = 10^(3 + 4 i/35), i = 0..35 (as NUM2).  Q_j for the three-term third-order Sol 2 ODE (sealed objects)."""
import sys, os, json, time, hashlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
assert hashlib.sha256(open(os.path.join(HERE, 'SEAL_NUM3.md'), 'rb').read()).hexdigest().startswith('b65ef9e8')
import mpmath as mp
from fractions import Fraction as Fr
from ops3 import sol2_three_term, fr2mp
kq = Fr(sys.argv[1]); pt = sys.argv[2]; dps = int(sys.argv[3]); i = int(sys.argv[4]); wrong = len(sys.argv) > 5 and sys.argv[5] == 'wrong'
mp.mp.dps = dps
PTS = {'A': (mp.mpf(1) / 3, mp.mpf(3) / 5), 'B': (mp.mpf(9) / 10, mp.mpf(1) / 5)}
l0, l1 = PTS[pt]; E = mp.mpf(10) ** (3 + mp.mpf(4) * i / 35)
tag = f'{str(kq).replace("/", "_")}_{pt}_{dps}_{i}' + ('_wrong' if wrong else '')
T0 = time.time(); op, th, F, rR = sol2_three_term(kq, l0, l1, E, wrong_order=wrong); n = fr2mp(F['n'])
lq, Q = op.Q(th, 2 * E ** (-1 / n))
os.makedirs(os.path.join(HERE, 'points'), exist_ok=True)
json.dump({'k': str(kq), 'point': pt, 'dps': dps, 'index': i, 'wrong_order': wrong, 'E': str(E), 'logabsQ': [str(mp.re(x)) for x in lq],
           'sign': [int(mp.sign(q)) for q in Q], 'thetas': [str(x) for x in th], 'seconds': time.time() - T0},
          open(os.path.join(HERE, 'points', f'pt_{tag}.json'), 'w'))
print(f'{tag}: E = {mp.nstr(E, 8)} log|Q| = {[mp.nstr(mp.re(x), 22) for x in lq]} signs {[int(mp.sign(q)) for q in Q]} ({time.time()-T0:.0f}s)', flush=True)
