"""NUM1 N1: M = 1, alpha = 0: numeric D(E) vs exact (2l+1) Gamma(l+1/2)/Gamma((2l+3-E)/4)."""
import sys, os, time, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
assert hashlib.sha256(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'SEAL_NUM1.md'), 'rb').read()).hexdigest().startswith('c1418ec1')
import mpmath as mp
from spec import *
LOG = open('run_n1.log', 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
for dps in (40,):
    mp.mp.dps = dps
    l = mp.mpf(1) / 3; worst = 0
    for e in (10, 100, 1000):
        T0 = time.time(); num = logD_neg(1, 0, l, e); ex = mp.log((2 * l + 1) * mp.gamma(l + mp.mpf(1) / 2) / mp.gamma((2 * l + 3 + e) / 4))
        say(f'dps {dps} M = 1 E = -{e}: logD num {mp.nstr(num, 30)}  exact {mp.nstr(ex, 30)}  diff {mp.nstr(num - ex, 5)} ({time.time()-T0:.1f}s)')
        worst = max(worst, abs(num - ex))
    for E in ('2.5', '7.25'):
        T0 = time.time(); num = logD_lin(1, 0, l, mp.mpf(E)); ex = (2 * l + 1) * mp.gamma(l + mp.mpf(1) / 2) / mp.gamma((2 * l + 3 - mp.mpf(E)) / 4)
        say(f'dps {dps} M = 1 E = +{E} (linear): D num {mp.nstr(num, 30)}  exact {mp.nstr(ex, 30)}  rel diff {mp.nstr((num - ex) / ex, 5)} ({time.time()-T0:.1f}s)')
