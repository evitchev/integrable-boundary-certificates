"""NUM1 N2: point A (M = 5, alpha = 1, l = 1/3): the lowest 3 zeros of D(E) (W[y, psi_+], inward) vs eigenvalues by outward shooting of psi_+ (bisection on
the sign of psi_+(x_far)).  dps 40."""
import sys, os, time, hashlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
assert hashlib.sha256(open(os.path.join(HERE, 'SEAL_NUM1.md'), 'rb').read()).hexdigest().startswith('c1418ec1')
import mpmath as mp
from spec import *
mp.mp.dps = 40
LOG = open(os.path.join(HERE, 'run_n2.log'), 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
M, al, l = 5, mp.mpf(1), mp.mpf(1) / 3
def bisect(f, a, b, tol):
    fa = f(a)
    while b - a > tol:
        m = (a + b) / 2; fm = f(m)
        if (fm > 0) == (fa > 0): a, fa = m, fm
        else: b = m
    return (a + b) / 2
Es = [mp.mpf(k) / 4 for k in range(0, 241)]   # v2: [0, 60] (v1 [0, 30] held only 2 levels; log kept)
T0 = time.time(); sh = [shoot(M, al, l, E) for E in Es]; say(f'shooting scan E in [0, 60] step 1/4 ({time.time()-T0:.0f}s)')
br = [(Es[i], Es[i + 1]) for i in range(len(Es) - 1) if (sh[i] > 0) != (sh[i + 1] > 0)][:3]
say(f'brackets: {[(float(a), float(b)) for a, b in br]}')
worst = 0
for a, b in br:
    Es_ = bisect(lambda E: shoot(M, al, l, E), a, b, mp.mpf(10) ** -30)
    Ed = bisect(lambda E: D_lin(M, al, l, E), a, b, mp.mpf(10) ** -30)
    d = abs(Es_ - Ed); worst = max(worst, d)
    say(f'eigenvalue: shooting {mp.nstr(Es_, 28)}  zero of D {mp.nstr(Ed, 28)}  |diff| {mp.nstr(d, 3)}')
say(f'N2: lowest {len(br)} levels agree to {int(-mp.log10(worst)) if worst else 99} digits (sealed >= 10): {len(br) == 3 and worst < mp.mpf(10) ** -10}')
