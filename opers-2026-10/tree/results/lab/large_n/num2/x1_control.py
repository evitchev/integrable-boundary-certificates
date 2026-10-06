"""NUM2 control X1 / X1': exact E = 0 Q-functions (Meijer G) vs the solver: three-term (t = 2 and t = 9/4) and sixth-order (t = 2); beta == c theta_G."""
import sys, os, time, hashlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
assert hashlib.sha256(open(os.path.join(HERE, 'SEAL_NUM2.md'), 'rb').read()).hexdigest().startswith('25ab6266')
import mpmath as mp
from ops import *
mp.mp.dps = 60
LOG = open(os.path.join(HERE, 'run_x1.log'), 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
PTS = {'A': (mp.mpf(1) / 3, mp.mpf(3) / 5), 'B': (mp.mpf(9) / 10, mp.mpf(1) / 5)}
ok = True
for k in (2, Fr(10, 3)):
    for nm, (l0, l1) in PTS.items():
        T0 = time.time(); op, th, F = three_term(k, l0, l1, 0)
        a, W = op.asym_coeffs(600)
        # beta = coefficient of x^0 in S
        i0 = int(op.gamma / op.delta); beta = a[i0]
        logQ, Q = op.Q(th, mp.mpf('0.5'))
        ex, beta_ex = exact_three_term_E0(k, l0, l1)
        ds = [abs(Q[j] - ex[j]) / abs(ex[j]) for j in range(4)]
        say(f'three-term k = {k} point {nm}: beta {mp.nstr(beta, 20)} vs c theta_G {mp.nstr(beta_ex, 20)} (diff {mp.nstr(abs(beta - beta_ex), 3)}); Q_j rel diffs {[mp.nstr(d, 3) for d in ds]} ({time.time()-T0:.0f}s)')
        ok &= max(ds) < mp.mpf(10) ** -30 and abs(beta - beta_ex) < mp.mpf(10) ** -40
for nm, (l0, l1) in PTS.items():
    T0 = time.time(); op, ex6 = sixth(l0, l1, 0)
    a, W = op.asym_coeffs(600); i0 = int(op.gamma / op.delta); beta = a[i0]
    logQ, Q = op.Q(ex6, mp.mpf('0.5'))
    ex, beta_ex = exact_sixth_E0(l0, l1)
    ds = [abs(Q[j] - ex[j]) / abs(ex[j]) for j in range(6)]
    say(f'sixth-order point {nm}: beta {mp.nstr(beta, 20)} vs 9 theta_G {mp.nstr(beta_ex, 20)}; Q_j rel diffs {[mp.nstr(d, 3) for d in ds]} ({time.time()-T0:.0f}s)')
    ok &= max(ds) < mp.mpf(10) ** -30 and abs(beta - beta_ex) < mp.mpf(10) ** -40
say(f'X1 / X1\' (sealed: >= 30 digits, beta exact): {ok}')
sys.exit(0 if ok else 1)
