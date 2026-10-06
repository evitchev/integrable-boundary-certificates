"""NUM2 control X0: the generic order-D solver at order 2 (Sol 1, theta-form: (theta^2 - theta - lam) y = x^2 (x^(2M) + alpha x^(M-1) - E) y)
(a) M = 1, alpha = 0: Q_(-l) = Gamma(l + 1/2)/Gamma((2l+3-E)/4) exactly (NUM1's D/(2l+1));  (b) M = 5, point A: log Q_(-l) + log(2l+1) == NUM1's archived log D."""
import sys, os, json, time, hashlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
assert hashlib.sha256(open(os.path.join(HERE, 'SEAL_NUM2.md'), 'rb').read()).hexdigest().startswith('25ab6266')
import mpmath as mp
from fractions import Fraction as Fr
from ode import Op
mp.mp.dps = 60
LOG = open(os.path.join(HERE, 'run_x0.log'), 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
def sol1(M, alpha, l, E):
    lam = l * (l + 1)
    return Op([-lam, -1, 1], [(Fr(2 * M + 2), [-1]), (Fr(M + 1), [-alpha]), (Fr(2), [mp.mpf(E)])])
l = mp.mpf(1) / 3; worst = mp.mpf(0)
for E in (-10, -100, -1000, 3, '7.5'):
    E = mp.mpf(E); T0 = time.time()
    op = sol1(1, mp.mpf(0), l, E); x0 = mp.mpf(2) / mp.sqrt(abs(E) + 1)
    logQ, Q = op.Q([-l, l + 1], x0)
    ex = mp.gamma(l + mp.mpf(1) / 2) / mp.gamma((2 * l + 3 - E) / 4)
    d = abs(Q[0] - ex) / abs(ex); worst = max(worst, d)
    say(f'(a) M = 1 E = {E}: Q_-l num {mp.nstr(Q[0], 25)} exact {mp.nstr(ex, 25)} rel diff {mp.nstr(d, 3)} ({time.time()-T0:.1f}s)')
say(f'(a) worst {mp.nstr(worst, 3)} -> {int(-mp.log10(worst))} digits (sealed >= 40): {worst < mp.mpf(10) ** -40}')
J = json.load(open(os.path.join(HERE, 'inputs', 'num1_dps60.json')))
worst = mp.mpf(0)
for i in (0, 17, 35):
    e = mp.mpf(J['grid'][i]); T0 = time.time()
    op = sol1(5, mp.mpf(1), l, -e); logQ, Q = op.Q([-l, l + 1], 2 / mp.sqrt(e))
    v = logQ[0] + mp.log(2 * l + 1); ref = mp.mpf(J['data']['A'][i]); d = abs(v - ref); worst = max(worst, d)
    say(f'(b) M = 5 alpha = 1 e = {mp.nstr(e, 6)}: log D generic {mp.nstr(v, 30)} NUM1 {mp.nstr(ref, 30)} diff {mp.nstr(d, 3)} ({time.time()-T0:.1f}s)')
say(f'(b) worst {mp.nstr(worst, 3)} (sealed <= 1e-40): {worst < mp.mpf(10) ** -40}')
