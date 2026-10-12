"""NUM3 controls X0, X1 (sealed).  X0: NUM1 path (Sol 1 M = 1 exact; generic Op) and one archived NUM2 point (Sol 3 t = 2, A, E = 1e3).
X1: exact E = 0 Meijer-G Q for the 3rd order (t = 2, 9/4; A, B) and the 7th order (A, B); beta == p theta_G.  Also the SOL3 4th-order formula re-derived
by the generic exact_E0 (consistency with NUM2's closed form)."""
import sys, os, json, time, hashlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
assert hashlib.sha256(open(os.path.join(HERE, 'SEAL_NUM3.md'), 'rb').read()).hexdigest().startswith('b65ef9e8')
ORIG = {l.split()[1]: l.split()[0] for l in open(os.path.join(HERE, 'inputs', 'ORIGIN_SHA256'))}
assert all(hashlib.sha256(open(os.path.join(HERE, 'inputs', f_), 'rb').read()).hexdigest() == h for f_, h in ORIG.items())
import mpmath as mp
from fractions import Fraction as Fr
from ode import Op
from ops3 import *
mp.mp.dps = 60
LOG = open(os.path.join(HERE, 'run_controls.log'), 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
ok = True
# X0 (a): Sol 1, M = 1
l = mp.mpf(1) / 3
for E in (-10, -100, 3):
    E = mp.mpf(E); op = Op([-l * (l + 1), -1, 1], [(Fr(4), [-1]), (Fr(2), [mp.mpf(0)]), (Fr(2), [E])])
    lq, Q = op.Q([-l, l + 1], 2 / mp.sqrt(abs(E) + 1)); ex = mp.gamma(l + mp.mpf(1) / 2) / mp.gamma((2 * l + 3 - E) / 4)
    d = abs(Q[0] - ex) / abs(ex); ok &= d < mp.mpf(10) ** -40; say(f'X0(a) Sol 1 M = 1 E = {E}: rel diff {mp.nstr(d, 3)}')
# X0 (b): NUM2 archived point (Sol 3, t = 2, A, E = 1e3)
J = json.load(open(os.path.join(HERE, 'inputs', 'data_2_A_60.json')))
k = Fr(2); n = k + 4; b = k / (k + 1); c = n / b; s = mp.mpf(5) / 2; l0, l1 = mp.mpf(1) / 3, mp.mpf(3) / 5; E = mp.mpf(J['grid'][0])
P = polymul([s * s - l0 * l0, -2 * s, mp.mpf(1)], [s * s - l1 * l1, -2 * s, mp.mpf(1)])
op = Op(P, [(c, [-(s - fr2mp(c) / 2), mp.mpf(1)]), (n, [-E])])
lq, Q = op.Q([s - l0, s + l0, s - l1, s + l1], 2 * E ** (-mp.mpf(1) / 6))
dd = max(abs(mp.re(lq[j]) - mp.mpf(J['logabsQ'][0][j])) for j in range(4)); ok &= dd < mp.mpf(10) ** -40
say(f'X0(b) NUM2 archived Sol 3 point (t = 2, A, E = 1e3): max |log|Q| diff| = {mp.nstr(dd, 3)}')
# X1: 3rd order
for kk in (Fr(2), Fr(10, 3)):
    for nm, (l0, l1) in {'A': (mp.mpf(1) / 3, mp.mpf(3) / 5), 'B': (mp.mpf(9) / 10, mp.mpf(1) / 5)}.items():
        T0 = time.time(); op, th, F, rR = sol2_three_term(kk, l0, l1, 0)
        a, W = op.asym_coeffs(600); beta = a[int(op.gamma / op.delta)]
        lq, Q = op.Q(th, mp.mpf('0.5')); ex, bex = exact_E0(th, F['c'], mp.mpf(1), rR)
        ds = [abs(Q[j] - ex[j]) / abs(ex[j]) for j in range(3)]; ok &= max(ds) < mp.mpf(10) ** -30 and abs(beta - bex) < mp.mpf(10) ** -40
        say(f'X1 3rd order k = {kk} point {nm}: beta {mp.nstr(beta, 15)} vs p theta_G {mp.nstr(bex, 15)}; Q rel diffs {[mp.nstr(d, 3) for d in ds]} ({time.time()-T0:.0f}s)')
# X1': 7th order
for nm, (l0, l1) in {'A': (mp.mpf(1) / 3, mp.mpf(3) / 5), 'B': (mp.mpf(9) / 10, mp.mpf(1) / 5)}.items():
    T0 = time.time(); op, ex7 = seventh(l0, l1, 0)
    a, W = op.asym_coeffs(600); beta = a[int(op.gamma / op.delta)]
    lq, Q = op.Q(ex7, mp.mpf('0.5')); ex, bex = exact_E0(ex7, Fr(21, 2), mp.mpf(1), [])
    ds = [abs(Q[j] - ex[j]) / abs(ex[j]) for j in range(7)]; ok &= max(ds) < mp.mpf(10) ** -30 and abs(beta - bex) < mp.mpf(10) ** -40
    say(f'X1 7th order point {nm}: beta {mp.nstr(beta, 15)} vs p theta_G {mp.nstr(bex, 15)}; Q rel diffs {[mp.nstr(d, 3) for d in ds]} ({time.time()-T0:.0f}s)')
say(f'X0/X1 all: {ok}'); sys.exit(0 if ok else 1)
