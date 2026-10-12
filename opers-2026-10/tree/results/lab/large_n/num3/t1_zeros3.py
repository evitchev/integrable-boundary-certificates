"""NUM3 (sealed T1): Q3_j(E)/Q7_j(E_Gamma = -E) == Q3_j(0)/Q7_j(0) (exact Meijer G) for j = s, s - l1, s + l1  <->  7th-order exponents 3, 3 - l1, 3 + l1,
at E in {-20, -5, 2, 50, 500}; the four l0-exponent Q7 reported.  Zeros of Q_(s - l1) at A on E in [-300, 0]: projection vs outward Wronskian."""
import sys, os, json, time, hashlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
assert hashlib.sha256(open(os.path.join(HERE, 'SEAL_NUM3.md'), 'rb').read()).hexdigest().startswith('b65ef9e8')
import mpmath as mp
from fractions import Fraction as Fr
from ops3 import *
mp.mp.dps = 60
LOG = open(os.path.join(HERE, 'run_t1_zeros3.log'), 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
PTS = {'A': (mp.mpf(1) / 3, mp.mpf(3) / 5), 'B': (mp.mpf(9) / 10, mp.mpf(1) / 5)}
MAP = {0: 0, 1: 5, 2: 6}      # three-term [s, s - l1, s + l1]  ->  seventh [3, ..., 3 - l1, 3 + l1]
OUT = {}; worst = mp.mpf(0)
for nm, (l0, l1) in PTS.items():
    op0, th0, F, rR = sol2_three_term(2, l0, l1, 0); e3, _ = exact_E0(th0, F['c'], mp.mpf(1), rR)
    op70, ex7 = seventh(l0, l1, 0); e7, _ = exact_E0(ex7, Fr(21, 2), mp.mpf(1), [])
    C = [e3[j] / e7[MAP[j]] for j in range(3)]; OUT[nm] = {}
    for E in (-20, -5, 2, 50, 500):
        E = mp.mpf(E); x0 = 2 / (abs(E) + 1) ** (mp.mpf(1) / 7); T0 = time.time()
        op3, th3, F, rR = sol2_three_term(2, l0, l1, E); lq3, Q3 = op3.Q(th3, x0)
        op7, ex7 = seventh(l0, l1, -E); lq7, Q7 = op7.Q(ex7, x0)
        dev = [abs(Q3[j] / Q7[MAP[j]] / C[j] - 1) for j in range(3)]; worst = max(worst, max(dev))
        OUT[nm][str(E)] = {'Q3': [str(q) for q in Q3], 'Q7': [str(q) for q in Q7], 'dev': [str(d) for d in dev]}
        say(f'T1 point {nm} E = {E}: |Q3_j/Q7_j/C_j - 1| = {[mp.nstr(d, 3) for d in dev]}; l0-exponent Q7 = {[mp.nstr(Q7[i], 10) for i in (1, 2, 3, 4)]} ({time.time()-T0:.0f}s)')
say(f'T1 (sealed >= 25 digits): worst {mp.nstr(worst, 3)} -> {worst < mp.mpf(10) ** -25}')
l0, l1 = PTS['A']; jd = 1
def Qd(E):
    op, th, F, rR = sol2_three_term(2, l0, l1, mp.mpf(E)); lq, Q = op.Q(th, mp.mpf('0.6')); return Q[jd]
def Wd(E, x1=mp.mpf(3) / 2):
    E = mp.mpf(E); op, th, F, rR = sol2_three_term(2, l0, l1, E); x0 = mp.mpf('0.6')
    xm = op.xmax(mp.mpf(2)); L, Wv = op.asym_eval(xm)
    Yy, _ = op.integrate(mp.log(xm), [Wv[m] / mp.factorial(m) for m in range(3)], mp.log(x1)); cols = [[Yy[m] * mp.factorial(m) for m in range(3)]]
    for i, t in enumerate(th):
        if i == jd: continue
        v = op.frob(t, x0); Yc, _ = op.integrate(mp.log(x0), [v[m] / mp.factorial(m) for m in range(3)], mp.log(x1)); cols.append([Yc[m] * mp.factorial(m) for m in range(3)])
    return mp.det(mp.matrix(cols).T)
Es = [-mp.mpf(5) * i for i in range(61)]; T0 = time.time(); vals = [Qd(E) for E in Es]
say(f'zeros scan Q_(s-l1) A, E in [-300, 0] step 5 ({time.time()-T0:.0f}s): {"".join("+" if v > 0 else "-" for v in vals)}')
br = [(Es[i + 1], Es[i]) for i in range(60) if (vals[i] > 0) != (vals[i + 1] > 0)]
def bis(f, a, b):
    fa = f(a)
    while b - a > mp.mpf(10) ** -25:
        m = (a + b) / 2; fm = f(m)
        if (fm > 0) == (fa > 0): a, fa = m, fm
        else: b = m
    return (a + b) / 2
zw = mp.mpf(0)
for a, b in br[:4]:
    z1, z2 = bis(Qd, a, b), bis(Wd, a, b); zw = max(zw, abs(z1 - z2)); say(f'zero: projected {mp.nstr(z1, 28)} Wronskian {mp.nstr(z2, 28)} |diff| {mp.nstr(abs(z1 - z2), 3)}')
say(f'zeros (sealed >= 10 digits): {zw < mp.mpf(10) ** -10}' if br else 'no real zeros on [-300, 0] at A')
json.dump(OUT, open(os.path.join(HERE, 't1_3.json'), 'w'), indent=1)
