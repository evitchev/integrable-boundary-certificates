"""NUM2 (3) Theorem 1 at spectral level, t = 2 (sealed T1): Q3_j(E) / Q6_j(E_Gamma = -E) == Q3_j(0)/Q6_j(0) (exact Meijer-G values) for the four
dynamical exponents, at E in {-20, -5, 2, 50, 500} (E = 0 is X1).  Also the frozen-exponent Q6 (theta = 1, 4).
(Zeros, sealed) Q_(s - l_max)(E) at point A on E in [-300, 0]: scan; zeros by bisection on the projected Q vs on the outward-Wronskian det
[chi_i (i != j) integrated OUTWARD from x0 to x1 = 3/2 ; y integrated inward to x1]."""
import sys, os, json, time, hashlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
assert hashlib.sha256(open(os.path.join(HERE, 'SEAL_NUM2.md'), 'rb').read()).hexdigest().startswith('25ab6266')
import mpmath as mp
from fractions import Fraction as Fr
from ops import *
mp.mp.dps = 60
LOG = open(os.path.join(HERE, 'run_t1_zeros.log'), 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
PTS = {'A': (mp.mpf(1) / 3, mp.mpf(3) / 5), 'B': (mp.mpf(9) / 10, mp.mpf(1) / 5)}
OUT = {}; worst = mp.mpf(0)
for nm, (l0, l1) in PTS.items():
    e3, _ = exact_three_term_E0(2, l0, l1); e6, _ = exact_sixth_E0(l0, l1)
    C = [e3[j] / e6[j + 2] for j in range(4)]
    OUT[nm] = {}
    for E in (-20, -5, 2, 50, 500):
        E = mp.mpf(E); x0 = 2 / (abs(E) + 1) ** (mp.mpf(1) / 6); T0 = time.time()
        op3, th3, F = three_term(2, l0, l1, E); lq3, Q3 = op3.Q(th3, x0)
        op6, ex6 = sixth(l0, l1, -E); lq6, Q6 = op6.Q(ex6, x0)
        rat = [Q3[j] / Q6[j + 2] for j in range(4)]; dev = [abs(rat[j] / C[j] - 1) for j in range(4)]
        worst = max(worst, max(dev))
        OUT[nm][str(E)] = {'Q3': [str(q) for q in Q3], 'Q6': [str(q) for q in Q6], 'ratio_over_C_minus_1': [str(d) for d in dev]}
        say(f'T1 point {nm} E = {E}: |Q3_j/Q6_j / C_j - 1| = {[mp.nstr(d, 3) for d in dev]}; frozen Q6 (theta = 1, 4) = {mp.nstr(Q6[0], 12)}, {mp.nstr(Q6[1], 12)} ({time.time()-T0:.0f}s)')
say(f'T1 (sealed >= 25 digits, ratio E-independent and equal to the exact E = 0 ratio): worst deviation {mp.nstr(worst, 3)} -> {worst < mp.mpf(10) ** -25}')
# zeros at point A
l0, l1 = PTS['A']; jdom = 2       # theta = s - l1 (l1 = 3/5 is l_max)
def Qdom(E):
    E = mp.mpf(E); op, th, F = three_term(2, l0, l1, E); lq, Q = op.Q(th, mp.mpf('0.6')); return Q[jdom]
def Wdet(E, x1=mp.mpf(3) / 2):
    E = mp.mpf(E); op, th, F = three_term(2, l0, l1, E); x0 = mp.mpf('0.6')
    xm = op.xmax(mp.mpf(2)); L, Wv = op.asym_eval(xm)
    Y1 = [Wv[m] / mp.factorial(m) for m in range(4)]; Yy, lsy = op.integrate(mp.log(xm), Y1, mp.log(x1))
    cols = [[Yy[m] * mp.factorial(m) for m in range(4)]]
    for i, t in enumerate(th):
        if i == jdom: continue
        v = op.frob(t, x0); Yc, lsc = op.integrate(mp.log(x0), [v[m] / mp.factorial(m) for m in range(4)], mp.log(x1))
        cols.append([Yc[m] * mp.factorial(m) for m in range(4)])
    return mp.det(mp.matrix(cols).T)
Es = [-mp.mpf(i) * 5 for i in range(0, 61)]
T0 = time.time(); vals = [Qdom(E) for E in Es]; say(f'zeros scan Q_(s-l1) on E in [-300, 0] step 5 ({time.time()-T0:.0f}s): signs {"".join("+" if v > 0 else "-" for v in vals)}')
br = [(Es[i + 1], Es[i]) for i in range(len(Es) - 1) if (vals[i] > 0) != (vals[i + 1] > 0)]
say(f'brackets: {[(float(a), float(b)) for a, b in br]}')
def bisect(f, a, b, tol):
    fa = f(a)
    while b - a > tol:
        m = (a + b) / 2; fm = f(m)
        if (fm > 0) == (fa > 0): a, fa = m, fm
        else: b = m
    return (a + b) / 2
zw = mp.mpf(0)
for a, b in br[:4]:
    z1 = bisect(Qdom, a, b, mp.mpf(10) ** -25); z2 = bisect(Wdet, a, b, mp.mpf(10) ** -25)
    zw = max(zw, abs(z1 - z2)); say(f'zero: projected Q {mp.nstr(z1, 28)}  outward-Wronskian det {mp.nstr(z2, 28)}  |diff| {mp.nstr(abs(z1 - z2), 3)}')
if br: say(f'zeros (sealed >= 10 digits): {zw < mp.mpf(10) ** -10}')
else: say('no real zeros of Q_(s-l1) on [-300, 0] at point A (sealed: then nothing to compare)')
json.dump(OUT, open(os.path.join(HERE, 't1.json'), 'w'), indent=1)
