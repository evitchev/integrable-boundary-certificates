"""NUM1 main (sealed SEAL_NUM1.md c1418ec1...).  Usage: num1.py DPS  (DPS = 60 main; 120 = precision doubling, point A only).
N1 (M = 1 exact); data grid at M = 5 for points A, B, C; E1 residual test; E2 blind fit; tamper M = 5.1; writes num1_<dps>.json."""
import sys, os, json, time, hashlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
assert hashlib.sha256(open(os.path.join(HERE, 'SEAL_NUM1.md'), 'rb').read()).hexdigest().startswith('c1418ec1')
import mpmath as mp
from spec import *
from rwkb import wkb, integral_coeff
DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 60
mp.mp.dps = DPS
LOG = open(os.path.join(HERE, f'run_num1_dps{DPS}.log'), 'w'); T00 = time.time()
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
OUT = {'dps': DPS}
M = 5; mu = mp.mpf(3) / 5
POINTS = {'A': (mp.mpf(1), mp.mpf(1) / 3), 'B': (mp.mpf(-3) / 2, mp.mpf(7) / 10), 'C': (mp.mpf(0), mp.mpf(1) / 3)}
if DPS != 60: POINTS = {'A': POINTS['A']}
# ---- N1
if DPS == 60:
    l = mp.mpf(1) / 3; worst = mp.mpf(0)
    for e in (10, 100, 1000, 10000):
        e = mp.mpf(e); num = logD_neg(1, 0, l, e); ex = mp.log((2 * l + 1) * mp.gamma(l + mp.mpf(1) / 2) / mp.gamma((2 * l + 3 + e) / 4))
        worst = max(worst, abs(num - ex)); say(f'N1 M = 1 E = -{e}: |num - exact| = {mp.nstr(abs(num - ex), 3)}')
    for E in ('2.5', '7.25', '-3.5'):
        num = D_lin(1, 0, l, mp.mpf(E)); ex = (2 * l + 1) * mp.gamma(l + mp.mpf(1) / 2) / mp.gamma((2 * l + 3 - mp.mpf(E)) / 4)
        rel = abs((num - ex) / ex); worst = max(worst, rel); say(f'N1 M = 1 E = {E} (linear): rel diff {mp.nstr(rel, 3)}')
    say(f'N1 worst: {mp.nstr(worst, 3)} -> digits {int(-mp.log10(worst))} (sealed >= 25): {worst < mp.mpf(10) ** -25}')
    OUT['N1_digits'] = int(-mp.log10(worst))
# ---- WKB coefficients
S = wkb(M, 24)
def Cpred(alpha, l, Mv=M):
    lam = l * (l + 1); Sx = S if Mv == M else wkb(Mv, 24)
    return {k: -integral_coeff(Sx[k], Mv, alpha, lam)[0] for k in range(25) if k != 1}     # log D = -sum int S_k
# ---- data grid
NP = 36
grid = [mp.mpf(10) ** (3 + mp.mpf(4) * i / (NP - 1)) for i in range(NP)]
DATA = {}
for nm, (al, l) in POINTS.items():
    T0 = time.time(); DATA[nm] = [logD_neg(M, al, l, e) for e in grid]
    say(f'data {nm}: {NP} points e in [1e3, 1e7] ({time.time()-T0:.0f}s)')
OUT['data'] = {nm: [str(v) for v in vs] for nm, vs in DATA.items()}; OUT['grid'] = [str(e) for e in grid]
def blind_fit(vals, K2, idx=None):
    """unknowns: e^mu, log e, 1, e^(mu(1-k)) for k = 2..K2; least squares at 3x working precision"""
    idx = idx or list(range(NP))
    with mp.workdps(3 * DPS):
        rows = []; rhs = []
        for i in idx:
            e = grid[i]; r = [e ** mu, mp.log(e), mp.mpf(1)] + [e ** (mu * (1 - k)) for k in range(2, K2 + 1)]
            rows.append(r); rhs.append(vals[i])
        A = mp.matrix(rows); b = mp.matrix(rhs)
        x, res = mp.qr_solve(A, b)
        out = {0: x[0], 'log': x[1], 'const': x[2]}
        for k in range(2, K2 + 1): out[k] = x[k + 1]
    return out
def digits(a, b):
    if b == 0: return mp.inf if a == 0 else -mp.log10(abs(a))
    d = abs(a - b) / abs(b); return mp.inf if d == 0 else -mp.log10(d)
RESULT = {}
for nm, (al, l) in POINTS.items():
    pred = Cpred(al, l)
    say(f'\n== point {nm}: alpha = {al}, l = {mp.nstr(l, 6)} (X = {mp.nstr(-al**2 / 120, 6)}, Y = {mp.nstr((4 - (l + mp.mpf(1) / 2) ** 2) / 6, 8)})')
    # E1 residual test
    for G in (4, 6, 8, 24):
        R = [DATA[nm][i] - sum(pred[k] * grid[i] ** (mu * (1 - k)) for k in pred if k <= G) for i in range(NP)]
        # remove grade 1 (c log e + d) by least squares
        with mp.workdps(3 * DPS):
            A = mp.matrix([[mp.log(e), 1] for e in grid]); x, _ = mp.qr_solve(A, mp.matrix(R))
        Rr = [R[i] - x[0] * mp.log(grid[i]) - x[1] for i in range(NP)]
        if G == 24:
            say(f'   E1 G = 24: grade-1 fit: log e coeff {mp.nstr(x[0], 20)}, const {mp.nstr(x[1], 20)}; max |residual| {mp.nstr(max(abs(r) for r in Rr), 3)}')
            OUT.setdefault('grade1', {})[nm] = (str(x[0]), str(x[1]))
        else:
            sl = (mp.log(abs(Rr[NP - 1])) - mp.log(abs(Rr[NP // 2]))) / (mp.log(grid[NP - 1]) - mp.log(grid[NP // 2]))
            say(f'   E1 G = {G}: residual at e = 1e7: {mp.nstr(Rr[-1], 4)}; log-slope {mp.nstr(sl, 6)} (expected -mu (G+1-1) = {mp.nstr(-mu * G, 4)}; if grade G+1 vanishes, -mu (G+1))')
    # E2 blind fit, two truncations for the error estimate
    f1 = blind_fit(DATA[nm], 30); f2 = blind_fit(DATA[nm], 26)
    res = {}
    for k in [0] + list(range(2, 9)):
        dg = digits(f1[k], pred[k]); stab = digits(f1[k], f2[k])
        res[k] = (str(f1[k]), str(pred[k]), float(min(dg, 999)), float(min(stab, 999)))
        say(f'   E2 grade {k} (e^{mp.nstr(mu*(1-k),3)}): extracted {mp.nstr(f1[k], 18)}  WKB {mp.nstr(pred[k], 18)}  agree {mp.nstr(dg, 4)} digits (fit stability {mp.nstr(stab, 4)})')
    RESULT[nm] = res
    if nm == 'A':
        pt = Cpred(al, l, Mv=mp.mpf(51) / 10) if False else None
OUT['E2'] = {nm: {str(k): v for k, v in r.items()} for nm, r in RESULT.items()}
# ---- tamper: WKB at M = 5.1 (point A), compared with the M = 5 extraction
from fractions import Fraction as Fr
St = wkb(Fr(51, 10), 8); al, l = POINTS['A']; lam = l * (l + 1)
for k in (2, 4, 6):
    ct = -integral_coeff(St[k], Fr(51, 10), al, lam)[0]
    say(f'TAMPER M = 5.1, point A grade {k}: WKB(5.1) {mp.nstr(ct, 12)} vs extracted {mp.nstr(mp.mpf(RESULT["A"][k][0]), 12)}: agree {mp.nstr(digits(mp.mpf(RESULT["A"][k][0]), ct), 3)} digits (must be << sealed)')
json.dump(OUT, open(os.path.join(HERE, f'num1_dps{DPS}.json'), 'w'), indent=1)
say(f'done ({time.time()-T00:.0f}s)')
