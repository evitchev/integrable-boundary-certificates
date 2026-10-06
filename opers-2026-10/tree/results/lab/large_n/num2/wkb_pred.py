"""NUM2 WKB predictions (absolute): coefficient of E^(power_K) in log Q = -int s_K du, s_K the decaying-branch theta-form WKB of the Gamma form
(archived gwkb.py, sol12), integrated by thwkb.integrate (Beta continuation).  Cross-check at k = 2: thwkb.wkb on the sixth-order polynomial
(growing branch c0 = +1, so s_K^dec = (-1)^(1-K) s_K^grow).  Tamper: l0 -> l0 + 1/10.  Writes wkb_pred.json."""
import sys, os, json, time, hashlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
assert hashlib.sha256(open(os.path.join(HERE, 'SEAL_NUM2.md'), 'rb').read()).hexdigest().startswith('25ab6266')
import mpmath as mp, sympy as sp
from fractions import Fraction as F
import gwkb, thwkb
mp.mp.dps = 60
LOG = open(os.path.join(HERE, 'run_wkb_pred.log'), 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
KMAX = 12
def sym_charges(k):
    k = F(k); n = k + 4; s = (n - 1) / 2; d = n / k; M = 1 / k
    blocks = [(1, {(1, 0): F(1)}, d), (1, {(1, 0): F(-1)}, d), (1, {(0, 1): F(1)}, d), (1, {(0, 1): F(-1)}, d), (k, {}, d)]
    T0 = time.time(); sv = gwkb.wkb(blocks, s, M, KMAX); bb = n * M
    out = {}
    for K in range(0, KMAX + 1):
        sk = sv[K]
        if not sk: out[K] = None; continue
        tot = {a + c for (a, c) in sk}
        if any(a.denominator == 1 and a <= 0 for (a, c) in sk) or any((-(a + c)).denominator == 1 and -(a + c) <= 0 for (a, c) in sk):
            out[K] = ('pole', tot); continue
        pw, pref, poly = thwkb.integrate(sk, bb)
        out[K] = (pw, pref, poly)
    say(f'gwkb k = {k}: s_0..s_{KMAX} in {time.time()-T0:.0f}s')
    return out, n
def evalK(entry, l0, l1):
    pw, pref, poly = entry
    v = sum(mp.mpf(c.numerator) / c.denominator * l0 ** i * l1 ** j for (i, j), c in poly.items())
    return mp.mpf(pw.numerator) / pw.denominator, -mp.mpf(sp.N(pref, 80)) * v
PTS = {'A': (mp.mpf(1) / 3, mp.mpf(3) / 5), 'B': (mp.mpf(9) / 10, mp.mpf(1) / 5), 'C': (mp.mpf(1) / 2, mp.mpf(1) / 4)}
RES = {}
for k in (F(2), F(10, 3)):
    ch, n = sym_charges(k)
    RES[str(k)] = {}
    for nm, (l0, l1) in PTS.items():
        row = {}
        for K, e in ch.items():
            if e is None: row[K] = ('zero', None); continue
            if e[0] == 'pole': row[K] = ('pole', str(e[1])); continue
            pw, val = evalK(e, l0, l1)
            row[K] = (str(pw), str(val))
            # tamper l0 + 1/10
        rowt = {K: (None if (e is None or e[0] == 'pole') else str(evalK(e, l0 + mp.mpf(1) / 10, l1)[1])) for K, e in ch.items()}
        RES[str(k)][nm] = {'pred': {str(K): v for K, v in row.items()}, 'tamper_l0': {str(K): v for K, v in rowt.items()}}
        say(f'k = {k} point {nm}: ' + '; '.join(f'K={K}: E^{v[0]} coef {mp.nstr(mp.mpf(v[1]), 18) if v[1] and v[0] not in ("zero", "pole") else v}' for K, v in row.items()))
# cross-check at k = 2 with thwkb on the sixth-order polynomial
l0s, l1s = sp.symbols('l0 l1')
lam = [{(0, 0): F(1)}, {(0, 0): F(4)}, {(0, 0): F(5, 2), (1, 0): F(1)}, {(0, 0): F(5, 2), (1, 0): F(-1)}, {(0, 0): F(5, 2), (0, 1): F(1)}, {(0, 0): F(5, 2), (0, 1): F(-1)}]
sv6 = thwkb.wkb(thwkb.charpoly_coeffs(lam, 6), 6, F(2), F(3), 8)
for nm, (l0, l1) in PTS.items():
    msg = []
    for K in (0, 2, 6, 8):          # K = 4 is the Gamma-pole grade at k = 2 (v1 crashed on it; log kept as run_wkb_pred_v1_crash.stdout)
        r = thwkb.integrate(sv6[K], F(3)); pw, pref, poly = r
        v = sum(mp.mpf(c.numerator) / c.denominator * l0 ** i * l1 ** j for (i, j), c in poly.items()) * mp.mpf(sp.N(pref, 80))
        dec = (-1) ** (1 - K) * v; pred = -dec
        g = mp.mpf(RES['2'][nm]['pred'][str(K)][1])
        msg.append(f'K={K}: thwkb {mp.nstr(pred, 15)} vs gwkb {mp.nstr(g, 15)} rel {mp.nstr(abs(pred - g) / abs(g), 3)}')
    say(f'cross-check k = 2 point {nm}: ' + '; '.join(msg))
json.dump(RES, open(os.path.join(HERE, 'wkb_pred.json'), 'w'), indent=1)
