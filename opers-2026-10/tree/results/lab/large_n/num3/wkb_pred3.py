"""NUM3 WKB predictions (absolute, -int s_K) from the archived Sol 2 Gamma engine gwkb (blocks (k, +-l0, delta), (1, +-l1, delta), (1, 0, delta),
s = k+1, M = 1/k, delta = n/k), K = 0..KMAX; Gamma-pole grades detected.  Pole grades: residue in k by the symmetric difference at k = 2 +- 1e-25, and the
predicted log E coefficient R (dp_K/dk - dp'_m/dk).  Tamper l0 + 1/10 values stored.  Writes wkb_pred3.json and PREDICTION_NUM3_logE.json."""
import sys, os, json, time, hashlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
assert hashlib.sha256(open(os.path.join(HERE, 'SEAL_NUM3.md'), 'rb').read()).hexdigest().startswith('b65ef9e8')
import mpmath as mp, sympy as sp
from fractions import Fraction as F
import gwkb, thwkb
mp.mp.dps = 80
LOG = open(os.path.join(HERE, 'run_wkb_pred3.log'), 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
KMAX = 16
PTS = {'A': (mp.mpf(1) / 3, mp.mpf(3) / 5), 'B': (mp.mpf(9) / 10, mp.mpf(1) / 5)}
def charges(k, KM):
    k = F(k); n = 2 * k + 3; s = k + 1; d = n / k; M = 1 / k
    blocks = [(k, {(1, 0): F(1)}, d), (k, {(1, 0): F(-1)}, d), (F(1), {(0, 1): F(1)}, d), (F(1), {(0, 1): F(-1)}, d), (F(1), {}, d)]
    sv = gwkb.wkb(blocks, s, M, KM); bb = n * M; out = {}
    for K in range(KM + 1):
        sk = sv[K]
        if not sk: out[K] = None; continue
        if any(a.denominator == 1 and a <= 0 for (a, c) in sk) or any((-(a + c)).denominator == 1 and -(a + c) <= 0 for (a, c) in sk):
            out[K] = ('pole', {a + c for (a, c) in sk}.pop()); continue
        out[K] = thwkb.integrate(sk, bb)
    return out, n
def ev(entry, l0, l1):
    pw, pref, poly = entry
    v = sum(mp.mpf(c.numerator) / c.denominator * l0 ** i * l1 ** j for (i, j), c in poly.items())
    return mp.mpf(pw.numerator) / pw.denominator, -mp.mpf(sp.N(pref, 100)) * v
RES = {}; poles = {}
for k in (F(2), F(10, 3)):
    T0 = time.time(); ch, n = charges(k, KMAX); say(f'gwkb Sol 2 k = {k}: s_0..s_{KMAX} ({time.time()-T0:.0f}s)')
    RES[str(k)] = {}
    for nm, (l0, l1) in PTS.items():
        row = {}; rowt = {}
        for K, e in ch.items():
            if e is None: row[K] = ('zero', None); rowt[K] = None; continue
            if e[0] == 'pole': row[K] = ('pole', str(e[1])); rowt[K] = None; poles.setdefault(str(k), set()).add(K); continue
            pw, v = ev(e, l0, l1); row[K] = (str(pw), str(v)); rowt[K] = str(ev(e, l0 + mp.mpf(1) / 10, l1)[1])
        RES[str(k)][nm] = {'pred': {str(K): v for K, v in row.items()}, 'tamper_l0': {str(K): v for K, v in rowt.items()}}
        say(f'k = {k} point {nm}: ' + '; '.join(f'K={K}: ' + (v[0] if v[0] in ('zero',) else ('POLE at E^' + v[1] if v[0] == 'pole' else f'E^{mp.nstr(mp.mpf(F(v[0]).numerator)/F(v[0]).denominator, 6)} {mp.nstr(mp.mpf(v[1]), 15)}')) for K, v in row.items()))
json.dump(RES, open(os.path.join(HERE, 'wkb_pred3.json'), 'w'), indent=1)
say(f'pole grades: { {k: sorted(v) for k, v in poles.items()} }')
# residues at k = 2 for the pole grades
PRED = {}
if '2' in poles:
    eps = F(1, 10 ** 25); KS = [K for K in sorted(poles['2']) if K > 1]   # v2: K = 1 (E^0 log E, a pole at every k) excluded; v1 crashed on it
    cp, _ = charges(F(2) + eps, max(KS)); cm, _ = charges(F(2) - eps, max(KS))
    for nm, (l0, l1) in PTS.items():
        PRED[nm] = {}
        for K in KS:
            vp = ev(cp[K], l0, l1)[1] if cp[K] and cp[K][0] != 'pole' else None
            vm = ev(cm[K], l0, l1)[1] if cm[K] and cm[K][0] != 'pole' else None
            e_ = mp.mpf(eps.numerator) / eps.denominator
            R = (vp - vm) * e_ / 2; fin = (vp + vm) / 2
            kk = mp.mpf(2); nn = 2 * kk + 3
            dpK = (1 - K) * ((nn) - 2 * (kk + 1)) / nn ** 2            # d/dk [(k+1)(1-K)/(2k+3)]
            pK = (kk + 1) * (1 - K) / nn; m = -pK * kk / (kk + 1)         # p'_m = -m(k+1)/k matched to p_K
            dpm = m / kk ** 2
            L = R * (dpK - dpm)
            PRED[nm][str(K)] = {'power': str(pK), 'm': str(m), 'residue': str(R), 'finite_sym': str(fin), 'logE_coefficient_pred': str(L)}
            say(f'point {nm} K = {K} (E^{mp.nstr(pK, 6)}, m = {mp.nstr(m, 6)}): residue R = {mp.nstr(R, 25)}; PREDICTED coefficient of E^{mp.nstr(pK, 6)} log E = {mp.nstr(L, 25)}')
json.dump(PRED, open(os.path.join(HERE, 'PREDICTION_NUM3_logE.json'), 'w'), indent=1)
