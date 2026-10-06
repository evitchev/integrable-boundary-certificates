"""NUM2 refined prediction (written AFTER the seal, BEFORE any data file was read): at t = 2 the WKB coefficient C_K(k) (= -int s_K, absolute) has a
POLE at k = 2 for K = 4, 7, 10 (powers E^(-3/2), E^(-3), E^(-9/2) = E^(-m(k+1)/k), m = 1, 2, 3).  If a non-WKB small-x term at E^(-m(k+1)/k) cancels the
pole, log Q contains R_K (dp_K/dk - dp'_m/dk) E^(p) log E with R_K = Res_(k=2) C_K(k), p_K = (k+1)(1-K)/(k+4), p'_m = -m(k+1)/k:
  K = 4: -R/2;  K = 7: -R;  K = 10: -(3/2) R.   R from the symmetric difference [C_K(2+eps) - C_K(2-eps)] eps/2 (eps = 1e-25)."""
import sys, os, json, time, hashlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
assert hashlib.sha256(open(os.path.join(HERE, 'SEAL_NUM2.md'), 'rb').read()).hexdigest().startswith('25ab6266')
import mpmath as mp, sympy as sp
from fractions import Fraction as F
import gwkb, thwkb
mp.mp.dps = 80
LOG = open(os.path.join(HERE, 'run_residue_pred.log'), 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
PTS = {'A': (mp.mpf(1) / 3, mp.mpf(3) / 5), 'B': (mp.mpf(9) / 10, mp.mpf(1) / 5), 'C': (mp.mpf(1) / 2, mp.mpf(1) / 4)}
def coefs(k, KS):
    k = F(k); n = k + 4; s = (n - 1) / 2; d = n / k; M = 1 / k
    blocks = [(1, {(1, 0): F(1)}, d), (1, {(1, 0): F(-1)}, d), (1, {(0, 1): F(1)}, d), (1, {(0, 1): F(-1)}, d), (k, {}, d)]
    sv = gwkb.wkb(blocks, s, M, max(KS)); out = {}
    for K in KS:
        r = thwkb.integrate(sv[K], n * M)
        out[K] = None if r is None else r
    return out
eps = F(1, 10 ** 25)
cp, cm = coefs(F(2) + eps, (4, 7, 10)), coefs(F(2) - eps, (4, 7, 10))
mult = {4: (1, -mp.mpf(1) / 2), 7: (2, -mp.mpf(1)), 10: (3, -mp.mpf(3) / 2)}
OUT = {}
for nm, (l0, l1) in PTS.items():
    OUT[nm] = {}
    for K in (4, 7, 10):
        vals = []
        for r in (cp[K], cm[K]):
            if r is None: vals.append(mp.mpf(0)); continue
            pw, pref, poly = r
            v = sum(mp.mpf(c.numerator) / c.denominator * l0 ** i * l1 ** j for (i, j), c in poly.items())
            vals.append(-mp.mpf(sp.N(pref, 100)) * v)
        e_ = mp.mpf(eps.numerator) / eps.denominator
        R = (vals[0] - vals[1]) * e_ / 2; fin = (vals[0] + vals[1]) / 2
        Lcoef = mult[K][1] * R
        OUT[nm][K] = {'residue': str(R), 'logE_coefficient_pred': str(Lcoef), 'power': str(-F(3, 2) * mult[K][0])}
        say(f'point {nm} K = {K} (E^{-1.5 * mult[K][0]}): C_K(2+-eps) = {mp.nstr(vals[0], 12)}, {mp.nstr(vals[1], 12)}; residue R = {mp.nstr(R, 25)}; symmetric finite part {mp.nstr(fin, 12)}; PREDICTED coefficient of E^({-1.5 * mult[K][0]}) log E: {mp.nstr(Lcoef, 25)}')
json.dump(OUT, open(os.path.join(HERE, 'PREDICTION_NUM2_logE.json'), 'w'), indent=1)
