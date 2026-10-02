"""NUM1 calibration of the WKB engine at M = 1, alpha = 0 against the EXACT asymptotic expansion of log D = log[(2l+1) Gamma(l+1/2)/Gamma((2l+3+e)/4)]:
sign/branch convention fixed here once (log D ~ sigma * sum_k int S_k du for k >= 2)."""
import sys, os, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
assert hashlib.sha256(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'SEAL_NUM1.md'), 'rb').read()).hexdigest().startswith('c1418ec1')
import mpmath as mp, sympy as sp
from rwkb import *
LOG = open('run_n1_wkb.log', 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
mp.mp.dps = 40
ls = sp.Symbol('l'); es = sp.Symbol('e', positive=True); w = sp.Symbol('w', positive=True)
z = (es + 2 * ls + 3) / 4
# -log Gamma(z) asymptotic, then expand in 1/e
st = -((z - sp.Rational(1, 2)) * sp.log(z) - z + sp.log(2 * sp.pi) / 2 + sum(sp.bernoulli(2 * n) / (2 * n * (2 * n - 1) * z**(2 * n - 1)) for n in range(1, 8)))
ser = sp.series(st.subs(es, 1 / w), w, 0, 7).removeO()
S = wkb(1, 8)
for lv in (sp.Rational(1, 3), sp.Rational(7, 10)):
    lam = mp.mpf(lv.p) / lv.q * (mp.mpf(lv.p) / lv.q + 1)
    for k in range(2, 8):
        C, poles = integral_coeff(S[k], 1, 0, lam)
        exact = sp.Rational(1) * ser.coeff(w, k - 1).subs(ls, lv)
        say(f'l = {lv} grade {k} (e^{1-k}): WKB int = {mp.nstr(C, 25)}  exact coeff = {mp.nstr(mp.mpf(sp.N(exact, 50)), 25)}  ratio = {mp.nstr(C / mp.mpf(sp.N(exact, 50)), 20) if exact != 0 else "exact 0"}  poles {poles}')
