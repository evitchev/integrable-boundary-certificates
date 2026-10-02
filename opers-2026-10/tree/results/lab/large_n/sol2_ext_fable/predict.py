"""SOL2F spin-11 prediction (BEFORE vev_profile_sol2_w12 is evaluated at t != 2): family S2(k), M = 1/k and dual, k = 4, 3, 10/3, 1.  Raw and monic
(monic only if the P_X^12 coefficient is nonzero)."""
import sys, time, json, hashlib
from fractions import Fraction as F
import sympy as sp
from common import *
LOG = open(os.path.join(HERE, 'run_predict.log'), 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
out = {'note': 'Sol 2, internal momenta (P_X, pi = P_phi - rho), VIR3-R dictionary; raw = WKB charge with C substituted; monic on P_X^12 when nonzero'}
for k in (F(4), F(3), F(10, 3), F(1)):
    tv = (3 * sp.Rational(k.numerator, k.denominator) + 2) / (sp.Rational(k.numerator, k.denominator) + 2); g = classG(2, tv); n = g['D']
    out[f'k={k},t={tv}'] = {}
    for M in (1 / k, -1 / (k + 1)):
        T0 = time.time(); Cv = 2 * n * (1 + M); Cs = sp.Rational(Cv.numerator, Cv.denominator)
        sv = gwkb.wkb(blocks_of(g, M), g['s'], M, 12); bb = n * M; sk = sv[12]
        if any(a.denominator == 1 and a <= 0 for (a, c) in sk): pw, poly, nonpole = thwkb.integrate_reg(sk, bb); poly = dict(poly)
        else:
            r = thwkb.integrate(sk, bb); poly = {k_: sp.Rational(v.numerator, v.denominator) for k_, v in r[2].items()} if r else {}
        e = sp.expand(sum(v * (Cs * PX**2 / g['aY'])**(i // 2) * (Cs * PI**2 / g['aX'])**(j // 2) for (i, j), v in poly.items()))
        top = sp.Poly(e, PX, PI).coeff_monomial(PX**12) if e != 0 else 0
        rec = {'raw': str(e), 'monic': str(sp.expand(e / top)) if top != 0 else None}
        out[f'k={k},t={tv}'][f'M={M}'] = rec
        say(f'k = {k} (t = {tv}, n = {n}) M = {M}: spin 11 {"IDENTICALLY 0" if e == 0 else ("P_X^12 coeff 0 (projective only)" if top == 0 else "monic ok")}, {len(sp.Poly(e, PX, PI).terms()) if e != 0 else 0} coefficients ({time.time()-T0:.0f}s)')
    a, b = out[f'k={k},t={tv}'].values()
    same = (a['monic'] == b['monic']) if a['monic'] else (a['raw'] == b['raw'] == '0')
    say(f'   M and dual agree: {same}')
json.dump(out, open(os.path.join(HERE, 'PREDICTION_SOL2F.json'), 'w'), indent=1)
h = hashlib.sha256(open(os.path.join(HERE, 'PREDICTION_SOL2F.json'), 'rb').read()).hexdigest()
open(os.path.join(HERE, 'PREDICTION_SOL2F.sha256'), 'w').write(h + '  PREDICTION_SOL2F.json\n'); say('PREDICTION sha256', h)
