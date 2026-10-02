"""SOL2F hold-out (after hash 48f1e50b was sent): spin 11 vs vev_profile_sol2_w12.json at t = 7/3, 11/5, 9/4, 5/3.  Monic where P_X^12 is nonzero, else
reported.  Negative control t + 1/100; scoring tamper after the hash guard."""
import os, sys, json, hashlib
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
assert hashlib.sha256(open(os.path.join(HERE, 'PREDICTION_SOL2F.json'), 'rb').read()).hexdigest().startswith('48f1e50b')
ORIG = {l.split()[1]: l.split()[0] for l in open(os.path.join(HERE, 'inputs_holdout', 'ORIGIN_SHA256'))}
assert all(hashlib.sha256(open(os.path.join(HERE, 'inputs_holdout', f_), 'rb').read()).hexdigest() == h for f_, h in ORIG.items())
LOG = open(os.path.join(HERE, 'run_holdout.log'), 'w')
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
PX, PI, t = sp.symbols('P_X pi t'); L = {'P_X': PX, 'pi': PI}
pred = json.load(open(os.path.join(HERE, 'PREDICTION_SOL2F.json')))
tab = json.load(open(os.path.join(HERE, 'inputs_holdout', 'vev_profile_sol2_w12.json')))['vev']
def cert_at(tv):
    for v in tab.values():
        if sp.denom(sp.cancel(sp.sympify(v, locals={'t': t}))).subs(t, tv) == 0: return None
    p2 = 2 * (tv - 1) / (tv + 1); rho2 = (p2 - 2)**2 / (4 * p2)
    return sp.expand(sum(sp.cancel(sp.sympify(v, locals={'t': t}).subs(t, tv)) * (-PX**2)**int(k.split(',')[0]) * (rho2 - PI**2)**int(k.split(',')[1]) for k, v in tab.items()))
allok = True
for key, d in pred.items():
    if key == 'note': continue
    tv = sp.sympify(key.split('t=')[1]); c = cert_at(tv)
    if c is None: say(f'{key}: certified table has a pole at t = {tv} -> not scored'); continue
    ctop = sp.Poly(c, PX, PI).coeff_monomial(PX**12); cm = sp.expand(c / ctop)
    for mk, rec in d.items():
        if rec['monic'] is None:
            say(f'{key} {mk}: WKB spin 11 = {rec["raw"][:40]} (sealed: identically 0 at n = 11); certified spin 11 regular with {len(sp.Poly(cm, PX, PI).terms())} coefficients -> MISSING SPIN'); continue
        pr = sp.sympify(rec['monic'], locals=L); nb = len([x for x in sp.Poly(sp.expand(pr - cm), PX, PI).coeffs() if x != 0])
        say(f'{key} {mk}: spin 11 mismatches {nb} / {len(sp.Poly(cm, PX, PI).terms())}'); allok &= nb == 0
        if mk.startswith('M=') and mk == list(d.keys())[0]:
            c2 = cert_at(tv + sp.Rational(1, 100)); c2 = sp.expand(c2 / sp.Poly(c2, PX, PI).coeff_monomial(PX**12))
            say(f'   negative control t + 1/100: {len([x for x in sp.Poly(sp.expand(pr - c2), PX, PI).coeffs() if x != 0])} mismatches (must be > 0)')
            say(f'   scoring tamper: {len([x for x in sp.Poly(sp.expand(pr - cm - PI**4 / 1000), PX, PI).coeffs() if x != 0])} mismatch (must be 1)')
say('HOLD-OUT (scored fibres):', 'PASS' if allok else 'FAIL'); sys.exit(0 if allok else 2)
