"""WUW2 (sealed 'if feasible'): full K_10 = ker ad I_3 at weight 10 mod d, at t = 7/5 (kappa = 3/2); GATE: certified vev_sol3_w10 at 7/5 in span eigen(K_10);
G_10 and the Sol 3 class span K_10?"""
import os, sys, json, hashlib, time
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
assert hashlib.sha256(open(os.path.join(HERE, 'SEAL_WUW2.md'), 'rb').read()).hexdigest().startswith('6ef71c92')
ORIG = {l.split()[1]: l.split()[0] for l in open(os.path.join(HERE, 'inputs', 'ORIGIN_SHA256'))}
assert all(hashlib.sha256(open(os.path.join(HERE, 'inputs', f_), 'rb').read()).hexdigest() == h for f_, h in ORIG.items())
src = open(os.path.join(HERE, 'wuw2.py')).read()
pre = src[:src.index('# ================= (0) engine control')].replace("open(os.path.join(HERE, 'run_wuw2.log'), 'w')", "open(os.path.join(HERE, 'run_k10_full.log'), 'w')")
exec(pre)
K.set_kappa(F(3, 2))
I3 = K.reps_mod_d(K.commutant(4, [{A: F(1, 3), B: F(1, 3)}, {A: F(2), B: F(0)}, {A: F(0), B: F(2)}]), 4)[0]
T1 = time.time(); kd13 = K.d_ech(13); gm = K.monomials(10)
vecs = [K.ope(I3, {g: F(1)}, poles={1}).get(1, {}) for g in gm]
say(f'{len(gm)} weight-10 monomials; images computed ({time.time()-T1:.0f}s)')
T1 = time.time(); rel = K.relations(vecs, kd13); say(f'kernel (before mod d): dim {len(rel)} ({time.time()-T1:.0f}s)')
fields = [{gm[j]: x[j] for j in range(len(gm)) if x[j] != 0} for x in rel]
kd10 = K.d_ech(10); reps = [f for f in fields if kd10.insert(f)]
say(f'dim K_10 mod d = {len(reps)}  (record / VIR5d: 2)')
W5 = K.current(2, A); W5Z = K.current(2, B); G10 = K.mul(W5, W5Z)
tv = sp.Rational(7, 5); p2 = sp.Rational(1, 3); rho2 = (p2 - 2)**2 / (4 * p2); Pa, Qa, tS = sp.symbols('P Q t')
d = json.load(open(os.path.join(HERE, 'inputs', 'vev_sol3_w10.json')))
cert = sp.expand(sum(sp.cancel(sp.sympify(v, locals={'t': tS}).subs(tS, tv)) * Pa**int(k_.split()[0][2:]) * Qa**int(k_.split()[1][2:]) for k_, v in d['coefficients'].items()))
cert = sp.expand(cert.subs({Pa: sp.I * PX}).subs({Qa: sp.sqrt(rho2 - PI**2)}))
def internal(e): return sp.expand(e.subs({PA: sp.sqrt(3) / 2 * (PX + PI), PB: sp.sqrt(3) / 2 * (PI - PX)}))
eK = [internal(eigen_raw(r)) for r in reps]
cs = sp.symbols(f'c0:{len(eK)}'); sol = sp.solve(sp.Poly(sp.expand(cert - sum(c * e for c, e in zip(cs, eK))), PX, PI).coeffs(), cs, dict=True)
say(f'GATE certified vev_sol3_w10 (t = 7/5) in span eigen(K_10): {bool(sol)}')
e_ = K.d_ech(10)
for r in reps: e_.insert(r)
say(f'G_10 in K_10 + d (cross-check of P3 with the full kernel): {not e_.reduce(G10)}')
eG = internal(eigen_raw(G10)); s2 = sp.solve(sp.Poly(sp.expand(cert - cs[0] * eG), PX, PI).coeffs(), [cs[0]], dict=True)
say(f'certified Sol 3 spin 9 itself proportional to eigen(G_10)? {bool(s2)} (expected False: the Sol 3 class is the symmetric one)')
sys.exit(0 if (len(reps) == 2 and sol) else 1)
