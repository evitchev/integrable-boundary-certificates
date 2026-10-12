"""WUW2 -- t = 7/5 (i = 2): V = U/s, Z = W/s (kappa = 3/2).  W^(5) = Res e^(2V) e^(-2V); G_10 = W_V^(5) W_Z^(5).  Sealed SEAL_WUW2.md 6ef71c92....
(0) engine control kappa = 1 reproduces WUW1; (1) W^(5) properties; (2) I_3 gate; P1 eigen(W^(5)); P2 bilinear space in K_10 + d; P3 G_10 in K_10, eigen ~ E_2;
P4 primariness: T_U + T_W yes; T_total (any G_10 + dH)?  Controls and tampers.  The full dim K_10 is a separate script (k10_full.py)."""
import os, sys, json, hashlib, time
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
assert hashlib.sha256(open(os.path.join(HERE, 'SEAL_WUW2.md'), 'rb').read()).hexdigest().startswith('6ef71c92')
import sympy as sp
import kwick as K
from kwick import A, B, F as _F
LOG = open(os.path.join(HERE, 'run_wuw2.log'), 'w'); T0 = time.time()
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
FAIL = []; RES = {}
def control(n_, ok):
    say(f'CONTROL {n_}: {ok}')
    if not ok: FAIL.append(n_)
def result(n_, ok):
    say(f'SEALED TEST {n_}: {ok}'); RES[n_] = bool(ok)
eps = sp.Symbol('eps')
Rser = sp.series(1 / (4 * sp.sinh(eps / 2)**2) - 1 / eps**2, eps, 0, 30).removeO()
Rd = [sp.Rational(sp.diff(Rser, eps, k).subs(eps, 0)) for k in range(28)]
PA, PB, PX, PI = sp.symbols('P_A P_B P_X pi')
def pairings(L):
    if not L: yield [], []; return
    first, rest = L[0], L[1:]
    for pr, un in pairings(rest): yield pr, [first] + un
    for i, o in enumerate(rest):
        if o[0] == first[0]:
            for pr, un in pairings(rest[:i] + rest[i + 1:]): yield [(first, o)] + pr, un
def eigen_raw(G):
    tot = 0; kap = sp.Rational(K.KAPPA.numerator, K.KAPPA.denominator)
    for m, c in G.items():
        for pr, un in pairings(list(m)):
            if any(q >= 1 for (_, q) in un): continue
            val = sp.Rational(c.numerator, c.denominator)
            for (a, b) in pr: val *= kap * (-1)**b[1] * Rd[a[1] + b[1]]
            for (cc, q) in un: val *= (PA if cc == A else PB)
            tot += val
    return sp.expand(tot)
def proportional(a, b, gens):
    if a == 0 or b == 0: return False
    pa, pb = sp.Poly(a, *gens).as_dict(), sp.Poly(b, *gens).as_dict(); ref = next(iter(pa))
    return pb.get(ref, 0) != 0 and all(sp.simplify(pa.get(m, 0) * pb[ref] - pb.get(m, 0) * pa[ref]) == 0 for m in set(pa) | set(pb))
def bilinear_space(I3, W1, kd):
    fA = [{m: F(1)} for m in K.monomials(W1, (A,))]; fB = [{m: F(1)} for m in K.monomials(W1, (B,))]
    prods = [K.mul(f, g) for f in fA for g in fB]
    vecs = [K.ope(I3, p_, poles={1}).get(1, {}) for p_ in prods]
    rel = K.relations(vecs, kd)
    return fA, fB, prods, rel
# ================= (0) engine control: kappa = 1 reproduces WUW1 at t = 5/3
K.set_kappa(1)
c4 = K.reps_mod_d(K.commutant(4, [{A: F(1, 2), B: F(1, 2)}, {A: F(2), B: F(0)}, {A: F(0), B: F(2)}]), 4)
WU1 = K.current(2, A)
_, _, _, rel1 = bilinear_space(c4[0], 3, K.d_ech(4 + 6 - 1))
ok0 = len(c4) == 1 and WU1 == {((A, 0),) * 3: F(4, 3), ((A, 0), (A, 1)): F(2), ((A, 2),): F(1, 3)} and len(rel1) == 1
if ok0:
    r = rel1[0]; ok0 = all(r[3 * i + j] * F(1) == r[0] * [F(1), F(3, 2), F(1, 4)][i] * [F(1), F(3, 2), F(1, 4)][j] for i in range(3) for j in range(3))
control('(0) kappa = 1 engine reproduces WUW1 (I_3 dim 1; W_U; unique bilinear (1,3/2,1/4)^2 in K_6)', ok0)
# ================= kappa = 3/2, t = 7/5
K.set_kappa(F(3, 2))
TU = {((A, 0), (A, 0)): F(1, 3), ((A, 1),): F(2, 3)}; TW = {((B, 0), (B, 0)): F(1, 3), ((B, 1),): F(2, 3)}
TT = {((A, 0), (A, 0)): F(1, 3), ((B, 0), (B, 0)): F(1, 3), ((A, 1),): F(-5, 6), ((B, 1),): F(-5, 6)}
control(f'c(T_U) = -7: 4th pole {K.ope(TU, TU).get(4)}', K.ope(TU, TU).get(4) == {(): F(-7, 2)})
control(f'c(T_total) = -23: 4th pole {K.ope(TT, TT).get(4)}', K.ope(TT, TT).get(4) == {(): F(-23, 2)})
for al, nm in ((F(2), 'e^(2V)'), (F(-2, 3), 'e^(-(2/3)V)')):
    o = K.ope(TU, {((A, 0),): al * K.KAPPA})   # weight check via h = kappa a^2/2 - Q' a: report only
W5 = K.current(2, A); W5Z = K.current(2, B)
say(f'(1) W^(5) = Res e^(2V) e^(-2V) = {W5}')
o = K.ope(TU, W5)
control('(1)(i) W^(5) T_U-primary of weight 5', all(not o.get(r) for r in o if r >= 3) and o.get(2, {}) == K.scale(W5, F(5)) and o.get(1, {}) == K.dpoly(W5))
control('(1)(ii) field level Res e^(-(2/3)V)(z) W^(5)(w) = 0', not K.field_screen(W5, F(-2, 3), A))
say(f'   (reported) field level Res e^(2V) W^(5) has {len(K.field_screen(W5, F(2), A))} terms')
d3T = K.dpoly(K.dpoly(K.dpoly(TU))); TdT = K.ope(TU, K.dpoly(TU)).get(0, {})
e_ = K.Ech(); e_.insert(d3T); e_.insert(TdT)
control('(1)(iii) W^(5) not in span{d^3 T_U, :T_U dT_U:}', bool(e_.reduce(W5)))
# I_3
c4 = K.reps_mod_d(K.commutant(4, [{A: F(1, 3), B: F(1, 3)}, {A: F(2), B: F(0)}, {A: F(0), B: F(2)}]), 4)
I3 = c4[0]
tv = sp.Rational(7, 5); p2 = sp.Rational(1, 3); rho2 = (p2 - 2)**2 / (4 * p2)
Pa, Qa = sp.symbols('P Q')
ROOT = os.path.join(HERE, '..', '..', '..', '..', 'code'); sys.path.insert(0, ROOT); _cwd = os.getcwd(); os.chdir(ROOT)  # KIT PATCH: repo code/ from __file__ (was ~/opus5-w8-clone/code; cyl_first_order.py, wkb.py, first_order/vev_w4_points.json byte-identical there)
from cyl_first_order import cyl, X as Xs, Y as Ys, t as tc
os.chdir(_cwd)
def internal(e): return sp.expand(e.subs({PA: sp.sqrt(3) / 2 * (PX + PI), PB: sp.sqrt(3) / 2 * (PI - PX)}))
cert3 = sp.expand(cyl(3, 4).subs(tc, tv).subs({Xs: Pa**2, Ys: Qa**2}).subs({Pa: sp.I * PX}).subs({Qa: sp.sqrt(rho2 - PI**2)}))
control('GATE I_3 (weight-4 commutant of {e^((V+Z)/3), e^(2V), e^(2Z)}) dim 1 and eigen proportional to cyl(3,4) at t = 7/5',
        len(c4) == 1 and proportional(internal(eigen_raw(I3)), cert3, (PX, PI)))
eW = sp.factor(eigen_raw(W5)); say(f'   eigen(W^(5)) = {eW}')
result('P1 eigen(W^(5)) proportional to P_V (P_V^2 - 1/4)(P_V^2 - 1)', proportional(eigen_raw(W5), sp.expand(PA * (PA**2 - sp.Rational(1, 4)) * (PA**2 - 1)), (PA,)))
# vacuity: shapes of one-boson weight-5 eigenvalues
shapes = [eigen_raw({m: F(1)}) for m in K.monomials(5, (A,))]
say(f'   vacuity: eigenvalues of the 7 one-boson weight-5 monomials span dim {sp.Matrix([[sp.Poly(s_, PA).coeff_monomial(PA**k) for k in (1, 3, 5)] for s_ in shapes]).rank()} (shape space a P^5 + b P^3 + c P)')
# P2
T1 = time.time(); kd13 = K.d_ech(13); say(f'   d-echelon at weight 13: {len(kd13.rows)} rows ({time.time()-T1:.0f}s)')
T1 = time.time(); fA, fB, prods, rel = bilinear_space(I3, 5, kd13)
kd10 = K.d_ech(10); nind = sum(1 for p_ in prods if kd10.insert(p_)) if False else None
say(f'   bilinears f_i(V) g_j(Z) (7 x 7 = 49) with Res_1 I_3 x . in d: relation space dim {len(rel)} ({time.time()-T1:.0f}s)')
wv = [W5.get(next(iter(f)), F(0)) for f in fA]
target = [wv[i] * wv[j] for i in range(7) for j in range(7)]
match = len(rel) == 1 and all(rel[0][k] * target[0] == target[k] * rel[0][0] for k in range(49)) if rel and target[0] != 0 else False
if len(rel) == 1 and target[0] == 0:
    k0 = next(k for k in range(49) if target[k] != 0); match = all(rel[0][k] * target[k0] == target[k] * rel[0][k0] for k in range(49))
result('P2 the bilinear space in K_10 + d is 1-dimensional and equals W^(5) (x) W^(5)', match)
# P3
G10 = K.mul(W5, W5Z)
inK = not kd13.reduce(K.ope(I3, G10, poles={1}).get(1, {}))
E2 = sp.expand(sp.prod([(PX + PI - j * sp.sqrt(p2)) * (PX - PI - j * sp.sqrt(p2)) for j in range(-2, 3)]))
eG = internal(eigen_raw(G10))
say(f'   eigen(G_10) = {sp.factor(eG)}'); say(f'   E_2 = {sp.factor(E2)}')
result('P3 G_10 in K_10 (mod d) and eigen(G_10) proportional to E_2', inK and proportional(E2, eG, (PX, PI)))
say(f'   (P3 parts: in K_10 {inK}; eigen ~ E_2 {proportional(E2, eG, (PX, PI))})')
control('antisymmetry: eigen(G_10) odd under P_X <-> pi', sp.expand(eG.subs({PX: PI, PI: PX}, simultaneous=True) + eG) == 0)
# P4 primariness
def poles_vec(T, G):
    o = K.ope(T, G); v = {}
    for r, f in o.items():
        if r >= 3:
            for m, c in f.items(): v[(r, m)] = v.get((r, m), 0) + c
    return {k: x for k, x in v.items() if x}
prim_UW = not poles_vec(K.add(TU, TW), G10) and K.ope(K.add(TU, TW), G10).get(2) == K.scale(G10, F(10))
result('P4a G_10 primary (weight 10) under T_U + T_W', prim_UW)
T1 = time.time(); eT = K.Ech()
for h in K.monomials(9): eT.insert(poles_vec(TT, K.dpoly({h: F(1)})))
resid = eT.reduce(poles_vec(TT, G10))
result('P4b NO representative G_10 + dH is T_total-primary (consistent with cyl_t2_primary)', bool(resid))
say(f'   ({time.time()-T1:.0f}s; T_total-primary representative exists: {not resid})')
# controls
Gs = K.add(K.mul(W5, W5), K.mul(W5Z, W5Z)); eS = internal(eigen_raw(Gs))
control('symmetric W_V^2 + W_Z^2: eigen NOT antisymmetric', sp.expand(eS.subs({PX: PI, PI: PX}, simultaneous=True) + eS) != 0)
say(f'   (reported) W_V^2 + W_Z^2 in K_10: {not kd13.reduce(K.ope(I3, Gs, poles={1}).get(1, {}))}')
tw = {((A, 0), (A, 0), (A, 2)): F(1, 10)}; twZ = {((B, 0), (B, 0), (B, 2)): F(1, 10)}
Gt = K.mul(K.add(W5, tw), K.add(W5Z, twZ))
control('TAMPER W + (dV)^2 d^3V/10: leaves K_10 and breaks P1 and P3',
        bool(kd13.reduce(K.ope(I3, Gt, poles={1}).get(1, {}))) and not proportional(eigen_raw(K.add(W5, tw)), sp.expand(PA * (PA**2 - sp.Rational(1, 4)) * (PA**2 - 1)), (PA,)) and not proportional(E2, internal(eigen_raw(Gt)), (PX, PI)))
Wa = K.schur(F(21, 10), A, 5); WaZ = K.schur(F(21, 10), B, 5); Ga = K.mul(Wa, WaZ)
control('TAMPER exponent 21/10 (S_5 shape): breaks P1 and leaves K_10',
        not proportional(eigen_raw(Wa), sp.expand(PA * (PA**2 - sp.Rational(1, 4)) * (PA**2 - 1)), (PA,)) and bool(kd13.reduce(K.ope(I3, Ga, poles={1}).get(1, {}))))
E2t = sp.expand(sp.prod([(PX + PI - j * sp.sqrt(p2 + sp.Rational(1, 1000))) * (PX - PI - j * sp.sqrt(p2 + sp.Rational(1, 1000))) for j in range(-2, 3)]))
control('SCORING TAMPER E_2 with p^2 + 1/1000 breaks P3', not proportional(E2t, eG, (PX, PI)))
json.dump({'W5': {str(k): str(v) for k, v in W5.items()}, 'I3': {str(k): str(v) for k, v in I3.items()}, 'eigen_G10': str(eG), 'E2': str(E2),
           'bilinear_relations': [[str(x) for x in r] for r in rel], 'results': RES}, open(os.path.join(HERE, 'wuw2.json'), 'w'), indent=1)
say(f'done ({time.time()-T0:.0f}s)')
if FAIL: say(f'CONTROLS FAILED {FAIL} -> exit 1'); sys.exit(1)
if not all(RES.values()): say('a sealed prediction failed -> exit 2 (a result)'); sys.exit(2)
sys.exit(0)
