"""q = 8 test, w2: P2 (the 225 bilinears), P3, controls C2-C6, P4b at weight 14 (the recorded predictions). Exit: 0 every sealed item passes
with the controls firing; 1 a control fails (or I_3 is not one-dimensional); 2 a sealed prediction fails (a result)."""
import os, sys, time
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import q8_lib as L
from q8_lib import K, A, B, say, PX, PI, PA, PB
import sympy as sp
say('w2_bilinear.py sha256', L.self_hash(__file__))
T0 = time.time(); FAIL = []; RES = {}
def control(name, ok):
    ok = bool(ok); say(f'CONTROL {name}: {ok}  ({time.time()-T0:.0f}s)')
    if not ok: FAIL.append(name)
def result(name, ok):
    ok = bool(ok); say(f'SEALED TEST {name}: {ok}  ({time.time()-T0:.0f}s)'); RES[name] = ok
K.set_kappa(2)
TV = {((A, 0), (A, 0)): F(1, 4), ((A, 1),): F(3, 4)}; TZ = {((B, 0), (B, 0)): F(1, 4), ((B, 1),): F(3, 4)}
TT = {((A, 0), (A, 0)): F(1, 4), ((B, 0), (B, 0)): F(1, 4), ((A, 1),): F(-7, 8), ((B, 1),): F(-7, 8)}
c4 = K.reps_mod_d(K.commutant(4, [{A: F(1, 4), B: F(1, 4)}, {A: F(2), B: F(0)}, {A: F(0), B: F(2)}]), 4)
control('I_3 dim 1 (recomputed)', len(c4) == 1); I3 = c4[0]
H7 = K.current(2, A); H7Z = K.current(2, B); G14 = K.mul(H7, H7Z)
P = PA; binom7 = sp.expand(sp.prod([(2 * P + 3 - k) for k in range(7)]) / sp.factorial(7))
E3 = sp.expand(sp.prod([(PX + PI - sp.Rational(j, 2)) * (PX - PI - sp.Rational(j, 2)) for j in range(-3, 4)]))
# ---------------- P2
fA = [{m: F(1)} for m in K.monomials(7, (A,))]; fB = [{m: F(1)} for m in K.monomials(7, (B,))]
prods = [K.mul(f, g) for f in fA for g in fB]; say(f'   bilinears f_i(V) g_j(Z): {len(fA)} x {len(fB)} = {len(prods)}')
t1 = time.time(); vecs = [L.euler_vec(K.ope(I3, p_, poles={1}).get(1, {})) for p_ in prods]
say(f'   images Res_1 I_3 x (f_i g_j) and their Euler images: {time.time()-t1:.0f}s; coordinates {len({k for v in vecs for k in v})}')
t1 = time.time(); rel, rk = L.nullspace_Q(vecs); rkp = L.rank_mod(vecs, 2147483647)
say(f'   relation space over Q: dim {len(rel)} (rank {rk}; rank mod 2^31-1 {rkp}) ({time.time()-t1:.0f}s)')
hv = [H7.get(next(iter(f)), F(0)) for f in fA]; target = [hv[i] * hv[j] for i in range(15) for j in range(15)]
k0 = next(k for k in range(225) if target[k] != 0)
match = len(rel) == 1 and all(rel[0][k] * target[k0] == target[k] * rel[0][k0] for k in range(225))
if rel:
    r0 = rel[0]; say('   relation vector / its first nonzero entry (15 x 15, row-major):', [str(x / r0[k0]) for x in r0])
    say('   h_i h_j / h_0 h_0 (H^(7) coefficients):               ', [str(x / target[k0]) for x in target])
result('P2 the 225-bilinear space in K_14 + d is ONE-dimensional and equals H^(7) (x) H^(7)', match)
# ---------------- P3
t1 = time.time(); br = K.ope(I3, G14, poles={1}).get(1, {}); inK = L.is_exact(br)
say(f'   Res_1 I_3 x G_14: {len(br)} terms; Euler images (V, Z): {len(L.euler(br, A))}, {len(L.euler(br, B))} terms -> in K_14 + d: {inK} ({time.time()-t1:.0f}s)')
t1 = time.time(); eGr = L.eigen_raw(G14); eG = L.internal(eGr, 2); say(f'   eigen(G_14) by the two-boson cylinder rule ({time.time()-t1:.0f}s) = {sp.factor(eG)}')
eGf = L.internal(L.eigen_fast(G14), 2)
control('eigen_fast(G_14) == eigen_raw(G_14)', sp.expand(eGf - eG) == 0)
control('eigen(G_14) = eigen(H_V) eigen(H_Z) (product of the one-boson rules)', sp.expand(eGr - L.eigen_raw(H7) * L.eigen_raw(H7Z)) == 0)
say('   -(64/99225) E_3 =', sp.factor(-sp.Rational(64, 99225) * E3))
c3 = sp.expand(eG + sp.Rational(64, 99225) * E3) == 0
result('P3 G_14 in K_14 (mod d) and eigen(G_14) = -(64/99225) E_3 EXACTLY', inK and c3)
say(f'   (P3 parts: in K_14 {inK}; constant exact {c3}; proportional {L.proportional(E3, eG, (PX, PI))})')
# ---------------- controls
control('C2 eigen(G_14) antisymmetric under P_X <-> pi', sp.expand(L.swap_sym(eG) + eG) == 0)
Gs = K.add(K.mul(H7, H7), K.mul(H7Z, H7Z)); eS = L.internal(L.eigen_fast(Gs), 2)
inKs = L.is_exact(K.ope(I3, Gs, poles={1}).get(1, {}))
control(f'C3 symmetric H_V^2 + H_Z^2: eigen symmetric and NOT antisymmetric, and NOT in K_14 + d (in K_14: {inKs})', sp.expand(L.swap_sym(eS) - eS) == 0 and sp.expand(L.swap_sym(eS) + eS) != 0 and not inKs)
tw = {((A, 0), (A, 0), (A, 1), (A, 2)): F(1, 10)}; twZ = {((B, 0), (B, 0), (B, 1), (B, 2)): F(1, 10)}
Ht = K.add(H7, tw); HtZ = K.add(H7Z, twZ); Gt = K.mul(Ht, HtZ)
leaves = not L.is_exact(K.ope(I3, Gt, poles={1}).get(1, {})); p1b = sp.expand(L.eigen_raw(Ht) - binom7) != 0
p3b = not L.proportional(E3, L.internal(L.eigen_fast(Gt), 2), (PX, PI))
control(f'C4 shape tamper H + (dV)^2 d^2V d^3V/10: leaves K_14 ({leaves}), breaks P1 ({p1b}) and P3 ({p3b})', leaves and p1b and p3b)
Wa = K.schur(F(21, 10), A, 7); WaZ = K.schur(F(21, 10), B, 7); Ga = K.mul(Wa, WaZ)
p1a = not L.proportional(L.eigen_raw(Wa), binom7, (PA,)); la = not L.is_exact(K.ope(I3, Ga, poles={1}).get(1, {}))
control(f'C5 exponent tamper S_7(21/10): breaks P1 ({p1a}) and leaves K_14 ({la})', p1a and la)
p2t = sp.Rational(1, 4) + sp.Rational(1, 1000)
E3t = sp.expand(sp.prod([(PX + PI - j * sp.sqrt(p2t)) * (PX - PI - j * sp.sqrt(p2t)) for j in range(-3, 4)]))
control('C6 scoring tamper E_3 with p^2 + 1/1000 is NOT proportional to eigen(G_14)', not L.proportional(E3t, eG, (PX, PI)))
# ---------------- P4b
t1 = time.time(); hm = K.monomials(13)
vecs = [L.poles_vec(TT, K.dpoly({h: F(1)})) for h in hm]; rhs = L.poles_vec(TT, G14)
say(f'   P4b system: {len(hm)} unknowns H, {len({k for v in vecs for k in v} | set(rhs))} equations (poles >= 3 of T_total); built in {time.time()-t1:.0f}s')
t1 = time.time(); okc, r0, r1 = L.consistent_Q(vecs, rhs)
say(f'   rank homogeneous {r0}, augmented {r1} ({time.time()-t1:.0f}s); a T_total-primary representative exists: {okc}')
result('P4b NO representative G_14 + dH (H weight 13) is T_total-primary', (not okc) and r0 == len(hm))
say(f'done ({time.time()-T0:.0f}s)')
if FAIL: say(f'CONTROLS FAILED {FAIL} -> exit 1'); sys.exit(1)
if not all(RES.values()): say('a sealed prediction failed -> exit 2 (a result)'); sys.exit(2)
say('w2: all sealed items pass with the controls firing'); sys.exit(0)
