"""q = 8 test, w1: the frame at kappa = 2 by OPE, GATE G1 (I_3), the other three screenings (reported), H^(7) and its properties, P1,
the vacuity statement, P4a, and the eigen_fast == eigen_raw control. Exit: 0 all pass; 1 a gate/control fails (INCONCLUSIVE); 2 P1 or P4a fails."""
import os, sys, time, random
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import q8_lib as L
from q8_lib import K, A, B, say, PX, PI, PA, PB
import sympy as sp
say('w1_gates.py sha256', L.self_hash(__file__))
T0 = time.time(); FAIL = []; RES = {}
def control(name, ok):
    ok = bool(ok); say(f'CONTROL {name}: {ok}  ({time.time()-T0:.0f}s)')
    if not ok: FAIL.append(name)
def result(name, ok):
    ok = bool(ok); say(f'SEALED TEST {name}: {ok}  ({time.time()-T0:.0f}s)'); RES[name] = ok
K.set_kappa(2)
TV = {((A, 0), (A, 0)): F(1, 4), ((A, 1),): F(3, 4)}; TZ = {((B, 0), (B, 0)): F(1, 4), ((B, 1),): F(3, 4)}
TT = {((A, 0), (A, 0)): F(1, 4), ((B, 0), (B, 0)): F(1, 4), ((A, 1),): F(-7, 8), ((B, 1),): F(-7, 8)}
control(f'c(T_V) = -25/2: 4th pole {K.ope(TV, TV).get(4)} == c/2', K.ope(TV, TV).get(4) == {(): F(-25, 4)})
control(f'c(T_total) = -139/4: 4th pole {K.ope(TT, TT).get(4)} == c/2', K.ope(TT, TT).get(4) == {(): F(-139, 8)})
control('T_V + T_Z - T_total = (13/8)(d^2V + d^2Z)', K.add(K.add(TV, TZ), TT, -1) == {((A, 1),): F(13, 8), ((B, 1),): F(13, 8)})
ws = {al: L.exp_pole(TV, {A: F(al), B: F(0)}, 2) for al in (F(2), F(-1, 2), F(-2))}
say('   weights of e^{alpha V} under T_V (pole-2 coefficient):', {str(k): v for k, v in ws.items()})
control('weights under T_V: e^{2V} 1, e^{-V/2} 1, e^{-2V} 7 (exponentials primary: pole 2 is a constant)', ws[F(2)] == {(): F(1)} and ws[F(-1, 2)] == {(): F(1)} and ws[F(-2)] == {(): F(7)})
wt = {al: L.exp_pole(TT, {A: F(al), B: F(al)}, 2) for al in (F(1, 4), F(-2))}
control('Virasoro screenings e^{(V+Z)/4}, e^{-2(V+Z)} have weight 1 under T_total', wt[F(1, 4)] == {(): F(1)} and wt[F(-2)] == {(): F(1)})
# GATE G1
c4 = K.reps_mod_d(K.commutant(4, [{A: F(1, 4), B: F(1, 4)}, {A: F(2), B: F(0)}, {A: F(0), B: F(2)}]), 4)
say(f'   weight-4 commutant of {{e^((V+Z)/4), e^(2V), e^(2Z)}}: dim {len(c4)}')
I3 = c4[0] if c4 else {}
say('   I_3 =', I3)
cert3 = L.read_cyl34('3', '9/7'); eI = L.internal(L.eigen_raw(I3), 2) if I3 else 0
say('   eigen(I_3) =', sp.factor(eI)); say("   record cyl(3,4) at 9/7 =", sp.factor(cert3))
g1 = len(c4) == 1 and L.proportional(eI, cert3, (PX, PI))
control('GATE G1: I_3 dim 1 and eigen(I_3) proportional to the record Sol-3 cyl(3,4) at t = 9/7', g1)
if g1:
    ratio = sp.cancel(sp.Poly(eI, PX, PI).coeff_monomial(PX**4) / sp.Poly(cert3, PX, PI).coeff_monomial(PX**4)); say('   ratio eigen(I_3)/cyl(3,4) =', ratio)
for nm, al in (('e^{-V/2}', {A: F(-1, 2), B: F(0)}), ('e^{-Z/2}', {A: F(0), B: F(-1, 2)}), ('e^{-2(V+Z)}', {A: F(-2), B: F(-2)})):
    say(f'   (reported) I_3 conserved at field level by {nm}: {L.field_conserved(I3, al)}')
say(f'   (reported) a generic exponent e^{{V/3}} does NOT conserve I_3: {not L.field_conserved(I3, {A: F(1, 3), B: F(0)})}')
# H^(7)
H7 = K.current(2, A); H7Z = K.current(2, B)
say(f'(1) H^(7) = Res e^(2V) e^(-2V) = S_7(2 d^kV/k!), {len(H7)} terms:', H7)
control('(1)(i) H^(7) T_V-primary of weight 7', L.is_primary(TV, H7, 7))
control('(1)(ii) field level Res e^(-V/2)(z) H^(7)(w) = 0', not K.field_screen(H7, F(-1, 2), A))
say(f'   (reported) field level Res e^(2V) H^(7) has {len(K.field_screen(H7, F(2), A))} terms')
vac = L.vir_vacuum_fields(TV, 7)
control('(1)(iii) H^(7) not in the weight-7 T_V vacuum-module span %s' % [p for p, _ in vac], not L.in_span([f_ for _, f_ in vac], H7) and len(vac) == 4)
# P1
eW = L.eigen_raw(H7); P = PA
binom7 = sp.expand(sp.prod([(2 * P + 3 - k) for k in range(7)]) / sp.factorial(7))
say('   eigen(H^(7)) =', sp.factor(eW)); say('   binom(2P+3, 7) =', sp.factor(binom7))
result('P1 eigen(H^(7)) = binom(2 P_V + 3, 7) = (8/315) P_V (P_V^2-1/4)(P_V^2-1)(P_V^2-9/4) EXACTLY', sp.expand(eW - binom7) == 0)
control('eigen_fast == eigen_raw on H^(7), on 10 random weight-7 one-boson fields and 5 random weight-6 two-boson fields', L.eigen_fast(H7) == eW and all(L.eigen_fast(f_) == L.eigen_raw(f_) for f_ in [{m: F(random.randint(-3, 3) or 1) for m in random.sample(K.monomials(7, (A,)), 4)} for _ in range(10)] + [{m: F(random.randint(-3, 3) or 1) for m in random.sample(K.monomials(6), 6)} for _ in range(5)]))
shapes = [L.eigen_raw({m: F(1)}) for m in K.monomials(7, (A,))]
rk = sp.Matrix([[sp.Poly(s_, PA).coeff_monomial(PA**k) for k in (1, 3, 5, 7)] for s_ in shapes]).rank()
say(f'   vacuity: the eigenvalues of the 15 one-boson weight-7 monomials span dim {rk} of the 4-dim odd shape space (predicted: 3)')
# P4a
G14 = K.mul(H7, H7Z)
result('P4a G_14 = H_V^(7) H_Z^(7) primary of weight 14 under T_V + T_Z', L.is_primary(K.add(TV, TZ), G14, 14))
say(f'done ({time.time()-T0:.0f}s)')
if FAIL: say(f'GATE/CONTROL FAILED {FAIL} -> exit 1 (INCONCLUSIVE if a gate)'); sys.exit(1)
if not all(RES.values()): say('a sealed prediction failed -> exit 2 (a result)'); sys.exit(2)
say('w1: all sealed items pass'); sys.exit(0)
