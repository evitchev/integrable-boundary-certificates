"""VIR1 lead addition 1: I_7 (weight 8) with QUANTUM products, own engine.
Post-seal extension (spin 7 was explicitly outside the seal).
Exit 0 iff I_7 lies in HV_8 + d at every tested fibre and the controls fire;
2 if some I_7 is outside; 1 on control failure."""
import sys, time, random
from fractions import Fraction as F
sys.dont_write_bytecode = True
import vir_lib as V

t0 = time.time()
LOOP0 = F(7, 3)
hv8 = V.hv_span(8, LOOP0)
f8, B8, i8 = V.quotient_functionals(8, hv8)
C8 = V.Commutant(8)
print('dim V_8 = %d, dim HV_8 span = %d, dim Q_8 = %d; complement of dV_7: %d; V_11/dV_10: %d  (%.1fs)'
      % (len(B8), len(hv8), len(f8), len(C8.comp), len(C8.proj), time.time() - t0), flush=True)
bad = ctlbad = 0
TS = [F(2), F(-2), F(7, 2), F(-7, 2), F(4), F(11, 5), F(-3, 7), F(9, 2), F(-13, 4), F(1, 2)]
SP = [F(-1, 3), F(1, 3), F(3), F(5), F(-5), F(5, 3)]
cl = {k: V.classical_span(8, k) for k in (6, 8)}
def member(vec, span):
    rows = [V.to_row(v, i8) for v in span] + [V.to_row(V.d_mono(m), i8) for m in V.basis(7)]
    r0 = V.rank(rows)
    return V.rank(rows + [vec]) == r0
for sol in (1, 2, 3):
    for t in TS + SP:
        N, s = V.curve_point(t)
        K = C8.kernel(V.P4(sol, N, s), N)
        coords = [V.apply(f8, v) for v in K]
        rk = V.rank(coords) if any(any(c) for c in coords) else 0
        tag = 'generic' if t in TS else 'SPECIAL'
        print('Sol %d t=%-5s (N=%s) %s: kernel dim %d = %d in HV+d (+) %d outside  (%.1fs)'
              % (sol, t, N, tag, len(K), len(K) - rk, rk, time.time() - t0), flush=True)
        if t in TS and (len(K) != 1 or rk):
            bad += 1
        if t == 2:
            c6, c8 = member(K[0], cl[6]), member(K[0], cl[8])
            print('      classical at t=2: I_7 in A<=6 %s, A<=8 %s (bcd table: no / YES)' % (c6, c8))
            if (c6, c8) != (False, True):
                ctlbad += 1
            pl = list(K[0]); pl[i8[((), ((4, 4),))]] += F(1, 1000)
            if not any(V.apply(f8, pl)):
                ctlbad += 1
            else:
                print('      plant I_7 + (1/1000)(44)_Y flagged')
rng = random.Random(11)
g = [F(rng.randint(-9, 9), rng.randint(1, 9)) for _ in B8]
nz = sum(1 for c in V.apply(f8, g) if c)
print('vacuity: generic V_8 vector has %d of %d coordinates nonzero' % (nz, len(f8)))
if nz != len(f8):
    ctlbad += 1
print('RESULT:', 'CONTROL FAILURE' if ctlbad else ('I_7 OUTSIDE HV at %d generic fibres' % bad if bad else
      'I_7 in <X, T_Y> + Im d at every generic fibre tested, all three solutions'))
sys.exit(1 if ctlbad else (2 if bad else 0))
