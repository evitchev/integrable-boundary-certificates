"""VIR1 main (Fable seat).  Seal: SEAL_VIR1.md (hash guard below).

Tests the sealed predictions
  P-spin3: the non-Virasoro (W_4) coordinate of I_3 vanishes identically in t;
  P-spin5: the three non-Virasoro coordinates of I_5 vanish identically in t;
  P-jump : the two-dimensional spin-5 kernels at t = 5/3 and t = 3 lie in HV + d
on Solutions 1, 2, 3, with controls C1..C5 of the seal, plus the lead's two
additions (quantum ladder at t = 2 through I_7; full kernel decompositions).

Exit 0 = P-spin3 and P-spin5 confirmed and every control fires;
     2 = a sealed prediction refuted;  1 = a control failed.
"""
import hashlib
import random
import sys
from fractions import Fraction as F

sys.dont_write_bytecode = True
import vir_lib as V

SEAL = '6f23e25680592942a8370a9417609c7c558323a4af0c32c4dfe246a33993cf23'
got = hashlib.sha256(open('SEAL_VIR1.md', 'rb').read()).hexdigest()
if got != SEAL:
    print('SEAL HASH MISMATCH', got)
    sys.exit(1)
print('seal hash ok', SEAL)

TAMPER = len(sys.argv) > 1 and sys.argv[1] == 'tamper'   # scoring tamper (after the guard)

ctl_fail = []
pred_fail = []


def ctl(cond, msg):
    print(('PASS  ' if cond else 'CONTROL-FAIL  ') + msg, flush=True)
    if not cond:
        ctl_fail.append(msg)


def pred(cond, msg):
    print(('HOLDS  ' if cond else 'REFUTED  ') + msg, flush=True)
    if not cond:
        pred_fail.append(msg)


LOOP0 = F(7, 3)     # any value: the Virasoro states are loop-independent (v0, check d)
f4, B4, i4 = V.quotient_functionals(4, V.hv_span(4, LOOP0))
f6, B6, i6 = V.quotient_functionals(6, V.hv_span(6, LOOP0))
if TAMPER:
    # scoring tamper: spoil the functionals (add an HV-direction component)
    f6 = [[a + (F(1, 1000) if j == i6[((1, 1, 1, 1, 1, 1), ())] else 0)
           for j, a in enumerate(lam)] for lam in f6]
    f4 = [[a + (F(1, 1000) if j == i4[((1, 1, 1, 1), ())] else 0)
           for j, a in enumerate(lam)] for lam in f4]
print('dim Q_4 = %d, dim Q_6 = %d' % (len(f4), len(f6)))
C6 = V.Commutant(6)
print('weight 6: complement of dV_5 has %d monomials; target quotient V_9/dV_8 has dim %d'
      % (len(C6.comp), len(C6.proj)))

SPECIAL = {F(1), F(-1), F(1, 3), F(-1, 3), F(3), F(-3), F(5), F(-5), F(5, 3), F(-5, 3), F(0)}
SAMPLES = []
for d, rng in ((2, range(-31, 32, 2)), (7, range(-45, 46)), (11, range(-40, 41))):
    for n in rng:
        t = F(n, d)
        if t not in SPECIAL and t not in SAMPLES and abs(t) != 1:
            SAMPLES.append(t)
random.Random(20261001).shuffle(SAMPLES)
SAMPLES = SAMPLES[:120]
print('%d generic sample fibres (fixed seed), both signs: %d negative'
      % (len(SAMPLES), sum(1 for t in SAMPLES if t < 0)))

# ------------------------------------------------------------ P-spin3
print('\n== P-spin3: W_4 coordinate of I_3 ==')
for sol in (1, 2, 3):
    vals = set()
    for t in SAMPLES:
        N, s = V.curve_point(t)
        vals.add(V.apply(f4, V.to_row(V.P4(sol, N, s), i4))[0])
    pred(vals == {F(0)}, 'Sol %d: W_4 coordinate of I_3 = 0 at all %d sample fibres (values seen: %s)'
         % (sol, len(SAMPLES), sorted(vals)))
# identically: the functional touches only the constant Y-only coefficients
supp = {B4[i] for i, c in enumerate(f4[0]) if c}
print('   functional support:', sorted(supp))
ctl(TAMPER or all(m[0] == () for m in supp),
    'W_4 functional is supported on Y-only monomials, whose I_3 coefficients are t-independent '
    '-> zero identically in t')

# ------------------------------------------------------------ P-spin5
print('\n== P-spin5: three non-Virasoro coordinates of I_5 ==')
xtop = i6[((1, 1, 1, 1, 1, 1), ())]
kern_store = {}
for sol in (1, 2, 3):
    dims = {}
    nonzero = 0
    xzero = 0
    for t in SAMPLES:
        N, s = V.curve_point(t)
        K = C6.kernel(V.P4(sol, N, s), N)
        dims[len(K)] = dims.get(len(K), 0) + 1
        for v in K:
            if any(V.apply(f6, v)):
                nonzero += 1
            if not v[xtop]:
                xzero += 1
        kern_store[(sol, t)] = K
    print('   Sol %d: kernel dimensions over samples %s; (dX)^6 coefficient zero in %d kernels'
          % (sol, dims, xzero))
    ctl(dims == {1: len(SAMPLES)}, 'Sol %d: kernel one-dimensional at all sample fibres' % sol)
    pred(nonzero == 0, 'Sol %d: all three non-Virasoro coordinates of I_5 vanish at all %d sample fibres'
         % (sol, len(SAMPLES)))
print('   degree bound: (t^2-1)^4 x the weight-6 system has polynomial entries of degree <= 8;')
print('   an adjugate kernel vector has degree <= 11 x 8 = 88 < %d sample fibres -> identically zero' % len(SAMPLES))

# ------------------------------------------------------------ special points
print('\n== special points (full kernels) ==')
POINTS = [('exact point N=28, Sol 1 t=-1/3', 1, F(-1, 3)),
          ('exact point N=28, Sol 2 t=+1/3', 2, F(1, 3)),
          ('anchor N=-2, Sol 3 t=3', 3, F(3)),
          ('N=-2, Sol 2 t=3', 2, F(3)),
          ('N=-2, Sol 1 t=3', 1, F(3)),
          ('N=0, Sol 3 t=5', 3, F(5)),
          ('N=0, Sol 2 t=-5', 2, F(-5)),
          ('jump N=-25/2, Sol 2 t=5/3', 2, F(5, 3)),
          ('jump N=-25/2, Sol 3 t=5/3', 3, F(5, 3)),
          ('N=-25/2 other sheet, Sol 1 t=5/3', 1, F(5, 3))]
special = {}
for name, sol, t in POINTS:
    N, s = V.curve_point(t)
    K = C6.kernel(V.P4(sol, N, s), N)
    coords = [V.apply(f6, v) for v in K]
    rk = V.rank(coords) if coords and any(any(c) for c in coords) else 0
    special[(sol, t)] = (len(K), rk)
    print('   %-36s (N,s)=(%s,%s): kernel dim %d = %d in HV+d (+) %d outside; I_3 W_4 coord %s'
          % (name, N, s, len(K), len(K) - rk, rk,
             V.apply(f4, V.to_row(V.P4(sol, N, s), i4))[0]))
# N = 1 (t -> infinity): the density is regular at (N, s) = (1, 0)
for sol in (1, 2, 3):
    K = C6.kernel(V.P4(sol, F(1), F(0)), F(1))
    coords = [V.apply(f6, v) for v in K]
    rk = V.rank(coords) if coords and any(any(c) for c in coords) else 0
    special[(sol, 'inf')] = (len(K), rk)
    print('   %-36s (N,s)=(1,0): kernel dim %d = %d in HV+d (+) %d outside'
          % ('N=1 (t=inf), Sol %d' % sol, len(K), len(K) - rk, rk))


def limit_class(sol, t0):
    """Limit of the generic one-dimensional class as t -> t0: evaluate at
    t0 + eps for shrinking exact eps and check which kernel line at t0 the
    normalised class approaches (exactly: project onto the kernel at t0)."""
    out = []
    for eps in (F(1, 10 ** 6), F(1, 10 ** 9)):
        N, s = V.curve_point(t0 + eps)
        K = C6.kernel(V.P4(sol, N, s), N)
        out.append((len(K), [V.apply(f6, v) for v in K]))
    return out


for name, sol, t in POINTS:
    lc = limit_class(sol, t)
    ok = all(d == 1 and not any(c[0]) for d, c in lc)
    print('   approach to %-30s kernel dims %s, non-Virasoro coords all zero: %s'
          % (name, [d for d, _ in lc], ok))

yes = lambda key: special[key][1] == 0 and special[key][0] >= 1
pred(yes((1, F(-1, 3))) and yes((2, F(1, 3))), 'N = 28: yes (whole spin-5 kernel in HV + d, both solutions)')
pred(yes((3, F(3))), 'N = -2: yes (Sol 3 anchor)')
pred(yes((3, F(5))) and yes((2, F(-5))), 'N = 0: yes')
pred(all(yes((sol, 'inf')) for sol in (1, 2, 3)), 'N = 1: yes')
pred(yes((2, F(5, 3))) and yes((3, F(5, 3))), 'N = -25/2: yes')
jump_ok = special[(2, F(5, 3))] == (2, 0) and special[(3, F(5, 3))] == (2, 0)
anch_ok = special[(3, F(3))] == (2, 0)
print(('HOLDS  ' if jump_ok else 'REFUTED (secondary)  ')
      + 'P-jump t=5/3: kernel dim 2, both directions in HV + d: %s' % (special[(2, F(5, 3))],))
print(('HOLDS  ' if anch_ok else 'REFUTED (secondary)  ')
      + 'P-jump t=3 (Sol 3): kernel dim 2, both directions in HV + d: %s' % (special[(3, F(3))],))

# ------------------------------------------------------------ C1: classical ladder
print('\n== C1: cyl_membership\'s CLASSICAL test, reproduced with my own span ==')


def member(vec, w, span):
    B = V.basis(w)
    idx = {m: i for i, m in enumerate(B)}
    rows = [V.to_row(v, idx) for v in span] + [V.to_row(V.d_mono(m), idx) for m in V.basis(w - 1)]
    r0 = V.rank(rows)
    return V.rank(rows + [vec]) == r0


cl4 = {k: V.classical_span(4, k) for k in (2, 4)}
cl6 = {k: V.classical_span(6, k) for k in (2, 4, 6)}
for sol in (1, 2, 3):
    N, s = V.curve_point(2)
    p4 = V.to_row(V.P4(sol, N, s), i4)
    i5 = C6.kernel(V.P4(sol, N, s), N)[0]
    row3 = [member(p4, 4, cl4[k]) for k in (2, 4)]
    row5 = [member(i5, 6, cl6[k]) for k in (2, 4, 6)]
    print('   Sol %d t=2 classical: I_3 in A<=2 %s, A<=4 %s;  I_5 in A<=2 %s, A<=4 %s, A<=6 %s'
          % (sol, *row3, *row5))
    ctl(row3 == [False, True] and row5 == [False, False, True],
        'Sol %d: bcd-note table reproduced at t = 2 (I_3: no/YES; I_5: no/no/YES)' % sol)
for sol, t in ((1, F(-1, 3)), (2, F(1, 3))):
    N, s = V.curve_point(t)
    p4 = V.to_row(V.P4(sol, N, s), i4)
    ctl(not member(p4, 4, cl4[2]),
        'the classical test ALSO rejects the exact point (Sol %d, t = %s): I_3 "not in A<=2"' % (sol, t))

# ------------------------------------------------------------ quantum ladder at t = 2
print('\n== lead addition 1: the ladder at t = 2 with QUANTUM products (HV = <X, T_Y>) ==')
hv4, hv6 = V.hv_span(4, LOOP0), V.hv_span(6, LOOP0)
for sol in (1, 2, 3):
    N, s = V.curve_point(2)
    p4 = V.to_row(V.P4(sol, N, s), i4)
    i5 = C6.kernel(V.P4(sol, N, s), N)[0]
    q3, q5 = member(p4, 4, hv4), member(i5, 6, hv6)
    print('   Sol %d t=2 quantum: I_3 in <X,T_Y> %s;  I_5 in <X,T_Y> %s' % (sol, q3, q5))
    pred(q3 and q5, 'Sol %d: I_3 and I_5 at t = 2 lie in the vertex subalgebra <X, T_Y> + Im d' % sol)

# ------------------------------------------------------------ C2: exact point (done above) / C3 tamper
print('\n== C3: loop-factor tamper N -> N + 1 at the exact point ==')
for sol, t in ((1, F(-1, 3)), (2, F(1, 3))):
    N, s = V.curve_point(t)
    K0 = C6.kernel(V.P4(sol, N, s), N)
    K1 = C6.kernel(V.P4(sol, N, s), N + 1)
    print('   Sol %d t=%s: kernel dim %d with loop N; %d with loop N+1' % (sol, t, len(K0), len(K1)))
    ctl(len(K0) == 1 and len(K1) == 0,
        'Sol %d: with the loop factor N+1 the spin-5 charge at the exact point disappears' % sol)

# ------------------------------------------------------------ C4 plants
print('\n== C4: plants ==')
N, s = V.curve_point(2)
for sol in (1, 2, 3):
    i5 = C6.kernel(V.P4(sol, N, s), N)[0]
    planted = list(i5)
    planted[i6[((), ((3, 3),))]] += F(1, 1000)
    c = V.apply(f6, planted)
    ctl(any(c), 'Sol %d: I_5 + (1/1000)(33)_Y is flagged: coords %s' % (sol, [str(x) for x in c]))
    p4 = dict(V.P4(sol, N, s))
    p4[((), ((2, 2),))] += F(1, 1000)
    c = V.apply(f4, V.to_row(p4, i4))
    ctl(any(c), 'Sol %d: I_3 with (22) -> -2 + 1/1000 is flagged: coord %s' % (sol, c[0]))
    planted = list(i5)
    planted[i6[((1, 1), ((2, 2),))]] += F(1, 1000)
    c = V.apply(f6, planted)
    ctl(any(c), 'Sol %d: I_5 + (1/1000)(dX)^2(22)_Y is flagged: coords %s' % (sol, [str(x) for x in c]))

# ------------------------------------------------------------ C5 vacuity
print('\n== C5: vacuity ==')
rng = random.Random(7)
g = [F(rng.randint(-9, 9), rng.randint(1, 9)) for _ in B6]
c = V.apply(f6, g)
ctl(all(c), 'a generic vector of V_6 has all three coordinates nonzero: %s' % [str(x) for x in c])
span_rows = [V.to_row(v, i6) for v in hv6] + [V.to_row(V.d_mono(m), i6) for m in V.basis(5)]
h = [F(0)] * len(B6)
for r in span_rows:
    k = F(rng.randint(-9, 9), rng.randint(1, 9))
    h = [a + k * b for a, b in zip(h, r)]
c = V.apply(f6, h)
ctl(not any(c), 'a generic element of HV_6 + dV_5 has all three coordinates zero')
ctl(V.rank(span_rows) == 17, 'rank(HV_6 + dV_5) = %d (expected 17 of 20)' % V.rank(span_rows))

print()
if ctl_fail:
    print('RESULT: CONTROL FAILURE (%d)' % len(ctl_fail))
    sys.exit(1)
if pred_fail:
    print('RESULT: SEALED PREDICTION REFUTED (%d)' % len(pred_fail))
    sys.exit(2)
print('RESULT: P-spin3 and P-spin5 CONFIRMED; all controls fired')
sys.exit(0)
