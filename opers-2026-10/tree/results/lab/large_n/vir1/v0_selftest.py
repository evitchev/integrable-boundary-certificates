"""VIR1 C0: engine self-test.  Exit 0 iff every check passes.

 a. my residue == repo residue_mixed_mono on all weight-4 x weight-6 and
    weight-4 x weight-4 basis pairs (symbolic in the loop power);
 b. my P4 == repo cyl_P4 at several curve points, all three solutions;
 c. Virasoro algebra of T_Y = (1/2)(11): [L_m, L_n] = (m-n) L_{m+n}
    + (c/12)(m^3-m) delta, c = loop value, on a set of test states;
 d. the Virasoro states are loop-independent in the monomial basis;
 e. (T_Y T_Y) = (1/4)(11)(11) + (1/2)(13)  (the sealed hand composite);
 f. counts: dim V_4 = 7, V_6 = 20; dim Q_4 = 1, dim Q_6 = 3.
"""
import sys
import os as _os_kit; _KITROOT = _os_kit.path.join(_os_kit.path.dirname(_os_kit.path.abspath(__file__)), '..', '..', '..', '..') + '/'  # KIT PATCH: repo root from __file__ (was an absolute path)
from fractions import Fraction as F

sys.dont_write_bytecode = True
sys.path.insert(0, _KITROOT + 'code')
import vir_lib as V

fails = []


def req(cond, msg):
    print(('PASS  ' if cond else 'FAIL  ') + msg, flush=True)
    if not cond:
        fails.append(msg)


# a, b: against the repo engine (read-only import, no bytecode written)
import mixed_engine as ME

req(sorted(V.basis(4)) == sorted(ME.gen_mixed_basis(4)), 'a0 weight-4 basis == repo (%d)' % len(V.basis(4)))
req(sorted(V.basis(6)) == sorted(ME.gen_mixed_basis(6)), 'a0 weight-6 basis == repo (%d)' % len(V.basis(6)))
bad = 0
npairs = 0
for p in V.basis(4):
    for q in V.basis(6) + V.basis(4):
        npairs += 1
        if V.ope(p, q, -1) != ME.residue_mixed_mono(p, q):
            bad += 1
req(bad == 0, 'a  residue == repo residue_mixed_mono on %d monomial pairs' % npairs)
bad = 0
for m in V.basis(5) + V.basis(8):
    if V.d_mono(m) != ME.d_mixed_mono(m):
        bad += 1
req(bad == 0, 'a2 d_mono == repo d_mixed_mono on weights 5 and 8')
bad = 0
for t in (2, F(-1, 3), F(1, 3), 3, F(5, 3), F(7, 2), -4, F(11, 5)):
    N, s = V.curve_point(t)
    for sol in (1, 2, 3):
        if V.P4(sol, N, s) != ME.cyl_P4(sol, N, s):
            bad += 1
req(bad == 0, 'b  P4 == repo cyl_P4 at 8 curve points x 3 solutions')

# c: Virasoro algebra
def add(u, v, cu=1, cv=1):
    out = {}
    for m, c in u.items():
        out[m] = out.get(m, 0) + cu * c
    for m, c in v.items():
        out[m] = out.get(m, 0) + cv * c
    return {m: c for m, c in out.items() if c}


bad = 0
nchk = 0
for loop in (F(7, 3), F(28), F(-2)):
    states = [dict(V.VAC)]
    states += [{m: F(1)} for m in V.y_monos(2) and [((), y) for y in V.y_monos(2)]]
    states += [{((), y): F(1)} for y in V.y_monos(4)]
    states += [{((), y): F(1)} for y in V.y_monos(3)]
    for st in states:
        for m in range(-3, 4):
            for n in range(-3, 4):
                if m == n:
                    continue
                lhs = add(V.L(m, V.L(n, st, loop), loop), V.L(n, V.L(m, st, loop), loop), 1, -1)
                rhs = {k: (m - n) * c for k, c in V.L(m + n, st, loop).items()}
                if m + n == 0:
                    rhs = add(rhs, st, 1, loop * F(m ** 3 - m, 12))
                nchk += 1
                if add(lhs, rhs, 1, -1):
                    bad += 1
req(bad == 0, 'c  Virasoro algebra with c = loop on %d (state, m, n) triples, loop in {7/3, 28, -2}' % nchk)

# d: loop independence of Virasoro states through weight 8
bad = 0
for w in range(2, 9):
    for p in V.vir_partitions(w):
        if V.vir_state(p, F(7, 3)) != V.vir_state(p, F(-11, 5)):
            bad += 1
req(bad == 0, 'd  Virasoro states L_{-n1}..L_{-nk}|0> loop-independent through weight 8')

# e: the composite
TT = dict(V.vir_state((2, 2), F(7, 3)))
expect = {((), ((1, 1), (1, 1))): F(1, 4), ((), ((1, 3),)): F(1, 2)}
req(TT == expect, 'e  (T_Y T_Y) = (1/4)(11)(11) + (1/2)(13): got %s' % TT)

# f: counts
f4, B4, i4 = V.quotient_functionals(4, V.hv_span(4))
f6, B6, i6 = V.quotient_functionals(6, V.hv_span(6))
req(len(B4) == 7 and len(B6) == 20, 'f  dim V_4 = %d, dim V_6 = %d' % (len(B4), len(B6)))
req(len(f4) == 1 and len(f6) == 3, 'f  dim Q_4 = %d, dim Q_6 = %d' % (len(f4), len(f6)))
print('functional on V_4 (basis order %s):' % B4)
for lam in f4:
    print('   ', {B4[i]: str(c) for i, c in enumerate(lam) if c})
print('functionals on V_6:')
for lam in f6:
    print('   ', {B6[i]: str(c) for i, c in enumerate(lam) if c})

print('SELFTEST', 'FAILED' if fails else 'PASSED')
sys.exit(1 if fails else 0)
