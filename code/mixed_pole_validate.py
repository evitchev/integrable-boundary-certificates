"""Validation of mixed_pole (M1 spec, three legs; fail-closed)."""
import random
import sys
from fractions import Fraction as F

from mixed_engine import gen_mixed_basis, residue_mixed_mono, canon_mixed, d_mixed_mono
from mixed_pole import ope_pole_mixed_mono, ope_pole_mixed_cur
from ope_poles import ope_pole_mono
from super_engine import canonicalize

fails = []
def check(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg, flush=True)
    if not cond:
        fails.append(msg)

# (a) k = 1 reproduces the residue on every monomial pair up to weight 4, and random weight-6 pairs
pairs = [(P, Q) for P in gen_mixed_basis(4) for Q in gen_mixed_basis(4)]
random.seed(7)
b6 = gen_mixed_basis(6)
pairs += [(random.choice(b6), random.choice(b6)) for _ in range(40)]
bad = 0
for P, Q in pairs:
    if ope_pole_mixed_mono(P, Q, 1) != residue_mixed_mono(P, Q):
        bad += 1
check(bad == 0, f"(a) k=1 == residue_mixed_mono on {len(pairs)} pairs (mismatches: {bad})")

# (b) T = 1/2 (dX)^2 + 1/2 (dY.dY): T(z)T(w) ~ (c/2)/(z-w)^4 + 2T/(z-w)^2 + dT/(z-w), c = N + 1
T = {((1, 1), ()): F(1, 2), ((), ((1, 1),)): F(1, 2)}
for Nval in (F(3), F(-24), F(7, 2)):
    p4 = ope_pole_mixed_cur(T, T, 4, Nval)
    p3 = ope_pole_mixed_cur(T, T, 3, Nval)
    p2 = ope_pole_mixed_cur(T, T, 2, Nval)
    ident = canon_mixed([], [])
    check(p4 == {ident: (Nval + 1) / 2}, f"(b) N={Nval}: pole-4 of (T,T) = c/2 with c = N+1 (got {p4})")
    check(p3 == {}, f"(b) N={Nval}: pole-3 of (T,T) vanishes")
    check(p2 == {m: 2 * c for m, c in T.items()}, f"(b) N={Nval}: pole-2 of (T,T) = 2T")
    p1 = ope_pole_mixed_cur(T, T, 1, Nval)
    dT = {}
    for m, c in T.items():
        for dm, dc in d_mixed_mono(m).items():
            dT[dm] = dT.get(dm, F(0)) + c * F(dc)
    dT = {m: c for m, c in dT.items() if c}
    check(p1 == dT, f"(b) N={Nval}: pole-1 of (T,T) = dT (Codex round 17: the marker said "
          f"dT-consistent without checking it)")

# (c) pure-X inputs agree with the single-boson ope_poles at every pole order
def to_super(xs):
    return canonicalize(tuple((0, m, 0) for m in xs))[0]
def from_super(mono):
    return canon_mixed([m for (_, m, _) in mono], [])
bad = 0; count = 0
for P in gen_mixed_basis(4):
    for Q in gen_mixed_basis(4):
        if P[1] or Q[1]:
            continue
        for k in range(0, 9):
            mixed = ope_pole_mixed_mono(P, Q, k)
            got = {m: d.get(0, F(0)) for m, d in mixed.items()}
            sup = ope_pole_mono(to_super(P[0]), to_super(Q[0]), k)
            exp = {}
            for sm, c in sup.items():
                exp[from_super(sm)] = exp.get(from_super(sm), 0) + F(c)
            exp = {m: c for m, c in exp.items() if c}
            got = {m: c for m, c in got.items() if c}
            count += 1
            if got != exp:
                bad += 1
check(bad == 0, f"(c) pure-X poles k=0..8 agree with ope_poles on {count} (pair, k) cases (mismatches: {bad})")
print()
if fails:
    print(f"MIXED_POLE VALIDATION: {len(fails)} FAILURE(S)"); sys.exit(1)
print("MIXED_POLE VALIDATED (k=1 == residue; (T,T) poles give c = N+1, 0, 2T, dT; pure-X == single-boson)")
