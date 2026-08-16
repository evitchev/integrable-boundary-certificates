"""Validation of the symbolic-N invariant engine.

W1: Res(T, M) = dM inside the invariant algebra, several N samples.
W2: EXACT agreement with the component engine at N = 2, 3: expand
    invariant monomials into components and compare residues.
W3: KdV [I_3, I_5] = 0 at N = 1..12 -> identity in N by degree bound.
W4: O(N) family [I_3, I_5] = 0 at 12 samples -> identity in N.
W5: A_2^(2): unique I_7 commuting with I_5 at a generic curve point
    (validates eq. (34) + the engine end-to-end).
"""

from fractions import Fraction as F
from itertools import product as iproduct

from invariant_engine import (canon, canon_pair, d_inv_mono, fam_current,
                              gen_inv_basis, residue_inv_cur,
                              genuine_kernel_inv, InvSector, spin_inv,
                              kdv_P4, kdv_P6, oN_P4, oN_P6, a2_P6,
                              curve_point)
import ope_engine as comp


def check(label, ok):
    print(("PASS  " if ok else "FAIL  ") + label)
    if not ok:
        raise SystemExit(1)


def d_inv_cur(cur):
    out = {}
    for m, c in cur.items():
        for m2, k in d_inv_mono(m).items():
            v = out.get(m2, 0) + c * k
            if v:
                out[m2] = v
            else:
                out.pop(m2, None)
    return out


# ------------------------------------------------------------------ W1
T = fam_current([(F(1, 2), [(1, 1)])])
ok = True
for M in [{canon(((1, 2), (1, 2))): F(1)},
          {canon(((2, 2), (1, 1))): F(3, 7)},
          {canon(((1, 3),)): F(1), canon(((2, 2),)): F(-2)}]:
    for Nval in [F(1), F(2), F(5), F(11, 3)]:
        lhs = residue_inv_cur(T, M, Nval)
        rhs = d_inv_cur(M)
        if lhs != rhs:
            ok = False
check("W1: Res(T, M) = dM in the invariant algebra", ok)

# ------------------------------------------------------------------ W2
def inv_to_comp(mono, N):
    cur = {(): F(1)}
    for a, b in mono:
        new = {}
        for m, c in cur.items():
            for mu in range(N):
                mm = tuple(sorted(m + ((a, mu), (b, mu))))
                new[mm] = new.get(mm, 0) + c
        cur = new
    return cur


def comp_of_inv_cur(cur, N):
    out = {}
    for mono, c in cur.items():
        for m, k in inv_to_comp(mono, N).items():
            v = out.get(m, 0) + c * k
            if v:
                out[m] = v
            else:
                out.pop(m, None)
    return out


ok = True
tests = [
    (((1, 1), (1, 1)), ((1, 1),)),
    (((1, 2), (1, 2)), ((1, 1), (1, 1))),
    (((2, 2),), ((1, 3),)),
    (((1, 2), (1, 2)), ((1, 2), (1, 2))),
    (((1, 1), (1, 1), (1, 1)), ((2, 2), (1, 1))),
    (((3, 3),), ((1, 1), (1, 1), (1, 1))),
]
for N in (2, 3):
    for P, Q in tests:
        P = canon(P)
        Q = canon(Q)
        got = comp_of_inv_cur(
            residue_inv_cur({P: F(1)}, {Q: F(1)}, F(N)), N)
        want = comp.residue_cur(inv_to_comp(P, N), inv_to_comp(Q, N))
        want = {m: c for m, c in want.items() if c}
        if got != want:
            ok = False
            print(f"    mismatch at N={N}, P={P}, Q={Q}")
check("W2: invariant engine == component engine at N = 2, 3", ok)

# ------------------------------------------------------------------ W3
ok = True
tgt = InvSector(9)
for N in range(1, 13):
    res = residue_inv_cur(kdv_P4(), kdv_P6(N), F(N))
    if not tgt.is_total_derivative(res):
        ok = False
check("W3: KdV [I_3, I_5] = 0 at N = 1..12  (=> identity in N)", ok)

# ------------------------------------------------------------------ W4
ok = True
for N in range(1, 13):
    res = residue_inv_cur(oN_P4(N), oN_P6(N), F(N))
    if not tgt.is_total_derivative(res):
        ok = False
check("W4: O(N) [I_3, I_5] = 0 at N = 1..12  (=> identity in N)", ok)

# ------------------------------------------------------------------ W5
Nval, sval = curve_point(F(7))   # N = 1/2, s = -7/2: generic curve point
P6 = a2_P6(Nval, sval)
sols = genuine_kernel_inv(P6, gen_inv_basis(8), Nval)
check(f"W5: A2(2) at curve point t=7 (N={Nval}): unique I_7 "
      f"(found {len(sols)})", len(sols) == 1)

print("\nAll invariant-engine validations passed.")
