"""Spin-content diagnosis of the three Z2 x O(N) cylindrical families
(hep-th/0404195 Sec. 3.3) at general N.

MV1: Res(T_total, M) = dM (translation) at several N.
MV2: mixed engine == component engine at N = 1, 2 (exact expansion).
MV3: each solution's I_3 admits a commuting I_5 (validates eqs. 53/55/57
     transcription -- a wrong I_3 has empty kernel).
Scan: genuine kernel dims of ad_{I_3} for sigma = 2..11 at rational
points of the double cover s^2 = (N-25)(N-1) (both sheets).

Interpretation: the set of spins supported diagnoses the exponents of
the underlying affine algebra (all odd <-> A_1^(1)-type paperclip
towers; 6n+-1 <-> A_2^(2); etc.).
"""

import sys
from fractions import Fraction as F
from multiprocessing import Pool

from mixed_engine import (canon_mixed, gen_mixed_basis, d_mixed_mono,
                          residue_mixed_mono, residue_mixed_cur,
                          MixedSector, genuine_kernel_mixed, cyl_P4, _M)
from invariant_engine import curve_point
import ope_engine as comp


def check(label, ok):
    print(("PASS  " if ok else "FAIL  ") + label)
    if not ok:
        raise SystemExit(1)


def mixed_to_comp(mono, N):
    """Expand a mixed monomial into component monomials:
    X -> field 0, Y -> fields 1..N."""
    xs, yp = mono
    cur = {tuple(sorted((o, 0) for o in xs)): F(1)}
    for a, b in yp:
        new = {}
        for m, c in cur.items():
            for mu in range(1, N + 1):
                mm = tuple(sorted(m + ((a, mu), (b, mu))))
                new[mm] = new.get(mm, 0) + c
        cur = new
    return cur


def comp_of_cur(cur, N):
    out = {}
    for mono, c in cur.items():
        for m, k in mixed_to_comp(mono, N).items():
            v = out.get(m, 0) + c * k
            if v:
                out[m] = v
            else:
                out.pop(m, None)
    return out


def validations():
    # MV1
    T = {_M((1, 1)): F(1, 2), _M((), [(1, 1)]): F(1, 2)}
    ok = True
    for M in [{_M((1, 2), [(1, 1)]): F(1)},
              {_M((), [(1, 2), (1, 1)]): F(3)},
              {_M((2, 2)): F(1), _M((1, 1, 1, 1)): F(-5)}]:
        for Nval in [F(1), F(3), F(9, 2)]:
            lhs = residue_mixed_cur(T, M, Nval)
            rhs = {}
            for m, c in M.items():
                for m2, k in d_mixed_mono(m).items():
                    rhs[m2] = rhs.get(m2, 0) + c * k
            rhs = {m: c for m, c in rhs.items() if c}
            if lhs != rhs:
                ok = False
    check("MV1: Res(T_total, M) = dM in the mixed algebra", ok)

    # MV2
    tests = [
        (_M((1, 1)), _M((1, 1))),
        (_M((1, 1)), _M((), [(1, 2), (1, 2)])),
        (_M((1, 2), [(1, 1)]), _M((1, 1), [(1, 1)])),
        (_M((), [(1, 1), (1, 1)]), _M((2, 2))),
        (_M((1, 1, 1, 1)), _M((1, 1), [(2, 2)])),
        (_M((), [(1, 2), (1, 2)]), _M((), [(1, 1), (1, 1)])),
    ]
    ok = True
    for N in (1, 2):
        for P, Q in tests:
            got = comp_of_cur(
                residue_mixed_cur({P: F(1)}, {Q: F(1)}, F(N)), N)
            want = comp.residue_cur(mixed_to_comp(P, N),
                                    mixed_to_comp(Q, N))
            want = {m: c for m, c in want.items() if c}
            if got != want:
                ok = False
                print(f"    mismatch N={N} P={P} Q={Q}")
    check("MV2: mixed engine == component engine at N = 1, 2", ok)

    # MV3
    for sol in (1, 2, 3):
        for t in (F(2), F(7)):
            Nval, sval = curve_point(t)
            gen = genuine_kernel_mixed(cyl_P4(sol, Nval, sval),
                                       gen_mixed_basis(6), Nval)
            check(f"MV3: Solution {sol} admits I_5 at t={t} "
                  f"(dim {len(gen)})", len(gen) >= 1)


POINTS = [("t=2 (N=-7)", curve_point(2)),
          ("t=3 (N=-2)", curve_point(3)),
          ("t=7 (N=1/2)", curve_point(7)),
          ("t=-2 (sheet 2)", curve_point(-2))]
SIGMAS = list(range(2, 12))


def scan_task(args):
    sol, sigma = args
    cand = gen_mixed_basis(sigma + 1)
    dims = []
    for _, (Nval, sval) in POINTS:
        gen = genuine_kernel_mixed(cyl_P4(sol, Nval, sval), cand, Nval)
        dims.append(len(gen))
    return sol, sigma, MixedSector(sigma + 1).q_dim(), dims


if __name__ == "__main__":
    validations()
    tasks = [(sol, s) for sol in (1, 2, 3) for s in SIGMAS]
    with Pool(16) as pool:
        results = pool.map(scan_task, tasks)
    table = {(sol, s): (q, dims) for sol, s, q, dims in results}
    print("\ngenuine ker ad_{I_3} dims; rows = parameter points")
    hdr = " ".join(f"{s:2d}" for s in SIGMAS)
    for sol in (1, 2, 3):
        print(f"\nSolution {sol}:            sigma | {hdr}")
        for i, (lbl, _) in enumerate(POINTS):
            row = " ".join(f"{table[(sol, s)][1][i]:2d}" for s in SIGMAS)
            print(f"  {lbl:<22}      | {row}")
    print("\n  dim Q^sigma (Z2xO(N) sector) | " +
          " ".join(f"{table[(1, s)][0]:2d}" for s in SIGMAS))
