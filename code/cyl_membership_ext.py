"""Ladder extension for the cylindrical membership computation:
I9 (weight 10) and I11 (weight 12) against the finite W(2,4,6,8)
candidate left open by cyl_membership.py.

Decisive questions (outcomes NOT asserted -- genuinely open):
  - I9  in A<=8 ?   (if yes for all solutions, the finite candidate
                     survives its first new test)
  - I11 in A<=8, A<=10 ?  (not-in A<=10 = the ladder continues to
                     its own weight; in A<=8 = finite candidate)
Own-weight membership is trivially true (the A<=w span at weight w
is the full sector) and is not tested.

Mechanical integrity fails closed: kernel dimensions must be 1, and
span controls at BOTH working weights must pass (weight 10: (33)(22)
positive, bare (55) negative for A<=8; weight 12, mirroring the
decisive I11 tests: (55)(11) in A<=10, (66) not in A<=10, (66) in
A<=12).  The success predicate also requires the mathematically
redundant I11-not-in-A<=8 as a regression detector on the
machinery.  Representatives are exact genuine kernels of
ad_I3 at the curve point t = 2; spans are t-independent, so each
non-membership is generic by semicontinuity."""
import sys
from fractions import Fraction as F

import collapsed
import mixed_engine
mixed_engine.residue_mixed_mono = collapsed.residue_mixed_mono_smart
from cyl_membership import membership, mult
from invariant_engine import curve_point
from mixed_engine import (cyl_P4, d_mixed_mono, gen_mixed_basis,
                          residue_mixed_cur, spin_mixed)
from sparse_linalg import genuine_kernel_sparse

A8 = [(1, 1), (2, 2), (3, 3), (4, 4)]
A10 = A8 + [(5, 5)]


def main():
    failures = []

    p11 = {((), ((1, 1),)): F(1)}
    p33 = {((), ((3, 3),)): F(1)}
    p22 = {((), ((2, 2),)): F(1)}
    p55 = {((), ((5, 5),)): F(1)}
    p66 = {((), ((6, 6),)): F(1)}
    A12 = A10 + [(6, 6)]
    controls = [
        ("(33)(22) in A<=8 (w10)", mult(p33, p22), 10, A8, True),
        ("(55) in A<=8 (w10)", p55, 10, A8, False),
        ("(55) in A<=10 (w10)", p55, 10, A10, True),
        # weight-12 controls mirroring the decisive I11 tests
        ("(55)(11) in A<=10 (w12)", mult(p55, p11), 12, A10, True),
        ("(66) in A<=10 (w12)", p66, 12, A10, False),
        ("(66) in A<=12 (w12)", p66, 12, A12, True),
    ]
    for name, cur, w, gens, expect in controls:
        got, _ = membership(cur, w, gens)
        print(f"control {'PASS' if got == expect else 'FAIL'}: {name} "
              f"-> {got} (expect {expect})", flush=True)
        if got != expect:
            failures.append(f"control {name}")

    N, s = curve_point(F(2))
    table = {}
    for sol in (1, 2, 3):
        P4 = cyl_P4(sol, N, s)
        for sigma in (9, 11):
            g = genuine_kernel_sparse(
                P4, gen_mixed_basis(sigma + 1), N, residue_mixed_cur,
                gen_mixed_basis, d_mixed_mono, spin_mixed)
            if len(g) != 1:
                print(f"FAILURES: kernel dim {len(g)} != 1 at "
                      f"sol {sol} sigma {sigma}", flush=True)
                print("CYL LADDER EXTENSION FAILED", flush=True)
                sys.exit(1)
            cur = g[0]
            w = sigma + 1
            tests = [("A<=8", A8)]
            if sigma == 11:
                tests.append(("A<=10", A10))
            for name, gens in tests:
                ok, _ = membership(cur, w, gens)
                table[(sol, sigma, name)] = ok
                print(f"sol {sol} I{sigma} (weight {w}) in {name}: "
                      f"{ok}", flush=True)

    print()
    if failures:
        print(f"FAILURES: {failures}", flush=True)
        print("CYL LADDER EXTENSION FAILED", flush=True)
        sys.exit(1)

    i9 = {table[(sol, 9, "A<=8")] for sol in (1, 2, 3)}
    i11_8 = {table[(sol, 11, "A<=8")] for sol in (1, 2, 3)}
    i11_10 = {table[(sol, 11, "A<=10")] for sol in (1, 2, 3)}
    # i11_8 is implied by i11_10 mathematically (A<=8 within A<=10)
    # but required here as a regression detector on the machinery
    if i9 == {False} and i11_8 == {False} and i11_10 == {False}:
        print("VERDICT: the ladder CONTINUES -- I9 needs spin-10 and "
              "I11 needs spin-12 singlet content; the finite "
              "W(2,4,6,8) candidate is excluded through the computed "
              "range.  (Finite-range statement: a still-larger finite "
              "closure is not excluded by any finite computation.)",
              flush=True)
    elif i9 == {True} and i11_8 == {True}:
        print("VERDICT: the towers CLOSE over the spin<=8 generators "
              "at I9 and I11 -- the finite W(2,4,6,8) candidate "
              "SURVIVES; an LV-style matching against a W(2,4,6,8) "
              "structure is well-posed.", flush=True)
    else:
        print(f"VERDICT: mixed pattern -- inspect the table: "
              f"I9 in A<=8: {sorted(i9)}, I11 in A<=8: "
              f"{sorted(i11_8)}, I11 in A<=10: {sorted(i11_10)}",
              flush=True)
    print("CYL LADDER EXTENSION COMPLETE", flush=True)


if __name__ == "__main__":
    main()
