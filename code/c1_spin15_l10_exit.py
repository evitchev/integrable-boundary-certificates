"""The paper's printed spin-15 prediction, tested in-stack.

Section 6.6 of the submitted article: the tower member at spin sigma
lies in L(1,0) exactly at the A2(2) exponent spins (sigma = +-1 mod 6;
verified 5,7,11,13) and exits at sigma = 3 mod 6 -- observed at 3 and
9, PREDICTED at 15.  Preregistered expectations: notes (commit
bff6630, before this computation).

This certificate: (1) computes ker ad_{I5} at charge spin 15 (weight
16, one boson, even sector) with the same machinery and normalization
as the certified sigma<=9 range and the spin-14 vacancy certificate --
expected dim 1, extending this project's own certified tower range to
sigma=15; (2) tests the member's first-order conservation under the
four momenta: EXPECTED conserved under +-1/sqrt2 and +-2sqrt2, NOT
conserved under +-sqrt2 (the L(1,0) exit -- the printed prediction)
nor under the sqrt3 control.  Fail-closed on every expectation; a
dim != 1 kernel or a CONSERVED +-sqrt2 row fails loudly (the latter
would contradict the submitted paper and demands escalation, not
silent acceptance).
"""
import sys
from fractions import Fraction as F

import sympy as sp

from families import L
from launch_provenance import capture_launch
from ope_engine import Sector, gen_basis
from screening_fit import screening_conditions
from solve_charges import genuine_kernel
from pathlib import Path

failures = []


def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok:
        failures.append(name)


def acc(pairs):
    out = {}
    for mono, c in pairs:
        out[mono] = out.get(mono, 0) + c
    return {m: c for m, c in out.items() if c}


A2_P6 = acc([(L((3, 0), (3, 0)), F(124, 96)),
             (L((2, 0), (2, 0), (1, 0), (1, 0)), F(-6)),
             (L((1, 0), (2, 0), (1, 0), (2, 0)), -F(-68, 8)),
             (L((1, 0),) * 6, F(1))])

if __name__ == "__main__":
    _ROOT = Path(__file__).resolve().parent.parent
    _LAUNCH = capture_launch(_ROOT, Path(__file__).resolve())

    ker3 = genuine_kernel(A2_P6, gen_basis(4, 1, True),
                          Sector(9, 1, True), Sector(4, 1, True))
    check("control: I'3 unique at sigma=3 (reproduces the certified "
          "low-spin computation)", len(ker3) == 1)

    w = 16
    ker15 = genuine_kernel(A2_P6, gen_basis(w, 1, True),
                          Sector(w + 5, 1, True), Sector(w, 1, True))
    check(f"tower multiplicity at sigma=15 is {len(ker15)} == 1 (extends "
          "this stack's certified range from sigma<=9; a value >= 2 "
          "would be a finding demanding separate scrutiny)",
          len(ker15) == 1)
    if len(ker15) != 1:
        print("\nSPIN-15 TEST ABORTED (kernel not one-dimensional)")
        sys.exit(1)
    P15 = ker15[0]

    r2 = sp.sqrt(2)
    momenta = [("+1/sqrt2", (r2 / 2,), True),
               ("-1/sqrt2", (-r2 / 2,), True),
               ("+2sqrt2", (2 * r2,), True),
               ("-2sqrt2", (-2 * r2,), True),
               ("+sqrt2 (L(1,0) row)", (r2,), False),
               ("-sqrt2 (L(1,0) row)", (-r2,), False),
               ("sqrt3 (control)", (sp.sqrt(3),), False)]
    results = {}
    for name, al, expect_conserved in momenta:
        conds = [c for c in screening_conditions(P15, al, 1)
                 if sp.simplify(c) != 0]
        conserved = not conds
        results[name] = conserved
        check(f"spin-15 member under alpha = {name}: "
              f"{'conserved' if conserved else 'NOT conserved'} "
              f"(expected {'conserved' if expect_conserved else 'not'})",
              conserved == expect_conserved)

    print()
    if failures:
        print(f"SPIN-15 L(1,0)-EXIT TEST FAILED ({len(failures)}) -- if a "
              "+-sqrt2 row is the failure, the submitted paper's printed "
              "prediction is contradicted: ESCALATE, do not suppress")
        sys.exit(1)
    print("SPIN-15 L(1,0) EXIT CERTIFIED (tower multiplicity one at "
          "sigma=15 in this stack; the member is conserved under "
          "{+-1/sqrt2, +-2sqrt2}, NOT under +-sqrt2 -- the paper's "
          "printed sigma=3-mod-6 exit prediction CONFIRMED at 15)")
