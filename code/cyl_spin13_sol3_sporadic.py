"""Sporadic spin-13 jump of cylindrical Solution 3 at t = 9/7,
(N,s) = (-143/4,-189/4) -- the first sporadic spin-13 point.

Found by the global weight-14 stratification lab (2026-08-27): the
1082x365 weight-14 matrix of Solution 3 has generic rank 364 over
Q(t), and the rank check kept exactly two determinant roots, t = 3
(the (-2,-9) decoupling, shared with Solution 2) and t = 9/7 --
nullity 2 there (jump +1 over the generic multiplicity 1).  Open-stack
confirmation the same day: cyl_spin13_kdim.py 3 -143/4 -189/4 gave
dim 2 in 777 s.  This certificate re-establishes the finding
fail-closed and adds the missing sheet coverage at the point:

  Solution 1: spin-13 kernel dim 1 (lab: no surviving root at 9/7)
  Solution 2: spin-13 kernel dim 1 (already pinned by
              cyl_spin13_prediction -- the point is spin-11-only on
              sheet 2; recomputed here so the trio is self-contained)
  Solution 3: spin-13 kernel dim 2 (THE FINDING)

Reading: the preregistered t_k = (2k+3)/(2k+1) ladder predicted
spin-13 action at t_4 = 11/9; that point is certified clean on all
sheets (cyl_spin13_prediction).  The sporadic point instead recurs at
t_3 = 9/7, which already carries the Solution-2 spin-11 jump -- the
same curve point jumping at different spins on different sheets.  The
observed jumping locus through spin 13 stays inside {3} u {t_k}, but
the spin assignment is NOT the naive ladder.  Three exact weight-14
kernels (929 densities each), ~13 min per kernel, in parallel;
expensive."""
import os
import sys
from fractions import Fraction as F
from multiprocessing import Pool

from mixed_engine import cyl_P4, gen_mixed_basis, genuine_kernel_mixed


def require(cond, msg):
    """Fail-closed guard that survives python -O (no bare assert)."""
    if not cond:
        print(f"FAIL  {msg}")
        sys.exit(1)

W = 14


def kdim(job):
    sol, N, s = job
    require(s * s == (N - 25) * (N - 1), f"({N},{s}) is not on the curve")
    return len(genuine_kernel_mixed(cyl_P4(sol, N, s), gen_mixed_basis(W), N))


def point(t):
    t = F(t)
    return (t * t - 25) / (t * t - 1), -24 * t / (t * t - 1)


PT = point("9/7")
require(PT == (F(-143, 4), F(-189, 4)), "point arithmetic at t = 9/7")

EXPECT = {                                        # (sol, N, s) -> expected dim
    (1,) + PT: 1,                                 # sheet-1 coverage (new)
    (2,) + PT: 1,                                 # spin-11-only on sheet 2
    (3,) + PT: 2,                                 # THE SPORADIC JUMP
}
require(len(EXPECT) == 3, f"job inventory {len(EXPECT)} != 3")

if __name__ == "__main__":
    jobs = list(EXPECT)
    with Pool(min(len(jobs), os.cpu_count() or 1)) as pool:
        dims = pool.map(kdim, jobs)
    failures = []
    for job, d in zip(jobs, dims):
        ok = d == EXPECT[job]
        tag = "FINDING" if job[0] == 3 else "control"
        print(f"{'PASS' if ok else 'FAIL'}  Solution {job[0]} at ({job[1]},{job[2]}) "
              f"[{tag}]: spin-13 kernel dim {d} (expected {EXPECT[job]})", flush=True)
        if not ok:
            failures.append(job)
    print()
    if failures:
        print(f"SPIN-13 SPORADIC POINT TEST FAILED TO MATCH PINS ({len(failures)})")
        sys.exit(1)
    print("SPIN-13 SPORADIC POINT CERTIFIED (Solution 3 spin-13 kernel dim 2 at "
          "t = 9/7, (N,s) = (-143/4,-189/4), against dim 1 on sheets 1 and 2 -- "
          "the first sporadic spin-13 jump; the t_k ladder's spin-13 slot t = 11/9 "
          "is certified clean by cyl_spin13_prediction, so the ladder in its naive "
          "form is refuted while the observed jumping locus through spin 13 stays "
          "inside {3} u {t_k}; same-point cross-sheet recurrence: t_3 jumps at "
          "spin 11 on sheet 2 and spin 13 on sheet 3)")
