"""Preregistered spin-15 (weight-16) probe of the two-ladder hypothesis H1
(results/cyl_spin15_preregistration.md, committed 43a2f65 BEFORE any
weight-16 cylindrical kernel was computed).

H1: sheet 3 jumps at t_k exactly at spin 4k+1, sheet 2 at t_k (odd k)
exactly at spin 3k+2, all jumps single-spin.  H1 predicted NOTHING at
spin 15.  Outcome (2026-08-27, eight kernels 2.3-2.9 h each): all eight
preregistered dimensions are 1 -- confirmed at the eight probed points (a null
confirmation: consistent with H1, against the echo reading; the
decisive double positive prediction sits at spin 17, weight 18).

  baselines t=2 on all sheets: 1,1,1 (generic spin-15 multiplicity
    probe; the symbolic statement awaits the weight-16 stratification)
  t3=9/7:  sheet 2 (spin-11-only) 1; sheet 3 (spin-13 sporadic point,
    cyl_spin13_sol3_sporadic) back to 1 -- single-spin like its
    predecessors
  t4=11/9: sheets 2,3 at 1 (H1: sheet 3 activates at 17, sheet 2 never
    at an odd k... k=4 is even)
  t2=7/5:  sheet 3 (spin-9-only) 1

Fail-closed: every dimension pinned to the preregistered expectation;
eight exact weight-16 kernels (2211 densities), in parallel; expensive."""
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

W = 16


def kdim(job):
    sol, N, s = job
    require(s * s == (N - 25) * (N - 1), f"({N},{s}) is not on the curve")
    return len(genuine_kernel_mixed(cyl_P4(sol, N, s), gen_mixed_basis(W), N))


def point(t):
    t = F(t)
    return (t * t - 25) / (t * t - 1), -24 * t / (t * t - 1)


T2, T3, T4 = point("7/5"), point("9/7"), point("11/9")
require(T2 == (F(-24), F(-35)), "t2 arithmetic")
require(T3 == (F(-143, 4), F(-189, 4)), "t3 arithmetic")
require(T4 == (F(-238, 5), F(-297, 5)), "t4 arithmetic")
require(point(2) == (F(-7), F(-16)), "baseline arithmetic")

EXPECT = {                                        # (sol, N, s) -> expected dim
    (1,) + point(2): 1, (2,) + point(2): 1, (3,) + point(2): 1,
    (2,) + T3: 1, (3,) + T3: 1,
    (2,) + T4: 1, (3,) + T4: 1,
    (3,) + T2: 1,
}
require(len(EXPECT) == 8, f"job inventory {len(EXPECT)} != 8")

if __name__ == "__main__":
    jobs = list(EXPECT)
    with Pool(min(len(jobs), os.cpu_count() or 1)) as pool:
        dims = pool.map(kdim, jobs)
    failures = []
    for job, d in zip(jobs, dims):
        ok = d == EXPECT[job]
        tag = "baseline" if job[1:] == point(2) else "H1 probe"
        print(f"{'PASS' if ok else 'FAIL'}  Solution {job[0]} at ({job[1]},{job[2]}) "
              f"[{tag}]: spin-15 kernel dim {d} (expected {EXPECT[job]})", flush=True)
        if not ok:
            failures.append(job)
    print()
    if failures:
        print(f"SPIN-15 LADDER PROBE FAILED TO MATCH PINS ({len(failures)})")
        sys.exit(1)
    print("SPIN-15 LADDER PROBE CERTIFIED (H1 two-ladder hypothesis predicted nothing "
          "at spin 15: all eight preregistered weight-16 kernel dimensions are 1 -- "
          "baselines on all sheets, t3 on sheets 2,3 (the spin-13 sporadic point is "
          "single-spin), t4 on sheets 2,3, t2 on sheet 3; null confirmation, the "
          "decisive double positive prediction is at spin 17)")
