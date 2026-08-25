"""Preregistered spin-13 test of the lab hypothesis t_k = (2k+3)/(2k+1)
(results/cyl_spin13_preregistration.md, committed f5d6521 BEFORE any
weight-14 kernel was computed).

Prediction: a multiplicity jump of some cylindrical Solution at t = 11/9,
(N,s) = (-238/5,-297/5), at spin 13.  Outcome (2026-08-23): NOT
CONFIRMED -- all three Solutions have spin-13 kernel dimension 1 there,
equal to the preregistered baseline (the t = 2 probe; the generic
spin-13 multiplicity itself is not yet established symbolically -- that
is the global weight-14 stratification's job).  Not a refutation (the sequence has
a gap at spin 7, so the point may appear at a higher spin).

The controls show the computation sees jumps at weight 14: (-2,-9)
decoupling gives dim 2 for Solutions 2 and 3; every single-spin point
is back to dim 1; the negative controls t = 5/4, 7/6 sit at the baseline.

Nineteen exact weight-14 kernels (929 densities each), ~13 min per
kernel, run in parallel (all cores); expensive.  Fail-closed: every
one of the 19 dimensions is pinned to the preregistered expectation
(prediction point: pinned to the observed non-jump)."""
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


PRED = point("11/9")
require(PRED == (F(-238, 5), F(-297, 5)), "predicted point arithmetic")
EXPECT = {}                                       # (sol, N, s) -> expected dim
for sol in (1, 2, 3):
    EXPECT[(sol,) + PRED] = 1                     # observed: no jump (pinned)
    EXPECT[(sol,) + point(2)] = 1                 # preregistered baseline (t = 2)
    EXPECT[(sol,) + point("5/4")] = 1             # negative control
    EXPECT[(sol,) + point("7/6")] = 1             # negative control
    EXPECT[(sol, F(-2), F(-9))] = 2 if sol > 1 else 1   # decoupling persists
EXPECT[(2, F(-143, 4), F(-189, 4))] = 1           # spin-11-only
EXPECT[(3, F(-24), F(-35))] = 1                   # spin-9-only
EXPECT[(2, F(-25, 2), F(-45, 2))] = 1             # spin-5-only
EXPECT[(3, F(-25, 2), F(-45, 2))] = 1
require(len(EXPECT) == 19, f"job inventory {len(EXPECT)} != 19")

if __name__ == "__main__":
    jobs = list(EXPECT)
    with Pool(min(len(jobs), os.cpu_count() or 1)) as pool:
        dims = pool.map(kdim, jobs)
    failures = []
    for job, d in zip(jobs, dims):
        ok = d == EXPECT[job]
        tag = "PREDICTION" if job[1:] == PRED else "control"
        print(f"{'PASS' if ok else 'FAIL'}  Solution {job[0]} at ({job[1]},{job[2]}) "
              f"[{tag}]: spin-13 kernel dim {d} (expected {EXPECT[job]})", flush=True)
        if not ok:
            failures.append(job)
    pred_dims = [dims[jobs.index((sol,) + PRED)] for sol in (1, 2, 3)]
    baseline = [dims[jobs.index((sol,) + point(2))] for sol in (1, 2, 3)]
    print()
    if failures:
        print(f"SPIN-13 PREDICTION TEST FAILED TO MATCH PINS ({len(failures)})")
        sys.exit(1)
    require(pred_dims == baseline == [1, 1, 1], f"final consistency {pred_dims} {baseline}")
    print("SPIN-13 PREDICTION TEST CERTIFIED (t = 11/9 NOT CONFIRMED: dims "
          f"{pred_dims} = preregistered baseline {baseline} (t = 2 probe) on all three sheets at spin 13; "
          "positive controls (-2,-9) dim 2 prove the weight-14 computation sees "
          "jumps; single-spin points and negative controls at the baseline; not a refutation; "
          "generic spin-13 multiplicity itself awaits the global weight-14 stratification)")
