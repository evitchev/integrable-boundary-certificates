"""Spin-17 (weight-18) double prediction of the two-ladder hypothesis,
preregistered in results/cyl_spin17_preregistration.md (8444876) BEFORE
any weight-18 kernel: sheet 3 jumps at t4 = 11/9 at spin 4*4+1 = 17 and
sheet 2 jumps at t5 = 13/11 -- a point that had never jumped at any
spin -- at spin 3*5+2 = 17; nine controls at 1.

Status (2026-08-29): fast-path outcome obtained -- all eleven
dimensions as predicted by lab/fast_kernel2.py (modular + CRT + exact
identities verified in this engine; lab record under
results/lab/fast_kernel2/); this pure-Python certificate is the
PREREGISTERED independent method, its eleven kernels (24-30 h each)
running on the worker nodes, ten of eleven in and agreeing at the time
of writing.  Per the pre-result addendum: the
predicted slots are strongly supported; the universal clauses of the
hypothesis are not established by eleven points (its single-spin
clause is already refuted by the t1 recurrence at spin 15).

Pins:   (3, -238/5,-297/5) = 2   P1      (2, -119/2,-143/2) = 2   P2
        (1|2|3, -7,-16) = 1 baselines     (2, t4) = 1   (3, t5) = 1
        (2, t3) = 1   (3, t3) = 1   (1, t4) = 1   (1, t5) = 1
Eleven exact weight-18 kernels (5106 densities, ~20 GB each);
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

W = 18
POOL = int(os.environ.get("IB_KERNEL_POOL", "5"))   # ~20 GB per kernel


def kdim(job):
    sol, N, s = job
    require(s * s == (N - 25) * (N - 1), f"({N},{s}) is not on the curve")
    return len(genuine_kernel_mixed(cyl_P4(sol, N, s), gen_mixed_basis(W), N))


def point(t):
    t = F(t)
    return (t * t - 25) / (t * t - 1), -24 * t / (t * t - 1)


T3, T4, T5 = point("9/7"), point("11/9"), point("13/11")
require(T3 == (F(-143, 4), F(-189, 4)), "t3 arithmetic")
require(T4 == (F(-238, 5), F(-297, 5)), "t4 arithmetic")
require(T5 == (F(-119, 2), F(-143, 2)), "t5 arithmetic")
require(point(2) == (F(-7), F(-16)), "baseline arithmetic")

EXPECT = {                                        # (sol, N, s) -> expected dim
    (3,) + T4: 2, (2,) + T5: 2,                   # P1, P2
    (1,) + point(2): 1, (2,) + point(2): 1, (3,) + point(2): 1,
    (2,) + T4: 1, (3,) + T5: 1,
    (2,) + T3: 1, (3,) + T3: 1,
    (1,) + T4: 1, (1,) + T5: 1,
}
require(len(EXPECT) == 11, f"job inventory {len(EXPECT)} != 11")

if __name__ == "__main__":
    jobs = list(EXPECT)
    with Pool(min(len(jobs), POOL, os.cpu_count() or 1)) as pool:
        dims = pool.map(kdim, jobs)
    failures = []
    for job, d in zip(jobs, dims):
        ok = d == EXPECT[job]
        tag = "P1" if job == (3,) + T4 else "P2" if job == (2,) + T5 else "control"
        print(f"{'PASS' if ok else 'FAIL'}  Solution {job[0]} at ({job[1]},{job[2]}) "
              f"[{tag}]: spin-17 kernel dim {d} (expected {EXPECT[job]})", flush=True)
        if not ok:
            failures.append(job)
    print()
    if failures:
        print(f"SPIN-17 DOUBLE PREDICTION FAILED TO MATCH PINS ({len(failures)})")
        sys.exit(1)
    print("SPIN-17 DOUBLE PREDICTION CERTIFIED (sheet 3 jumps at t4 = 11/9 and "
          "sheet 2 at t5 = 13/11 at spin 17, both dim 2, exactly as preregistered; "
          "the nine controls at dim 1; the predicted slots of the two ladders are "
          "supported, their universal clauses are not thereby established)")
