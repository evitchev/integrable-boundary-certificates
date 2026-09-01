"""Odd-q controls (preregistered results/cyl_oddq_controls_preregistration.md).

The (1,q) identity reads the certified jumping locus through spin 17
as {beta^2 = 1/q, q EVEN}.  This certificate pins the odd-q points
t = 2, 3/2, 4/3 (q = 3, 5, 7) at spins 9, 11, 13 on Solutions 2 and 3
-- 18 exact kernels, every dimension expected 1 (generic).  A
dimension > 1 is a FINDING (a jump at odd q), reported, not re-pinned.
Marker: ODD-Q CONTROLS CERTIFIED."""
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


def point(t):
    t = F(t)
    return (t * t - 25) / (t * t - 1), -24 * t / (t * t - 1)


POINTS = {"2": (F(-7), F(-16)), "3/2": (F(-91, 5), F(-144, 5)), "4/3": (F(-209, 7), F(-288, 7))}
for tt, ns in POINTS.items():
    require(point(tt) == ns, f"point arithmetic at t = {tt}")
    require(ns[1] ** 2 == (ns[0] - 25) * (ns[0] - 1), f"({ns[0]},{ns[1]}) on the curve")

JOBS = [(w, sol, tt) for w in (10, 12, 14) for sol in (2, 3) for tt in POINTS]
require(len(JOBS) == 18, f"job inventory {len(JOBS)} != 18")


def kdim(job):
    w, sol, tt = job
    N, s = POINTS[tt]
    return len(genuine_kernel_mixed(cyl_P4(sol, N, s), gen_mixed_basis(w), N))


Q_OF = {"2": 3, "3/2": 5, "4/3": 7}

if __name__ == "__main__":
    with Pool(min(9, os.cpu_count() or 2)) as pool:
        dims = pool.map(kdim, JOBS)
    bad = []
    for job, d in zip(JOBS, dims):
        w, sol, tt = job
        if d != 1:
            bad.append((job, d))
        print(f"{'PASS' if d == 1 else 'FAIL'}  Solution {sol}, t = {tt} "
              f"(q = {Q_OF[tt]}, beta^2 = 1/{Q_OF[tt]}), spin {w - 1}: "
              f"kernel dim {d} (expected 1, generic)", flush=True)
    if bad:
        print(f"ODD-Q CONTROLS: {len(bad)} FINDING(S) -- a jump at odd q, "
              f"refuting the even-q clause: {bad}")
        sys.exit(1)
    print("ODD-Q CONTROLS CERTIFIED (18 exact kernels: t = 2, 3/2, 4/3 "
          "(q = 3, 5, 7) at spins 9, 11, 13 on Solutions 2 and 3, every "
          "dimension 1 -- the odd-q points are generic exactly where their "
          "even-q neighbours jump)")
