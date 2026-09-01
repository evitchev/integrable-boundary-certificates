"""Spin-15 (weight-16) stratification of the cylindrical curve: the exact
kernel dimensions behind the complete weight-16 residue-scan
stratification (2026-08-28; pre-run expectations bd4aa59).

Residue scans over GF(10007) (every non-pole residue) and GF(10009):
sheet 1 no rank drop; sheets 2 and 3 drop exactly at t = 3 (the (-2,-9)
all-spin decoupling) and at t = 5/3, (N,s) = (-25/2,-45/2) -- the point
where Solutions 2 and 3 merge, which had jumped at spin 5 and returned
to dim 1 at spins 7, 9, 11, 13.  Exact open-stack kernels: dim 2 there
on both merged sheets (13257 s / 12422 s), independently reproduced by
the modular+CRT+exact-verification path.  This certificate pins the
exact dimensions the stratification rests on:

    (2, -25/2,-45/2) = 2   (3, -25/2,-45/2) = 2    THE RECURRENCE
    (2,  -2,  -9  ) = 2   (3,  -2,  -9  ) = 2    all-spin decoupling
    (1, -25/2,-45/2) = 1   (1,  -2,  -9  ) = 1    sheet-1 controls at these two points

(the four lower-row pins were preregistered from the scan before their
exact kernels finished).  Baselines t = 2 on all sheets = 1 and the
ladder-probe points are pinned by cyl_spin15_ladder_probe.  SCOPE
(Codex round 5): this certifies exact dimensions at the two lifted
rational points; the residue scans over GF(10007) (every non-pole
residue) and GF(10009) at the lifted point detect no other reduction,
which is evidence for -- not a characteristic-zero proof of -- the
absence of further spin-15 points (an algebraic point, or a rational
one with unsuitable reduction, could escape a single prime).  Reading:
first recurrence at one (sheet, point): t1 = 5/3 jumps at spins 5 and
15; the single-spin clause of the two-ladder hypothesis is refuted, its
ladder clauses and the spin-17 double prediction are unchanged.
Six exact weight-16 kernels (2211 densities), 2.5-3.7 h each, in
parallel; expensive."""
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


T1 = point("5/3")
require(T1 == (F(-25, 2), F(-45, 2)), "t1 arithmetic")
require(point(3) == (F(-2), F(-9)), "decoupling point arithmetic")
DEC = (F(-2), F(-9))

EXPECT = {                                        # (sol, N, s) -> expected dim
    (2,) + T1: 2, (3,) + T1: 2,                   # the recurrence
    (2,) + DEC: 2, (3,) + DEC: 2,                 # all-spin decoupling
    (1,) + T1: 1, (1,) + DEC: 1,                  # sheet 1 controls
}
require(len(EXPECT) == 6, f"job inventory {len(EXPECT)} != 6")

if __name__ == "__main__":
    jobs = list(EXPECT)
    with Pool(min(len(jobs), os.cpu_count() or 1)) as pool:
        dims = pool.map(kdim, jobs)
    failures = []
    for job, d in zip(jobs, dims):
        ok = d == EXPECT[job]
        tag = ("RECURRENCE" if job[1:] == T1 and job[0] > 1
               else "decoupling" if job[1:] == DEC and job[0] > 1 else "control")
        print(f"{'PASS' if ok else 'FAIL'}  Solution {job[0]} at ({job[1]},{job[2]}) "
              f"[{tag}]: spin-15 kernel dim {d} (expected {EXPECT[job]})", flush=True)
        if not ok:
            failures.append(job)
    print()
    if failures:
        print(f"SPIN-15 STRATIFICATION FAILED TO MATCH PINS ({len(failures)})")
        sys.exit(1)
    print("SPIN-15 STRATIFICATION CERTIFIED (exact spin-15 kernel dimensions at "
          "the two lifted rational points: (-25/2,-45/2) dim 2 on both merged "
          "sheets, dim 1 on sheet 1 -- the recurrence; (-2,-9) dim 2 on sheets 2,3, "
          "dim 1 on sheet 1 -- the decoupling.  These are the exact anchors of the "
          "weight-16 residue-scan stratification; the complete GF(10007) scans "
          "detecting no other reduction are evidence for, not a proof of, the "
          "absence of further characteristic-zero points)")
