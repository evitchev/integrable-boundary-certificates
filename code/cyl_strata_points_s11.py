"""Open-stack certification of the lab's spin-11 stratification of the
cylindrical curve (wolfram/lab/cyl_strata.wls at weight 12; see
cyl_strata_points.py for spins 5-9 and results/cyl_strata_lab.md).

Lab result at spin 11 (generic multiplicity 1):
  Solution 1: no jump anywhere.
  Solution 2: (-2,-9) and (-143/4,-189/4)  [NEW, t = 9/7].
  Solution 3: (-2,-9) only.
Pinned here (all weight-12 kernels of ad_{I_3}, ~1 min each, hence the
expensive tier):
 1. Solution 2 at (-143/4,-189/4): kernel dims [1,1,1,2] at spins
    5,7,9,11 -- a spin-11-only jump; Solutions 1 and 3 have dim 1 there.
 2. Mechanism: the enlarged kernel has no e^{sqrt2 X}-screened direction
    (nor sqrt3) and no pure-X / pure-Y element -- a pure multiplicity
    jump, like (-24,-35) for Solution 3 at spin 9.
 3. (-24,-35), Solution 3: dim 1 at spin 11 (its jump is confined to
    spin 9).
 4. (-2,-9): dims 2 for Solutions 2 and 3 at spin 11 (the decoupling
    persists), 1 for Solution 1.
Finite-range (spins <= 11).  Fail-closed.
"""
import sys
from fractions import Fraction as F
from pathlib import Path

from mixed_engine import cyl_P4, gen_mixed_basis, genuine_kernel_mixed

_jp = Path(__file__).resolve().parent / "cyl_jump_point_tower.py"
_src = _jp.read_text()
_ns = {"__file__": str(_jp)}
exec(_src[:_src.index("N, s = F(-25, 2), F(-45, 2)")], _ns)
screened_dim = _ns["screened_dim"]
mixed_projection_ranks = _ns["mixed_projection_ranks"]

failures = []


def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok:
        failures.append(name)


def kdim(sol, N, s, w):
    return genuine_kernel_mixed(cyl_P4(sol, N, s), gen_mixed_basis(w), N)


N, s = F(-143, 4), F(-189, 4)
assert s * s == (N - 25) * (N - 1)
dims2 = [len(kdim(2, N, s, w)) for w in (6, 8, 10, 12)]
check(f"T1: Solution 2 at (-143/4,-189/4): kernel dims at spins 5,7,9,11 = {dims2} "
      "(spin-11-only jump)", dims2 == [1, 1, 1, 2])
for sol in (1, 3):
    d = len(kdim(sol, N, s, 12))
    check(f"T1: Solution {sol} at (-143/4,-189/4): spin-11 kernel dim {d} (no jump)",
          d == 1)
ker12 = kdim(2, N, s, 12)
check("T2: mechanism -- no e^{sqrt2 X}-screened direction (dim 0), none at sqrt3",
      screened_dim(ker12, 12, 2) == 0 and screened_dim(ker12, 12, 3) == 0)
rx, ry = mixed_projection_ranks(ker12, 12)
check(f"T2: ... and X-/Y-projection ranks ({rx},{ry}) of 2: no pure-X or pure-Y "
      "COMBINATION (not a decoupling; basis-independent)", (rx, ry) == (2, 2))
d = len(kdim(3, F(-24), F(-35), 12))
check(f"T3: Solution 3 at (-24,-35): spin-11 kernel dim {d} (jump confined to spin 9)",
      d == 1)
for sol, want in ((1, 1), (2, 2), (3, 2)):
    d = len(kdim(sol, F(-2), F(-9), 12))
    check(f"T4: Solution {sol} at (-2,-9): spin-11 kernel dim {d} (expected {want})",
          d == want)

print()
if failures:
    print(f"SPIN-11 STRATIFICATION CHECKS FAILED ({len(failures)}):")
    for f in failures:
        print(" -", f)
    sys.exit(1)
print("SPIN-11 STRATIFICATION POINTS CERTIFIED ((-143/4,-189/4): Solution 2 "
      "spin-11-only pure multiplicity jump; (-24,-35) and (-2,-9) as stated; "
      "finite-range)")
