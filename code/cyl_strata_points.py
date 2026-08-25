"""Open-stack certification of the multiplicity stratification found by
the Mathematica lab (wolfram/lab/cyl_strata.wls, spins 5, 7, 9).

The lab computes, exactly and for the whole curve s^2 = (N-25)(N-1),
the points where the charge-class kernel of ad_{I_3} (the multiplicity)
jumps above its generic value 1, for each cylindrical family.  That
completeness statement is a LAB result (typed as such in its JSON).
This certificate pins, in the open stack, the POSITIVE content: the
kernel dimensions at every lab-reported jump point, the control that
the other solutions do NOT jump there, the N = 1 endpoint (unreachable
by the t-parametrization), and the mechanism at the new point.

Lab stratification, spins 5 / 7 / 9 (generic multiplicity 1):
  Solution 1: no jump anywhere (finite t != +-1) at spins 5, 7, 9.
  Solution 2: (-2,-9) at 5, 7, 9;  (-25/2,-45/2) at 5 only.
  Solution 3: (-2,-9) at 5, 7, 9;  (-25/2,-45/2) at 5 only;
              (-24,-35) at 9 only  [NEW, found by the lab].
  N = 1 endpoint: multiplicity 1 for all three solutions at 5, 7, 9.
Mechanisms: (-2,-9) decoupling (kernel = X (+) Y at every odd spin);
(-25/2,-45/2) self-dual embedding with one extra mixed class at spin 5
(cyl_jump_point_tower.py); (-24,-35) a pure multiplicity jump -- no
e^{sqrt2 X}-screened direction, no pure-X/pure-Y element, and I_3 not
X-screened there at all (asserted below).

Fail-closed.  Finite-range (spins <= 9 here; spin 11 in a follow-up).
"""
import sys
from fractions import Fraction as F
from pathlib import Path

import sympy as sp

from mixed_engine import cyl_P4, gen_mixed_basis, genuine_kernel_mixed

_jp = Path(__file__).resolve().parent / "cyl_jump_point_tower.py"
_src = _jp.read_text()
_ns = {"__file__": str(_jp)}
exec(_src[:_src.index("N, s = F(-25, 2), F(-45, 2)")], _ns)
screened_dim = _ns["screened_dim"]
mixed_projection_ranks = _ns["mixed_projection_ranks"]
_sc = Path(__file__).resolve().parent / "cyl_hidden_kdv_scan.py"
_src2 = _sc.read_text()
_ns2 = {}
exec(_src2[:_src2.index("k1 = kdv_mixed_class(1)")], _ns2)
exec(_src2[_src2.index("def apply_d(cur):"):_src2.index("a = sp.Symbol('a')\n# controls")], _ns2)
x_screening_conditions, nonzero_roots = _ns2["x_screening_conditions"], _ns2["nonzero_roots"]

failures = []


def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok:
        failures.append(name)


def dims(sol, N, s, weights=(6, 8, 10)):
    P = cyl_P4(sol, N, s)
    return [len(genuine_kernel_mixed(P, gen_mixed_basis(w), N)) for w in weights]


# (point, expected dims at spins 5,7,9 per solution)
TABLE = {
    (F(-2), F(-9)):          {1: [1, 1, 1], 2: [2, 2, 2], 3: [2, 2, 2]},
    (F(-25, 2), F(-45, 2)):  {1: [1, 1, 1], 2: [2, 1, 1], 3: [2, 1, 1]},
    (F(-24), F(-35)):        {1: [1, 1, 1], 2: [1, 1, 1], 3: [1, 1, 2]},
    # N=1 endpoint: FORMAL invariant algebra (no finite-N Gram relations
    # imposed); the physical two-boson component kernels were independently
    # confirmed to be 1,1,1 at spins 5,7,9 (Codex audit 2026-08-23).
    (F(1), F(0)):            {1: [1, 1, 1], 2: [1, 1, 1], 3: [1, 1, 1]},
    (F(28), F(9)):           {1: [1, 1, 1], 2: [1, 1, 1], 3: [1, 1, 1]},   # embedding point: no jump
}
for (N, s), row in TABLE.items():
    assert s * s == (N - 25) * (N - 1)
    for sol, want in row.items():
        got = dims(sol, N, s)
        check(f"S1: Solution {sol} at ({N},{s}): kernel dims at spins 5,7,9 = "
              f"{got} (lab: {want})", got == want)

# mechanism at the new point (-24,-35), Solution 3
N, s = F(-24), F(-35)
P3 = cyl_P4(3, N, s)
ker10 = genuine_kernel_mixed(P3, gen_mixed_basis(10), N)
check("S2: (-24,-35) Solution 3: spin-9 kernel dim 2 with NO e^{sqrt2 X}-screened "
      "direction (dim 0) and none at sqrt3",
      len(ker10) == 2 and screened_dim(ker10, 10, 2) == 0
      and screened_dim(ker10, 10, 3) == 0)
rx, ry = mixed_projection_ranks(ker10, 10)
check(f"S2: ... and X-/Y-projection ranks ({rx},{ry}) of 2: no pure-X or pure-Y "
      "COMBINATION (not a decoupling; basis-independent)", (rx, ry) == (2, 2))
a = sp.Symbol('a')
check("S2: ... and I_3 itself has no X-screening at this point (not on the "
      "screening locus): a pure multiplicity jump",
      nonzero_roots(x_screening_conditions(P3, a, 4), a) == [])

print()
if failures:
    print(f"STRATIFICATION-POINT CHECKS FAILED ({len(failures)}):")
    for f in failures:
        print(" -", f)
    sys.exit(1)
print("STRATIFICATION POINTS CERTIFIED (lab jump points of spins 5-9 confirmed in "
      "the open stack: (-2,-9), (-25/2,-45/2), (-24,-35) [new, Solution 3, spin 9, "
      "pure multiplicity jump]; N=1 endpoint and N=28 multiplicity one)")
