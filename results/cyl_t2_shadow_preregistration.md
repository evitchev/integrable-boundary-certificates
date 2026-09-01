# The t_2 shadow test: preregistration (2026-08-31, before the run)

Committed BEFORE the computation.  Template: the certified
cyl_jump_point_tower (t_1 up the tower).  Point: t_2 = 7/5, (N, s) =
(-24, -35) = c_{1,6}, Solution 3 (the sheet that jumps at spin 9).
The curve's own screening relation c(a) - 1 = -3(a^2-2)^2/a^2 gives,
at c = -24, exactly a^2 in {1/3, 12} (a = sqrt3/3, 2 sqrt3): the
c_{1,6} pair, product a_+ a_- = 2 as at the other identified points.

Mechanism reading under test (mechanism note secs. 4b, 4c): the
generic tower at t_2 lies on a screened c_{1,6} line, and the extra
spin-9 class -- the discriminant's (1-6x) factor -- is the SHADOW: the
unique direction NOT screened by the c_{1,6} pair.

## Pins (H_S)

Solution 3 at (-24, -35), weights w = 6, 8, 10, 12 (spins 5, 7, 9, 11):
1. ad_{I_3} kernel dimensions [1, 1, 2, 1] (the certified record: the
   jump is confined to spin 9 in this range).
2. Screened-subspace dimensions under a^2 = 1/3: [1, 1, 1, 1]; under
   a^2 = 12: [1, 1, 1, 1] -- the screened line runs through the tower,
   and at spin 9 it is a PROPER line inside the 2-dim kernel (the
   shadow statement).
3. Negative controls: screened dims under a^2 = 2 and a^2 = 3 are
   [0, 0, 0, 0] (those are the c = 1 and c = -2-side values, wrong
   point).
4. At w = 10 the two screened lines (a^2 = 1/3 and 12) COINCIDE (same
   1-dim subspace: the dual pair screens the same tower member), and
   the extra class is genuinely mixed (no pure-X or pure-Y combination
   in the kernel; projection ranks 2, 2 as in the t_1 J3 check).

Fail-closed; any deviation is a FINDING (in particular, a screened dim
of 2 at w = 10 would refute the shadow reading -- the extra class would
be INSIDE the screened structure, pointing at a bigger algebra, not a
quotient shadow).  Expected runtime: minutes (weight <= 12, one point,
one sheet; Q(sqrt3) exact).

## Refutation addendum (2026-08-31, after the development smoke, BEFORE the registered run)

The smoke run REFUTED pins H_S 2-4 as preregistered: at t_2 the
screened dimensions are 0 at every weight under BOTH Kac-pair values
a^2 = 1/3 and 12 (and the joint dim is 0); only H_S 1 (kernel dims
[1, 1, 2, 1]) and the mixedness check passed.  Worse: the refutation
was FORESEEABLE -- results/cyl_strata_lab.md (certified) already
states that the spin-9 jump at (-24, -35) has "no screening
structure": no e^{sqrt2 X}- or sqrt3-screened direction, and I_3 is
not X-screened at this point at all.  This preregistration failed to
check the certified mechanism classification of the point it targeted;
the a^2 in {1/3, 12} derivation transplanted the t_1 embedding
mechanism to a point certified to be of the third (non-screening)
kind.  Process rule adopted: a preregistration about a specific point
MUST cite and check the certified record for that point first.

The registered certificate therefore pins the OBSERVED result: the
no-screening-structure statement EXTENDS to the c_{1,6} Kac pair
(a^2 = 1/3, 12 screen nothing at weights 6-12, individually and
jointly), alongside the reconfirmed kernel dims and mixedness.
Consequence for the mechanism note: the Phi_{1,3}-shadow reading has
NO exponential-screening operationalization at t_2; its proper test is
the singular-vector / quotient construction (mode-action extraction
beyond the residue-only engine -- the "small extension" of
[private-note]), which is real work, not a
cheap check.
