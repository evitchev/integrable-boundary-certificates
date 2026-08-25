# Preregistration: the spin-13 test of the t_k = (2k+3)/(2k+1) hypothesis

Written and committed BEFORE any weight-14 (spin-13) kernel of the
cylindrical Solutions was computed.  The commit hash of this file is
the timestamp; the computation (code/cyl_spin13_prediction.py) is
launched after the commit.

## Hypothesis (LAB HYPOTHESIS from results/cyl_strata_lab.md)

The sporadic single-spin multiplicity jumps of the cylindrical families
(hep-th/0404195 Sec. 3.3, curve s^2 = (N-25)(N-1), t-parametrization
N = (t^2-25)/(t^2-1), s = -24t/(t^2-1)) sit at

    t_k = (2k+3)/(2k+1):  5/3 (spin 5), 7/5 (spin 9), 9/7 (spin 11), ...

## Prediction under test

    t = 11/9,  (N, s) = (-238/5, -297/5)    [s^2 = (N-25)(N-1) verified]

- Spin: 13 (weight-14 ad_{I_3} kernel, class coordinates mod total
  derivatives, exact rationals, open Python stack -- NOT Mathematica).
- Sheet: Solution 3 is the alternation GUESS (Sol 2+3 merger at 5/3,
  Sol 3 at 7/5, Sol 2 at 9/7).  A jump on ANY Solution counts.

## Acceptance rules (fixed now)

- CONFIRMED: the spin-13 kernel dimension at (-238/5,-297/5) exceeds
  the generic spin-13 multiplicity for at least one Solution.
- NOT CONFIRMED (not refuted): no Solution jumps there at spin 13.
  The known sequence has a gap at spin 7, so absence at spin 13 does
  not exclude t = 11/9 appearing at a higher spin; a refutation needs
  a structural argument, not a search.
- A jump anywhere else at spin 13 is outside this test; it belongs to
  the global weight-14 stratification (Mathematica, next step either
  way).

## Controls (all three Solutions unless stated; expectations fixed now)

- Generic multiplicity at spin 13: t = 2, (N,s) = (-7,-16), not a known
  jump point.  Expected dim 1 (generic multiplicity was 1 at every
  spin 5-11).  This value defines "generic" for the acceptance rule.
- Negative controls on the same branch, off the sequence:
  t = 5/4, (N,s) = (-125/3,-160/3);  t = 7/6, (N,s) = (-851/13,-1008/13).
  Expected dim = generic (1) for every Solution.
- Positive controls (known points, continuation to spin 13):
  (-2,-9): Solutions 2,3 expected dim 2 (decoupling persists, dims
  were 2,2,2,2 at spins 5-11); Solution 1 expected 1.
  (-143/4,-189/4) Solution 2: expected 1 (jump confined to spin 11).
  (-24,-35) Solution 3: expected 1 (jump confined to spin 9).
  (-25/2,-45/2) Solutions 2,3: expected 1 (jump confined to spin 5).

## What happens next

- Hit: the certificate pins the dimension, the mechanism analysis
  (screened directions, X/Y projection ranks) follows, and the global
  weight-14 stratification tests uniqueness and completeness at spin 13.
- Miss: the hypothesis stays open; the global weight-14 run decides
  whether spin 13 has any new sporadic point at all and where.
