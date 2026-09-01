# Spin-17 (weight-18) double prediction: preregistration

Committed BEFORE any weight-18 cylindrical kernel is computed
(2026-08-27, late evening).  Follows the spin-13 sporadic point
(cyl_spin13_sol3_sporadic) and the spin-15 null confirmation
(cyl_spin15_ladder_probe, preregistered 43a2f65).

## Hypothesis H1 (two arithmetic ladders), unchanged from 43a2f65

    sheet 3 jumps at t_k exactly at spin 4k+1        (5, 9, 13, 17, ...)
    sheet 2 jumps at t_k, odd k, exactly at spin 3k+2 (5, 11, 17, ...)

Certified fits: k=1 spin 5 (sheets merge at t1), k=2 sheet-3 spin 9
and no sheet-2 point (3k+2 even: the spin-7 gap), k=3 sheet-2 spin 11
and sheet-3 spin 13, nothing at spin 15 (eight probes).  All jumps
single-spin, +1 over the generic multiplicity 1; sheet 1 never.

## The double positive prediction at spin 17

    P1: sheet 3 at t4 = 11/9, (N,s) = (-238/5,-297/5):   dim 2
    P2: sheet 2 at t5 = 13/11, (N,s) = (-119/2,-143/2):  dim 2

This is the first time H1 predicts something that has never been
seen: a jump at a curve point (t5) where no jump has ever been
certified at any spin, and the return of t4 -- clean through spin
15 on every sheet -- on sheet 3 exactly.  Outcomes:

- both P1 and P2 dim 2 (everything else at 1): H1 confirmed by a
  double positive prediction -- the jumping-loci paper's centerpiece.
- exactly one: H1 refuted as stated; the surviving ladder is kept and
  the other sheet's rule reopened.
- neither: H1 refuted; the locus {3} u {t_k} through spin 15 stands
  but its spin assignment is not arithmetic.
- any dim 2 at a control below: a new sporadic point (certified
  before being called a result) AND a refutation of "sheet s jumps
  only per its ladder".

## The eleven preregistered jobs (weight-18 kernels, 5106 densities)

    baselines: (1, t=2) (2, t=2) (3, t=2)              expect 1,1,1
    P1:        (3, -238/5,-297/5)                      expect 2
    P2:        (2, -119/2,-143/2)                      expect 2
    controls:  (2, -238/5,-297/5)  t4 sheet 2          expect 1 (k=4 even)
               (3, -119/2,-143/2)  t5 sheet 3          expect 1 (4*5+1 = 21)
               (2, -143/4,-189/4)  t3 sheet 2          expect 1 (spin-11-only)
               (3, -143/4,-189/4)  t3 sheet 3          expect 1 (spin-13-only)
               (1, -238/5,-297/5)  t4 sheet 1          expect 1 (sheet 1 never)
               (1, -119/2,-143/2)  t5 sheet 1          expect 1 (sheet 1 never)

Generic spin-17 multiplicity is not established symbolically (that
is a future weight-18 stratification's job); the t=2 baselines
calibrate it pointwise as before.

## Mechanics (first cluster campaign)

Tool: code/cyl_spin17_kdim.py (faithful sibling of the weight-14/16
tools, diff-checked), committed with this note.  Cost: measured
scaling 929 -> 2211 densities took 13 min -> 2.3-2.9 h (exponent
~2.8), so ~27 h per kernel at 5106; single-threaded, ~20 GB RAM.
Jobs are dispatched to the worker nodes by code/cluster_dispatch.sh
(hostnames live in an untracked hosts file, never in the repo):
workers compute raw dims from a git-HEAD snapshot of code/ and a
byte-identical copy of the python312 stack (sympy 1.14.0); the
artifact of record and the registry pin are produced on the main
machine.  Prerequisite, executed before dispatch: every worker must
reproduce the certified spin-13 sporadic point (Solution 3 at
(-143/4,-189/4) -> dim 2) from the same snapshot -- cross-node
determinism check, results recorded in the program note.

## Addendum (2026-08-28, 17:30, BEFORE any spin-17 result is in): outcome wording

"Both P1 and P2 dim 2 ... H1 confirmed" above overstates: two positive
predictions would strongly support H1's predicted slots at spin 17; they
cannot confirm H1's universal clauses ("exactly", "single-spin", "sheet 1
never"), which remain conjectures tested one rung at a time.  The
sheet-2 weight-16 candidate at t1 = 5/3 (pending exact kernels) already
puts the single-spin clause in question.  Predictions and job list are
unchanged.
