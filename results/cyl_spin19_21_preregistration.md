# Spins 19 and 21 (weights 20 and 22): preregistration

Committed BEFORE any weight-20 or weight-22 kernel, matrix scan, or
fast-path run (2026-08-30).  Follows the spin-17 double prediction
(results/cyl_spin17_preregistration.md, confirmed by two independent
exact methods and by residue scans under two primes on all sheets).

## Hypothesis H1 (two arithmetic ladders), with the spin-15 amendment

    sheet 3 jumps at t_k at spin 4k+1        (5, 9, 13, 17, 21, ...)
    sheet 2 jumps at t_k, odd k, at spin 3k+2 (5, 11, 17, 23, ...)

Amendment (certified 2026-08-28): jumps are NOT all single-spin -- t_1
(where the sheets merge) jumps at spins 5 AND 15.  Whether other points
recur, and by what rule, is OPEN.  [Corrected 2026-08-30, Codex round
11: the original sentence said no recurrence had been seen for any
point other than t_1; that is false -- t_3 = 9/7 recurs ACROSS sheets
(spin 11 on sheet 2, spin 13 on sheet 3).  Correct statement: no point
other than t_1 has jumped twice on the SAME sheet at spins 5..17.]

## Predictions

Spin 19 (weight 20): NOTHING on any sheet beyond the all-spin
decoupling at t = 3 -- neither ladder has a slot at 19 (4k+1 = 19 and
3k+2 = 19 have no integer solution).  A null test.

Spin 21 (weight 22): t_5 = 13/11, (N,s) = (-119/2,-143/2), jumps on
SHEET 3 (4*5+1 = 21); sheet 2 has no slot (3k+2 = 21 gives k = 19/3);
sheet 1 nothing.  Note t_5 already jumped on sheet 2 at spin 17: this
would be the SECOND nonmerged point jumping on both sheets at different
spins, after t_3 (11 on sheet 2, 13 on sheet 3) -- that cross-sheet
pattern repeated one rung up.  [Corrected 2026-08-30: the original
said "first", contradicting the t_3 record cited in the same sentence.]

Open watch-list at both spins: recurrences of t_1 = 5/3 (spins 5, 15;
next at 25 if the spacing is 10, at 45 if the pattern is 3x), t_2 =
7/5, t_3 = 9/7, t_4 = 11/9; any drop there is a candidate, second
prime, then exact kernel by both methods.

## Method

Weight-20 matrices built on the worker nodes (wolfram/lab/
cyl_strata_build.py 20, one sheet per node); residue scans over two
primes with lab/modp_scan.py (GPU kernels if available by then); exact
dimensions at every drop and at the predicted points by
lab/fast_kernel2.py (witness records) and by the pure-Python engine
(distributed attested mode when built); certificates registered before
their runs with these pins.  Scope: single-prime scans are evidence,
not completeness proofs; predictions concern the named slots and do
not establish H1's universal clauses.

## Recurrence clause (added 2026-08-30, before any weight-26 data)

Certified so far: t_1 = 5/3 (merged sheets 2 = 3) jumps at spins 5 AND
15; no other point has jumped twice on the same sheet through spin 17.
Three competing hypotheses, stated now so that the prediction predates
the data whichever weight first reaches it:

  R1 (period 10):     t_1 recurs at spin 25 (weight 26), then 35, ...
  R2 (geometric x3):  t_1 recurs at spin 45 (weight 46) -- out of reach.
  R3 (no rule):       t_1 does not recur at 25.

Universality sub-clause of R1 (R1u): if the period-10 rule is a
property of every jump and not of the merged point alone, then t_2
recurs on sheet 3 at spin 19, t_3 on sheet 2 at spin 21 and on sheet 3
at spin 23, t_4 on sheet 3 at spin 27, t_5 on sheet 2 at spin 27.
H1 predicts NOTHING at spin 19, so the spin-19 null test already
decides R1u: a t_2 drop at spin 19 on sheet 3 (weight-20 scans, both
primes, then exact kernel by both methods) would refute H1's silence
and support R1u; its absence refutes R1u and leaves R1/R2/R3 for t_1
open until weight 26.  Decisive weights: 20 (R1u), 26 (R1 vs R3),
46 (R2).  Nothing here changes the spin-19/21 predictions above.
