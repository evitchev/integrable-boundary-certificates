# Odd-q controls: preregistration (2026-08-30, before the run)

Context: the certified identity (cyl_coulomb_gas_parametrization)
reads every rational t as beta^2 = (t-1)/(t+1); the certified jumping
locus through spin 17 is {beta^2 = 1/q, q even, 2 <= q <= 12}.  The
odd-q points are equally good rational points of the curve.  Their
silence so far rests on the weight-14 symbolic stratification, the
weight-16/18 exhaustive residue scans (two primes), and the t = 2
baselines at spins 13 and 15.  This certificate anchors the "even q
only" statement with EXACT kernel dimensions at the odd-q points, at
exactly the spins where their even-q neighbours jump.

Points: t = 2   (q = 3,  (N,s) = (-7, -16)),
        t = 3/2 (q = 5,  (N,s) = (-91/5, -144/5)),
        t = 4/3 (q = 7,  (N,s) = (-209/7, -288/7)).
Sheets: Solutions 2 and 3.  Spins: 9, 11, 13 (weights 10, 12, 14) --
the spins certified to jump at q = 6, 8 (t_2 = 7/5, t_3 = 9/7).
18 exact kernels (mixed_engine, the spin-13 sporadic machinery).

EXPECTED (H-odd): every one of the 18 kernel dimensions equals 1 (the
generic multiplicity).  Any dimension > 1 is a FINDING: a jump at odd
q, refuting the even-q clause of the (1,q) reading and flagged as
such, never re-pinned silently.  Fail-closed; certificate registered
before its run; expensive with a predeclared runner record.
