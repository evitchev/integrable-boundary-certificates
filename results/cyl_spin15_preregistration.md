# Spin-15 (weight-16) targeted stratification probe: preregistration

Committed BEFORE any weight-16 cylindrical kernel is computed
(2026-08-27, same day the weight-14 global stratification found the
sporadic spin-13 point).

## Certified record this builds on

- t1=5/3, (-25/2,-45/2): Solutions 2 and 3 MERGE there (same I_3);
  kernel dims 2,1,1,1 at spins 5,7,9,11 -- jump at spin 5 only,
  attributable to both sheets (c_{1,p}-tower certificate); dim 1 both
  sheets at spin 13 (cyl_spin13_prediction).
- t2=7/5, (-24,-35): sheet-3 jump at spin 9 (spin-9 stratification);
  the certified spin-11 stratification lists no jump there; dim 1 at
  spin 13 (cyl_spin13_prediction).  No sheet-2 jump is known at any
  spin (in particular the global weight-8 scan found nothing at
  spin 7 anywhere).
- t3=9/7, (-143/4,-189/4): sheet-2 jump at spin 11 (spin-11
  stratification), back to dim 1 at spin 13 (cyl_spin13_prediction);
  sheet-3 jump at spin 13 (TODAY: weight-14 global stratification +
  cyl_spin13_sol3_sporadic, dims 1/1/2 across sheets); sheet 1 dim 1.
- t4=11/9, (-238/5,-297/5): dim 1 on all sheets at spin 13
  (cyl_spin13_prediction -- the naive-ladder refutation).
- Solution 1: no jump anywhere through spin 13; generic multiplicity
  1 at every odd spin <= 13, symbolically at weight 14.

## Hypothesis H1: two arithmetic ladders

    sheet 3 jumps at t_k exactly at spin 4k+1        (5, 9, 13, ...)
    sheet 2 jumps at t_k, odd k, exactly at spin 3k+2 (5, 11, 17, ...)

Fits every certified point: k=1 gives spin 5 on both formulas (the
sheets literally merge at t1); k=2 gives sheet-3 spin 9 and NO
sheet-2 point (3k+2=8 is even -- the spin-7 gap is explained, not
patched); k=3 gives sheet-2 spin 11 and sheet-3 spin 13 (today's
finding is the ladder's continuation, not an anomaly).  All jumps
single-spin (+1 over the generic multiplicity 1).

H1 PREDICTS AT SPIN 15 (weight 16): NOTHING.  All eight probe
kernels below have dimension 1.  The next activations sit at spin 17
(weight 18): t4=11/9 on sheet 3 (4*4+1) AND t5=13/11 on sheet 2
(3*5+2) -- a double prediction at the next rung.

Alternative readings: an "echo" pairing (sheet-2 spin s -> sheet-3
spin s+2, from t3's 11->13) already fails at t2 (no sheet-2 spin-7
partner exists); kept only as a foil.  Agnostic reading: any dim-2
below is a new sporadic point regardless of hypotheses.

## The eight preregistered jobs (weight-16 kernels, 2211 densities)

    baselines:  (1, t=2)  (2, t=2)  (3, t=2)      expect 1,1,1
    t3=9/7:     (2, -143/4,-189/4)                 expect 1 (spin-11-only)
                (3, -143/4,-189/4)                 expect 1 (single-spin at 13)
    t4=11/9:    (2, -238/5,-297/5)                 expect 1 (H1: activates never at 15)
                (3, -238/5,-297/5)                 expect 1 (H1: activates at 17)
    t2=7/5:     (3, -24,-35)                       expect 1 (spin-9-only)

Any expectation failure is (a) a refutation of H1 at that slot and
(b) a finding: certified fail-closed (open stack) before being called
a result, per house protocol.  Generic spin-15 multiplicity is NOT
yet established symbolically (that is the weight-16 global
stratification's job); the t=2 baselines calibrate it pointwise
exactly as in the spin-13 preregistration.

Tool: code/cyl_spin15_kdim.py (weight-16 sibling of
cyl_spin13_kdim.py), committed with this note before any run.
