# VIR5a (Fable seat, 2026-10-01) -- the generic-t operator for Solution 3: all sealed cells pass

Seal `SEAL_VIR5a.md` 2efea474... (20:15:41Z).  Prediction `PREDICTION_VIR5a.json` 0c12a7e9...
(20:17:42Z), registered and snapshotted by the lead before any comparison.  The spin-11 and
spin-13 tables were already open to me (VIR4d), so this run is "hash-registered, not blind".

## Verdict

**For every fibre tested -- in all four regions of t -- the certified Solution-3 vacuum
eigenvalues of spins 1, 3, ..., 13 equal, coefficient by coefficient, the WKB coefficients of
ONE operator with no free parameter:**

    H(theta) psi = x^n ( x^{n/k} - E ) psi,          theta = x d/dx,   n = k + 4,   s = (n-1)/2,
    H(theta) = [ (s-theta)^2 - l_0^2 ] [ (s-theta)^2 - l_1^2 ] (n/k)^k Gamma(u + (k+1)/2) / Gamma(u - (k-1)/2),
    u = k (s - theta)/n,    (l_0, l_1) = n sqrt((k+1)/k) (PX, pi),    k = p^2/(1 - p^2),  t = (3k+2)/(k+2).

In Mellin space it is the three-term difference equation
    H(nu) psihat(nu) = psihat(nu - n - n/k) - E psihat(nu - n).

| cell | k | t | n | mismatches / coefficients (spins 1-13) | sealed confidence |
|---|---|---|---|---|---|
| A | 10/3 | 9/4 | 22/3 | 0 / 119 | 85% |
| A | 5/3 | 21/11 | 17/3 | 0 / 119 | 85% |
| B | 2/3 | 3/2 | 14/3 | 0 / 119 | 65% |
| B | 1/2 | 7/5 | 9/2 | 0 / 98; spin 9 not produced (9 = 2n), as sealed | 65% |
| C | -6 | 4 | -2 | 0 / 119 (regression: cc's item 297 already had this fibre) | 50% |
| C | -18/5 | 11/2 | 2/5 | 0 / 119 | 50% |
| C | -6/5 | -2 | 14/5 | 0 / 119 | 50% |
| D | -2/5 | 1/2 | 18/5 | 0 / 119 | 40% |
| control | 4 | 7/3 | 8 | 0 / 119; spins 11, 13 string-identical to the registered VIR4d prediction | |

Per spin the counts are 3, 6, 10, 15, 21, 28, 36 monomials; 0 mismatches in every (fibre, spin).
Tampers at t = 9/4: Gamma ratio replaced by u^k: 84 of 119 differ (already the spin-1 constant,
-19/72 against -1/6); M = 1/k + 1/10: 84 of 119 differ.  `e3_compare.py` exit 0.

## What was done

- `mellin_wkb.py`: WKB for an Euler-type operator with a general symbol, through the exact
  conjugation identity e^{-S} H(theta) e^{S} . 1 = sum_m Ht^(m)(T) C_m and the large-T expansion
  of the symbol; the Gamma ratio enters through c_j(k) = k(k-1)...(k-2j+1) [t^(2j)]
  (t/(2 sinh(t/2)))^(k+1).  Exact arithmetic; k is a rational number, nothing is interpolated.
- Validation before sealing (`e0_validate.py`, exit 0): at integer k = 2, 3, 4 the engine
  reproduces the verified chain engine at spins 1-9, including the missing spin 7 at k = 3.
- Nine fibres computed to order 14 in parallel, merged (`e2_merge.py`), registered, compared.

## Two points the lead asked about

- k = 1/2 (t = 7/5), spin 9.  The operator gives no charge of spin 9 (the WKB coefficient is
  identically zero), as sealed: 9 = 2n.  The record has a two-dimensional spin-9 kernel at
  t = 7/5 (cyl_t2_primary).  The same pattern as k = 1, 3, 5 (spin n) and as the record's
  grade-0 rank-drop locus a_X = 2i/(2k'-1): with a_X = 2/n that locus is "spin = i n".
  So the rule is: the operator has no charge at spins that are integer multiples of n, and
  those are exactly the record's degenerate (fibre, spin) pairs I know of.  What supplies the
  charge there is open (VIR5 (b)).
- k = -6, odd-order flags.  The "False" flags in the registered file were an artefact of my
  Beta-ratio formula (a Pochhammer in (i-1)/n vanishes for odd i at n = -2).  Direct test
  (`e4_oddflag.py`, exit 0): the odd-order WKB terms are exact derivatives at k = -6 (and at
  10/3, 1/2), so there is no even-spin charge.

## Errors and caveats

- The engine's first version lost every derivative term (a truncation applied in the wrong
  order); found by the pre-seal validation; log kept.
- No a-priori degree bound in k (said in the seal).  The statement tested is equality at
  eight fibres, 931 coefficients, plus the integer fibres of items 294-295.
- Vacuum eigenvalues only; spins <= 13; Solution 3 only.  The WKB series is formal: for
  n < 0 (t > 3) I have not examined which asymptotic regime of E it describes, only that the
  coefficients agree.  No spectral problem has been solved numerically.
- The remark that (t/(2 sinh(t/2)))^(k+1) may be the record's "x/sinh x law" is unchecked
  (follow-up agreed with the lead).

## Files

`SEAL_VIR5a.md`, `PREDICTION_VIR5a.json` (+ .time), `prediction_k*.json`, `mellin_wkb.py`,
`e0_validate.py`, `e0_debug.py`, `e1_predict.py`, `e2_merge.py`, `e3_compare.py`,
`e4_oddflag.py`, logs, `wkb_lib.py`, `cft_F.pkl`, `SHA256SUMS`.
