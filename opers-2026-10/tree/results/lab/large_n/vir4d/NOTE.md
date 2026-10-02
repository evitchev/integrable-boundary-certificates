# VIR4d (Fable seat, 2026-10-01) -- blind spin-11 and spin-13 hold-out for the family C_k: all eight cells pass

Protocol seal `SEAL_VIR4d.md` ab1bd1d9... (19:56:10Z).  Prediction `PREDICTION_VIR4d.json`
b8446a9f... (20:04:33Z), registered and snapshotted by the lead BEFORE any spin >= 11 Solution-3
table was opened by me.  Tables opened only after the lead's "REGISTERED".

## Result

| k | t | order | spin 11: mismatches / monomials | spin 13: mismatches / monomials |
|---|---|---|---|---|
| 3 | 11/5 | 7 | 0 / 28 | 0 / 36 |
| 4 | 7/3 | 8 | 0 / 28 | 0 / 36 |
| 6 | 5/2 | 10 | 0 / 28 | 0 / 36 |
| 8 | 13/5 | 12 | 0 / 28 | 0 / 36 |

256 exact rational coefficients, nothing fitted.  Each polynomial is symmetric under PX <-> pi,
so the independent numbers are 16 (spin 11) and 20 (spin 13) per cell, including the leading 1:
15 and 19 after normalisation, 136 in all.  Example (k = 4, spin 11, constant term):
6029497493/21112000000 on both sides.

Controls: tamper M = 1/4 + 1/10 at k = 4: 21 of 28 and 28 of 36 coefficients differ (12 of 16 and
16 of 20 independent).  Labelled NON-blind cell k = 2 (t = 2): 0 of 28 and 0 of 36, agreeing with
cc's registered result.  Screen: every table coefficient finite, PX^(s+1) coefficient nonzero,
at all fibres.  `d3_compare.py` exit 0.

The lead's independent comparison of the registered snapshot: 0 mismatches in all eight cells;
negative controls at t_k + 1/100 give 26 of 28.

## Errors and deviations

- The seal's count was wrong: it said "27 further coefficients with i >= j" at degree 12 and 35
  at degree 14.  The totals are 28 and 36 monomials in all; with i >= j there are 16 and 20.
  The pass criterion (every coefficient agrees) is unaffected.
- The four fibres were computed in parallel processes by the same script and merged
  (`d2_merge.py`); the merged file is the registered one.
- k = 3 at spin 13 was my addition to the lead's list.
- The table format ("a,b" -> coefficient of (P_art^2)^a (Q_art^2)^b) was read only after
  registration; the non-blind k = 2 cell confirms the conversion.

## Files

`SEAL_VIR4d.md`, `PREDICTION_VIR4d.json` (+ .time), `prediction_partial_{3,4,6,8,2}.json`,
`prediction_tamper_k4.json`, `d1_predict.py`, `d2_merge.py`, `d3_compare.py`,
`d4_tamper_predict.py`, logs, `wkb_lib.py`, `SHA256SUMS`.
