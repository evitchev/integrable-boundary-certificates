# NUMTAB -- the spectral-evidence table with stability tests (for the journal version)

cc (Opus 5.5 seat), 2026-10-10. Lead's commission (brief 8a7dc736...), answering ChatGPT's review sec. 4.2 / referee 3.2.
- Seal: `SEAL_NUMTAB.md` 79f073bf52bedd026445c98760503b3730a0d45f9a9b6ded4a358c3058972688 (15:09:26), sent before any refit.
- ARCHIVED data only: inputs/ hashed from results/lab/large_n/{num1,num2,num3}. No ODE integration was recomputed. WKB predictions come from the archived engines / stored predictions.

## Result
- **167 rows** (D1 NUM1 Sol 1 M = 5, A/B/C + ±alpha averages; D2 NUM2 Sol 3 t = 2, A/B/C; D3 NUM2 t = 9/4, A/B; D4 NUM3 Sol 2 t = 2 extended grid, A/B; D5 NUM3 t = 9/4, A/B).
  - 117 rows carry a prediction.
  - **No row is in tension: every predicted row agrees within 10 x its error estimate.**
- **Error estimate epsilon** = the largest spread over:
  - (a) two shifted windows (drop the lowest / highest quarter of the energies; refit at the largest admissible K2);
  - (b) two truncations (K2 vs K2 - 4);
  - (c) 60- vs 120-digit archived data (D1 A, D2 A, D4 A on [1e4, 1e7]; no 120-digit archive elsewhere).
  - It is dominated by the W+ window (dropping the highest energies). Precision spreads are small: ≤ 8e-22 on every determined predicted row (largest: D1 K8). The largest overall is 7e-13 (D4 E^(-15/2), undetermined), then 5e-16 (D4 E^-6 log E, an unconstrained R row). None of them sets the error estimate.
  - The error-estimate digits e are typically 5-12 below the agreement digits d. So the error estimate is CONSERVATIVE, and the agreement digits overstate what a window-robust statement supports. **Report e, not d, as the claim.**
- Point-A summary (d / e):

  | | grade 0 | spin 1 | spin 3 | spin 5 |
  |---|---|---|---|---|
  | D1 | 52.6 / 39.6 | 42.4 / 30.6 | 35.1 / 24.4 | e^-3 resonant post hoc: 28.4 / 18.9 |
  | D2 | 32.1 / 21.8 | 23.7 / 14.4 | E^-3/2 log E, registered: 16.9 / 8.6 | 11.4 / 4.1 |
  | D3 | 30.8 / 25.8 | 21.3 / 17.3 | 13.9 / 11.0 | 8.0 / 6.1 |
  | D4 | 34.1 / 26.5 | 25.3 / 18.7 | 18.1 / 12.5 | 12.4 / 7.6 |

- **UNDETERMINED (e < 0, error estimate exceeds the value), 9 rows:**
  - D2 E^(-9/2) log E (registered -(3/2) R_10; the fit agrees to 3-4 digits but the window spread is larger);
  - D4 E^(-9/2), E^(-15/2);
  - D5 E^(-39/10).
  - Also D4 E^-6 log E: its prediction is 0 and the fit is unconstrained (|fit| = 2e13 with a larger epsilon). It is not a test.
- **UNRESOLVED resonant finite parts** (no prediction):
  - D2 E^(-3/2) (e = 8.1) and E^-3 (e = 1.7);
  - D4 E^-3 (e = 2.9).

## Classes
- W: WKB predicted (-int s_K; odd grades predicted 0).
- R: non-WKB registered before the data.
- S: non-WKB power sealed, coefficient free.
- H: non-WKB post hoc.
- F: fitted free (grade 1).
- U: resonant finite part, unresolved.
- *: undetermined.

## Files for the journal
- `table.tex` (booktabs, one page, point A of D1-D4; compiled check in texcheck/t.pdf, no overfull boxes).
- `table.csv` / `table_final.json` (all 167 rows with the four spreads).
- Caption (5 lines, also in table.tex):
  > Spectral evidence: large-E fits of log D or log|Q| (dominant Frobenius exponent) for the three operators, point A (all points in the supplementary table).
  > Class: W = WKB coefficient (predicted); R = non-WKB term registered before the data; S = non-WKB power sealed, coefficient fitted; H = post hoc; F = fitted free; U = resonant finite part, not identified; * = undetermined.
  > d = agreement digits; for a zero prediction the absolute fit is shown.
  > e = error-estimate digits, from the largest spread over two shifted fit windows, two basis truncations and, where archived, 60- vs 120-digit data.
  > All predicted rows agree within 10 epsilon.
- **Method note (one paragraph):** For each operator and momentum point we computed log D(-e) (second order) or log|Q_j(E)| for the Frobenius exponent dominant at x = 0 (third and fourth order) at 36 energies, 1e3 ≤ E ≤ 1e7 (Solution 2: 39 energies, 1e4 ≤ E ≤ 1e9), with 60-digit arithmetic (120 digits at one point per operator), and fitted the large-E expansion by least squares at three times the working precision. The basis contains the WKB lattice E^((k+1)(1-K)/n), log E, a constant, and the non-WKB powers E^(-m(k+1)/k) (with log E at resonant powers), with no WKB input to the fit. Each fitted coefficient is compared with its prediction where one exists: the absolute WKB coefficient -int s_K, or a registered non-WKB value. The error estimate is the largest change of the coefficient under (a) removing the lowest or highest quarter of the energies, (b) reducing the basis by four powers, and (c) replacing 60- by 120-digit data. We report both the agreement digits d and the error-estimate digits e. Terms introduced after the initial fits (H), terms whose power but not coefficient was fixed in advance (S), and resonant finite parts we cannot yet identify (U) are marked separately.

## Disclosures (post-processing rules, stated here; not in the seal)
- A WKB prediction with |pred| < 1e-40 is a round-off ZERO (odd grades at alpha = 0 and in the alpha-even average) and is treated as pred = 0.
- "Undetermined" = e < 0 for a row with a non-zero or no prediction.
- `make_table.py` v1 recomputed d from 30-digit strings, which printed 'exact' above ~30 digits (a display bug, caught on the rendered page). It now keeps numtab.py's full-precision d and flags.
- `numtab.py` first stopped on an ORIGIN_SHA256 that listed itself; the file was regenerated without the self-entry.
- The window refits use a smaller K2 (18-21) than the base fits, so (a) mixes window and truncation effects. That is intended (conservative) and stated.
- The dominant-at-0 Q is used: it is j-dependent between points (e.g. NUM2 A uses s - l1, B uses s - l0). The other Q_j are not tabulated here; the earlier fits show the same agreement.

## Files
- Code: `numtab.py`, `make_table.py`.
- Logs: `run_numtab.log/stdout`, `run_make_table.log`.
- Outputs: `table.json`, `table_final.json`, `table.csv`, `table.tex`, `texcheck/`.
- Inputs: `inputs/` (79 files, ORIGIN_SHA256).
