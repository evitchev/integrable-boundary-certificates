# NUM3 -- numerical spectral check of the third-order Sol 2 operator

cc (Opus 5.5 seat), 2026-10-07. Lead's commission at the operator's request (briefs NUM3 26366974..., COMMON dbcb573e...).
- Seal: `SEAL_NUM3.md` b65ef9e83c1fc55be3842b492bf93b3ed185b9ea505ee7122150a5af602e990d (03:17:18 EDT); snapshotted by the lead at 07:17:30Z.
- Refinement registered before any data: `REGISTRATION_NUM3.txt` 114f105fb562f02ccd3f3012a862dfff8260a72837ca48edc4c2a661cdd8d685 (04:18:47 EDT), which records:
  - PREDICTION_NUM3_logE.json 048c36e6...;
  - wkb_pred3.json 4fa7008b...;
  - "no points/ or data_* present".
  - It was registered LOCALLY during the lead switch (the old lead asked for no messages), sent to the new lead before any data job, and snapshotted at 08:19:20Z.
- Conventions as in NUM2: K is the WKB order, spin = K - 1; E is the three-term energy; E -> +infinity.

## Verdict
1. **Theorem 1 at SPECTRAL level, t = 2 (T1, sealed): PASS.**
   - The Mellin multiplier's zeros are exactly the four extra l0-exponents (5/4 ± l0, 19/4 ± l0) of the seventh-order Gamma form.
   - The three three-term exponents s, s ± l1 correspond one-to-one to the seventh-order 3, 3 ± l1, and Q3_j(E)/Q7_j(E_Gamma = -E) equals the exact E = 0 Meijer-G ratio.
   - Both points, E = -20 … 500: worst deviation 8.2e-37 (at E = 2), typically 1e-55 (sealed ≥ 25 digits).
2. **P1 (sealed): the WKB charges are the asymptotic coefficients of log Q_j for all three j, at both fibres, with zero parameters.** This holds once the sealed non-WKB terms are in the basis, and on the post-hoc extended grid (point 3).
   - Final fits (E ∈ [1e4, 1e9], sealed bases), in digits:

     | | grade 0 | spin 1 | spin 3 | spin 5 | K = 8 | odd K |
     |---|---|---|---|---|---|---|
     | t = 2 | 33.8-34.1 | 25.0-25.3 | 17.9-18.1 | 12.1-12.4 | (collision cell) | ≤ 5e-8 |
     | t = 9/4 | 31.8-32.5 | 22.9-23.6 | 16.0-16.7 | 10.7-11.4 | 6.3-7.0 | ≤ 3e-5 |

   - Sealed targets: grade 0 ≥ 15, spin 1 ≥ 10, spin 3 ≥ 7, spin 5 ≥ 4.
   - Fit stability (two truncations): 28-30 / 19-22 / 13-15 / 7.7-9.6 digits.
3. **The sealed grid [1e3, 1e7] was TRUNCATION-LIMITED for Sol 2** (lattice spacing E^(-3/7) and E^(-13/29); fast coefficient growth).
   - Sealed-basis fits there: grade 0 at ~20 digits, spin 1 at 12-13, spin 3 at 6-7, spin 5 at ~2. **Spins 3 and 5 are below the sealed levels on the sealed grid.**
   - Post-hoc diagnostic (`diag_window.py`): dropping E < 1e4 improves every grade by ~8-10 digits. Hence the **POST-HOC GRID EXTENSION** E = 10^(7 + 2i/12), i = 1..12 (`num3_point_ext.py`), with the sealed bases unchanged, giving the table above.
4. **The non-WKB family E^(-m(k+1)/k) is confirmed, with the sealed collision arithmetic:**
   - **t = 2:** the Gamma-pole grades of the WKB engine are exactly K = 8 and 15 (m = 2, 4), as sealed.
     - **The sealed coefficient -(9/14) R_8 of E^-3 log E was refined before any data to R_8 = 0**: the residue is ~1e-47 at eps = 1e-25, a 0·∞ cell. So no log term is predicted.
     - Fits WITH and WITHOUT the log terms are equally good (residual 2e-36 vs 4e-35; same grade digits). The fitted log coefficient is ~2e-5 with no stability, i.e. undetermined, consistent with 0.
     - Odd m: E^(-3/2) is a separate term. N_(-3/2) at point A for (s, s-l1, s+l1) = 6.0629455, 8.3740312, 4.7515900; at B 7.3886716, 8.1194338, 6.7785876. Stable to ~13 digits.
   - **E^-3 collision cell (post hoc; contains the non-WKB m = 2 term):**
     - fitted E^-3 coefficient minus lim_(k->2) C_8: A: -257.0, +109.6, +122.4; B: -38.0, +2.4, +7.4.
     - The j-dependent non-WKB contamination is as expected at a collision. K = 8 is not a WKB test at t = 2.
   - **t = 9/4:** no collision (the first is at K = 30). N1 (E^(-13/10)) is A: -39.997601, -38.091397, -41.555227; B: -43.664195, -43.026251, -44.261182, stable to 13 digits. The m = 2, 3 terms are in the sealed basis.

## Controls
| control | result |
|---|---|
| X0: NUM1 path (Sol 1, M = 1 exact D) | 57-59 digits |
| X0: archived NUM2 point (Sol 3, t = 2, A, E = 1e3) | max |Δ log|Q|| = 5e-59 |
| X1: exact E = 0 Meijer G, 3rd order (t = 2, 9/4; A, B) | 55-56 digits; beta = p theta_G exactly |
| X1': exact E = 0 Meijer G, 7th order (A, B) | 53-57 digits; beta exact |
| WKB: gwkb pole grades | t = 2: K = 1, 8, 15 (as sealed); t = 9/4: K = 1 only; odd K = 0 |
| Tamper (a) WKB l0 + 1/10 | 1.4-2.8 digits: fires |
| Tamper (b) wrong ordering x^c L(T) | spin 1 -0.3, spin 3 -0.6, spin 5 -0.5 digits; odd K ≠ 0: fires |
| Precision doubling (t = 2, A, dps 120, sealed grid) | raw 6e-56; fit coefficients agree to 30-49 digits |
| Zeros (Q_(s-l1), A, E ∈ [-300, 0]) | one zero, -38.5234778416685267237408972: projection = outward Wronskian (1e-25) |

## Registered vs post hoc
- **Sealed:** the objects; P1; the non-WKB family with the collision arithmetic (even m at K = 1 + 7m/2); the fit bases; T1 with the one-to-one exponent correspondence; controls.
- **Registered before data:** R_8 = R_15 = 0 (refinement of the sealed -(9/14) R_8).
- **POST HOC:**
  - the truncation diagnosis and the grid extension to 1e9;
  - the no-log-term refit;
  - the E^-3 collision cell vs lim C_8.

## Disclosures
- `wkb_pred3.py` v1 crashed in the residue step: the trivial K = 1 pole (E^0 log E) was included. It was fixed to exclude K = 1 and the log kept (`run_wkb_pred3_v1_crash.stdout`).
- The fits on the sealed grid miss the sealed spin-3/5 levels (point 3). The verdict rests on the post-hoc extended grid.
- Memory: one job measured at 20.7 MB RSS, 45 MB VSZ. Runs used `ulimit -v 120000`, 14 concurrent, on the head only (free memory checked). Runner PIDs are in pid_*.txt. No outer timeout on the data runs.
- The cluster nodes dsk002-005 do not share ~ (checked) and were not used.

## Files
- Code: `ode.py`, `ops3.py`, `gwkb.py`, `thwkb.py`, `controls.py`, `wkb_pred3.py`, `num3_point.py`, `num3_point_ext.py`, `assemble3.py`, `num3_fit.py`, `num3_fit_ext.py`, `t1_zeros3.py`, `diag_window.py`, `diag_k8.py`, `tasks*.txt`.
- Logs: `run_*.log/stdout`. Data: `points/`, `points_ext/`, `data_*.json`.
- Results: `fit_*.json`, `fitext_*.json`, `t1_3.json`, `wkb_pred3.json`, `PREDICTION_NUM3_logE.json`, `REGISTRATION_NUM3.txt`.
- Inputs: `inputs/` (num1_dps60.json, NUM2 data_2_A_60.json; hashed).
