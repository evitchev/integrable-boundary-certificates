# NUM2 -- numerical spectral check of the fourth-order Sol 3 operator

cc (Opus 5.5 seat), 2026-10-04/05. Lead's commission at the operator's request (brief dfbc1440...).
- Seal: `SEAL_NUM2.md` 25ab6266cc9db62737e803213e21bc6795ed67a5265629af5a423837a04c9074 (2026-10-04 23:23:13 EDT). Snapshotted by the lead at 03:23:29Z.
- Refined prediction: `PREDICTION_NUM2_logE.json` 01174712d94ac43e11e253932c21c5bf7692dd68c896b4eb33ed46b0167a17d7 (23:50:32), registered before any t = 2 data existed.
  - At that time no data_*.json existed. Two t = 9/4 run logs already held raw log|Q| at E = 1000 and 1301. The t = 2 logs were empty or not yet started (lead's check).
- Every hash was sent only after its command returned.

## Conventions
- **K is the WKB order** (s_K), and **spin = K - 1** for even K. The power is E^((k+1)(1-K)/n).
- At t = 2: E^(-1/2) is K = 2 (spin 1), E^(-3/2) is K = 4 (spin 3), E^(-5/2) is K = 6 (spin 5). (The paper's index is spin 2K - 1; mine is not.)
- E is the three-term energy. The WKB regime is E -> +infinity, i.e. E_Gamma = -E -> -infinity.

## Verdict
1. **Theorem 1 holds at the SPECTRAL level at t = 2 (T1, sealed): PASS at ~56 digits.**
   - For both points, E = -20 … 500 and all four dynamical exponents: Q3_j(E) / Q6_j(E_Gamma = -E) = Q3_j(0)/Q6_j(0) (exact Meijer-G values), with deviation ≤ 2.7e-56 (sealed: ≥ 25 digits).
   - The sixth-order form's two frozen-exponent Q (theta = 1, 4) have no counterpart; they are listed in t1.json.
2. **The WKB charges ARE the asymptotic coefficients of log Q_j for every j**, with zero parameters, once non-WKB terms are included.
   - All non-WKB terms belong to ONE family, at E^(-m(k+1)/k) = E^(-m c/n), m = 1, 2, 3. They come from the small-x region at m-th order in the x^(c/2) T x^(c/2) coupling (flagged in the seal for m = 1).
   - **t = 2 (k = 2): the family COLLIDES with WKB grades K = 4, 7, 10.** There the absolute WKB coefficient has a Gamma pole. log Q carries **E^(-3/2) log E**, with coefficient -R/2, where R is the residue in k of the WKB coefficient.
     - **This was registered before any t = 2 data. Fitted: 16-17 digits at all 12 Q_j** (3 points x 4 exponents), against fit stability of 12-13.
     - The sealed basis without log E is INCOHERENT at t = 2 (grades 6-14 digits, odd grades non-zero). **So the t = 2 agreement rests on the refined, registered prediction.**
     - With it:
       - grade 0: 31-32 digits;
       - spin 1: 23 digits;
       - spin 5 (K = 6): 10.6-11.5 digits;
       - K = 8: 6.4-7.3 digits;
       - odd K: 0 (~1e-20).
     - The E^(-9/2) log E coefficient agrees to 3.4-4.1 digits (fit-limited); E^(-3) log E is predicted 0 and fitted 3e-8..1e-7 (stability ~0).
   - **t = 9/4 (k = 10/3): the family is OFF the lattice.**
     - Sealed F1 (lattice only): FAILS.
     - Sealed F2 (+ E^(-13/10)): grade 0 at 21 digits, spin 1 at 11.6, but spin 3 only 4.5 and spin 5 FAILS. **It falls short of the sealed levels at spins 3 and 5.**
     - **Post-hoc F3 (+ E^(-2.6), E^(-3.9), the m = 2, 3 members of the same family):** grade 0 at 30.8 digits, spin 1 at 21.2, spin 3 at 13.8-13.9, spin 5 at 7.9-8.0, K = 8 at 2.5-2.6 (fit-limited), odd K ≤ 2e-17.
     - **F4** (WKB-dressed copies E^(-13m/10 - 13j/22)): the dressing coefficients come out ≈ 0 and unstable, so F3's structure is the right one.
3. **Which object:** every Q_j (exponents s ± l0, s ± l1) carries the same WKB series. It has the same coefficients at the WKB powers for all four, as predicted from Weyl symmetry (P1).
   - The non-WKB coefficients depend on j.
   - **P2 (sealed hypothesis: the non-WKB terms are odd under l -> -l): REFUTED** at t = 9/4. The E^(-13/10) coefficients at point A are -0.5788 (s-l0), -0.9237 (s+l0), -0.4897 (s-l1), -1.0917 (s+l1); they are neither odd nor equal between the two Weyl pairs.
   - At t = 2 the log E coefficient is common to all four Q_j, and the finite E^(-3/2) part has a Weyl-odd j-dependence (both pair averages agree: 1.8920280410761 at A). That pattern belongs to the resonance; it is not a general parity rule.

## Controls (all fire or pass as required)
| control | result |
|---|---|
| X0: generic order-D solver at order 2 vs NUM1 (Sol 1) | M = 1 exact D: 55 digits; M = 5 archived log D: ≤ 8.2e-53 (sealed ≥ 40) |
| X1/X1': exact E = 0 Meijer G, 4th order (t = 2, 9/4) and 6th order | 55-57 digits for every Q_j; beta = c theta_G exactly (sealed ≥ 30) |
| WKB engines: gwkb (Gamma, archived sol12) vs thwkb (6th-order polynomial) | absolute coefficients agree to ~60 digits (branch sign as derived) |
| Tamper (a): WKB l0 + 1/10 | 1.3-2.2 digits at spin 1, worse beyond: fires |
| Tamper (b): numerics with the wrong ordering x^c T | spin 1: -0.2 digits, spin 5: 0.2, odd K = 2.13: fires. (The E^(-3/2) log E coefficient is unchanged: the log term is insensitive to the ordering) |
| Precision doubling (t = 2, A, dps 60 vs 120) | raw log|Q| agree to 7e-55; extracted coefficients agree to 27-52 digits |
| Zeros (Q_(s-l1), A, E ∈ [-300, 0]) | two zeros, -18.666905482576714873965815 and -118.293359232304173588752319: projection = outward-Wronskian at the 1e-25 bisection tolerance (sealed ≥ 10) |

## Registered vs post hoc
- Sealed: the objects; P1 (WKB at non-resonant grades); the flagged non-WKB family at E^(-m c/n) (m = 1 named, with the t = 2 collision); P2; T1; X0/X1; tampers; zeros; F1/F2 at t = 9/4.
- Registered after the seal, before any t = 2 data: the E^p log E coefficients (PREDICTION_NUM2_logE.json).
- POST HOC: F3/F4 at t = 9/4 (`diag_fibre2.py`); the Weyl-pair observations; the ordering-insensitivity of the log term.

## Disclosures
- **The timeout kill (my error).** The first t = 9/4 runs and the dps-120 run were killed by my own `timeout 20000`: launched 23:46:32, killed at 05:19:52, after their last writes at 04:43 and 04:56, with no traceback and no OOM in dmesg. They wrote JSON only at the end.
  - Rerun with per-energy checkpoints (`num2_point.py`, `points/`), one process per energy, 14 concurrent. Peak measured with /usr/bin/time -v: RSS 20.9 MB, VSZ 45 MB. Final cap `ulimit -v 120000` (120 MB per job; 1.7 GB worst case), per the lead's memory constraint.
  - A first relaunch with a 16 GB ulimit was stopped after ~13 min (PIDs read from the saved runner PID with ps --ppid). Logs in `v1_timeout/`.
- `wkb_pred.py` v1 crashed in its thwkb cross-check on the K = 4 pole, before writing JSON. It was fixed to skip pole grades (log `run_wkb_pred_v1_crash.stdout`).
- `ode.py` v1 (slow Taylor step) kept as `ode_v1_slow.py`. v2 only precomputes the integer tables. X0 and the first X1 run used v1 (log `run_x1_v1engine.*`). X1 was RE-RUN with v2 (`run_x1.*`) and passes identically (55-57 digits). X0 was not re-run with v2.
- Fit stability is the agreement between two basis truncations (K2 = 26 vs 22; 24 vs 20 for F3).
- Not done: a closed-form prediction of the non-WKB coefficients (first-order perturbation of the generalised-Bessel region); the finite E^(-3/2) part at t = 2; the even grades of the sixth-order Q6 at large E (only T1).

## Files
- Code: `ode.py`, `ops.py`, `gwkb.py`, `thwkb.py`, `x0_control.py`, `x1_control.py`, `wkb_pred.py`, `residue_pred.py`, `num2_data.py`, `num2_point.py`, `assemble.py`, `num2_fit.py`, `t1_zeros.py`, `diag_fibre2.py`, `tasks*.txt`.
- Logs: `run_*.log/stdout`.
- Data: `data_*.json`, `points/`. Results: `fit_*.json`, `t1.json`, `wkb_pred.json`, `PREDICTION_NUM2_logE.json`.
- Inputs: `inputs/` (num1_dps60.json, hashed).
