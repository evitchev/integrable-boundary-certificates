# SEAL NUMTAB -- the spectral-evidence table with stability tests (lead's commission; brief 8a7dc736...)

cc (Opus 5.5 seat), 2026-10-10. Sealed before any refit.
- Inputs: ARCHIVED data only, copied read-only and hashed from results/lab/large_n/{num1,num2,num3}. No ODE integration is recomputed.
- WKB predictions are re-evaluated from the archived engines (rwkb for NUM1; the stored wkb_pred*.json for NUM2/3). That is symbolic/Beta evaluation, not integration.

## Datasets and BASE models (the final fits of record)
- **D1 NUM1** (Sol 1, M = 5): points A, B, C; 36 energies in [1e3, 1e7].
  - Base: e^mu, log e, 1, e^(mu(1-K)), K = 2..30 (mu = 3/5).
  - Also the ±alpha doublet average at A, B (from avgpm.json).
- **D2 NUM2** (Sol 3, t = 2): A, B, C; [1e3, 1e7].
  - Base: lattice E^((1-K)/2), K = 2..26 + log E + 1 + E^p log E (p = -3/2, -3, -9/2), the registered refinement.
- **D3 NUM2** (Sol 3, t = 9/4): A, B; [1e3, 1e7].
  - Base: lattice + log E + 1 + E^(-13m/10), m = 1, 2, 3 ("F3"), K2 = 24.
- **D4 NUM3** (Sol 2, t = 2): A, B; the extended grid E ∈ [1e4, 1e9] (27 + 12 = 39 points).
  - Base: lattice E^(3(1-K)/7), K = 2..24 + log E + 1 + E^(-3/2), E^(-9/2), E^(-15/2) + E^(-3) log E, E^(-6) log E (the sealed basis).
- **D5 NUM3** (Sol 2, t = 9/4, extra, not in the brief): A, B; extended grid; lattice + E^(-13m/10), m = 1..3; K2 = 24.
- Every fit: least squares (mp.qr_solve) at 3x the data's working precision, on log|Q| of the j = dominant-at-0 Q (NUM2/3; D1: log D). Other j are reported in JSON only.

## Row status classes (one row per fitted power)
- **WKB-predicted:** grades K = 0, 2..8 at non-resonant powers. Prediction = -int s_K; odd K are predicted 0.
- **Non-WKB, registered before data:** NUM2 t = 2 log E coefficients (-R/2 etc.); NUM3 t = 2 R_8 = R_15 = 0 (i.e. no E^-3/E^-6 log E). The predicted value is shown.
- **Non-WKB, position sealed, coefficient fitted free:** NUM2 t = 9/4 E^(-13/10) (m = 1, sealed F2); NUM3 E^(-3/2) (t = 2) and E^(-13/10) (t = 9/4).
- **Non-WKB, post hoc:**
  - NUM1 grade 6 at alpha ≠ 0 (lim WKB + alpha Mellin(6), derived after the sealed failure, confirmed at M = 51/10 by a registered prediction);
  - NUM2 t = 9/4 E^(-2.6), E^(-3.9);
  - NUM3 t = 2 E^(-9/2), E^(-15/2) positions.
- **Fitted free:** log e / log E (grade 1) and the constant.
- **Resonant finite part, UNRESOLVED:** NUM2 t = 2 E^(-3/2) and E^(-3) finite parts; NUM3 t = 2 E^-3 (lim C_8 + non-WKB m = 2).
  - Shown without a prediction; lim C_K is listed for reference where it exists.

## Agreement digits
- d = -log10(|fit - pred|/|pred|). For pred = 0: the absolute |fit|.

## STABILITY tests (sealed)
- **(a) Fit window.** Two shifted windows of the base grid:
  - W-: drop the lowest-energy quarter of the points (floor(NP/4));
  - W+: drop the highest-energy quarter.
  - Each is refitted with the base model at the largest K2 ≤ base K2 such that #unknowns ≤ NP_window - 4.
- **(b) Truncation.** The base window at the base K2 and at K2 - 4.
- **(c) Working precision.** The archived 120-digit runs, refitted with the identical model and window:
  - D1 A;
  - D2 A;
  - D4 A, on the common subset E ∈ [1e4, 1e7] only. 120-digit data do not exist above 1e7, so the D4 precision test uses that subset at the K2 given by the window rule; stated as such.
  - Everywhere else: "no 120-digit archive".
- **ERROR ESTIMATE** epsilon(row) = max of the spreads |base - W-|, |base - W+|, |base(K2) - base(K2 - 4)|, |dps60 - dps120| (where available).
  - Reported as e = -log10(epsilon/|value|) next to d.
  - Consistency flag: "consistent" if |fit - pred| ≤ 10 epsilon, else "TENSION". Rows without a prediction: e only.

## Deliverables
- table.json (all datasets, points, rows, the four spreads);
- table.tex: booktabs, one page; point A of D1-D4, all rows, d and e columns, status symbols. Other points are in the JSON;
- a 5-line caption and a one-paragraph method note;
- NOTE.md + SHA256SUMS.
