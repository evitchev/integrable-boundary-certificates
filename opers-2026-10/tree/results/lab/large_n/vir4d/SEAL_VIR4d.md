# SEAL_VIR4d -- Fable seat.  Blind spin-11 and spin-13 hold-out for the family C_k.
Written before computing any order k >= 11 of any C_k equation and before opening any
Solution-3 data of spin >= 11 (I have opened none, ever: not vev_profile_sol3_w12.json, not
...w14.json, nor any other spin >= 11 Sol-3 table or log line).

## Object
C_k (VIR4c, seal 34cfbf14...): t_k = (3k+2)/(k+2), n = k + 4, M = 1/k, s = (n-1)/2, theta = x d/dx,
   prod_{e in g} (theta - e) psi = x^n (x^{n/k} - E) psi,   sign (-1)^n as in DDMST (3.18),
   g = {s +- l_0, s +- l_1} U {s + (n/k) j : j = -(k-1)/2, ..., (k-1)/2},
   (l_0, l_1) = n sqrt((k+1)/k) (PX, pi).

## Prediction file (written BEFORE any table is opened)
PREDICTION_VIR4d.json: for each (k, spin) below, the polynomial R_(spin+1) of my WKB engine
(wkb_lib.py, unchanged from VIR4) at M = 1/k, with l_i^2 -> n^2 ((k+1)/k) (PX^2, pi^2), divided
by its PX^(spin+1) coefficient ("monic in PX"), as a dictionary "i,j" -> exact rational string
for the coefficient of PX^i pi^j.  Also recorded per entry: the leading coefficient before
normalisation (must be nonzero), and that the neighbouring odd orders R_11, R_13 vanish.
Fibres and spins: k = 4 (t = 7/3), k = 6 (t = 5/2), k = 8 (t = 13/5): spins 11 and 13;
k = 3 (t = 11/5): spin 11 (lead's list) and spin 13 (my addition; 13 is not a multiple of 7).
Its sha256 goes to the lead; I open no table until the lead answers "registered".

## Comparison (after registration)
Tables: results/lab/anchor11/vev_profile_sol3_w12.json and results/lab/anchor13/...w14.json,
the 'vev' field, evaluated at t_k.  Conversion: P_art^2 = -PX^2, Q_art^2 = rho^2 - pi^2,
rho = p/2 - 1/p, p^2 = 2(t-1)/(t+1); then divide by the PX^(spin+1) coefficient.
Screen: every table coefficient finite at t_k and the PX^(spin+1) coefficient nonzero; a fibre
failing the screen is reported as unscreened.
Count: a monic even symmetric polynomial of degree 12 has 27 further coefficients with i >= j
(28 monomials with i >= j including the leading one); degree 14: 35 further (36 in all).
PASS(k, spin) = every coefficient agrees exactly (and the table is symmetric under PX <-> pi).

## Prediction
All eight (k, spin) cells PASS: confidence 85%.

## Controls
- the engine at these orders reproduces BLZ: R_12 and R_14 of the BLZ equation are checked
  against the k-recurrence only indirectly; so the control is instead the registered t = 2
  result of cc (independent engine) -- I recompute k = 2 at spins 11, 13 AFTER registration and
  compare with the tables as a ninth and tenth cell, labelled non-blind (the lead told me they match).
- tamper: after the comparison, the same script run with M = 1/k + 1/10 at k = 4 must FAIL.
Exit codes of the comparison script: 0 = all blind cells pass and the tamper fails;
2 = some blind cell fails; 1 = screen/control failure.
