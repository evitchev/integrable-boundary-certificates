# SEAL SOL2F -- the Sol 2 operator at generic t: a zero-parameter family (lead's commission)

cc (Opus 5.5 seat), 2026-10-01. Sealed BEFORE testing any fibre other than t = 2 (SOL12, item 301 pending).
- Fable's stages are not opened. The VIR5d proposal is known only through the lead's message.

## Derivation (from the record's exponents, not fitted)
- With k = 2(t-1)/(3-t) (t = (3k+2)/(k+2)), Sol 2's a_X = -2(t-3)/(t+5) = 2/(2k+3) and a_Y = 4(t-1)/(t+5) = 2k/(2k+3).
- So SOL12's class G with beta = 1 has D = n = 2/a_X = 2k+3, alpha = a_Y D/2 = k, gamma = D(1 - a_X - a_Y) = 1.
  - This is exactly Fable's Sol 2 top family Lambda = T(T^2 - Y)(T^2 - X)^k (item-294-type crossed curve).

## THE SEALED FAMILY S2(k) (zero parameters)
- Equation: **G(s - theta) psi = x^n (x^(n/k) - E) psi**, with n = 2k+3, s = (n-1)/2 = k+1, M = 1/k, delta = nM = (2k+3)/k.
- Symbol: **G(sigma) = delta^k R_k((sigma - l_0)/delta) delta^k R_k((sigma + l_0)/delta) (sigma - l_1)(sigma + l_1) sigma**, with R_k(u) = Gamma(u + (k+1)/2)/Gamma(u - (k-1)/2).
- Momenta: C = 2n(1+M) = 2(2k+3)(k+1)/k; l_0^2 = C P_X^2/a_Y = C (2k+3)/(2k) P_X^2; l_1^2 = C pi^2/a_X = C (2k+3)/2 pi^2.
- For integer k >= 1 it is the POLYNOMIAL ODE of order 2k+3 with exponents:
  - X: s ± l_0 + delta(j - (k-1)/2), j = 0..k-1 (two arithmetic k-strings of step delta, centred at s ± l_0);
  - Y: s ± l_1;
  - frozen: s.
  - Sum = n(n-1)/2.
  - The k = 2 member is SOL12's t = 2 operator, reproduced not refitted.
- DUAL (also sealed): M' = -1/(k+1), with the same C rule C = 2n(1+M') (at k = 2: -1/3, 28/3 = SOL12's second solution).
  - Prediction: the dual gives identical charges at every fibre.

## Fibres and SEALED PREDICTIONS (monic, every coefficient, nothing fitted)
- Certified side: VIR3-R dictionary; spin 1; cyl(2, 4); vev_sol2_w6/8/10.
- Pre-registered fallback if a table has a pole: k ± 1/100 is NOT used; the fibre is dropped and reported.

| fibre | k, n | missing spins (n Z) | sealed prediction |
|---|---|---|---|
| t = 7/3 | k = 4, n = 11 | spin 11 | spins 1-9 EQUAL; spin 11 is IDENTICALLY 0 in the WKB, so it MISSES the regular certified spin 11 (if regular) |
| t = 11/5 | k = 3, n = 9 | spin 9 | spins 1-7 EQUAL; spin 9 identically 0; spin 11 EQUAL (hold-out) |
| t = 9/4 | k = 10/3, n = 29/3 | no odd integer spin in n Z | spins 1-9 EQUAL; spin 11 EQUAL (hold-out); genuine Gamma symbol (X block R_(10/3)) |
| t = 5/3 | k = 1, n = 5 | spin 5 | see below |

- t = 5/3 extra fibre: the S2(1) symbol R_1 R_1 sigma with M = 1, C = 20 coincides IDENTICALLY with Sol 3's k = 1 member. So it predicts Sol 3's charges.
  - It can equal Sol 2 only where the record's Sol 2 = Sol 3 at t = 5/3 (W >= 8, i.e. spins 7, 9; item-270 annotation).
  - Prediction: spins 7, 9 EQUAL; spins 1, 3 EQUAL iff Sol 2 = Sol 3 there (checked and reported); spin 5 missing.
- **Priors:** k = 4, 3: 65%; k = 10/3: 55%; dual identical: 85%.
- If a fibre fails: report the first failing spin and coefficient. No refit, except ONE disclosed diagnostic: refit (M, C) at that fibre from spins 1-3, as in SOL12.

## Spin-11 hold-out
- Compute spin 11 at k = 3 and 10/3 (and k = 4, where it is predicted to be 0), write PREDICTION_SOL2F.json, hash it, and send the hash BEFORE evaluating vev_profile_sol2_w12.json at t ≠ 2.
- Disclosure: that table was opened in SOL12 and evaluated at t = 2 only (and at t = 2.01 as a negative control). Its symbolic form was printed nowhere. Hence the hold-out is hash-registered but not strictly blind.

## Controls
- Engine controls re-run (BLZ, Sol 3 k = 2, -6, 10/3).
- Regression: S2(2) reproduces SOL12 at t = 2, spins 1-9.
- TAMPER: the X-string step delta -> delta + 1/10 (integer k) must fail at spins 3-9.
- Scoring tamper after the hash guard.
