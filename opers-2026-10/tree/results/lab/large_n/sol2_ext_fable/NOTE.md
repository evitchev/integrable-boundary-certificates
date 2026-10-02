# SOL2F -- the Sol 2 operator at generic t: a zero-parameter family

cc (Opus 5.5 seat), 2026-10-01. Lead's commission.
- Seal: `SEAL_SOL2F.md` 80274804e222a595756ca2002e442f2f13f039dd399d5ee3c1aba529e2478c4d (17:13:24), sent before any new fibre.
- Prediction: `PREDICTION_SOL2F.json` 48f1e50bfde295b6205d7c680aa7c3d8f8ab7d9d89c3a9cd5af466b7a0b0e3df (17:18:52). The lead registered it; the table was evaluated only after the lead's go-ahead.
- Fable's stages were not opened.
- Exit codes: `t_engine.py` 0; `sol2f.py` 0; `diag_k3_spin3.py` 0 and `diag_k3_projective.py` 0 (post-hoc, labelled); `predict.py` 0; `holdout.py` 0.

## The family S2(k) (sealed; zero parameters)
- Equation: **G(s - theta) psi = x^n (x^(n/k) - E) psi**, with n = 2k+3, s = k+1, M = 1/k, delta = n/k.
- Symbol: **G(sigma) = delta^(2k) R_k((sigma - l_0)/delta) R_k((sigma + l_0)/delta) (sigma^2 - l_1^2) sigma**, with R_k(u) = Gamma(u+(k+1)/2)/Gamma(u-(k-1)/2).
- Momenta: C = 2n(1+M), l_0^2 = C(2k+3)/(2k) P_X^2, l_1^2 = C(2k+3)/2 pi^2, with k = 2(t-1)/(3-t).
- Derived, not fitted:
  - Sol 2's record exponents give a_X = 2/(2k+3) and a_Y = 2k/(2k+3), so SOL12's class G (beta = 1) has (alpha, beta, gamma) = (k, 1, 1) and D = 2k+3. This is Fable's T(T^2-Y)(T^2-X)^k.
  - M and C are SOL12's t = 2 fit read as M = 1/k and C = 2n(1+M).
- Integer k: a polynomial ODE of order 2k+3. Its exponents are:
  - X: two k-strings s ± l_0 + delta(j - (k-1)/2);
  - Y: the pair s ± l_1;
  - frozen: s.
- Dual (sealed): M' = -1/(k+1) with the same C rule. It gives IDENTICAL charges at every fibre and spin tested.

## Results (both M and dual; every coefficient; nothing fitted)
| k (t), n | spins 1-9 | spin 11 (hash-registered) | missing spin (n Z) |
|---|---|---|---|
| 2 (2), 7 -- regression = SOL12 | 1, 3, 5, 9 EXACT | 0/28 (SOL12) | 7: WKB identically 0 |
| 4 (7/3), 11 | 1, 3, 5, 7, 9 ALL EXACT | WKB identically 0 as sealed; certified regular (28 terms) | 11 |
| 3 (11/5), 9 | 1 EXACT; 3 see below; 5, 7 EXACT | **0/28** | 9: WKB identically 0 |
| 10/3 (9/4), 29/3 | 1-9 ALL EXACT (genuine Gamma X block) | **0/28** | none (n not an integer) |
| 1 (5/3), 5 | 1, 3, 7, 9 EXACT | **0/28** | 5: WKB identically 0 |

- Hold-out controls: t + 1/100 gives 27/28 mismatches at each scored fibre; a scoring tamper gives exactly 1. The lead's independent projective comparison gives the same (0/28; 27/28).
- Tamper delta + 1/10 (k = 2, 3, 4): fails at every non-missing spin 3-9.

## The k = 3 spin-3 top-zero (as the lead asked)
- **Sealed score:** the sealed monic scoring printed **'DIFF (1 coeffs)' at spin 3**. That is the sealed output, recorded as such.
- **Cause (post-hoc, `diag_k3_spin3.py`):**
  - The record's cyl(2,4), normalised X^2 = 1, has the factor (5t - 11) in the denominator of EVERY other coefficient. So at t = 11/5 the true spin-3 charge has a vanishing P_X^4 coefficient.
  - The WKB's P_X^4 coefficient also vanishes there.
  - Monic normalisation divides by zero on both sides (sympy produced zoo/nan). The '1 coefficient' is an artefact of that division, not a mismatch.
- **Projective comparison (post-hoc, labelled, `diag_k3_projective.py`):**
  - Certified side: the regular limit lim_(t -> 11/5) (5t - 11) cyl(2,4), giving 72/5 P_X^2 pi^2 - 6/5 P_X^2 - 72/5 pi^4 + 6 pi^2 - 2/5.
  - The WKB charge is exactly 30 times it (PROPORTIONAL in all 5 monomials, both M and dual). A pi^2/1000 tamper breaks proportionality.
- The same (5t - 11) is the family point t_3 of Sol 3 (VIR4c-R's P^4Q^4 zero).
- Pre-registration gap (disclosed): the seal's drop-the-fibre rule screened only the w6-w10 tables for poles, not cyl(2,4).

## What 'Sol 2 = Sol 3 at k = 1 except spin 5' means for the meeting at t = 5/3
- At k = 1, S2(1) and Sol 3's k = 1 member are THE SAME OPERATOR, identically, as written:
  - G = sigma(sigma^2 - l_0^2)(sigma^2 - l_1^2), M = 1, C = 20 (order 5, potential x^5(x^5 - E)).
  - The only difference is which exponent carries the R_k block: R_1 = identity block, so the X block (k = 1) and Sol 3's frozen block (k = 1) coincide.
- Its charges equal the certified Sol 2 AND the certified Sol 3 at spins 1, 3, 7, 9, 11.
  - Checked here: certified Sol 2 == certified Sol 3 at spins 1, 3, 7, 9 (and both equal the operator at 11).
  - They differ ONLY at spin 5 = n, where the operator has NO charge (WKB identically 0).
- **Reading (observation, not proved):** at t = 5/3 the Sol 2 and Sol 3 curves meet in ONE oper.
  - The record's "fourth coincidence" (Sol 2 = Sol 3 at t = 5/3, W >= 8; here also spins 1, 3) is that shared oper.
  - The two solutions differ only in the spin-n charge, which the oper does not supply. That is item 298's missing-spin locus, where the record's spin-n kernel is two-dimensional, so each solution can pick a different spin-5 element.
  - That the two certified spin-5 charges are the two kernel directions was NOT computed here.

## Disclosures
- The w12 table was opened in SOL12 (at t = 2 and 2.01 only), so this hold-out is hash-registered, not strictly blind.
- Engine, common code and inputs were copied from SOL12; only the seal guard changed. The engine controls were re-run (PASS).
- The two diagnostics are post-hoc.

## Files
- Code: `gwkb.py`, `thwkb.py`, `common.py`, `t_engine.py`, `sol2f.py`, `diag_k3_spin3.py`, `diag_k3_projective.py`, `predict.py`, `holdout.py`.
- Logs: `run_*.log/stdout`. JSON: `sol2f.json`, `PREDICTION_SOL2F.json` (+ .sha256).
- Inputs: `inputs/`, `inputs_holdout/` (hashed).
