# NUM1 -- numerical spectral check of the Sol 1 oper: the formal WKB IS the asymptotic expansion of the spectral determinant, up to one exactly identified alpha-odd non-local term

cc (Opus 5.5 seat), 2026-10-01/02. Lead's commission, at the operator's request.
- Seal: `SEAL_NUM1.md` c1418ec1351aa2672f2521a390da7ff017cff2d0a84c87735ad2908211e9d526 (23:03:11). The hash was sent after the command returned.
- Post-hoc prediction: `PREDICTION_extra_M51.txt` 0f2896c45b30279f4452df37e3c9440e9d987093d7313c150eddad57362e779a (23:58:09), written and sent before the M = 51/10 fit printed.
- Nothing was fitted to the charge tables.

## Object
- -psi'' + (x^(2M) + alpha x^(M-1) + l(l+1)/x^2) psi = E psi (item 306), at t = 2: M = 5, mu = 3/5.
- D(E) = W[y, psi_+], with y the subdominant solution at +infinity, normalised exactly by y ~ x^(-(M+alpha)/2) exp(-x^(M+1)/(M+1)), and psi_+ = x^(l+1)(1 + ...).
- Points (all real):
  - A: alpha = 1, l = 1/3 (X = -1/120, Y = 119/216);
  - B: alpha = -3/2, l = 7/10;
  - C: alpha = 0, l = 1/3.

## Methods (own code)
- `spec.py` (mpmath):
  - y at x_max from the large-x asymptotic Riccati series, with convergence judged on block maxima;
  - Riccati q' = q^2 - V - e integrated INWARD (stable) by own Taylor steps with step-doubling error control;
  - psi_+ from its Frobenius series at x0 = 2/sqrt(e);
  - a linear (y, y') integrator for E > 0 and for outward shooting;
  - rational M supported (exponents in units 1/den(M)).
- `rwkb.py`: a new theta-form Riccati WKB, S^2 + S' - S - lam = x^(2M+2) + alpha x^(M+1) + e x^2.
  - Graded in alpha, d/du and lam; Beta continuation.
  - Sign fixed ONCE: log D = -sum_k int S_k du, from the M = 1 exact case. There, every even grade has ratio exactly -1 (`n1_wkb.py`).
- Development defects, all fixed before any M = 5 comparison (v1 code kept in `v1_slow/`):
  - a zero-term stopping bug in the asymptotic series;
  - a Taylor radius heuristic fooled by stiff noise, giving 10^4 steps;
  - x_max too large;
  - integer-M indexing.

## Sealed results
| test | sealed | result |
|---|---|---|
| N1: M = 1 exact D = (2l+1) Gamma(l+1/2)/Gamma((2l+3-E)/4), E = -10..-10^4 and +2.5, +7.25, -3.5 | >= 25 digits | **52 digits PASS** |
| N2: lowest 3 levels at A, zeros of D vs outward shooting | >= 10 digits | **PASS**: 7.227874944457305085949920, 23.04108840264275887428250, 45.56839039043296518791373, identical to the 1e-30 bisection tolerance |
| E2 blind fit (36 energies 1e3..1e7, dps 60; no WKB input) vs WKB | grade 0 >= 15, spin 1 >= 10, spin 3 >= 7, spin 5 >= 4, odd 3/5 >= 6/4 | grades 0, 2, 3, 4, 5, 7, 8 at A and B: **51-52, 41-42, 37-38, 34-35, 30-31, 24-25, 21-23 digits PASS**. **Grade 6 (e^-3, spin 5) at A and B: < 1 digit -- FAIL.** C (alpha = 0): even grades 0-8 at 21-51 digits, including grade 6 at 27 |
| odd grades at C vanish | yes | PASS: fit values ~1e-23..1e-39 = noise; WKB 1e-55..1e-60 |
| E1 residual scaling | e^(-mu G) | **FAIL**: the residual sits at ~1e-10, set entirely by the grade-6 discrepancy |
| N3 precision doubling (dps 120, A) | unchanged within error | **PASS**: every non-resonant grade improves by ~8 digits (spin 1: 50 digits); the grade-6 discrepancy is identical |
| TAMPER: WKB at M = 5.1 | must miss | PASS: 1.2 / 1.3 / 0.7 digits at grades 2 / 4 / 6 |

## The grade-6 failure, diagnosed (post-hoc, labelled)
1. **Resonance.** At integer M, the WKB grade M+1 sits at the integer power e^(-(M+1)/2).
   - The WKB coefficient there is DISCONTINUOUS in M at alpha ≠ 0 (a 0 × infinity term). The value at M = 5 exactly is not the M -> 5 limit (`diag_grade6.py`).
   - At alpha = 0 it is continuous.
2. **A non-WKB, alpha-odd small-x term.** To first order in alpha, the Bessel region x ~ e^(-1/2) contributes
   **+alpha · Mellin(M+1) · e^(-(M+1)/2)**, with Mellin(s) = ∫_0^∞ t^(s-1) I_nu(t) K_nu(t) dt = Gamma(s/2) Gamma(s/2+nu) Gamma(1/2-s/2) / (4 sqrt(pi) Gamma(1+nu-s/2)), nu = l + 1/2.
   - At M = 5 this is -(4/15) alpha nu(nu^2-1)(nu^2-4).
   - **M = 5 check (`diag_resonance.py`):** lim_(M->5) C_6^WKB + alpha Mellin(6) = the extracted e^-3 coefficient to **28.4 digits (A), 26.7 (B), 27.2 (C)**, i.e. at the fit-stability level.
   - **Generic M = 51/10 (`diag_generic_M.py`; prediction hashed BEFORE the fit printed):**
     - a fit with the WKB family only FAILS (grade 6: 1.3 digits; grade 7/8 garbage);
     - adding e^(-(M+1)/2) restores grades 2-8 at 10-30 digits;
     - its coefficient is **-0.30947207843779164 vs predicted -0.30947207843781443: 12 digits** (fit stability 12).
3. **Doublet average (lead's added test, `avgpm.py`):** with -alpha run at the same l, the ±alpha AVERAGE of every extracted grade equals the alpha-even WKB at M = 5 exactly:
   - A: grade 6 at **28.5 digits**, all other even grades 22-53;
   - B: grade 6 at **27.1 digits**;
   - averaged odd grades are noise-level zero.
   - Both the discontinuity and the Mellin term are alpha-odd and cancel. So **the doublet-averaged determinant has NO anomalous grade**, consistent with item 309 (the cylindrical charges are the alpha-even part).

## Reading
- The formal WKB charges of the Sol 1 oper ARE the asymptotic coefficients of a genuine spectral determinant: every grade to 21-52 digits, with zero parameters. The only exception is the alpha-odd sector at e^(-(M+1)/2).
- There, D also carries a non-local alpha-linear term, the Bessel Mellin transform above. At integer M it collides with grade M+1.
- The ±alpha average (the object the charge tables describe) is anomaly-free to 27-28 digits.
- Grade 1 (fitted): log e coefficient -7/15 at (A, +alpha) and -11/30 at -alpha; at B, -21/40 and -27/40. They are alpha-odd shifted around -(M+1)/(4M)-type values. Not analysed further.

## Not done
- Part (2) (Sol 3 sixth-order): not run. TQ1 and BR1 are queued behind this.
- No higher-order (alpha^3) non-local terms were derived. The fits show none at the achieved precision.

## Files
- Code: `spec.py`, `rwkb.py`, `n1_check.py`, `n1_wkb.py`, `num1.py`, `n2_shoot.py`, `diag_grade6.py`, `diag_resonance.py`, `diag_generic_M.py`, `avgpm.py`.
- Logs: `run_*.log/stdout`, including `run_n2_v1_range30.*`, the first N2 scan to E = 30, which held only 2 levels.
- JSON: `num1_dps60.json`, `num1_dps120.json`, `avgpm.json`. Old code: `v1_slow/`.
