# SEAL NUM1 -- numerical spectral check of the Sol 1 oper (lead's commission, operator's request)

cc (Opus 5.5 seat), 2026-10-01. Sealed BEFORE any computation. Nothing is fitted to the charge tables. The hash is sent only after the command returns.

## Object
- -psi'' + (x^(2M) + alpha x^(M-1) + l(l+1)/x^2) psi = E psi (item 306), at t = 2: **M = 5**, mu = (M+1)/(2M) = 3/5.
- Points, all real (alpha^2 = -4M(M+1)X = -120 X, (l+1/2)^2 = (M-1)^2/4 - (M+1)Y = 4 - 6Y):
  - A: alpha = 1, l = 1/3 (X = -1/120, Y = 119/216);
  - B: alpha = -3/2, l = 7/10 (X = -3/160, Y = 59/150);
  - C: alpha = 0, l = 1/3 (control: no alpha).
- **D(E) := W[y, psi_+]**, with:
  - y the subdominant solution at +infinity, normalised exactly by y ~ x^(-(M+alpha)/2) exp(-x^(M+1)/(M+1)) (all decaying corrections included);
  - psi_+ = x^(l+1)(1 + O(x^2)) (the l-branch).

## Numerical method (mpmath, own code)
- At x_max with x_max^(2M) >= 10^4 e: log y and q = -y'/y from the large-x asymptotic Riccati series.
  - q = sum_j a_j x^(M-j), from q^2 - q' = V, with every term down to 10^(-dps).
- Inward integration of the Riccati equation q' = q^2 - V(x) - e (stable inward) together with log y, from x_max to x0 = 2/sqrt(e), by mpmath odefun (Taylor). For E = -e:
  - log D = log y(x0) + log(psi_+'(x0) + q(x0) psi_+(x0));
  - psi_+ from its Frobenius series (entire).
- **Checks of the numerics:**
  - (N1) M = 1, alpha = 0: D(E) = (2l+1) Gamma(l + 1/2)/Gamma((2l+3-E)/4) EXACTLY (derived from Kummer U in the same normalisation). Numeric vs exact at several E (both signs), at least 25 digits.
  - (N2) zeros of D(E) (E > 0) at point A vs direct eigenvalues by shooting (Dirichlet-at-infinity bisection on psi_+ integrated outward), lowest 3 levels, at least 10 digits.
  - (N3) precision doubling (dps 40 -> 80): extracted coefficients unchanged beyond the stated error.

## WKB side (own engine, new: theta-form Riccati)
- In u = log x: S^2 + S' - S - lambda = x^(2M+2) + alpha x^(M+1) + e x^2, with lambda = l(l+1).
- Graded expansion: alpha, the derivative and the -S term have grade 1; lambda has grade 2.
- Terms x^a W^b with W = x^(2M) + e, integrated by the Beta continuation (1/(2M)) e^(b + a/(2M)) B(a/(2M), -b - a/(2M)).
- This is the scheme that reproduced BLZ exactly in VIR3-R/VIR4-R; it is recalibrated here against N1's exact D at M = 1.
- Grade k <-> e^(mu(1-k)):
  - even k = 2j: the spin-(2j-1) charges;
  - odd k >= 3: the alpha-odd "even-spin" charges;
  - k = 0: the leading e^mu;
  - k = 1: e^0 and log e (Beta pole): fitted, not predicted.
- The overall sign/branch convention is fixed ONCE by N1 (M = 1). No per-spin constants.
- Not used: suzuki_wkb.py. The engine is new, and it is cross-checked against item 306's WKB statement only through this comparison.

## Extraction and comparison
- (E1) Residual test: R_G(e) = log D_num(-e) - sum_(k<=G, k != 1) WKB_k(e) - (c log e + d), with c, d fitted.
  - It must scale like e^(-mu G) for G = 4, 6, 8, over e in [10^2, 10^5] (geometric grid, about 12 points).
- (E2) Blind joint fit, no WKB input: log D_num(-e) fitted with free coefficients of e^mu, log e, 1, e^(-mu), ..., e^(-11 mu) (exact solve or least squares on about 16-20 points).
  - Extracted grades 0, 2, 3, 4, 5, 6 compared with WKB.

## SEALED PREDICTIONS (agreement level)
- N1: at least 25 digits. N2: at least 10 digits.
- E1: residual scaling holds, with exponents within 0.05 of -mu G.
- E2, at points A and B: grade 0 (e^mu) at least 15 digits; grade 2 (spin 1) at least 10; grade 4 (spin 3) at least 7; grade 6 (spin 5) at least 4; odd grades 3, 5 at least 6 / 4 digits.
- Prior that all hold: 60%. Most likely failure mode: numerical precision at grade 6, not a WKB mismatch.

## CONTROLS
- WRONG M: the WKB side with M = 5.1 (numerics at M = 5) must miss grades 2, 4 by orders above the stated error.
- Precision doubling: N3.
- alpha = 0 (point C): the odd grades vanish numerically to the error level.

## Part (2) (Sol 3, t = 2, sixth order): only if (1) passes and time allows
- Method: the generalised subdominant solution of the order-6 equation (the solution decaying fastest along the positive axis, normalised by its asymptotics), and D(E) = the coefficient of the s - l0 exponent branch at x -> 0, obtained by matching to the six Frobenius solutions at x0 (a 6x6 Wronskian-ratio solve).
- Compare spins 1, 3 at the item-294 parameters.
- No prediction sealed for (2) beyond 'same scheme, same agreement'; it would be sealed separately if run.
