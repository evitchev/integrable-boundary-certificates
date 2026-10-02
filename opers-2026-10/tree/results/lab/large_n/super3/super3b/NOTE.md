# SUPER3b (POST-HOC DIAGNOSTIC; prompted by the lead's observation). Opus 5 (2), 2026-10-01. NOT sealed: this corrects SUPER3 (item 293).
## The error in SUPER3
- run_super3.py used alpha^2 = X + tA and lam = r Y + tL, which FIXES the common momentum scale.
- The seal's sentence "only the ratio sY/sX, tX and tY matter after X^k normalisation" is WRONG. Monic normalisation removes one overall factor per spin. A common scale s weights the degree-d layer by s^(d-k), so s is a genuine dictionary parameter.
- SUPER3's "lower layers FAIL" is therefore a DICTIONARY ARTIFACT. SUPER3's top-layer and parity results are unaffected.
## Result
- run_super3b.py: with alpha^2 = s X + tA, lam = s r Y + tL, the 4 unknowns are fitted on spin 3 (XY, X, Y, 1), or on spin 5 where no spin-3 point exists. At every one of 10 fibres, ONE of the two roots gives a full MATCH at spins 3, 5, 7, 9. That root always has tA = 0 and s < 0.
- Closed form read off the fits, then checked with ZERO parameters (closed_dict.py, exit 0):
      alpha^2 = -4 M (M+1) X,   (l + 1/2)^2 = (M-1)^2/4 - (M+1) Y,   M = (t+3)/(t-1),
  in Suzuki's ODE  -psi'' + (x^(2M) + alpha x^(M-1) + l(l+1)/x^2) psi = E psi  (the l-term is our extension of Suzuki's eq. (1)).
  - **ALL 882 non-normalised certified coefficients at spins 3, 5, 7, 9 MATCH, at 19 fibres.**
  - The fibres are 2, 4, 6, 7, 8, -4, 5/2, 9/2, the exact point -1/3, 13/3, 17/5, 25/3, 12, 14, -7/5, 23/7, 3/7, -17/4 and 1e6.
  - 10 of them were never used in any fit: 13/3, 17/5, 25/3, 12, 14, -7/5, 23/7, 3/7, -17/4, and 8 for the closed form.
- The record's Sol 1 charges are, through spin 9 at every tested fibre (including the exact point and the paperclip regime t = 1e6), the WKB local IMs of this scalar ODE. It is consistent with Fable's VIR7 at loss 1 through K = 20 (not compared here).
- Caveats:
  - Spin 3 is untested at the fibres without a certified spin-3 point.
  - The Gamma-pole fibres (2, 4, 6, 7, -4) are evaluated as limits (stageSUPER3/limit_check_t7.log: linear convergence).
  - Spins beyond 9: no tables were used.
