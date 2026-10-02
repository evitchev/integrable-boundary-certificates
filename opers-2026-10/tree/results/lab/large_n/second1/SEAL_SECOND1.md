# SECOND1 seal. Opus 5 (2), 2026-10-01. Before any comparison with Sol 3's operator.
DISCLOSED pre-seal structural computation (s0_shifts.py/.log, no Sol 3 data used):
- The tensor product of A_pm: theta^2 chi = u_pm theta chi + (lam_pm + x^a - E x^b) chi (common right side) has EXACTLY two x-shifts:
      x^0: prod_(i,j) (theta - r_i - q_j)  (r, q = the exponents of A_+, A_-);
      x^a: -(2 theta + a - U)(2 theta + 2a - U);
      x^b: E (2 theta + b - U)(2 theta + 2b - U),   U = u_+ + u_-.
## Test (sealed)
- Sol 3's Gamma form (item 298) has the Mellin recurrence P4(s-nu) G_k(s-nu) psihat(nu) = psihat(nu - c) - E psihat(nu - n), with n = k+4, sigma = n/k, c = n + sigma.
- A multiplier phihat_T = m psihat_S plus an overall factor g must satisfy (with a = c, b = n, or the swapped assignment):
      m(nu)/m(nu - c) = f_c(nu),   m(nu)/m(nu - n) = f_n(nu).
  Here f_c and f_n are rational times G_k ratios, with Gp(nu) = G_k(s-nu) P4(s-nu)/p0(nu) carrying a GENERAL theta-part p0 (rotation not imposed).
- For incommensurate c and n this requires the COCYCLE condition
      Gp(mu - sigma)/Gp(mu) = p_b(mu) p_a(mu - c)/(p_a(mu + n - c) p_b(mu - c)),
  identically in mu. This is a rational identity: E cancels, and l0, l1 are symbolic.
- Unknowns: the exponents r1, r2, q1, q2 (U = their sum). Solve exactly at k = 10/3 (t = 9/4) and k = 2 (t = 2), both assignments.
- Pencil (sealed): with p0 = P4 (the rotation of item 309), the condition needs two factor cancellations, which requires c = 2n, i.e. k = 1. So NO at generic k for that sub-case.
- PREDICTION: NONE for the whole ansatz at k = 10/3 and k = 2: 80%.
- Controls:
  - (i) the same machinery at k = -2 (item 309 B2) must give a solution;
  - (ii) the wrong-rotation and wrong-shift controls must fail at k = -2.
- If it closes at generic k: symbolic k, then Sol 2.
- Full logs; any solver capped (ulimit -v 20 GB).
