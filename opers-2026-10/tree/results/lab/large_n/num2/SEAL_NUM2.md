# SEAL NUM2 -- numerical spectral check of the fourth-order Sol 3 operator (lead's commission; brief dfbc1440...)

cc (Opus 5.5 seat), 2026-10-04. Sealed BEFORE any computation. The hash is sent only after the command returns.
- Read before sealing: the brief; register items 294, 298, 307, 311, 317 (as already read); FILE NAMES of results/lab/large_n/vir8 and sol12 (no content opened today).
- Engines: my own. The new numerics extend NUM1's spec.py to arbitrary order. The WKB side is the archived Gamma engine gwkb.py (sol12 archive, my code), with thwkb.py as a cross-check at t = 2.

## 1. The objects (sealed definitions)
- Three-term ODE (item 307): [P(T) - x^(c/2) T x^(c/2) - E x^n] phi = 0.
  - P(T) = (T^2 - l0^2)(T^2 - l1^2), T = s - theta, theta = x d/dx, s = (n-1)/2, n = k+4, c = n/b, b = k/(k+1).
  - Since x^(c/2) T x^(c/2) = x^c (T - c/2), the Mellin recurrence is P(s - nu) phihat(nu) = (s - nu + c/2) phihat(nu - c) + E phihat(nu - n).
- Frobenius basis at 0: chi_j = x^(theta_j)(1 + O(x^min(c,n))), with theta_j ∈ {s - l0, s + l0, s - l1, s + l1}.
- Subdominant solution: as x -> infinity, theta^4 ~ -x^c theta gives S^3 = -x^c, so y ~ x^beta exp(-(3/c) x^(c/3)).
  - NORMALISATION as in NUM1: y = exp(int S du) with S the full large-x asymptotic series. Every decaying term is included and there is no additive constant; beta is the coefficient of the x^0 term of S.
- **Q-functions: y = sum_j Q_j(E) chi_j.** Each Q_j is entire in E.
- **Which object matches the WKB charges, and why:**
  - The WKB integral log y(x) = int S du gives, as x -> 0, the coefficient of the Frobenius solution that DOMINATES at 0: Q_(s - l_max).
  - By analyticity in (l0, l1), every Q_j is a Weyl image of one function: l0 -> -l0, l1 -> -l1, l0 <-> l1.
  - The formal WKB series is even in l0, l1 and symmetric. Hence
    **(P1) for each of the four j, log Q_j(E) as E -> +infinity has, at every WKB power E^((k+1)(1-K)/n), K ≥ 2, the coefficient -∫ s_K du of the decaying-branch theta-form WKB of the Gamma form** (gwkb.py; sign as NUM1: log D = -sum ∫ S_K).
  - Grades K ≥ 2 only. K = 0 (E^((k+1)/n)) is also compared. K = 1 (E^0 and log E) is fitted, not compared.
  - The regime is E -> +infinity in the three-term convention, i.e. E_Gamma = -E -> -infinity.
  - Odd K: zero, at the precision reached.
- **Possible non-WKB terms (stated in advance, NOT hidden in a fit):** at first order in the coupling x^(c/2) T x^(c/2), the small-x (generalised-Bessel) region can contribute at E^(-c/n) = E^(-(k+1)/k), and at multiples of it.
  - At t = 2 this is E^(-3/2), which coincides with WKB grade K = 4 (spin 3), as NUM1's resonance did.
  - At t = 9/4 it is E^(-13/10), off the WKB lattice E^(-(13/22)(K-1)).
  - **(P2)** If such a term exists, its coefficient is ODD under l -> -l of the dominant pair, so it cancels in log Q_(s-l) + log Q_(s+l). Hypothesis, prior 45%.

## 2. Computation (fibre 1: t = 2, k = 2, n = 6, c = 9, s = 5/2, sigma = n/k = 3)
- Momentum points (real exponents, chosen off every resonance theta_i - theta_j ∈ 3Z):
  - A: (l0, l1) = (1/3, 3/5);
  - B: (9/10, 1/5);
  - C: (1/2, 1/4).
- Numerics (mpmath, dps 60, precision doubling to 120 at A):
  - large-x asymptotic series of S (generic order-D recursion, block-max convergence) at x_max;
  - INWARD Taylor integration of the linear ODE in u = log x with step doubling and log-rescaling (stable: y grows fastest inward);
  - Frobenius series of the chi_j at x0 = 2 E^(-1/n);
  - 4 x 4 projection on (theta^m y)(x0), m = 0..3.
- Blind fit (as NUM1):
  - 36 energies E ∈ [1e3, 1e7];
  - basis E^((k+1)/n), log E, 1, E^((k+1)(1-K)/n) for K = 2..K2 (K2 = 30 and 26 for stability). At t = 2 this is every half-integer power, including the odd-K integer powers.
- **Sealed agreement levels (P1), at non-resonant grades:**
  - grade 0 at least 15 digits; spin 1 at least 10; spin 5 at least 4; odd K: zero to the fit error;
  - spin 3 at least 7 digits **unless resonant**. If a non-WKB term appears at E^(-3/2), I report it the NUM1 way: the M -> integer limit analogue (k -> 2), its size, and its l-parity.
- Prior that P1 holds at all non-resonant grades: 65%.

## 3. Spectral-level test of Theorem 1 at t = 2 (sealed relation)
- Sixth-order Gamma form (items 294/298): (theta-1)(theta-4) prod_± (theta - 5/2 ∓ l0)(theta - 5/2 ∓ l1) psi = x^6 (x^3 - E_Gamma) psi, with E_Gamma = -E.
  - Its subdominant solution: S^6 = x^9, so psi ~ x^(beta') exp(-(2/3) x^(3/2)), normalised the same way.
  - Its Q-functions Q6_j sit on its six Frobenius exponents {1, 4, 5/2 ± l0, 5/2 ± l1}.
- **PREDICTION (T1):** for each of the four DYNAMICAL exponents j,
  **Q3_j(E) / Q6_j(E_Gamma = -E) = C_j, independent of E, and C_j = Q3_j(0)/Q6_j(0)**, both values being known exactly (control X1).
  - Tested at E ∈ {-20, -5, 0, 2, 50, 500} (both signs), to at least 25 digits. Prior 70%.
  - (Failure would mean H1-H3 / boundary matching fail: the subdominant solution is not mapped to the subdominant solution.)
- The two frozen-exponent Q6 (theta = 1, 4) have no counterpart; reported.

## 4. Second fibre: t = 9/4 (k = 10/3, n = 22/3, c = 143/15, sigma = 11/5 > 0)
- Points A and B (as above).
- Fits: F1 with the WKB lattice only; F2 with F1 + E^(-13/10).
- The decisive question: does a non-WKB term exist off the lattice? Reported as observed.

## Controls (each can fail)
- **X0 (engine, NUM1 path):** the generic order-D solver at order 2 must reproduce:
  - Sol 1 at M = 1: the exact D = (2l+1) Gamma(l+1/2)/Gamma((2l+3-E)/4), to at least 40 digits;
  - M = 5: NUM1's archived log D values (num1_dps60.json, three grid points), to at least 40 digits.
- **X1 (exact, fourth order):** at E = 0 the three-term equation is two-term. Its subdominant solution is c · G^(4,0)_(1,4)(c^-3 x^c | b1; a_1..a_4), with a_i = theta_i/c and b1 = (s + c/2)/c, so
  Q3_i(0) = prod_(j≠i) Gamma(a_j - a_i)/Gamma(b1 - a_i) · c^(-3 a_i) · sqrt3/(2 pi) · c^(3 theta_G), where theta_G = (1/3)(-1 + Σ a - b1).
  - This uses the standard G^(q,0)_(p,q) asymptotics (2 pi)^((sigma_G - 1)/2) sigma_G^(-1/2) exp(-sigma_G z^(1/sigma_G)) z^(theta_G), sigma_G = q - p.
  - Numerics must match at all four j to at least 30 digits, and beta must equal c theta_G exactly.
  - X1' is the same for the sixth order (G^(6,0)_(0,6)(9^-6 x^9)).
  - (If X1 fails, I check the asymptotic constant of the G-function first and disclose.)
- **Tampers (must break P1):**
  - (a) WKB side l0 -> l0 + 1/10;
  - (b) numerics with the WRONG ORDERING x^c T (the recurrence factor s - nu + c instead of s - nu + c/2).
- **Precision doubling:** extracted coefficients unchanged within the fit stability.
- **Zeros:** scan Q_(s-l_max)(E) for real zeros on E ∈ [-300, 0] at A. If any exist, compare zeros of the projected Q with zeros of the Wronskian det[chi_i (i ≠ j) integrated OUTWARD to x1 = 1.5 ; y] (a different numerical path), to at least 10 digits.
