# SEAL NUM3 -- numerical spectral check of the third-order Sol 2 operator (lead's commission; briefs 26366974..., dbcb573e...)

cc (Opus 5.5 seat), 2026-10-07. Sealed BEFORE any computation. The hash is sent only after the command returns.
- Read: both briefs; register items 301/303/307/325 (already known); my NUM2 stage (my own code).
- Conventions as in NUM2: K = WKB order (s_K), spin = K - 1; E = three-term energy; WKB regime E -> +infinity (E_Gamma = -E -> -infinity).

## Objects
- Three-term ODE (item 307), multiplied by -1 so that P is monic in theta:
  [-T(T^2 - l1^2) + x^c ((T - c/2)^2 - l0^2) + E x^n] phi = 0, T = s - theta, s = k+1, n = 2k+3, c = n/b, b = k/(k+1).
  (x^(c/2) L(T) x^(c/2) = x^c L(T - c/2).)
- Frobenius exponents theta ∈ {s, s - l1, s + l1}.
- Large x: S^3 + x^c S^2 = 0, so S = -x^c and y ~ x^beta exp(-x^c/c). Normalised as in NUM1/NUM2 (full asymptotic series of S, no additive constant).
- **Q_j: y = sum_j Q_j chi_j**, with chi_j = x^(theta_j)(1 + ...).
- Fibres:
  - t = 2: k = 2, n = 7, c = 21/2, s = 3, sigma = n/k = 7/2;
  - t = 9/4: k = 10/3, n = 29/3, c = 377/30, s = 13/3, sigma = 29/10.
- Points (operator exponents, real): A (l0, l1) = (1/3, 3/5); B (9/10, 1/5).

## Predictions (sealed)
- **P1 (WKB):** for each j, the coefficient of log Q_j at E^((k+1)(1-K)/n) (K = 0 and K ≥ 2, non-resonant) equals -∫ s_K of the archived Sol 2 Gamma engine.
  - Engine: gwkb, blocks (k, ±l0, δ), (1, ±l1, δ), (1, 0, δ), M = 1/k, as in sol2f. Odd K: zero.
  - Digits: grade 0 ≥ 15; spin 1 ≥ 10; spin 3 ≥ 7; spin 5 ≥ 4. Prior 70%.
- **Non-WKB family** (as in NUM2): E^(-m(k+1)/k) = E^(-m c/n), from the small-x region at m-th order in the coupling.
  - **t = 2:** WKB powers are E^(-3(K-1)/7) and the family is E^(-3m/2). They COLLIDE only for even m: m = 2 at K = 8 (E^-3, spin 7), m = 4 at K = 15.
    - m = 1 (E^(-3/2)) and m = 3 (E^(-9/2)) are OFF the lattice: separate terms.
    - At K = 8 I predict a Gamma pole in the absolute WKB coefficient. Then log Q carries E^(-3) log E with coefficient R_8 (dp_8/dk - dp'_2/dk) = **-(9/14) R_8**, where R_8 = Res_(k=2) of the WKB coefficient and p_K = (k+1)(1-K)/(2k+3), p'_m = -m(k+1)/k.
    - The numerical R_8 (and any other pole the engine shows) will be REGISTERED before any data is read.
  - **t = 9/4:** WKB E^(-13(K-1)/29), family E^(-13m/10). No collision before K = 30.
- **Sealed fit bases** (NUM2's lesson, no post-hoc additions needed):
  - t = 2: WKB lattice E^(3(1-K)/7) for K = 0, 2..K2, log E, 1, plus E^(-3/2), E^(-9/2), E^(-15/2) (odd m), plus E^(-3) log E and E^(-6) log E (even m at the collisions).
  - t = 9/4: lattice + log E + 1 + E^(-13m/10), m = 1, 2, 3.
  - Also reported: the fit WITHOUT the non-WKB terms.

## Theorem 1 at spectral level, t = 2 (T1, sealed)
- Seventh-order Gamma form (item 301): -(theta-3) prod_± (theta - 5/4 ∓ l0)(theta - 19/4 ∓ l0)(theta - 3 ∓ l1) psi = x^7 (x^(7/2) - E_Gamma) psi, E_Gamma = -E.
  - Its Frobenius exponents: {3, 3 ± l1} ∪ {5/4 ± l0, 19/4 ± l0}.
- **Correspondence:** the Mellin multiplier m(nu) = prod_(r = ±l0) sigma^(nu/sigma)/Gamma(gamma_r - nu/sigma) VANISHES at nu = s - r + sigma(1-k)/2 + sigma j = 5/4 ∓ l0 (j = 0) and 19/4 ∓ l0 (j = 1). Those are exactly the four l0-exponents, so phi = m psi kills those components.
  - Prediction: the three-term exponents s, s ± l1 = 3, 3 ± l1 correspond one-to-one to the seventh-order exponents 3, 3 ± l1.
  - **Q3_j(E) / Q7_j(E_Gamma = -E) = C_j = Q3_j(0)/Q7_j(0)** (exact Meijer G), constant in E, for the three j. Tested at E ∈ {-20, -5, 2, 50, 500}, to ≥ 25 digits. Prior 75%.
- The four l0-exponent Q7 have no counterpart; reported.

## Controls
- X0: the NUM1 path (Sol 1, M = 1, exact Gamma D) and one archived NUM2 Sol 3 point (data_2_A_60, E = 1e3): reproduced to ≥ 40 digits.
- X1: exact E = 0 Q (Meijer G^(D,0)_(dR,D)) for the 3rd order (both fibres) and the 7th order, generic formula
  Q_i(0) = prod_(j≠i) Gamma(a_j - a_i) / prod_m Gamma(B_m - a_i) · C^(a_i) / [(2 pi)^((sG-1)/2) sG^(-1/2) C^theta_G],
  with a_i = theta_i/p, B_m = 1 + r_m/p (r_m the roots of R), C = rho (-1)^(dR-D+1) p^(dR-D), sG = D - dR, theta_G = [(1-sG)/2 + Σa - ΣB]/sG.
  ≥ 30 digits; beta = p theta_G exactly.
- Tampers: WKB l0 + 1/10; the WRONG ORDERING x^c L(T) (no c/2 shift) in the numerics. Each must break P1.
- Precision doubling (t = 2, A, dps 120): raw values and extracted coefficients unchanged within fit stability.
- Zeros: Q_(s - l1)(E) on E ∈ [-300, 0] at A. If real zeros occur: projection vs outward Wronskian, ≥ 10 digits.
- Cluster discipline: per-energy checkpoints; RSS/VSZ measured before fan-out; ulimit -v ~2-3 x VSZ; no timeout under 6 h; saved PIDs.
