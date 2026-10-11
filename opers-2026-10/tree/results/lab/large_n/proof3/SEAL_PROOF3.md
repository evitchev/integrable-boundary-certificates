# PROOF3 seal: Theorem 1 at the spectral level. Opus 5 (2), 2026-10-07. Written before any computation.
Read before sealing (definitions only):
- the briefs PROOF3 35faa7f4 and COMMON dbcb573e;
- NUM2 NOTE.md, ode.py header, ops.py (three_term, sixth, exact_three_term_E0, exact_sixth_E0), the t1.json keys;
- register items 311, 321, 325.
A pencil remark made while reading ops.py is part of this seal (see C below), not a computed result.
## Objects (as in NUM2)
- Three-term ODE: L3 phi = 0, L3 = P(T) - x^(c/2) L(T) x^(c/2) - E x^n, T = s - theta, L monic. Gamma form: L6 psi = 0, L6 = P(T) prod_r G_a(T - r) - x^c - E x^n (the SAME E; NUM2's E_Gamma = -E).
- Q-functions:
  - y_sub is the solution decaying along x -> +infinity, NORMALISED by log y_sub = sum over e != 0 of a_e x^e / e + beta log x (no constant term). This is NUM2's asym() convention.
  - chi_j = x^(theta_j)(1 + ...) are the Frobenius solutions at the DYNAMICAL exponents theta_j (the roots of P(s - theta)).
  - y_sub = sum_j Q_j chi_j + (frozen components for the Gamma form).
## Statement to be proved (THEOREM S, finite-order case)
- Assume:
  - a = 1/M a positive integer, so the Gamma form is an ODE of order deg P + a deg L;
  - sigma = n/a > 0;
  - deg P > deg L, L monic;
  - H1: generic momenta (no Frobenius resonances; dynamical exponents off the frozen lattice).
- Then, for every E:
      Q^ODE_j(E) = [ m(theta_j) / kappa ] Q^Gamma_j(E),   m(nu) = prod_r sigma^(nu/sigma) / Gamma(gamma_r - nu/sigma),   gamma_r = (s - r)/sigma + (1-a)/2,
  with kappa != 0 INDEPENDENT of E and of j.
- Hence Q^ODE_j(E)/Q^Gamma_j(E) = Q^ODE_j(0)/Q^Gamma_j(0): NUM2's T1.
- Proof plan:
  (S1) K, the operator with multiplier m (K x^nu = m(nu) x^nu), is realised by iterated Hankel integrals:
           K psi (x) = prod_r [ (1/2 pi i) Int_H dt e^t t^(-gamma_r) ] psi( x (sigma t)^(1/sigma) ).
       It CONVERGES on Gamma-form solutions because their growth order rho = c/n = 1 + M satisfies rho < sigma/deg L <=> deg L < deg P (H2, automatic). It needs sigma > 0 and analyticity in |arg x| < pi deg L/sigma (H3, automatic for ODE solutions on the universal cover).
  (S2) INTERTWINING (from PROOF1's identities I1-I3 plus integration by parts on the Hankel contour): L3 K = x^n K x^(-n) L6. So K maps Gamma-form solutions to ODE solutions.
  (S3) K chi^Gamma_j = m(theta_j) chi^ODE_j (termwise on convergent Frobenius series). K kills the frozen Frobenius solutions, because the frozen exponents plus their shift lattice are exactly zeros of m (finite-order case).
  (S4) K y^Gamma_sub = kappa y^ODE_sub. K y_sub solves L3 and decays like exp(-q x^(c/n') ...) by Laplace's method, so it is the unique subdominant solution up to the constant kappa. kappa is E-independent because E enters both normalised asymptotic expansions only at negative powers of x (beta does not depend on E).
- C (pencil, from reading ops.py; to be machine-checked): at E = 0, NUM2's exact Meijer-G Q's give Q3_j(0)/Q6_j(0) ∝ 3^(theta_j/3)/Gamma(1/3 - theta_j/3) = m(theta_j) by Gauss's multiplication formula (t = 2: sigma = 3, gamma = 1/3, frozen theta = 1, 4 = zeros of m).
## Machine checks (sealed)
- M1: the intertwining identity on monomials and the Hankel integration-by-parts identity, symbolic.
- M2: frozen exponents = zeros of m for integer a (symbolic in a = 1..4, deg L = 1, 2).
- M3 (the CONSTANT check vs NUM2): (Q3_j/Q6_j) * Gamma(gamma - theta_j/sigma) * sigma^(-theta_j/sigma) is the SAME for all four j, at both points and at every archived E. Prediction: agreement to the data's precision (~55 digits), with one common value 1/kappa per point, E-independent.
  - Also: 1/kappa from the exact E = 0 formulas in closed form.
  - TAMPER: replace m by 1/Gamma(gamma + 1/7 - theta/sigma); it must fail.
- M4: kappa from the Laplace method (an independent route), if it can be done exactly; otherwise NOT DONE.
## General sigma > 0 and non-integer a
- Claim, to be argued: the same proof goes through for ANY Gamma-form solution in the class S (analytic in the sector, growth order < sigma/deg L, decaying on R+, small-x expansion over the dynamical exponents plus the frozen lattice).
- The OBSTRUCTION expected: for non-integer a the Gamma form is not a finite-order ODE. Existence and uniqueness of a subdominant solution in S, i.e. a definition of Q^Gamma, is not provided by ODE theory.
- sigma < 0: the Hankel representation diverges (obstruction).
## Lean (for the lead to relay)
- Candidates: M1's algebraic identities and M2's zero-set identity (finite Gamma-product identities).
- NOT the analytic lemmas (convergence, Laplace asymptotics, uniqueness of the subdominant solution).
