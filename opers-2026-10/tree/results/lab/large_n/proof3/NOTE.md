# PROOF3: Theorem 1 at the spectral level. Opus 5 (2), 2026-10-07.
- Seal: SEAL_PROOF3.md, sha256 2803c8e07551b9a78a465cf24c3dc3ae20cd3d1e9eb2f6c5bfe2918cec986f96. It was sent to the lead before any computation.
- Inputs: copies of NUM2's t1.json, ops.py and ode.py. Their hashes are in INPUTS.SHA256.
- All runs were on the head node, under ulimit -v 8 GB, with full stdout logs.

## Result
- THEOREM S is PROVED for Sol 3 at t = 2 (k = a = 2, fourth order vs sixth order). The proof uses the Hankel route S1-S4 sealed in advance.
  - The algebraic steps are machine-checked symbolically (M1, M2). The analytic steps are argued below.
  - The constant is exact. For all j and all E,
        Q3_j(E) / Q6_j(E_Gamma = -E) = m(theta_j) / kappa,   with
        m(nu) = 3^(nu/3) / Gamma(1/3 - nu/3)   and   kappa = 3^(1/3) / sqrt(3 pi) = 3^(-1/6) pi^(-1/2).
  - Cross-check against NUM2's archived Q's (M3):
    - j-spread ≤ 1.9e-56;
    - E-drift ≤ 1.2e-56 over E = -20 ... 500;
    - this holds at both NUM2 points;
    - the tamper (gamma + 1/7) spreads by 0.36 and 0.78.
  - kappa was derived independently by Laplace's method (M4). It agrees with M3's number to 28+ digits.
  - HONESTY NOTE: the Laplace derivation was done after M3 had produced the number. It has no free parameter.
- General sigma > 0:
  - The same proof holds verbatim for every integer a ≥ 1 with deg L = 1 and deg P ≥ 2.
  - Non-integer a is OBSTRUCTED only at the DEFINITION of Q^Gamma (details below).
- sigma < 0: the Hankel route is OBSTRUCTED. A candidate path is given.

## Proof (t = 2; general integer a in brackets)
- Objects as in the seal:
  - L3 = P(s - theta) - x^(9/2)(s - theta)x^(9/2) - E x^6, with s = 5/2;
  - L6 = P(s - theta)(theta - 1)(theta - 4) - x^9 - E x^6.
- S1, the Hankel kernel:
  - Define K psi(x) = (1/2 pi i) Int_H e^t t^(-1/3) psi(x (3t)^(1/3)) dt.
  - On monomials, K x^nu = m(nu) x^nu, by Hankel's formula for 1/Gamma.
    - M1a checks this by quadrature at three nu, including a complex one: relative error ≤ 2e-41.
  - Convergence condition:
    - A Gamma-form solution has growth order rho_Gamma = c/(deg P + a deg L). Here this is 9/6 = 3/2.
    - The Hankel integral converges iff rho_Gamma < sigma.
    - With c = n(1+a)/a and sigma = n/a, this is ⇔ 1 + a < deg P + a deg L. For deg L = 1 that is ⇔ deg P ≥ 2, so it is automatic here.
  - CORRECTION TO THE SEAL: the seal wrote this as "rho = c/n < sigma/deg L ⇔ deg L < deg P". The deg-L = 1 statement above is the one proved.
    - For deg L ≥ 2 the iterated kernel enlarges the growth after each factor. Convergence is then a CONDITION to check, not automatic.
  - Analyticity: ODE solutions continue to the universal cover of C*, so H3 holds.
- S2, intertwining L3 K = x^6 K x^(-6) L6:
  - On monomials this reduces to two Gamma recurrences:
    - m(mu)/m(mu - 6) = (mu - 1)(mu - 4);
    - m(mu + 3)/m(mu) = -(2 + mu).
  - M1c checks it symbolically in (l0, l1, E, mu). Its tamper fails as it should.
  - For general psi, theta under K becomes the integration-by-parts relation K_(g-1) psi = K_g[(g - 1 - theta/sigma) psi]. The boundary terms vanish on H because e^t decays on both rays.
    - M1b checks it symbolically on monomials, and by quadrature on the non-monomial psi = x^0.4 exp(-x^2) (relative error 1.3e-41).
- S3, Frobenius mapping:
  - Each Frobenius series has exponents in theta_j + 3N (lattice of n = 6 and c = 9). It is entire of order 3/2 times x^theta_j, so termwise application of K is justified by dominated convergence (same bound as S1).
  - Hence K chi6_j = m(theta_j) chi3_j for the dynamical exponents. The image solves L3 (S2) and has a unique Frobenius normalisation (H1).
  - K kills the frozen components (theta = 1, 4 + lattice): every exponent there is a zero of m.
    - M2 checks this symbolically for a = 1..4, per root factor of m. deg L = 2 follows because m is a product over roots; it was not run separately.
    - The tamper (frozen + sigma/2) fails, as it should, for all a.
- S4, subdominance and kappa:
  - The decaying branch of L6 is psi ~ y^beta6 exp(-(2/3) y^(3/2)). It is recessive for |arg y| < pi/3, and its asymptotic expansion is valid for |arg y| < pi.
  - On H, y = x (3t)^(1/3) has |arg y| ≤ pi/3. So the expansion holds on the whole contour.
  - Laplace's method:
    - f(t) = t - A sqrt t, with A = (2/sqrt3) x^(3/2);
    - saddle t* = x^3/3, with f(t*) = -x^3/3;
    - y* = x^2, and Gaussian factor A/(2 sqrt pi).
  - Result: K psi = kappa x^(2 beta6 + 1/2) e^(-x^3/3) (1 + O(x^-3)), with kappa = 3^(1/3)/sqrt(3 pi).
  - This is L3's subdominant branch:
    - beta3 = 2 beta6 + 1/2 = -2, checked with NUM2's own asymptotic solver at E = 0 and E = 37 (M4);
    - L3's other solutions are either exponentially growing or power-law, so the subdominant solution is unique up to scale.
  - kappa is E-independent:
    - E enters both normalised log-expansions only at negative powers;
    - beta3 and beta6 do not depend on E (M4);
    - so no constant can be generated.
- Conclusion: K y6_sub = kappa y3_sub. Expanding both sides in Frobenius bases gives Q3_j = m(theta_j) Q6_j / kappa.
- Pencil remark C (Gauss multiplication) is consistent:
  - at E = 0, M3 evaluates NUM2's exact Meijer-G formulas and gets 1/kappa = sqrt(3 pi)/3^(1/3) = 2.12860533089645112119295775073;
  - the j-spread is ≤ 1.5e-61.

## What is proved vs assumed
- PROVED by machine (symbolic or exact, with failing tampers):
  - the multiplier identities;
  - intertwining on monomials;
  - the frozen-zero lattice (a = 1..4);
  - the constant 1/kappa = sqrt(3 pi) 3^(-1/3), against NUM2 to ~1e-56.
- PROVED on paper (standard analysis, not machine-checked):
  - Hankel convergence;
  - termwise application on Frobenius series;
  - the Laplace asymptotics with the sector bookkeeping above;
  - uniqueness of the subdominant solution.
- ASSUMED: H1 (no resonances; dynamical exponents off the frozen lattice). Both NUM2 points satisfy it.

## General sigma > 0, and the obstruction
- Integer a, deg L = 1, deg P ≥ 2: S1-S4 go through unchanged.
  - kappa is a product of sigma-, sqrt(pi)- and power factors from one Gaussian per Hankel factor.
  - Not computed beyond t = 2: NOT DONE, with this formula sketch as a witness.
- deg L ≥ 2: the iterated-Hankel convergence is a condition to verify case by case (see the S1 correction).
- Non-integer a:
  - The Gamma form is not a finite-order ODE. Its frozen exponents form an infinite set, namely exactly the zeros of m.
  - Steps S1-S4 apply verbatim to any solution in class S (analytic in the sector, growth < sigma, small-x expansion over dynamical exponents plus frozen lattice, decaying on R+).
  - What is missing is EXISTENCE AND UNIQUENESS of the subdominant class-S solution, i.e. a definition of Q^Gamma. ODE theory does not supply it.
  - This is the obstruction. The ODE-side Q's need no such input.
- sigma < 0: (sigma t)^(1/sigma) sends t → ∞ to y → 0, and the Hankel integral diverges on Gamma-form solutions.
  - Candidate path (unchecked): realise 1/m by the Euler integral Gamma(z) = Int_0^∞ e^-t t^(z-1) dt. That gives the inverse map, from ODE solutions to Gamma-form solutions.

## Lean candidates (for the lead to relay to cognitive)
- (i) The Gamma-recurrence identities behind S2: m(mu)/m(mu - n) = prod G_a, the c-shift identity, and PROOF1's I1-I3, as finite Gamma-product identities.
- (ii) The frozen-zero lattice statement (M2): gamma_r - nu/sigma ∈ Z≤0 on the frozen lattice, for all integer a. This is pure arithmetic.
- (iii) Gauss's multiplication formula at m = 3. Mathlib may already have it; check.
- (iv) Optionally, the Laplace constant as an identity of Gaussian integrals.
- NOT candidates (heavy analysis): Hankel convergence, termwise interchange, Laplace asymptotics with sectors, uniqueness of the subdominant solution.

## Files
- SEAL_PROOF3.md, INPUTS.SHA256, t1.json, ops.py, ode.py (inputs);
- m3_constant.py + m3_constant.log (M3 PASS);
- m124_checks.py + m124_checks.log (M1/M2/M4 PASS, tampers fail);
- m124_checks_v0_branchbug.log: the first run, which FAILED M1a/M1b because of a principal-branch error in my quadrature (u e^(-i pi) landed on the wrong sheet). Fixed by parametrising by log t. Kept for the record.
- Shared memory: no changes in PROOF3.
