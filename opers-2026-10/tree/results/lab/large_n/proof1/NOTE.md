# PROOF1: proofs identically in t. Opus 5 (2), 2026-10-01. Seal b0cbfcf7 (sent before computing).
Exact sympy; no solver runs; full logs. Inputs (read-only copies, INPUTS.SHA256): Codex's projection_loss1.json for Sols 1-3 (audit codex_wardsol12b_27403fc) and Fable's g1_lib.py (vir7).
Labels: PROVED = a complete argument whose every algebraic step is machine-checked symbolically; CONDITIONAL = proved under stated hypotheses; VERIFIED = finite checks; NOT DONE.

## (1) Equivalence theorem: PROVED at the formal Mellin level (general P, L, real a); CONDITIONAL at the solution / spectral level; the WKB-level "same B x polynomial" NOT proved in general
Setting: theta = x d/dx, T = s - theta, L(T) = prod_r (T - r) monic, M = 1/a, sigma = n/a, c = n + sigma = n(1 + M), and
    G_a(T) = sigma^a Gamma(T/sigma + (1+a)/2) / Gamma(T/sigma + (1-a)/2)      (item 298/306's sigma^a R_a(T/sigma); a polynomial of degree a for integer a).
THEOREM 1 (formal).
- Let phi(x) = int phihat(nu) x^nu dnu and psi(x) = int psihat(nu) x^nu dnu, with
      phihat(nu) = m(nu) psihat(nu),   m(nu) = prod_r sigma^(nu/sigma) / Gamma(gamma_r - nu/sigma),   gamma_r = (s - r)/sigma + (1 - a)/2.
- Then the Mellin recurrence of  [P(T) - x^(c/2) L(T) x^(c/2) - E x^n] phi = 0  is, after division by -m(nu - n), exactly the Mellin recurrence of
      P(T) prod_r G_a(T - r) psi = x^n (x^(n/a) + E) psi.
- The ODE's -E x^n corresponds to the Gamma form's +E: this fixes the sign convention.
PROOF (p1_equivalence.py, exit 0; every identity symbolic in a, sigma, s, r, nu, Gamma kept symbolic):
- (I0) T x^(c/2) = x^(c/2)(T - c/2), so x^(c/2) L(T) x^(c/2) x^nu = L(s - nu - c/2) x^(nu + c), and the ODE recurrence is
       P(s-nu) phihat(nu) - L(s-nu+c/2) phihat(nu-c) - E phihat(nu-n) = 0.
- (I1) m(mu)/m(mu - sigma) = prod_r (sigma gamma_r - mu): Gamma(z+1) = z Gamma(z), applied per root.
- (I2) m(nu)/m(nu - n) = prod_r G_a(s - nu - r), using n = a sigma. For real a this is a Gamma identity, not a finite product. It needs the 1/Gamma(gamma_r - nu/sigma) form; the alternative Gamma(nu/sigma + ...) form differs by a sigma-periodic factor, by reflection.
- (I3) L(s - nu + c/2) m(nu - c)/m(nu - n) = 1.
- (I4) The ODE recurrence divided by -m(nu - n) equals P G psihat(nu) - psihat(nu - c) - E psihat(nu - n) for symbolic P; deg L = 1, 2, 3 checked.
- For general deg L the identities factor root by root: each is a product of the single-root identity.
- Negative controls fire: the non-symmetric ordering x^c L(T) (VIR8's tamper) breaks (I3); gamma_r + 1/3 breaks (I2).
REMARKS.
- m is ENTIRE in nu, since 1/Gamma is entire, so phihat = m psihat introduces no poles. The inverse psihat = phihat/m has poles at the zeros of m, nu = sigma(gamma_r + j), j = 0, 1, 2, ...
- In x-space m is a Mellin multiplier: an Erdelyi-Kober / Meijer-G type fractional operator, non-local. This is the same mechanism as RED2 B2.
COROLLARY (CONDITIONAL: spectral / solution level).
- Hypotheses:
  - (H1) psihat is analytic in a vertical strip S containing the contour and its shifts by n and c;
  - (H2) m(nu) psihat(nu) x^nu decays on vertical lines fast enough to shift contours (1/Gamma grows like exp(pi |Im nu|/(2 sigma)) per root, so psihat must decay faster);
  - (H3) no pole of psihat lies between the contour and its shifts.
- Under H1-H3 the map psi -> phi sends solutions of the Gamma form to solutions of the ODE. The multiplier does not depend on E.
- So the Frobenius exponents and connection data are related by E-independent constants. The spectral determinants then coincide up to an E-independent factor, provided the boundary-value solutions (regular at x = 0, subdominant at infinity) map to each other.
- Hence any charges defined as asymptotic coefficients of log D(E) COINCIDE. H1-H3 and the boundary-condition matching are NOT proved for the record's solutions.
NOT PROVED in general: Fable's formal-WKB statement "both forms give cB x B(A0, B0+1) with the same polynomial" for all P, L and real a.
- It is VERIFIED by Fable (item 307) at k = 2 and 10/3, orders 2..8, and in the second-order case.
- The present theorem gives the mechanism: the forms differ by log m(nu) - log m(nu - n), whose Stirling expansion is exactly the Gamma symbol's correction. But it is not a formal-WKB proof.

## (2) Rank-drop vanishing at all layers: NOT DONE
- What is available: Fable's proof for the top and loss 1 (item 304; Z0(a), Z0(b)). No all-order closed form of the Gamma-symbol WKB exists, so no all-loss proof.
- Reduction (argued, not machine-proved). In the three-term ODE form the charge is cB x B(A0, B0+1). Normalise it by the Pochhammer contour, i.e. multiply by (1 - e^(2 pi i A0))(1 - e^(2 pi i B0)): the result is entire in A0. At A0 = -m (m a non-negative integer) it equals a nonzero constant times Res_(u=0)[u^(A0-1)(1-u)^B0 w W(w)], and it is proportional to the limit of cB.
  - So the all-loss rank-drop statement is EQUIVALENT to the vanishing of that residue at every WKB order, for 2i <= K - 1.
  - A degree lemma would suffice: the loss-m generating polynomial at integer nu has z-degree <= 4nu - 2m. It holds at the top (plain degree) and at loss 1 (Fable), and is not proved beyond.
- Witness of the gap: no all-order formula is available, and I did not derive one in this stage.

## (3) G1 == the record's Ward-derived loss-1 layer, IDENTICALLY in (a, K, t), for ALL THREE solutions: PROVED (p3_g1_ward.py, exit 0)
- Object: R(a, K, t) = [X^a Y^(K-1-a)] / [X^a Y^(K-a)] of the record-normalised charge. Codex's projection_over_c_a, derived from the Ward grades (items 180/185), is an explicit rational function in Q(a, K, t).
- G1 side: every G1 term is [z^N] of prod_e (1 - e z)^(nu alpha_e - d_e) times a monomial in (z, l).
  - A paired family (1 - l z)^(beta - d+)(1 + l z)^(beta - d-) = (1 - l^2 z^2)^(beta - D) times a FINITE polynomial. A single family (Sol 1's l0) is one binomial.
  - So the coefficient of l0^(2i) l1^(2b) is a finite sum of binomials whose arguments differ from the top's by INTEGERS. Each ratio binom(x+dx, m+dm)/binom(x, m) = rf(x+1, dx)/rf(m+1, dm)/rf(x-m+1, dx-dm) is rational in (i, K, t).
  - Negative lower indices give 0 automatically (the rising factorial carries the zero), so the identity holds at every 0 <= a <= K - 1, not only generically.
- Record conversion (Fable's to_record, as used in item 306): l0^2 = -CX X, l1^2 = CY (rho^2 - Y), rho^2 = 2/(t^2 - 1). The pi^2 -> rho^2 - Y shift adds -(K - a) rho^2 to R.
- Operator data:
  - Sol 1: phi = (1 - l0 z)^a_s (1 - l1^2 z^2); one string of length a_s = (t-1)/2 at l0; n = (t+3)/2; CY = n^2(M+1)/2.
  - Sol 2: phi = (1 - l1^2 z^2)(1 - l0^2 z^2)^k; strings of length k at +-l0; n = 2k+3; CY = 2n(1+M)/a_X.
  - Sol 3: phi = (1 - l0^2 z^2)(1 - l1^2 z^2); a string of length k at 0; n = k+4; CY = 2n(1+M)/a_X.
  - In all cases k = 2(t-1)/(3-t) and M = 1/(string length).
- RESULT: sympy.cancel(R_G1 - R_Codex) = 0 in Q(a, K, t) for Sols 1, 2 and 3.
- Controls (p3_tamper.py, exit 0): CY x 11/10, and string length + 1/10, each break the identity for all three solutions (6/6 fire).
- Scope:
  - The inputs are TRUSTED: Codex's Ward-derived laws and the operator dictionaries of items 298/303/306.
  - The G1 formula itself (Fable's first-quantum-order formula for a Gamma symbol) is used as defined. Its own derivation from the operator is Fable's (item 306) and is not re-proved here.
  - What is PROVED: the operator side's G1 layer equals the record's Ward loss-1 layer identically in (a, K, t). Item 306's check at K <= 40 thereby becomes an identity.
## Errors and disclosures
- p1_equivalence v0: sympy's gammasimp does not reduce Gamma(z+1)/Gamma(z) for composite z, so the identities first printed False. Fixed with expand_func; the identities themselves were unchanged.
- p3 v0: Python's (-1)**negative_int made floats, so Sols 1-2 first showed float residues. Fixed with integer parity; Sol 3 was exact from the start.
- p3_tamper: a dead third-control stub was removed before the stored run.
- No part of (2) is claimed beyond the reduction.
