# PROOF2: rank-drop vanishing at losses 2-3, and the formal-WKB equivalence. Opus 5 (2), 2026-10-02. Seal 1102bec9 (sent before computing).
Exact sympy; all runs capped (ulimit -v 16 GB); full logs.
Inputs (read-only copies; INPUTS.SHA256, VIR8COPY.SHA256): Codex's Sol 3 Ward laws, losses 1-3 (loss 4 exists only as stencils, so it is not included), and Fable's VIR8 engines (wkb3.py, mellin_gen.py, mellin_wkb.py, a_equiv.py).

## (1) Rank-drop vanishing of the RECORD's Sol 3 charges at losses 2 and 3: PROVED (loss 1 reproduced as a control)
- Statement: on nu = (2K-1)/n = i (i >= 1 an integer, 2i <= K - 1), every coefficient X^a Y^(K-m-a) (0 <= a <= K - m) of the record's loss-m layer vanishes, for m = 1, 2, 3.
- Object: coef_m(a) = L_m(a, K, t) c_a, where L_m is Codex's Ward-derived law (exact in Q(a, K, t)) and c_a ∝ binom(nu, a) binom(nu, K - a) (BL top at x = -nu; item 304 Z0(a)).
- Proof (p1_factor.log, p1_proof.log, p1_proof_b.log, all exit 0):
  (i) With t = (10nu - 6K + 3)/(2nu - 2K + 1), every denominator factor of L_m is LINEAR in nu:
      (2K - 4nu - 1), (2K - 3nu - 1), (nu - j0) for j0 < m, and (-K + a + nu + j) for j = 1..m.
      - (2K - 4nu - 1) never vanishes at integer nu (parity).
      - (2K - 3nu - 1) = 0 on the locus needs i <= -1.
      - So poles on the locus come only from (nu - j0) (at i = j0) and (-K + a + nu + j) (at a = K - i - j). All are simple.
  (ii) The top's zero order at nu = i is z = [a >= i+1] + [K - a >= i+1] >= 1 everywhere, since K >= 2i + 1. Outside the pole families coef = 0.
  (iii) Pole families with z = 2 against one pole vanish automatically. The z = 1 sub-families, and the coincidence families (both poles, z = 2), are enumerated with the free parameter (K or i) SYMBOLIC. In every one, ord_eps(N) - ord_eps(D) + z = 1 at nu = i + eps, i.e. coef -> 0.
       - The z = 1 sub-families: loss 2 has 5, loss 3 has 11, loss 1 has 1.
       - The coincidence families: loss 2 has 2, loss 3 has 5.
  (iv) The finitely many parameter values where the leading denominator coefficient vanishes (small K, or i = 1, 2) are covered by an EXACT check of every cell with i = 1..6, K = 2i+1..2i+16, 0 <= a <= K - m (1488 / 1392 / 1296 cells for m = 1/2/3). There are 0 nonzero limits; the binomials are exact polynomials in eps there.
- Control: off the locus (2i = K) the layer does NOT vanish (6 / 9 / 10 nonzero coefficients).
- Scope:
  - A THEOREM about the record's (Ward-derived, certified) charges at losses 1-3 on the Sol 3 rank-drop locus.
  - The operator (Gamma-form) side at loss >= 2 equals the record only as VERIFIED (item 298, spins <= 13). So "the operator's polynomial factor vanishes at losses 2-3" is VERIFIED, with PROVED record-side vanishing.
  - Loss >= 4: NOT DONE (no law file). Sols 1 and 2: not attempted.
## (2) Formal-WKB equivalence (item 311(1)): PROVED at the formal level for all P, L (deg P > deg L) and real a, generic exponents; VERIFIED beyond Fable's cases
- Argument:
  (a) PROOF1: phi = m(theta) psi as Mellin multipliers, with m(nu) = prod_r sigma^(nu/sigma)/Gamma(gamma_r - nu/sigma), E-INDEPENDENT. If psi solves the Gamma form, phi solves the three-term ODE (formal Mellin level).
  (b) With psi = e^S and phi = m(theta) e^S, the conjugated log-derivatives satisfy theta log phi - theta log psi = theta log M_full, M_full := e^(-S) m(theta) e^S. EXACTLY a total theta-derivative.
  (c) Stirling: log Gamma(z) = (z - 1/2) log z - z + const + sum_r B_2r /(2r(2r-1) z^(2r-1)). The log z terms sit only at the two leading orders (tau_0^1, tau_0^0). At every WKB order i >= 2 (spin >= 1), (log M_full)_i is a pure power tau_0^(1-i) times a Laurent polynomial in w (the uniformisation of wkb3.py).
      Since tau_0^(1-i) = u^(A0)(1-u)^(B0+1), with Fable's exponents A0 = (1-i)/sigma and B0 + 1 = (i-1)c/(n sigma), the difference of the charge integrands is d/du[u^(A0) (1-u)^(B0+1) R(w)], with R a Laurent polynomial and w = dP - dL u > 0 on [0, 1].
  (d) The continued-Beta functional is defined by analytic continuation in (A0, B0) from Re A0, Re(B0+1) > 0, where the boundary terms vanish. It therefore annihilates such weighted total derivatives IDENTICALLY in the exponents (the same lemma as 38(el)).
  (e) Hence J_ODE(i) = J_Gamma(i) as meromorphic functions of the exponents. Writing J_ODE = cB M_0 + cG M_(-1), J_Gamma = R M_0, and using that M_(-1)/M_0 (a 2F1 ratio) is not rational in the exponents, we get cG = 0 and cB = R. This is Fable's "same cB x B(A0, B0+1)".
- Gaps (stated):
  - The conjugation identity for the NON-polynomial multiplier is used FORMALLY (asymptotic Stirling series). Convergence is not addressed.
  - 'Not rational' for M_(-1)/M_0 is used generically.
  - At Beta-pole cells (A0 a non-positive integer) the statement holds as an identity of limits, not of values.
- VERIFIED extension (Fable's engines, copied; p2_aequiv_*.log, p2_general_case*.log, all exit 0):
  - Sol 3 at k = 3/2, 5/3, -6/5, 7/2 and Sol 2 at k = 3/2, 5/2, -6/5, orders <= 8: AGREE. A k = 2 regression agrees.
  - GENERAL operators P = (T^2 - l1^2) prod (T - r), L = (T - l0) prod (T - q), with (dP, dL) = (3, 1), (3, 2), (4, 2), non-integer k, rational extra roots: AGREE at all non-degenerate orders.
  - Tampers (x^c L(T) ordering): DISAGREE (they fire).
  - At k = 3/2 (sigma = 3) the cells i = 4, 7 have A0 = -1, -2 (Beta poles). There the engines' values differ, but at k = 3/2 + 1e-6 ALL orders agree (p2_general_case1_near.log). The pole-cell discrepancy is an evaluation-at-the-pole artefact, consistent with (e) as an identity of meromorphic functions.
## (3) H1-H3 for Sol 1: NOT DONE (not attempted; it was conditional on time).
## Disclosures
- p2_general v0 accepted an invalid input (deg P = deg L; wkb3 requires dP > dL) and a degenerate fibre. The invalid case was replaced, and the degenerate cells are now labelled and excluded from the AGREE count, but reported.
- The shell quoting of the first p2_general loop failed before any run.
- p1_proof.py was refactored (a main() guard) so the supplement could import it. It was rerun; the stored logs are from the final versions.
