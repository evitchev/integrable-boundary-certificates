# PROOF1 seal: toward proofs identically in t. Opus 5 (2), 2026-10-01. Written before any computation.
Read before sealing: ansatz note items 180, 185, 294, 298, 304, 306, 307; Fable's archived NOTEs vir5d, vir7, vir8 (results/lab/large_n/).
Nothing computed. Rules: exact sympy; solver memory capped (ulimit -v 20 GB); full logs; each statement labelled PROVED (a complete argument, with every algebraic step machine-checked symbolically) or VERIFIED (finite checks).

## (1) Equivalence theorem: the three-term ODE <=> the Gamma-symbol operator
Plan, at the FORMAL MELLIN level:
- With T = s - theta and x^(c/2) L(T) x^(c/2) = x^c L(T - c/2), the ODE's Mellin recurrence is
      P(s - nu) phihat(nu) - L(s - nu + c/2) phihat(nu - c) - E phihat(nu - n) = 0.
- Set phihat = m psihat, with m defined by the first-order difference equation m(mu)/m(mu - sigma) = prod_r (affine in mu), step sigma = c - n = n/a.
  This gives the Gamma-symbol recurrence, with prod_r G_a(T - r) = m(nu)/m(nu - n), a Gamma ratio of "length a" (n = a sigma).
  This is the same mechanism as RED2 B2.
- Machine-check every identity symbolically in a, with the Gamma functions kept symbolic (gammasimp / functional equations), for general P and L (roots symbolic).
- State the branch/contour conditions precisely: the Mellin contour, a strip of analyticity, decay of m psihat, and the absence of pole crossing. These are stated as HYPOTHESES under which the solution spaces correspond, not proved for the record's solutions.
- WKB level (both charges = the same Beta prefactor times the same polynomial):
  - Plan: show that the multiplier m contributes to the Mellin-space log-derivative only through log m(nu) - log m(nu - n), whose asymptotic (Stirling) expansion is exactly the Gamma symbol's correction. The charges of the two forms then coincide as formal series.
  - Prior: formal Mellin equivalence PROVED 85%; WKB-level statement PROVED 35% (otherwise VERIFIED at further k, plus a proof sketch).
## (2) Rank-drop vanishing at all layers
- Plan: show that the loss-m layer of the Gamma form is [w^(K-m)] of Sum_j c_j(nu) phi(w)^(nu-j) P_(m,j)(w), with falling factorials c_j(nu) vanishing for j > nu (nu integer) and deg P_(m,j) <= 2j - m.
  - Then at nu = i integer the generating function is a polynomial of degree <= 2nu - m, and [w^(K-m)] = 0 iff K > 2nu, i.e. 2i <= K - 1, at EVERY loss m.
  - This requires a structural lemma on the all-order Mellin-WKB of the Gamma symbol (z-degree counting).
- Prior: PROVED 25%; otherwise the lemma VERIFIED to high order plus a reduction statement.
## (3) G1 == the record's Ward loss-1 amplitudes, identically in (K, t)
- Plan: write both sides as finite sums over a (the X-power) of Pochhammer products in (a, K, t). Codex's K_d is a coefficient of Appell F1, and G1 is coefficient extraction from powers of phi.
- Reduce the identity, per monomial X^a Y^(K-1-a), to a rational-function identity in (a, K, t) after dividing by a common Pochhammer term. Machine-check it symbolically, for each of the three solutions.
- Inputs: Codex's published loss-1 amplitudes (results/lab/audits/codex_wardsol12*, read-only) and Fable's g1_lib.py definitions (read-only).
- Prior: PROVED for all three solutions 45%.
## Order of work
(1), then (3), then (2). Each part is reported as PROVED, VERIFIED or NOT DONE, with a witness.
