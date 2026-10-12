# SEAL VIR2-R -- independent replication of Fable's VIR2 (two-boson screening system for the cylindrical charges)

cc (Opus 5.5 seat), 2026-10-01. Lead's commission. Sealed BEFORE any computation; the hash goes to the lead first.
- Fable's VIR2 (~/fable-work/vir2, NOTE a6005095, seal 9ac0c556) is NOT opened. I use only the claim as stated in the brief.
- NO code is shared with Fable or with mixed_engine. mixed_engine (w8 clone, read-only) is used ONLY to produce the record's charge representatives (inputs).

## Claim under test (from the brief)
- **Setup:** two bosons X, phi; dX dX ~ dphi dphi ~ +1/(z-w)^2. T_phi = (1/2)(dphi)^2 + rho d^2 phi, with p^2 = 2(t-1)/(t+1), rho = p/2 - 1/p, c = 1 - 12 rho^2 = N(t).
- **Screening pairs** V_+- = exp(+-a X + beta phi), conserved against I_3, I_5, I_7:

| solution | beta | a^2 |
|---|---|---|
| Sol 1 | -p/2 | 1 - p^2/4 |
| Sol 2 | -p | 1 - p^2 |
| Sol 3 (pair A) | 1/p | 1/p^2 |
| Sol 3 (pair B) | -p | p^2 |

- **Determination:** the commutant of {exp(p phi), V_+, V_-} among ALL two-boson densities, mod d, is 1-dimensional at weights 4, 6 and equals I_3, I_5. At weight 8 it equals I_7.

## Method (own code; stated first)
- **(M1) Inputs:** the record's I_3 = cyl_P4 and I_5, I_7 (genuine kernels at weights 6, 8), at curve_point(t) = (N, s), from mixed_engine.
- **(M2) O(N) pair-language Wick engine with FORMAL N** (new, `onwick.py`):
  - builds on my WINF1/VIR1-R Wick code (`wick.py`);
  - each Y pair (a, b) gets a fresh dummy label;
  - Y-Y contractions unify labels (X only with X);
  - remaining letters are regrouped by label class: 2 letters give a new pair; 0 letters (a closed cycle) give a factor N.
  - CONTROL (e): it must agree with my explicit-component Wick code at integer N = 3, 4 (mod-p evaluation) on :T_Y T_Y: and :T_Y :T_Y T_Y::.
- **(M3) Decomposition:** at each fibre (rational N, s), solve EXACTLY over Q in pair-monomial coordinates (independent for non-integer N):
  - I_k = sum_j x_j (X monomial × Vir word in T_Y) + d(V_(W-1)).
  - Vir words are nested normal-ordered products of d^n T_Y (one per partition into parts >= 2).
  - Report the residual (must be 0, i.e. VIR1 holds at that fibre) and the nullspace dimension (the ambiguity of x).
- **(M4) Translation:** I_k^(X,phi) = sum_j x_j (X monomial × the SAME Vir word in T_phi). T_phi uses rho as above, built with my explicit Wick code (letters x, f).
- **(M5) Conservation test** for a density G and V = exp(alpha_x X + alpha_f phi):
  - the OPE G(z)V(w) by contractions (d^p a_c(z) with exp: alpha_c (-1)^p p! (z-w)^-(p+1)), with G's remaining letters Taylor-expanded at w;
  - R = the residue (letter weight W-1);
  - conserved iff R = dQ + (alpha_x a_x + alpha_f a_f) Q for some Q of weight W-2 (an exact linear solve).
  - Arithmetic is in Q(sqrt(d)) (a^2 = d rational, a possibly irrational), with an own exact field class.
- **(M6) Commutant** at weights 4, 6, 8: G = general two-boson polynomial of weight W (all monomials in letters x, f).
  - Joint linear conditions from exp(p phi), V_+, V_-, with unknown Q's.
  - dim(commutant mod d) = dim(solution projection to G) - dim(d(weight W-1)).
  - Compare with I_k^(X,phi) (membership, and equality as 1-dimensional spaces).

## Fibres (rational p; t = (2 + p^2)/(2 - p^2))
- **Sol 1:** p = 8/5 (t = -57/7), 6/5 (t = 43/7), 40/29 (t = 1641/41).
- **Sol 2:** p = 4/5 (t = 33/17), 3/5 (t = 59/41), 8/17 (t = 321/257).
- **Sol 3:** p = 1/2 (t = 9/7), 2/3 (t = 11/7), 3/2 (t = -17).
- The commutant (M6) runs at the FIRST fibre of each solution, weights 4, 6, 8. I_7 (weight 8) is tested at all fibres if time allows; otherwise at the first, stated.

## SEALED PREDICTION (= the claim)
- (P1) Every listed pair is conserved against I_3, I_5, I_7 at every fibre of its solution. VIR1 holds at every fibre: decomposition residual 0.
- (P2) The commutant of {exp(p phi), V_+, V_-} mod d has dimension 1 at weights 4, 6, 8 (Sol 3: with pair A), spanned by I_3, I_5, I_7 respectively.

## Controls
- **(a)** a^2 -> a^2 + 1/10 must FAIL conservation (the I_5 test).
- **(b)** rho -> -rho in the translation must FAIL (V_+- against I_5).
- **(c)** the Virasoro screenings exp(p phi) and exp(-(2/p) phi) must PASS for every translated charge.
- **(d)** generic (a^2, beta) = (2/7, 3/11) must FAIL.
- **(e)** the formal-N engine must agree with the explicit-component engine at N = 3, 4.
- **(f)** sanity: a density known to commute (d of anything) must pass every test.
