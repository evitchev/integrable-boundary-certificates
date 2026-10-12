# SEAL WUW1 -- density-level test of the VIR5d reading E_1 = <W_U W_W> at t = 5/3 (lead's commission)

cc (Opus 5.5 seat), 2026-10-01. Sealed BEFORE any computation; hash to the lead first.
- Read before the seal: ~/fable-work/vir5d/NOTE.md only (the statement, as allowed). No Fable code was opened.
- Also read (own workspace): my stageVIR2R/VIR3R code, for conventions.

## Setting (t = 5/3, i = 1)
- p^2 = 1/2, rho = p/2 - 1/p = -3/(2 sqrt2), N = -25/2.
- Bosons: X and phi with d d ~ +1/(z-w)^2, phi phi ~ +log (my VIR2-R conventions).
- U = (phi + X)/sqrt2, W = (phi - X)/sqrt2 (orthonormal). In (U, W) the Sol 3 screenings at t = 5/3 have RATIONAL exponents:
  - pair A (beta = 1/p, a = 1/p): e^(2U), e^(2W);
  - pair B (beta = -p, a = p): e^(-U), e^(-W);
  - Virasoro: e^(p phi) = e^((U+W)/2), and e^(-(2/p) phi) = e^(-2(U+W)).
- Screening tests (residue mod D) do not involve the background charge.
- With T_U := (1/2)(dU)^2 + (1/2) d^2 U, both e^(2U) and e^(-U) have weight 1: the c_U = 1 - 12/4 = -2 Coulomb gas (alpha_+ = 2, alpha_- = -1), likewise W. This is the VIR2 Vir_U x Vir_W.
  - Note: T_U + T_W differs from the record's T by a multiple of d^2 phi. That does not affect zero modes.

## (1) The weight-3 objects (no printed source invoked; everything computed and verified)
- **W_U := Res_(z -> w) e^(2U)(z) e^(-2U)(w)**: the long screening applied to e^(-2U), whose weight under T_U is 2 + 1 = 3.
  - Explicitly: the coefficient of (z-w)^3 in exp(2 sum_(k>=1) d^k U(w) (z-w)^k/k!), a weight-3 polynomial in dU, d^2U, d^3U. Likewise W_W.
- Verified, not assumed:
  - (i) W_U is T_U-primary of weight 3 (OPE poles 5..3);
  - (ii) at FIELD level it is annihilated by the short screening (Res e^(-U)(z) W_U(w) = 0);
  - (iii) it is not a T_U-descendant (not in span{dT_U}).
- This is the Coulomb-gas weight-3 generator of the c = -2 singlet algebra in this realisation; the name is a description of verified properties.
- Its relation to the symplectic-fermion / triplet currents is NOT claimed or computed.

## (2) Tests
- **I_3:** the weight-4 commutant (mod d) of {e^((U+W)/2), e^(2U), e^(2W)} (VIR2's screening set).
  - Expected dim 1. GATE: its cylinder vacuum eigenvalue equals the record's cyl(3,4) at t = 5/3 (VIR3-R recipe and dictionary).
- **K_6:** the kernel of ad I_3 at weight 6 mod d, i.e. {G : Res_1 I_3(z) G(w) in d(...)}.
  - Expected dim 2 (VIR5d). GATES: it contains the weight-6 screening commutant, whose eigenvalue is vev_sol3_w6 at 5/3; and the certified Sol 2 spin-5 charge (vev_sol2_w6 at 5/3) lies in the span of the K_6 eigenvalues.
- **G_6 := W_U W_W** (U, W independent, so the normal-ordered product is the plain product). Sealed predictions:
  - P1: eigen(W_U) ∝ P_U (P_U^2 - 1/4), with P_U = (P_X + pi)/sqrt2 (equally spaced roots, spacing 1/2 = alpha_-/2). Prior 70%.
  - P2: G_6 lies in K_6 (mod d) and is NOT in span{I_5} + d. Prior 55%.
  - P3: eigen(G_6) = const × E_1, E_1 = prod_(j=-1..1)(P_X + pi - j p)(P_X - pi - j p). Prior 60% (implied by P1, given factorisation).
  - P4: G_6 commutes (mod d) with pair B {e^(-U), e^(-W)} and the Virasoro screening? Reported, not sealed.
- If P2 fails: report the residual. No rescue.

## (3) Controls
- Antisymmetry under P_X <-> pi (P_U -> P_U, P_W -> -P_W): eigen(G_6) antisymmetric. Symmetric combination W_U^2 + W_W^2: its eigenvalue must be symmetric, i.e. NOT antisymmetric, and whether it lies in K_6 is reported.
- TAMPER 1: G_6 + (dU)^6/100 must leave K_6 and break P3.
- TAMPER 2: W_U built with exponent 2 + 1/10 (Res e^(aU) e^(-aU), a = 21/10) must break P1 and P3.
- Scoring tamper on E_1 (p^2 + 1/1000) must break P3.
- Engine control: my wick.py and residue code from VIR2-R/VIR3-R (validated). The cylinder recipe is re-gated by the I_3 and I_5 gates above.
