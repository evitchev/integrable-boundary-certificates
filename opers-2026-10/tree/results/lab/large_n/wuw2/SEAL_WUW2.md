# SEAL WUW2 -- i = 2 at t = 7/5: is the antisymmetric spin-9 class W_U^(5) W_W^(5)? (lead's commission)

cc (Opus 5.5 seat), 2026-10-01. Sealed BEFORE any computation. The hash is sent only after the command returns.
- Read before the seal: the record's notes on cyl_t2_primary (mode-action note M2; post-publication programme 2026-09-01), and the head of code/cyl_t2_primary.py (docstring).
- No Fable code was opened.

## Setting (t = 7/5, i = 2)
- p^2 = 1/3, rho = p/2 - 1/p = -5/(2 sqrt3), N = -24. U, W = (phi ± X)/sqrt2.
- Pair A: e^(sqrt6 U), e^(sqrt6 W). Pair B: e^(-s U), e^(-s W), with s = sqrt(2/3). Virasoro: e^((s/2)(U+W)), e^(-3s(U+W)).
- **Rescaled bosons V = U/s, Z = W/s:** d V d V ~ kappa/(z-w)^2 with kappa = 3/2. Then every exponent is RATIONAL:
  - pair A e^(2V), e^(2Z); pair B e^(-(2/3)V), e^(-(2/3)Z); Virasoro e^((V+Z)/3), e^(-2(V+Z)).
- Engine: WUW1's code with every contraction multiplied by kappa (OPE pairs, exponential residues, field-level screenings, cylinder pairings).
  - Engine control: at kappa = 1 it must reproduce WUW1's I_3, K_6 and W_U W_W results.
- T_U = (1/3)(dV)^2 + (2/3) d^2 V (c_U = 1 - 12 s^2 = -7 = c_(1,3)); likewise T_W.
- The record's T_total = (1/3)((dV)^2 + (dZ)^2) - (5/6)(d^2V + d^2Z) = T_U + T_W - (3/2)(d^2V + d^2Z).

## Construction (explicit)
- **W_U^(5) := Res_(z->w) e^(2V)(z) e^(-2V)(w) = S_5(2 d^k V/k!)**, the coefficient of x^5 in exp(2 sum_k d^kV x^k/k!). Here (z-w)^(-a^2 kappa) = (z-w)^(-6), and e^(-2V) has T_U-weight 5. Likewise W_W^(5).
- Checks (sealed as controls):
  - T_U-primary of weight 5;
  - annihilated at field level by the short screening e^(-(2/3)V);
  - NOT in the weight-5 span of the Virasoro vacuum module {d^3 T_U, T_U dT_U}.
- G_10 := W_U^(5) W_W^(5) (independent bosons: plain product).

## Tests and SEALED PREDICTIONS
- **I_3** = the weight-4 commutant of {e^((V+Z)/3), e^(2V), e^(2Z)}, dim 1. GATE: its eigenvalue is proportional to cyl(3,4) at t = 7/5.
- **P1:** eigen(W_U^(5)) ∝ P_V (P_V^2 - 1/4)(P_V^2 - 1), with P_V = (sqrt3/2)(P_X + pi): roots j/2, j = -2..2, which is E_2's root string in V-momentum. Prior 70%.
- **VACUITY STATEMENT (in advance):**
  - One-boson weight-5 fields form a 7-dimensional space (partitions of 5). Their vacuum eigenvalues are odd quintics a P^5 + b P^3 + c P, a 3-dimensional shape space. So matching E_2's specific shape is NOT automatic at i = 2, but it is not unique to W^(5) either.
  - The discriminating test is the bilinear one. Of the 7 × 7 = 49 bilinears sum M_ij f_i(V) g_j(Z), I compute how many lie in K_10 + d.
- **P2:** that space is ONE-dimensional and equals W_U^(5) ⊗ W_W^(5). Prior 55%.
- **P3:** G_10 lies in K_10 = ker ad I_3 at weight 10 (mod d), and eigen(G_10) = const × E_2 with E_2 = prod_(j=-2..2)(P_X + pi - j p)(P_X - pi - j p), p^2 = 1/3. Prior 60%.
- The full dim K_10 (record and VIR5d: 2) is computed if feasible, and otherwise stated as not computed. GATE if computed: the certified vev_sol3_w10 at 7/5 lies in span eigen(K_10).
- **P4 (reconciliation with cyl_t2_primary, "0 primary classes at 7/5 under T_total, T_X, T_Y"):**
  - G_10 IS primary under T_U + T_W (by construction).
  - It is NOT primary under T_total, which differs by -(3/2)(d^2V + d^2Z). Prediction: no representative G_10 + dH (H weight 9) is T_total-primary.
  - So the two statements do not conflict: primariness is frame-dependent, and the record tested the X/Y and total frames, not the U/W Coulomb-gas frame.
  - If G_10 + dH IS T_total-primary for some H, that is a CONFLICT with the record and will be reported as such.

## Controls
- Antisymmetry of eigen(G_10) under P_X <-> pi.
- Symmetric W_U^2 + W_W^2: not antisymmetric; K_10 membership reported.
- TAMPERS that change SHAPE, not scale:
  - W_U^(5) + (dV)^2 d^3V/10 (and the same for W): must leave K_10 and break P1/P3;
  - the exponent tamper a = 2 -> 2 + 1/10, i.e. S_5 with the shifted a (shape changes at weight 5): must break P1 and leave K_10;
  - scoring tamper on E_2 (p^2 + 1/1000).
