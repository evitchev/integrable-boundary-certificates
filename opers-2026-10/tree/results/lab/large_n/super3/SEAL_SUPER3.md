# SUPER3 seal. Opus 5 (2), 2026-10-01. Step A (reasoning) and the Step B predictions, before any WKB computation.
Read before sealing: the item-252 text (ansatz note l.21114-21130) and my SECTOR1 NOTE; ansatz note 38(ah)(x) and 38(el) on kappa-odd sources;
Suzuki quant-ph/0003066 text layer (QLIT2 copy) eqs (1)-(3) and sec. 2; Tanabe 2604.14899 and Ito-Zhu 2206.08024 passages that cite Babenko-Smirnov;
the header and coefficients of results/lab/vev/vev_sol1_w6.json (format only, to fix the normalisation). No WKB was computed.

## Step A: verdict (ii). Suzuki's ODE lies OUTSIDE S1's scope.
A1. S1's class (item 252, SECTOR1 NOTE): scalar opers psi'' = [kappa^2 Lambda + X v_X + Y v_Y + V_0] psi.
    The momenta enter ONLY at kappa^0 and LINEARLY in the record's (X, Y).
    Then the top at spin 2k-1 is binom(1/2,k) int (X v_X + Y v_Y)^k Lambda^(1/2-k) dz.
    This is ONE moment family of the single measure v_Y^k Lambda^(-s/2) dz, and moment uniqueness gives the pillow.
    CURVE1c and SFI1' close the kappa^0 V_0 classes. The kappa-odd results in the record (38(ah)(x), 38(el)) are PERTURBATIVE, at the exact point and first or second order in an analytic source.
A2. Suzuki's eq. (1): -psi'' + (x^(2M) + eps alpha x^(M-1)) psi = E psi.
    The printed ODE has NO l(l+1)/x^2 term; the l-term is our extension, flagged as such.
    With x = E^(1/(2M)) y and kappa = E^((M+1)/(2M)) it becomes psi_yy = [kappa^2 (y^(2M) - 1) + kappa alpha y^(M-1) + lambda y^(-2)] psi.
    So alpha enters at order kappa^1: a kappa-ODD term, LINEAR in alpha. The local IMs depend on alpha^2, so alpha is a SQUARE ROOT of a record-type momentum.
A3. THE GAP IN S1's SCOPE, named precisely.
    - S1 excludes nothing with a kappa^1 term whose coefficient is linear in sqrt(momentum).
    - At spin 2k-1, the top-degree terms are (alpha U)^(2j) (lambda W)^(k-j) Lambda^(1/2-k-j), for j = 0..k, with U = y^(M-1) and W = y^(-2).
    - The top is therefore a SUM of moment families, not the single moment family that S1's uniqueness argument needs. S1 does not apply.
    - The record's other kappa-odd results do not cover this either: they are perturbative in an analytic source at the exact point (38(ah)(x) explicitly leaves "sqrt(h) kappa Q" open).
    - So the class "kappa^2 Lambda + kappa sqrt(X) U + V" at generic fibres is UNTESTED in the record. This is a finding about the record's scope.
A4. Pencil (sealed, to be verified by machine): the extended Suzuki top EQUALS the record's BL top law.
    - The multinomial Gamma(3/2-k-j) cancels against the Beta's, giving top(j) proportional to binom(k,j) (-s/(2M))_j / (1/2)_j (alpha^2/4)^j lambda^(k-j).
    - The BL top binom(k,a) (x_Y)_a (x_X)_(k-a) X^a Y^(k-a) has this form iff x_X = -(k - 1/2), i.e. a_X = 1, and a_Y = 1/M, with X <-> alpha^2 and Y <-> lambda (up to a scale).
    - Sol 1 has a_X = 1 identically (n = eps - 1; VIR2: B = -1). So the match requires M = 1/a_Y = (t+3)/(t-1).
    - Then M = -2 at the exact point t = -1/3; M = 3 at t = 3; M -> 1 in the paperclip limit t -> infinity.
A5. Structural corollary (pencil): Q is invariant under (kappa, alpha) -> (-kappa, -alpha), so:
    - the odd-spin charges are EVEN in alpha (the X-even cylindrical family);
    - the even-spin charges exist only for alpha != 0 and are ODD in alpha, i.e. odd in sqrt(X): "X-odd".
    This is the pattern VIR2 found for Sol 1: a charge at every weight, X-odd at odd weight.
A6. Literature.
    - Babenko-Smirnov, "Suzuki equations and integrals of motion for supersymmetric CFT", Nucl. Phys. B 924 (2017) 406, arXiv:1706.03349. The number is consistent in two printed bibliographies (Tanabe [28] and Ito-Zhu [15]). I have not fetched the paper.
    - Per those printed passages, it concerns the N = 1 SCFT (the SU(2)_2 coset, K = 2), NOT N = 2.
    - Tanabe 2604.14899 eqs (4.56)-(4.58) print the N = 1 NS-sector cylinder IM eigenvalues against C(2)^(2) WKB periods. These depend on one highest weight only, so they cannot be compared with the two-momentum Sol 1.
    - Kojima and Dorey-Tateo N = 2 work: NOT located locally, and not fetched (no-arXiv rule).
    - No printed N = 2 local-IM vacuum eigenvalues are located.
## Step B plan and sealed predictions
- Charges: WKB of the extended Suzuki ODE (Riccati, all hbar orders), with integrals by continued Beta on (0, turning point), int x^e Lambda^f dx -> (-1)^(f-1/2) B((e+1)/(2M), f+1)/(2M), exact Fractions at rational M.
- Dictionary: no certified spin-1 table exists, so the dictionary is fixed at the LOWEST CERTIFIED SPIN, spin 3. It is the affine map alpha^2 = sX X + tX, lambda = sY Y + tY. Only the ratio sY/sX, tX and tY matter after X^k normalisation; they are fixed by the XY, X and Y coefficients of spin 3.
- M = (t+3)/(t-1) is NOT fitted.
- PREDICTIONS:
  - P1: the spin-3 Y^2 coefficient (top, from A4) matches.
  - P2: the spin-3 constant matches. UNDECIDED, a genuine test.
  - P3: all 9 non-normalised spin-5 coefficients match. UNDECIDED.
  - P4: as hold-outs, spins 7 and 9 (w8, w10 tables).
  - P5 (A5): the alpha-odd part of each odd-spin charge integrates to zero.
  - I make NO prediction that P2/P3 pass. The prior is weak: S1's exclusion does not reach here, but nothing forces a match of the lower layers.
- Controls:
  - the exact point t = -1/3 (M = -2);
  - the paperclip limit (large t, M -> 1), approached numerically;
  - tamper M -> M + 1/7, which must fail P1;
  - t-generic fibres 2, 4, 6, 7, -2, -4, 5/2, plus 9/2 if the degeneracy screen allows.
