# Independent checks for the two Vitchev drafts (Aug 2026)

Pure Python 3 (+ sympy for the symbolic pieces). `ope.py` is a small exact free-boson
OPE engine (multi-species Wick contractions with <d^a X_i d^b X_j> = delta_ij (-1)^(a-1)(a+b-1)!/(z-w)^(a+b),
Euler-operator test for total derivatives, twisted-derivative test for first-order conservation under
e^{alpha.X}, exact linear algebra over Q, Q(sqrt2) (`qsqrt2.py`) or sympy). Conventions coincide with
the repo's `ope_engine.py` / `screening_fit.py`.

Scripts (each prints its own verdicts; runtimes on one core):

- `01_engine_sanity.py` (1 s): T commutes with everything, weights, Table 1 KdV dims 0,2,0,3,1,5,2,8,5,
  ker ad_{I3}(c=1) dims 0,1,0,1,0,1,0,1,0.
- `02_kdv_background_charge_identification.py` (2 s): T_Q = J^2/2 + Q J'; :T_Q^2: = 1/4[J^4 + (4Q^2-2)J'^2] mod d;
  the coefficient A(Q) in I5(Q) = :T_Q^3: + A (dT_Q)^2 fixed by [I3(Q),I5(Q)]=0; class coordinates of I5(Q)
  in (J^6, J^2J'^2, J''^2); at Q^2 = 9/8 (c = -25/2) they equal (1/8)(1, 5/2, 31/24) = paper-2 P6, and equal
  the c=1 A2(2) charge :T_0^3: + 41/16 (dT_0)^2; [I'_3,P6]=0, [I3^{KdV,c=1},P6] != 0.
- `03_tower_equals_kdv_at_c_minus_25_over_2.py` (35 s): ker ad_{I5}, ker ad_{I'_3} (KdV c=-25/2),
  ker ad_{I3^{KdV,c=1}} through charge spin 13; class-by-class equality of the first two. Writes towers.pkl.
- `04_vertex_conditions_and_kernels.py` (15 s): condition polynomial for J^4 + kappa J'^2 under e^{alpha phi}
  = -2 alpha (alpha^4 - (kappa+6) alpha^2 + 4); its factorizations at kappa = 5/2, -2, a^2-6+4/a^2;
  conservation table of the tower members under 1/sqrt2, 2sqrt2, sqrt2 (spins 3..13); single and joint
  screening kernels through spin 11; L(1,0) counts 1,0,1,0,2,0,3,1,4,2,7.
- `05_a22_joint_kernel_at_c1.py` (30 s): joint kernel of the c=1 A2(2) pair {sqrt2, -1/sqrt2} in the even sector:
  1,0,0,0,1,0,1,0,0,0,1,0,1 at spins 1..13.
- `06_spin14_15_extension.py` (50 s): spin 14 vacant; at spin 15 ker ad_{I5} = ker ad_{I'_3} (dim 1) and
  I'_15 is conserved under e^{phi/sqrt2} but not under e^{sqrt2 phi}.
- `07_paperclip_two_boson_checks.py` (25 s) and `08_paperclip_P10_and_spin9_10.py` (40 s): need the repo's
  results/paperclip_P8.json and paperclip_P10.json (set REPO_RESULTS=/path/to/results). At the rational points
  n = -18/25 and n = -50/169 (screening vectors (±3/5, ±4/5), (±5/13, ±12/13)): joint screening kernel dims
  1,0,1,0,1,0,1 at spins 1..7; P8, P10 conserved under all four screenings and equal to the joint-kernel classes;
  [I3,I5] = [I3,P8] = [I5,P8] = [I3,P10] = [I5,P10] = [P8,P10] = 0; ker ad_{I3} dims 0,1,0,1,0,1,0,1,0 at spins 2..10;
  dim Q_sigma = 0,5,1,11,7,25,21,55,57.
