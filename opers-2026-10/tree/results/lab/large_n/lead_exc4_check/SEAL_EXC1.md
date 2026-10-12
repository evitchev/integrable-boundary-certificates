# SEAL_EXC1 -- Fable seat, 2026-10-02.  Excited states of the new opers: level 1, Solution 1 first.
Commissioned by the lead (queued behind RED1).  Written before any computation of this stage: no excited-state matrix has
been computed by this seat, and no modified-potential WKB has been run.

## 0. Disclosure
Read: ansatz note 38(k)-(l) (level-shift law; level-1 singlet block = two states, irreducible quadratic; branch locus
D_q at the q-fibres) and lab/excited_states.py (docstring and code).  The data quoted there are for the q-fibres
(t = 5/3, 7/5, 9/7), not for the fibres I use.
Hand derivations made before sealing (algebra only):
(i) Sol 1's vacuum oper in Schrodinger form (item 306; my VIR7 dictionary):
      -psi'' + [ x^(2M) + alpha x^(M-1) + (lam^2 - 1/4)/x^2 ] psi = E psi,   M = (t+3)/(t-1),  lam^2 = (M+1) pi^2,  alpha^2 = 4M(M+1) PX^2.
(ii) For psi'' = [W(x) + 2/(x - x0)^2] psi with W regular at x0, the monodromy at x0 is trivial iff W'(x0) = 0.
(iii) With xi = x^(M+1) and dV_z = -2 d_x^2 log(xi - z):  x^2 dV_z = 2(M+1) [ 1 + sum_{j>=1} (j(M+1)+1) (z/xi)^j ]  at large xi, and the
      regular part of dV_z at a root x0 of xi = z has derivative  -M(2-M)/(2 x0^3).
(iv) ONE point: the condition is the quadratic  2M z^2 + alpha(M-1) z - 2 lam^2 + (M-1)^2/2 = 0.
(v) TWO points z1, z2: for k = 1, 2, with r = z_k/(z_k - z_k') and u = M+1,
      2M z_k^2 + alpha(M-1) z_k - 2 lam^2 + (M-1)^2/2 - 2 [ u^3 (r - r^2)(1 - 2r) - 3 u^2 (r - r^2) + 2 u r ] = 0.

## 1. Data side (the record's tool, exact)
D1  Level-1 O(N-1)-singlet block (two states: the X oscillator and the Y oscillator along the momentum) of I_3 and I_5 for
    SOLUTION 1 and SOLUTION 3 at t = 2 (N = -7) and t = 9/4: the 2x2 matrices, their commutator, and the characteristic
    polynomials in (PX, pi) (bare P = PX, bare Q^2 = pi^2 - rho^2).  I_5 densities from my VIR1 kernel (ad I_3 mod d).
    Also the transverse block (one state) as a check of the level-shift law Q^2 -> Q^2 - 2 (in my variables: the vacuum
    polynomial at pi^2 -> pi^2 + 2).  In the two-boson language the singlet block is the level-1 subspace of the Fock module;
    the transverse state is a new Virasoro primary of weight h_Y + 1 (not a descendant of the Vir_N primary).
    PREDICTION: [I_3, I_5] = 0 on every block; singlet block irreducible quadratic at generic momenta for both solutions;
    transverse block = level-shift law.  90%.

## 2. Oper side, Solution 1 (second order; the classical BLZ construction adapted to the alpha term)
ANSATZ (sealed): level-1 potentials  V = V_vac - 2 d_x^2 [ log(xi - z1) + log(xi - z2) ],  xi = x^(M+1), trivial monodromy at
both points (system (v)).  Charges from the Mellin-WKB engine extended by potential terms (validated first on the vacuum:
the alpha-term engine must reproduce VIR7's S1 charges).
O1  One point (system (iv)): what does a single apparent singularity describe?  PREDICTION: I_1 = Delta + 1/2, i.e. NOT a
    state of the level-1 block; I expect it to be a vacuum oper with shifted momenta (reported, with the shift if found).  60%.
O2  Two points: the system (v) has, up to the exchange z1 <-> z2, solutions whose I_1 is Delta + 1 and whose I_3 values are
    the two roots of the characteristic polynomial of D1's singlet block at the same (PX, pi), at t = 2 and t = 9/4,
    tested at three rational momentum points each (exact algebraic numbers / resultants; no floating point).  45%.
    If the system has more solutions than two, the extra ones are reported with their I_1, I_3.
O3  HOLD-OUT: for the solutions that pass O2, the I_5 eigenvalues are predicted from the oper BEFORE the I_5 level-1 matrix
    is computed (the I_5 part of D1 is run only after the prediction file is hashed and sent to the lead).  If O2 passes: 80%.
If O2 FAILS: reported as a finding (which invariant fails, at which point); then, labelled post-hoc, the smallest
modification that works, if any (e.g. points in the x^(2M+2) plane, or different residue structure).

## 3. Solution 3 (fourth-order three-term ODE)
Only the data (D1) in this stage.  The trivial-monodromy construction for the fourth-order operator needs the local
form of an apparent singularity of a class-U operator, which I have not derived; I will not guess it.  If Sol 1 works,
I propose it as EXC2 with the Sol-3 level-1 data of D1 as its hash-registered target.

## 4. Controls
Vacuum engine control (alpha-term engine == VIR7 Gamma-form S1 charges, spins 1..7, and even spins).  BLZ control: at
alpha = 0 the symmetric two-point solution z2 = -z1 must exist and give the BLZ level-1 value.  Tamper: coefficient 2 of the
double pole replaced by 2 + 1/10 in the WKB (i.e. the series (iii) scaled) must break O2.
Exit codes: d1_data.py 0 = predictions hold; o_sol1.py 0 = O2 passes and controls behave, 2 = O2 fails, 1 = control misbehaves.
