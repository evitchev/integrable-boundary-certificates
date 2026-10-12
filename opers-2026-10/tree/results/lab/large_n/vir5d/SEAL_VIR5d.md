# SEAL_VIR5d -- Fable seat, 2026-10-01.  The spins in nZ: what the operator does there, and the second kernel class.
Follow-up (ii) to VIR5a, commissioned by the lead ("GO, seal first"; Z0 the most valuable piece; for Z2
use VIR1's finding that at t = 5/3 the spin-5 kernel is spanned by the Sol 2 and Sol 3 limit classes).
EXPLORATORY, labelled as such.  Written before any computation of this stage.

Notation as in VIR5c: K = spin index (spin 2K-1), k = fibre parameter, n = k + 4, nu = (2K-1)/n,
phi = (1 - X w)(1 - Y w) = 1 - e1 w + e2 w^2 in operator variables, D = w d/dw, psi_1 = n phi - 2 D phi.
"Lattice point (i, K)": nu = i a positive integer; "drop point": in addition 2i <= K - 1.

## 0. Disclosure
Read before sealing: memory files grade-0-derives-the-symbol, grade0-rank-drop-locus, sheet1-...; the
ansatz note 38(dw) item 3 and 38(dx) (drop sets {5/3}, {11/5}, {7/5, 17/7}, {13/7, 23/9} at K = 3..6;
multiplicity law min(2(K-i), K+1); record's chain "kernel = BL weights (empirical) => signed minors
proportional to the weights (theorem) => drop condition i + j <= K-1 (derived)"); the docstring of
code/cyl_t2_primary.py (kernel dims 2 at t = 7/5 and 1 at t = 11/5, 7/3, -7/5, weight 10); my own VIR1
NOTE (t = 5/3: spin-5 kernel = span of the Sol 2 and Sol 3 limit classes).
Derived by hand before sealing (algebra only, no data): the identity and F1' of sec. 1.
Known from VIR4/VIR5a engine runs: the operator's charge vanishes identically (all layers) at
(k, K) = (1, 3), (3, 4), (5, 5), (1/2, 5).  NOT done: any engine run at the other lattice points; any
deformed-operator run; any kernel computation at t = 11/5, 17/7, 7/5 by this seat; any look at the
second kernel class's eigenvalue anywhere (at t = 5/3 I know only that Sol 2's table differs from
Sol 3's in 7 of 10 coefficients, from SHEET1).

## 1. Z0 -- the vanishing locus of the operator IS the record's diagonal rank-drop locus
(a) PROOF, top layer.  By Lagrange inversion the operator's un-normalised top at spin index K is
    [w^K] phi^nu = (-1)^K sum_a binom(nu, a) binom(nu, K-a) X^a Y^(K-a).  binom(nu, a) = 0 iff nu is an
    integer with 0 <= nu < a.  Since binom(nu, 0) = 1, all K+1 coefficients vanish iff nu = i is an
    integer with i < K and there is no a with K-i <= a <= i, i.e. iff 2i <= K-1.  The number of vanishing
    coefficients at nu = i is min(K+1, 2(K-i)): the record's multiplicity law (38(dx)).  The coefficients
    are the record's BL weights before normalisation ((-nu)_a (-nu)_(K-a) up to a!(K-a)!), so by the
    record's kernel/minor theorem the locus is its diagonal rank-drop locus i + j <= K-1, i = j.
    What this does NOT prove: the record's one empirical link (grade-0 kernel = BL weights).
(b) PROOF, loss-1 layer.  Identity (d = K-1, nu n = 2K-1 = 2d+1):
        (D - d)( phi^nu D log psi_1 ) = phi^nu ( (1/2) D log psi_1 + D^2 log psi_1 ) - (nu/2) phi^(nu-1) D psi_1,
    because nu D log phi = (d + 1/2) - (nu/2) psi_1/phi.  Since [w^d](D - d)(...) = 0, formula F1 becomes
        loss-1 = (nu/6) [w^(K-1)] phi^(nu-1) { 1 + (K - 3/2 - nu) e1 w + (4 nu - 2K + 2) e2 w^2 }.          (F1')
    At nu = i the bracketed series is a polynomial of degree 2i in w whose top coefficient
    (4i - 2K + 2) e2^i vanishes at K = 2i+1; hence the loss-1 layer vanishes identically iff 2i <= K-1:
    the same locus.  F1' also shows why the record needs exactly three shift amplitudes at loss 1.
    To be machine-checked: F1' == F1 (K = 1..12, three fibres) and F1' == Codex's amplitudes (as S3).
(c) CHECKED, not proved: all layers.  Engine (mellin_wkb.py 60795426..., unchanged): the full
    un-normalised R_(2K) is identically zero at the drop points
        K = 6: t = 23/9 (k = 7, i = 1), t = 13/7 (k = 3/2, i = 2);
        K = 7: t = 29/11 (k = 9, i = 1), 19/9 (k = 5/2, i = 2), 9/7 (k = 1/3, i = 3);
        K = 9: t = 21/11 (k = 5/3, i = 3);
    and is NOT zero at the lattice non-drop point K = 4, i = 2 (t = 1/3, k = -1/2), where only the
    X^2 Y^2 top coefficient survives (the record: rank stays 4 there).
    PREDICTION: all six drop points give R = 0 in every layer (90%); t = 1/3 gives R != 0 with top
    proportional to X^2 Y^2 (95%).
(d) Scope remark, part of the claim: for the STRIPPED operator the loss-1 layer at K = 2i+1 does NOT
    vanish (constant n-1 instead of 4 leaves (n-5) e2^i); so at k = 1/i the Gamma ratio is what makes
    the operator's charge vanish.  Check at k = 1/2, K = 5 with the stripped engine: R != 0 (90%).

## 2. Z1 -- the certified Solution-3 charge at a drop point is the limit of the operator's normalised charge
At a drop point the operator at k_0 has no charge; the monic polynomial of the family at k_0 + h has a
limit as h -> 0 (the ratio of k-derivatives; a logarithmic residue).  TEST (exact rational arithmetic at
k = k_0 + h, h = 10^-8 and 10^-16, no floating point): max over coefficients of |monic(k_0 + h) -
certified(t_0)| is O(h): finite, ratio of the two values within a factor 2 of 10^8.  Points: (5/3, spin 5),
(11/5, 7), (17/7, 9), (7/5, 9).  PREDICTION: PASS at all four (90%).  This is a limit demonstration, not an
exact identity; stated as such.

## 3. Z2 -- the second kernel class (the lead's handle), hypothesis H-A
CFT side (my VIR1 engine for the kernel of ad I_3 mod d, VIR2 conversion to two bosons, VIR3 vacuum rule):
at each drop point compute the kernel at weight 2K and the vacuum eigenvalues of its classes as
polynomials in (PX, pi).  Points: (5/3, K = 3), (11/5, K = 4); (17/7, K = 5) and (7/5, K = 5) if the
weight-10 pipeline runs in reasonable time (stated either way).  Write S = the Sol-3 limit charge
(certified table), E = a second independent class.  At t = 5/3 E is checked against the Sol-2 table.
Operator side: the ANTISYMMETRIC power deformation of the same operator, everything else fixed,
    H_delta(T) = T^n (1 - l0^2/T^2)^(1+delta) (1 - l1^2/T^2)^(1-delta) sum_j c_j(k) (n/(kT))^(2j),
whose WKB coefficient R_(2K)(delta) vanishes at delta = 0 at a drop point; A := dR_(2K)/d delta at 0
(engine with coefficients in Q[l0, l1, e]/(e^2); l^2 = c2 (PX^2, pi^2)).  A is odd under PX <-> pi; its top
is the difference of the two corner blocks of the limit top.
H-A: A lies in the certified kernel, i.e. rank{S, E, A} = 2 (exactly, as coefficient vectors), at every
     drop point tested.  Then the kernel is {S, A}: the two classes are the symmetric and the
     antisymmetric logarithmic derivative of ONE operator, and "the charge the operator cannot see" is
     not another family's -- at t = 5/3 the Sol-2 charge would be S + cA.
PREDICTIONS.  (i) kernel dimension 2 at all four points: 85%.  (ii) the kernel is invariant under the
swap PX <-> pi (contains an antisymmetric element): 60%.  (iii) H-A: 40%.  If (ii) holds and H-A fails I
report where A departs from the antisymmetric class (which loss).  If (ii) fails, the second class is
not a Weyl-odd object and H-A is dead as posed.
Controls: plant = the symmetric deformation with both powers 1 + delta is NOT expected in the kernel
(reported, not gated); negative = A from the STRIPPED operator must fail if H-A passes; scoring tamper =
one coefficient of E + 1/1000 must change rank 2 to 3 (after the hash guard).

## 4. Not in this stage (proposal only, nothing computed)
In the two-power family T^n (1 - X/T^2)^alpha (1 - Y/T^2)^beta the record's three top laws are
Sol 3: (alpha, beta, n) = (1, 1, k+4);  Sol 2: (k, 1, 2k+3), i.e. Lambda = T (T^2 - Y)(T^2 - X)^k;
Sol 1: beta = 1, n = 2, alpha = k/(3k+4).  For integer k the Sol-2 symbol is a polynomial of degree 2k+3
with a k-fold X-factor -- a candidate for a string/level-k DDMST-type equation.  I propose it as the
next stage (VIR6) and will not compute on it without a GO.

## 5. Deliverables and exit codes
z0a_f1prime.py (0 = F1' == F1 and == Ward amplitudes); z0c_engine.py (0 = predictions of (c), (d) hold;
2 otherwise); z1_limit.py (0 = PASS at all four; 2 otherwise); z2_kernel.py, z2_deformed.py, z2_compare.py
(0 = H-A passes and controls fire; 2 = H-A fails as sealed; 1 = control failure).  NOTE.md, logs, SHA256SUMS.
