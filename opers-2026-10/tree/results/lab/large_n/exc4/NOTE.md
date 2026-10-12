# EXC4 (Fable secondary seat, 2026-10-07) -- (A) why one apparent singularity is half a level for orders 3-4; (B) level two for Solution 3

Seal: `SEAL_EXC4.md` 10949976... (2026-10-07T07:20:18Z, snapshotted by the lead at 07:24Z before any computation).
Prediction: `PREDICTION_EXC4.json` cff1aead... (11:19:51Z), registered by the lead at 11:20:35Z; the lead computed and compared the
level-two blocks first (match at all three points); this seat opened them at 11:24:08Z.
STATUS: A1 one point HOLDS; A1 two points UNDECIDED (resource); A2 NOT TESTABLE as sealed (object mismatch); B1 count = 5 as
expected; B2 (post hoc) IDENTICAL; B3 (blind) IDENTICAL at all three hold-out points.

## A. The half-level rule

**A.1 Rule and derivation (sealed).**  In a class-U normal form P(vartheta) + xi L(vartheta) - E xi^beta (P monic of degree m, the
non-E middle term with xi^1), an apparent singularity with the standard exponent set S_m = {-1, 1, ..., m-2, m} has its
vartheta^(m-2) simple-pole residue forced to -beta by the E-linear no-logarithm condition at the top exponent.  That residue
adds -beta to the vartheta^(m-2) coefficient of the effective polynomial at xi -> infinity, i.e. +beta to the sum of the squared
exponents, in effect to a1^2 (the Y exponent), and the monic I_1 = P_X^2 + pi^2 + const shifts by
        Delta I_1 = beta / kappa,   kappa = a1^2 / pi^2.
With the universal dictionary (l1 = n pi/p, c = n/beta, a1 = beta pi/p, kappa = beta^2/p^2) this is the sealed formula
        Delta I_1 = p^2 / beta,   standard points per level = 2 beta / p^2
(one level = +2 in the monic I_1): Sol 1's three-term form (beta = p^2/2) one point = one level; Sols 2, 3 (beta = p^2) one point =
half a level, and the sigma pairing supplies exactly the two.  In terms of the Gamma-string length a (b = a/(1+a)): points per
level = 2a(k+1)/(k(1+a)), which is 2 for a = k (Sols 2, 3) and 1 for a = k/(k+2) (Sol 1).  The residue mechanism was observed at
m = 2, 3, 4 (EXC1-3) and is confirmed at m = 6 below.

**A.2 Prediction A1, one point (sixth-order polynomial form of Sol 3 at t = 2): HOLDS.**  Q6(vartheta) - xi - E xi^b with
Q6 = (vartheta^2 - A)(vartheta^2 - B)(vartheta^2 - 1/36) (the frozen T-exponents +-3/2 are +-3/(2c) = +-1/6 in xi; the seal's "9/4 up
to the centring" meant this), L = -1, b = 2/3, xi = x^c with the same c = 9 as the three-term form.  T0: the Gamma-form engine on the
vacuum Q6 reproduces the certified I_3 and I_1 = P^2 + Q^2 + 1/2 (`f_six_t2_a.log`).  One point with S_6 = {-1,1,2,3,4,6}: 5 indicial
conditions, 16 no-log equations in the remaining 16 unknowns, zero-dimensional, TWO solutions (z = 7/5832, r32 quadratic);
the E-linear condition at the top exponent is r05 + 2 r14 z + 11 r23 z^2 + 27 r32 z^3 + 72 r41 z^4 + 372 z^4 = 0, which with the
indicial values gives r41 = -2/3 = -b; the monic I_1 is vacuum - 3 r41/2 = vacuum + 1 (`f_six_charges_t2_a.log`): HALF a level,
as sealed (prior 80%).  Both solutions have I_3 = 41/64, the same number as one of the two one-point values of the three-term
form in EXC2 (105/64 and 41/64), and not a root of the level-one data quadratic (1.39 +- 1.73 i at (1/2,1/3)).

**A.3 Prediction A1, two points: UNDECIDED (resource).**  The two-point class in the sixth-order form, residues free (the adjoint
invariance was to be read off afterwards), has 32 equations in 33 unknowns (with q = 1/(z1 - z2)); Singular reached the 16 GB cap
after 3 h 03 ("halt 14", `f_six_two_t2_a.log`).  Not retried with the self-adjoint reduction (V_3 = 2 vartheta V_4, V_1 = vartheta V_2
- vartheta^3 V_4, halving the unknowns) for lack of time; the sealed 50% prior stands untested.  So in the sixth-order form it is
established that one point is half a level, not what carries level one.

**A.4 Prediction A2 (the dual Solution 2 form): NOT TESTABLE AS SEALED -- object mismatch, disclosed.**  The seal assumed that
item 303's dual (M' = -1/(k+1), identical charges) is the class-U three-term form with the same P, L and beta' = 1/b.  It is not:
SOL2F's dual is a Gamma SYMBOL with strings of the same length k but spacing delta' = nM', whereas Theorem 1 ties the spacing to
the length (sigma = n/a); the three-term form (P, L, beta' = 1/b) is Mellin-equivalent to strings of length -(k+1), a different
operator.  Its vacuum through the three-term engine at M' = -1/3 (unvalidated at negative M) does not reproduce the certified
Sol 2 vacuum with either momentum dictionary (`f_dual_sol2_t2_v0_same_l.log`, `f_dual_sol2_t2.log`).  The sealed number 4/9 is
therefore neither confirmed nor refuted.  The exercise showed that the formula p^2/beta rests on a1^2 = (beta^2/p^2) pi^2; the
form Delta I_1 = beta/kappa (A.1) is the one to keep.

## B. Level two for Solution 3: four apparent singularities, five states

**B.1 Class and count (sealed).**  EXC2's normal form (vartheta^2 - A)(vartheta^2 - B) + xi(vartheta + 1/2) - E xi^b with FOUR points at
xi = +-z1, +-z3, each of type {-1,1,2,4}, the sigma relations of EXC2 inside each pair (s21 = r21 = -b, s22 = 4z, s13 = 8z^2,
s04 = 8z^3, s11 = -r11, s12 = r12 - 2(b+4)z; s01, s02 eliminated by the point's own conditions).  After the indicial conditions and
the linear eliminations (r02, r01, s02, s01 and the same for the second pair) the system has 20 equations in SIX unknowns
(r11, r12, z, u11, u12, z3).  Solved with the inverse-position variables 1/z, 1/z3, 1/(z - z3), 1/(z + z3) as polynomial variables
(so that nothing is cleared to high degree): msolve (grevlex) then Singular fglm (lex).  Result at every point tested:
zero-dimensional, vdim 40 ordered = FIVE unordered configurations (the 8 relabellings of {+-z1, +-z3}), I_1 = vacuum + 4 on all
(level two), as the sealed expectation (five level-two singlet states besides the level-shift state; prior 45%).
One-pair regression (same code, 'pair' mode): vdim 4, I_1 = vacuum + 2, and the I_3 eliminant equals EXC2's level-one data
quadratic at bare (P,Q) = (1/2,1/3) (`chk_pair.log`).  Ordering tamper xi(vartheta + 1): NO solutions (`h_fourq_t2_P1-2_Q1-3_tamper.log`).

**B.2 Post hoc (t = 2, bare (P,Q) = (1/2,1/3)).**  The level-two singlet block (6 states; EXC1's tool `d2_level2`) has characteristic
polynomial (e - level-shift) x quintic, the linear root being the vacuum at Q^2 -> Q^2 + 4 (EXC1's law, 24433/1296 here); the
monic I_3 eliminant over the five four-point solutions is IDENTICAL to the quintic (`chk_four_gb.log`).

**B.3 Blind hold-out: IDENTICAL at all three points.**  Registered quintics at t = 9/4 (1/2,1/3), t = 2 (1/3,1/5), t = 9/4 (2/5,1/4)
(and I_1 = vacuum + 4).  The lead computed the blocks and compared first: identical at all three; this seat's comparison
(`u_compare_exc4.py real`, exit 0): 0 of 6 coefficients differ at each point, level-shift root present; tamper (+1/1000 on the
data quintic) fires at each point.

**B.4 Not done.**  I_5 at level two (needs the first order in 1/xi: for Sol 3 it enters only at I_5, so the sums of the first-order
residues over the four points would do -- not computed); closed form of the two-pair system; level two for Sol 2; other fibres.

## Parameter count against conditions (sealed vs found)
Sixth-order one point: 21 unknowns, 5 indicial, 16 equations: 2 solutions.  Four points (Sol 3): per pair after sigma and the
indicial conditions 7 unknowns (r11, r12, r01, r02, s01, s02, z), 20 conditions for the four points; linear eliminations leave 6
unknowns and 20 equations; 5 unordered solutions.  Nothing fitted at any step.

## Errors, deviations, disclosures
- A2: object mismatch in the seal (A.4); no rescue fitted.  A1 two points: undecided at the memory cap (A.3).
- The direct route (sympy eliminations with factor, denominators cleared) produced equations of degree up to 26 in six unknowns;
  its Singular runs (t = 2 post hoc 2 h 22, t = 9/4 (2/5,1/4) 1 h 53, ordering tamper 1 h 10) and an msolve run on the same system
  (1 h 49) were killed by saved PID when the inverse-variable route (`h_four_q.py`) solved all four points in minutes
  (`pids_killed.txt`).  The Singular `std` cross-checks of the inverse-variable systems were also killed (unfinished after ~55 min);
  the record is msolve's basis + Singular's fglm, cross-validated by the one-pair regression and by the data at four points.
- `h_four_q.py` v1 generated the systems without the elimination loop (a slicing error when deriving it from `h_four_sol3.py`) and
  its Singular run was killed; v2 (with the loop) first failed on Singular's parsing of `x^2/4` (fixed by clearing denominators,
  as in EXC3); the tamper result of that run is the recorded one.
- `f_six_sol3.py` crashed after its Singular run when building the engine symbol (generator list); the one-point charges are in
  `f_six_charges.py` with the generic-symbol copy `mellin_gen6.py`; its first run hit the same `x^2/4` parsing point.
- `f_six_two.py` was launched once with that parsing defect and killed by saved PID; `pids_killed.txt` records two mistaken kills
  of this seat's own shell (pgrep -f matched its command line); no computation was affected.
- The lead changed session twice during the stage; messages were resent to the current address as instructed.
- Exit codes: f_six_sol3 1 (post-result crash; Singular output complete); f_six_charges 0; f_dual_sol2 2 (A.4); chk_pair equal;
  h_four_q pair/four/tamper: vdim 4 / 40 / 0; gb_elim 0 x 4; chk_four_gb 0; d2_level2_exc4 0 x 4; u_compare_exc4 real 0, tamper 0 (fires).

## Files
`SEAL_EXC4.md` (+ .sha256, .time), `PREDICTION_EXC4.json` (+ .sha256, .time), `HOLDOUT_OPENED.time`, `frob.py`, `wkb3x2.py`, `vir_lib.py`,
`mellin_gen.py`, `mellin_gen6.py`, `d1_data_exc4.py`, `d2_level2_exc4.py`, `f_six_sol3.py`, `f_six_charges.py`, `f_six_two.py`, `f_dual_sol2.py`
(+ v0), `h_four_sol3.py`, `h_four_fast.py`, `h_four_q.py`, `gb_elim.py`, `chk_pair.py`, `chk_four.py`, `chk_four_gb.py`, `u_predict_exc4.py`,
`u_compare_exc4.py`, `d2_sol3_*.json`, `h_fourq_*` (sing/out/ms/msgb), `gb_elim_*`, logs, `pids.txt`, `pids_killed.txt`, `SHA256SUMS`.
