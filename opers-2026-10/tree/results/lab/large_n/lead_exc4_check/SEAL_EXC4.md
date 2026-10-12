# SEAL_EXC4 -- Fable (secondary) seat, 2026-10-07.  (A) Why one apparent singularity is half a level for orders 3-4; (B) level two.
Commissioned by the lead (EXC4_brief_2026-10-07.md 2c5639ec..., COMMON_cluster_rules dbcb573e...).  Written before any computation of
this stage.  Known to me: EXC1-3 (items 313, 320, 324) and their archives; nothing of level two for Sols 2, 3 has been computed by anyone.

## A. The rule (hand derivation; to be tested on cases NOT yet computed)
Setting: a class-U normal form  L_0 = P(vartheta) + xi L(vartheta) - E xi^beta,  P monic of degree m, the non-E middle term with exponent
exactly 1 in xi, the universal dictionary l1 = n pi/p and xi = x^c/(const), c = n/beta.  The Y exponent at xi = 0 is then
a1 = l1/c = beta pi/p, i.e.  a1^2 = (beta^2/p^2) pi^2  (Sol 3: beta = p^2, a1^2 = p^2 pi^2; Sol 2: beta = p^2, a1^2 = p^2 pi^2;
Sol 1 three-term form: beta = p^2/2, a1^2 = p^2 pi^2/4).  The monic I_1 is P_X^2 + pi^2 + const in every case.
Apparent singularity at xi = z with the standard exponent set S_m = {-1, 1, 2, ..., m-2, m} (sum m(m-1)/2), deformation
xi sum_(j<=m-2) R_j(vartheta) vartheta^j with Fuchs pole orders.  Observed in EXC1 (m = 2: r01 = -beta), EXC2 (m = 4: r21 = -beta),
EXC3 (m = 3: r11 = -beta): the E-LINEAR no-logarithm condition at the top exponent forces the simple-pole residue of the
vartheta^(m-2) term to be -beta.  Reason: at first order in E only -E xi^beta = -E z^beta (1 + beta delta/z + ...) enters, so the
condition is linear in the residues with the coefficient beta multiplying the E z^beta term and the vartheta^(m-2) residue the only
deformation term at the same order; the two must cancel.  (To be re-derived by machine as part of the tests, not assumed.)
That residue adds -beta to the vartheta^(m-2) coefficient of the effective polynomial at xi -> infinity, whose momentum content
is -(sum of the squared exponents) = -a1^2 - (a0^2 if the X momentum is an exponent): so a1^2 -> a1^2 + beta, and since
a1^2 = (beta^2/p^2) pi^2 the monic I_1 shifts by
        Delta I_1 (one point) = beta / (beta^2/p^2) = p^2 / beta,        level fraction per point = p^2 / (2 beta)
(one level = +2 in the monic I_1, the convention of EXC1-3).  Sol 1: beta = p^2/2 -> Delta I_1 = 2 = ONE level.  Sols 2, 3:
beta = p^2 -> Delta I_1 = 1 = HALF a level, so level one needs two points, and the sigma pairing supplies exactly two.
RULE (sealed): Delta I_1 per point = p^2/beta, independent of the order m and of the fibre; the number of standard points per level
is 2 beta/p^2.
PREDICTIONS (not yet checked):
  A1. Sol 3 at t = 2 in its sixth-order polynomial form (k = 2: P = Q_6 = (vartheta^2 - A)(vartheta^2 - B)(vartheta^2 - 9/4) up to the
      centring, L = -1 constant, beta = p^2 = 2/3, xi = x^c with the SAME c as the three-term form): one point with S_6 = {-1,1,2,3,4,6}
      has Delta I_1 = 1 (half a level).  In this form the formal adjoint ALONE is a symmetry (P and L are even), so it maps a point
      at z to a point at z (exponents 5 - e, S_6 invariant) and does NOT pair z with -z: the level-one states are predicted to be
      self-adjoint TWO-point solutions at unrelated positions (z1, z2), not +-z pairs, with the I_3 values of EXC2's level-one
      block.  Priors: one-point solutions exist 70%; Delta I_1 = 1 given existence 80%; the two-point class reproduces the level-one
      I_3 block at (1/2, 1/3)*: 50%.   (* momenta as in EXC2: a0 = 1/2, a1 = 1/3.)
  A2. The DUAL Sol 2 form (item 303: M' = -1/(k+1), identical charges), normal form with beta' = 1/b = 1/p^2 (> 1: the E term grows
      faster than the middle term).  One point with S_3 = {-1,1,3}: Delta I_1 = p^2/beta' = p^4 = b^2 = 4/9 at t = 2, a NON-integer
      multiple of half a level; in particular no configuration of standard points gives a level-one state in that form.
      Priors: one-point solutions exist 45%; Delta I_1 = 4/9 given existence 65%.
  A3 (consistency, already known): Sol 1 three-term form, one point {-1,2}: Delta I_1 = 2 (EXC2 T0 reproduced EXC1's quadratic).
If A1 or A2 refutes the rule, the rule is reported refuted with the computed Delta I_1 and no rescue is fitted.

## B. Level two for Solution 3
Class (sealed): FOUR apparent singularities at xi = +-z1, +-z3, each of type {-1,1,2,4}, in the EXC2 normal form
(vartheta^2 - A)(vartheta^2 - B) + xi(vartheta + 1/2) - E xi^b, the configuration invariant under sigma (adjoint composed with xi -> -xi):
the sigma relations of EXC2 within each pair (s21 = r21 = -b, s22 = +4z, s13 = 8z^2, s04 = 8z^3, s11 = -r11, s12 = r12 - 2(b+4)z; s01, s02
left free and eliminated by the point's own linear conditions).  Count: per pair, after the indicial conditions and the linear
eliminations, (r11, r12, z): SIX unknowns against the two pairs' no-log conditions (3 independent per pair in EXC2).
Expected: zero-dimensional; FIVE sigma-invariant solutions = the five level-two singlet states that are not the level-shift state
(EXC1's level-2 count for two bosons; the sixth state is the vacuum at pi^2 -> pi^2 + 4).  I_1 = Delta + 2 on all of them.
Tests: B1 existence and count at t = 2, bare momentum (P, Q) = (1/2, 1/3) [A = b P^2, B = b (Q^2 + rho^2)]; B2 (post hoc) the
I_3 eliminant over the solutions equals the degree-5 factor of the characteristic polynomial of the 6 x 6 level-two singlet block
at that point (d2_level2 for Sol 3, which I compute: one fibre, post hoc); B3 (hold-out, registered): the same at t = 9/4 and at
further points, PREDICTION_EXC4.json hashed to the lead BEFORE those blocks exist; the lead computes and compares first.
I_3 at level two needs only the effective quartic (xi has weight 3 for Sol 3), computed from the sums of residues over the
four points -- no engine extension.  Priors: B1 zero-dimensional with 5 solutions 45% (other counts 35%, none 20%); B2 given 5: 60%;
B3 given B2: 85%.
Controls: EXC2 regression (one pair with the same code = EXC2's quadratic in u^2); ordering tamper xi(vartheta + 1) on the four-point
class (must destroy the solutions or the I_3 match); data tamper on the hold-out; exponent tamper.
Exact arithmetic; Singular for eliminations; runs capped (ulimit -v 16 GB; no outer timeout below 6 h on long jobs; PIDs saved);
full logs; repository read-only; workspace ~/fable-work/exc4.
