/-
  The three-term TQ relation of the two sector Q-functions at the c = -2 anchor.

  SPECIFICATION ITEM CERTIFIED: item A3 of the oper specification note
  ("The anchor t = 3 (c = -2)"), the clause

      both towers satisfy  Q(E+4) + a_S(E) Q(E-4) = T(E) Q(E),
      a_S(E) = [(1 - E/2)^2 - beta_S^2]/4,  T(E) = -(E/2),

  established in the ansatz note sec. 38(bh) (batch 56) from the Gamma
  recurrence and checked numerically there to 1e-38 at sixteen points, both
  signs of beta.  Here it is a kernel-checked algebraic identity.

  WHAT THIS FILE DOES NOT COVER.  It is a statement about the ARGUMENT
  ARITHMETIC only, and deliberately so:

  * it does not prove that D_S(E) = Gamma(1+beta_S)/Gamma((1+beta_S)/2 - E/4)
    is the spectral determinant of anything, nor that the anchor's determinant
    is D_X D_Y;
  * it does not touch the Gamma function at all.  The Gamma recurrence
    Gamma(z+1) = z Gamma(z) is used ON PAPER to reduce the ratios
    Q(E+4)/Q(E) and Q(E-4)/Q(E) to x-1 and 1/x, where x = (1+beta)/2 - E/4;
    what is formalized is the identity that remains after that reduction.
    `x_shift_up` and `x_shift_down` certify only that the shifts E -> E +- 4
    move x by exactly -+1, which is the step the reduction needs;
  * it says nothing about beta_X^2 = -4X or 1 - 4Y, i.e. nothing about the
    dictionary from the certified charges to the momenta;
  * it is over the rationals.  The physics statement is over the complex
    numbers; nothing here depends on the field beyond being of characteristic
    zero, but no such generality is claimed;
  * it has NO bearing on the anchor's DEFORMATION.  In particular it must not
    be read as supporting the "Y-shift" reading of stages Q-V, which the
    record refutes (sec. 38(bw), batch 71); this file is about the undeformed
    anchor only.

  Negative controls are included in the file (`tq_wrong_T_fails`,
  `tq_wrong_a_fails`): a perturbed T or a perturbed a must break the identity,
  and is proved to.  The repository's code/lean_check.py supplies a further
  negative control of its own (a false proof must be rejected).

  Contributed by the cognitive agent, 2026-09-18, under the standing
  arrangement recorded in sec. 38(bg).
-/
import Mathlib.Tactic

namespace IntegrableBoundary.AnchorTQ

/-- `x beta E` is the argument of the Gamma in the denominator of
    `D_S(E) = Gamma(1+beta)/Gamma(x)`. -/
def x (beta E : ℚ) : ℚ := (1 + beta) / 2 - E / 4

/-- The coefficient `a_S(E)` of the three-term relation: linear in the
    sector Casimir `beta^2`. -/
def a (beta E : ℚ) : ℚ := ((1 - E / 2) ^ 2 - beta ^ 2) / 4

/-- `T(E)`, carrying no sector dependence. -/
def T (E : ℚ) : ℚ := -E / 2

/-- The step `E -> E + 4` lowers `x` by exactly one.  This is what makes
    `Gamma(x)/Gamma(x-1) = x - 1` the right reduction of `Q(E+4)/Q(E)`. -/
theorem x_shift_up (beta E : ℚ) : x beta (E + 4) = x beta E - 1 := by
  unfold x; ring

/-- The step `E -> E - 4` raises `x` by exactly one, giving
    `Q(E-4)/Q(E) = Gamma(x)/Gamma(x+1) = 1/x`. -/
theorem x_shift_down (beta E : ℚ) : x beta (E - 4) = x beta E + 1 := by
  unfold x; ring

/-- The polynomial core of the relation, with no hypotheses: after clearing
    the single denominator, the three-term relation IS this identity. -/
theorem tq_core (beta E : ℚ) :
    (x beta E) ^ 2 + (x beta E) * (E / 2 - 1) + a beta E = 0 := by
  unfold x a; ring

/-- The three-term relation itself, in the form it takes after division by
    `Q(E)`: `(x - 1) + a/x = T`.  The hypothesis `x ≠ 0` is exactly the
    condition under which `Q(E-4)/Q(E) = 1/x` is meaningful. -/
theorem tq_three_term (beta E : ℚ) (hx : x beta E ≠ 0) :
    (x beta E - 1) + a beta E / x beta E = T E := by
  -- Eliminate `a` with the core identity; what remains is a ring identity.
  have ha : a beta E = -(x beta E) ^ 2 - (x beta E) * (E / 2 - 1) := by
    linear_combination tq_core beta E
  rw [ha]
  unfold T
  field_simp
  ring

/-- NEGATIVE CONTROL.  Perturbing `T` by one breaks the relation: the claim
    fails already at `beta = 3`, `E = 2`, where the true value is `-1`. -/
theorem tq_wrong_T_fails :
    ¬ (∀ beta E : ℚ, x beta E ≠ 0 →
        (x beta E - 1) + a beta E / x beta E = T E + 1) := by
  intro h
  have hx : x 3 2 ≠ 0 := by unfold x; norm_num
  have := h 3 2 hx
  unfold x a T at this
  norm_num at this

/-- NEGATIVE CONTROL.  Perturbing `a` by one breaks the core identity at
    EVERY point, not merely at a witness. -/
theorem tq_wrong_a_fails (beta E : ℚ) :
    (x beta E) ^ 2 + (x beta E) * (E / 2 - 1) + (a beta E + 1) ≠ 0 := by
  have h := tq_core beta E
  intro hbad
  rw [show (x beta E) ^ 2 + (x beta E) * (E / 2 - 1) + (a beta E + 1)
        = ((x beta E) ^ 2 + (x beta E) * (E / 2 - 1) + a beta E) + 1 by ring] at hbad
  rw [h] at hbad
  norm_num at hbad

end IntegrableBoundary.AnchorTQ
