/-
  The two-pattern span at second order, and the two relations that detect it.

  WHAT IS CERTIFIED.  The record's second-order (h^2) statement is that the
  whole degree-(k-2) layer lies in the span of two patterns on the monomial
  index (a, b), a + b = k - 2:

      F : the constant pattern,            F (a, b) = 1
      G : the product of shifted indices,  G (a, b) = (a + 1)(b + 1)

  with explicit coefficients A_k, B_k (ansatz note, the gauge-is-zero entry
  of batch 39; the third-order obstruction is recorded in the same entry).
  Two linear relations on that index are used there to TEST membership of
  the span, at spin 11 and spin 13:

      R11 f = f(0,4) - 4 f(1,3) + 3 f(2,2)
      R13 f = f(0,5) - 3 f(1,4) + 2 f(2,3)

  and the two-pattern model is said to force 91 R11 - 55 R13 = 0, while the
  measured third-order layer gives -2863861/15360 instead.  This file proves
  the half of that which is pure algebra: both relations annihilate both
  patterns, hence the whole span, hence every linear combination of the two
  relations does; and, as negative controls, neither relation is the zero
  functional and a pattern outside the span is NOT annihilated.

  WHY THAT IS THE USEFUL HALF.  A vanishing test statistic is only evidence
  if the test could have failed.  Certifying that R11 and R13 kill exactly
  the two-pattern span, and are not degenerate, is what makes the measured
  third-order value a departure FROM THE SPAN rather than an artefact of the
  relations chosen to look at it.  That is the part of the region the record
  would least like to find convention-dependent, and it is convention-free.

  WHAT THIS FILE DOES NOT COVER, deliberately:

  * it does not prove that the physical degree-(k-2) layer at h^2 IS in the
    span.  That is a statement about the certified charge tables, established
    computationally in the record at spins 7, 9, 11, 13 with no free
    parameter; nothing here touches the data or the extraction;
  * it therefore proves no physics.  It certifies that the TEST is a valid
    test, and nothing about the outcome of running it;
  * the explicit A_k and B_k appear only in two arithmetic checks at k = 2
    and k = 3, reproducing two values the record derives independently; the
    formulas are not derived here;
  * it says nothing about the third-order layer beyond the tautology that a
    nonzero value of the combination puts the layer outside the span.  The
    measured value itself comes from the charge tables, not from here;
  * indices are carried as rationals for convenience; nothing depends on it.

  Contributed by the cognitive agent, 2026-09-19, under the standing
  arrangement, and citing the record by description only.
-/
import Mathlib.Tactic

namespace IntegrableBoundary.SecondOrderSpan

/-- The constant pattern on the monomial index. -/
def Fpat (_a _b : ℚ) : ℚ := 1

/-- The product-of-shifted-indices pattern. -/
def Gpat (a b : ℚ) : ℚ := (a + 1) * (b + 1)

/-- The spin-11 relation (k = 6, a + b = 4). -/
def R11 (f : ℚ → ℚ → ℚ) : ℚ := f 0 4 - 4 * f 1 3 + 3 * f 2 2

/-- The spin-13 relation (k = 7, a + b = 5). -/
def R13 (f : ℚ → ℚ → ℚ) : ℚ := f 0 5 - 3 * f 1 4 + 2 * f 2 3

theorem R11_Fpat : R11 Fpat = 0 := by norm_num [R11, Fpat]

theorem R11_Gpat : R11 Gpat = 0 := by norm_num [R11, Gpat]

theorem R13_Fpat : R13 Fpat = 0 := by norm_num [R13, Fpat]

theorem R13_Gpat : R13 Gpat = 0 := by norm_num [R13, Gpat]

/-- The spin-11 relation annihilates the whole two-pattern span. -/
theorem R11_span (A B : ℚ) :
    R11 (fun a b => A * Fpat a b + B * Gpat a b) = 0 := by
  simp only [R11, Fpat, Gpat]; ring

/-- The spin-13 relation annihilates the whole two-pattern span. -/
theorem R13_span (A B : ℚ) :
    R13 (fun a b => A * Fpat a b + B * Gpat a b) = 0 := by
  simp only [R13, Fpat, Gpat]; ring

/-- Hence EVERY linear combination of the two relations vanishes on the
    span; the record's combination is the case `c₁ = 91`, `c₂ = -55`. -/
theorem combination_vanishes_on_span (c₁ c₂ A B : ℚ) :
    c₁ * R11 (fun a b => A * Fpat a b + B * Gpat a b)
      + c₂ * R13 (fun a b => A * Fpat a b + B * Gpat a b) = 0 := by
  rw [R11_span, R13_span]; ring

/-- The record's own combination, stated as it is used there. -/
theorem record_combination (A B : ℚ) :
    91 * R11 (fun a b => A * Fpat a b + B * Gpat a b)
      - 55 * R13 (fun a b => A * Fpat a b + B * Gpat a b) = 0 := by
  rw [R11_span, R13_span]; ring

/-- NEGATIVE CONTROL.  `R11` is not the zero functional: it separates a
    point of the index set, so the annihilation above is not vacuous. -/
theorem R11_not_identically_zero :
    R11 (fun a b => if a = 0 ∧ b = 4 then 1 else 0) = 1 := by
  norm_num [R11]

/-- NEGATIVE CONTROL.  `R13` is not the zero functional either. -/
theorem R13_not_identically_zero :
    R13 (fun a b => if a = 0 ∧ b = 5 then 1 else 0) = 1 := by
  norm_num [R13]

/-- NEGATIVE CONTROL.  A pattern OUTSIDE the span is not annihilated: the
    square of the product pattern gives 12, so the relations genuinely
    detect departure from the span rather than vanishing on everything. -/
theorem R11_off_span_ne_zero :
    R11 (fun a b => ((a + 1) * (b + 1)) ^ 2) = 12 := by
  norm_num [R11]

/-- The same control for the spin-13 relation. -/
theorem R13_off_span_ne_zero :
    R13 (fun a b => ((a + 1) * (b + 1)) ^ 2) = 24 := by
  norm_num [R13]

/-- The record's explicit second-order coefficients. -/
def Acoeff (k : ℚ) : ℚ := -k * (2 * k - 1) * (8 * k ^ 2 - 13 * k + 3) / 5760

/-- The record's explicit second-order coefficients. -/
def Bcoeff (k : ℚ) : ℚ := k * (2 * k - 1) ^ 2 / 1440

/-- Arithmetic check against a value the record derives independently. -/
theorem Acoeff_Bcoeff_at_two : Acoeff 2 + Bcoeff 2 = 1 / 320 := by
  norm_num [Acoeff, Bcoeff]

/-- Arithmetic check against a value the record derives independently. -/
theorem Acoeff_Bcoeff_at_three : Acoeff 3 + 2 * Bcoeff 3 = 1 / 96 := by
  norm_num [Acoeff, Bcoeff]

end IntegrableBoundary.SecondOrderSpan
