/-
  The classical jumping locus as the pole lattice of the pillow exponents
  (ansatz note sec. 38(w), 2026-09-07).  Statement-first formalization:
  each sheet's double-lattice point t(s, j) puts the momentum-carrying
  exponents of the pillow-type leading symbol on the WKB Beta lattice.
  Pure field identities over ℚ with explicit non-degeneracy hypotheses;
  the physics (which exponent carries which momentum, and the lattice
  itself) is the input of secs. 26, 33 and 38(o).
-/
import Mathlib.Tactic

namespace IntegrableBoundary.LatticeExponents

/-- Solution 3 at t = (3s-10j)/(s-2j): the exponent -(t-3)/(t-5) equals -2j/s. -/
theorem sol3_exponent_on_lattice (s j : ℚ) (hs : s ≠ 0) (h1 : s - 2 * j ≠ 0) :
    -((3 * s - 10 * j) / (s - 2 * j) - 3) / ((3 * s - 10 * j) / (s - 2 * j) - 5) = -2 * j / s := by
  have h1' : s - j * 2 ≠ 0 := by rwa [mul_comm]
  have e1 : (3 * s - 10 * j) / (s - 2 * j) - 3 = -4 * j / (s - 2 * j) := by field_simp; ring
  have e2 : (3 * s - 10 * j) / (s - 2 * j) - 5 = -2 * s / (s - 2 * j) := by field_simp; ring
  rw [e1, e2]
  field_simp
  ring

/-- Solution 1 at t = (6j+s)/(s-2j): the Y-point exponent -(t-1)/(t+3) equals -2j/s
    (the X-end exponent is frozen at -1, not a statement about t). -/
theorem sol1_exponent_on_lattice (s j : ℚ) (hs : s ≠ 0) (h1 : s - 2 * j ≠ 0) :
    -((6 * j + s) / (s - 2 * j) - 1) / ((6 * j + s) / (s - 2 * j) + 3) = -2 * j / s := by
  have h1' : s - j * 2 ≠ 0 := by rwa [mul_comm]
  have e1 : (6 * j + s) / (s - 2 * j) - 1 = 8 * j / (s - 2 * j) := by field_simp; ring
  have e2 : (6 * j + s) / (s - 2 * j) + 3 = 4 * s / (s - 2 * j) := by field_simp; ring
  rw [e1, e2]
  field_simp
  ring

/-- Solution 2 at t = (5j+2s)/(2s-j): the Y-point exponent -4(t-1)/(t+5) equals -2j/s. -/
theorem sol2_Y_exponent_on_lattice (s j : ℚ) (hs : s ≠ 0) (h1 : 2 * s - j ≠ 0) :
    -4 * ((5 * j + 2 * s) / (2 * s - j) - 1) / ((5 * j + 2 * s) / (2 * s - j) + 5) = -2 * j / s := by
  have e1 : (5 * j + 2 * s) / (2 * s - j) - 1 = 6 * j / (2 * s - j) := by field_simp; ring
  have e2 : (5 * j + 2 * s) / (2 * s - j) + 5 = 12 * s / (2 * s - j) := by field_simp; ring
  rw [e1, e2]
  field_simp
  ring

/-- Solution 2 at the same point: the X-end exponent 2(t-3)/(t+5) equals -2j'/s with j' = (s-2j)/3. -/
theorem sol2_X_exponent_on_lattice (s j : ℚ) (hs : s ≠ 0) (h1 : 2 * s - j ≠ 0) :
    2 * ((5 * j + 2 * s) / (2 * s - j) - 3) / ((5 * j + 2 * s) / (2 * s - j) + 5) = -2 * ((s - 2 * j) / 3) / s := by
  have e1 : (5 * j + 2 * s) / (2 * s - j) - 3 = (8 * j - 4 * s) / (2 * s - j) := by field_simp; ring
  have e2 : (5 * j + 2 * s) / (2 * s - j) + 5 = 12 * s / (2 * s - j) := by field_simp; ring
  rw [e1, e2]
  field_simp
  ring

/-- The coupling at the Solution-3 lattice point: (t-1)/(t+1) = (s-4j)/(2(s-3j)). -/
theorem sol3_coupling_at_lattice (s j : ℚ) (h1 : s - 2 * j ≠ 0) (h2 : s - 3 * j ≠ 0) :
    ((3 * s - 10 * j) / (s - 2 * j) - 1) / ((3 * s - 10 * j) / (s - 2 * j) + 1) = (s - 4 * j) / (2 * (s - 3 * j)) := by
  have h1' : s - j * 2 ≠ 0 := by rwa [mul_comm]
  have h2' : s - j * 3 ≠ 0 := by rwa [mul_comm]
  have e1 : (3 * s - 10 * j) / (s - 2 * j) - 1 = (2 * s - 8 * j) / (s - 2 * j) := by field_simp; ring
  have e2 : (3 * s - 10 * j) / (s - 2 * j) + 1 = (4 * s - 12 * j) / (s - 2 * j) := by field_simp; ring
  have h3 : (4 * s - 12 * j) ≠ 0 := by intro h; exact h2 (by linarith)
  have h4 : (2 * (s - 3 * j)) ≠ 0 := by intro h; exact h2 (by linarith)
  rw [e1, e2, div_div_div_cancel_right₀ h1]
  rw [div_eq_div_iff h3 h4]
  ring

/-- Sheet 1 never lifts: beta^2 = 1/(2k) forces s = 2j(4k-1), an even integer. -/
theorem sol1_lifting_spin_even (j k : ℤ) : Even (2 * j * (4 * k - 1)) :=
  ⟨j * (4 * k - 1), by ring⟩

end IntegrableBoundary.LatticeExponents
