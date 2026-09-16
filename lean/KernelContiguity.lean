/-
  The three-term contiguity relation among the shifted top-layer kernels
  (client ansatz note sec. 38(az) item 2; lab/kernel_contiguity.py checked it
  symbolically for k = 2..9).  Statement-first formalization, then the proof:
  with the rising factorial (x)_n and the shifted kernel

    w(k, i; dX, dY) = C(k, i) (xY + dY)_i (xX + dX)_(k-i) r^(k-i) / (xY + dY)_k ,

  for every k and every i <= k,

    xX * w(1, 0) + (xY + k) * w(0, 1) = (xX + xY + k) * w(0, 0).

  The proof is the two elementary shift identities x (x+1)_n = (x)_n (x+n),
  applied once on the X side and twice on the Y side.  Pure field identities
  over Q with explicit non-degeneracy hypotheses; the physics (which kernel,
  which shifts) is the input.
-/
import Mathlib.Tactic
import Mathlib.RingTheory.Polynomial.Pochhammer

namespace IntegrableBoundary.KernelContiguity

open Polynomial

/-- The rising factorial (x)_n = x (x+1) ... (x+n-1) as a rational function of x. -/
noncomputable def rf (x : ℚ) (n : ℕ) : ℚ := (ascPochhammer ℚ n).eval x

lemma rf_succ_left (x : ℚ) (n : ℕ) : rf x (n + 1) = x * rf (x + 1) n := by
  unfold rf
  rw [ascPochhammer_succ_left]
  simp [Polynomial.eval_mul, Polynomial.eval_comp]

lemma rf_succ_right (x : ℚ) (n : ℕ) : rf x (n + 1) = rf x n * (x + n) := by
  unfold rf
  rw [ascPochhammer_succ_eval]

/-- The shift identity x (x+1)_n = (x)_n (x+n). -/
lemma rf_shift (x : ℚ) (n : ℕ) : x * rf (x + 1) n = rf x n * (x + n) := by
  rw [← rf_succ_left, rf_succ_right]

/-- The shifted two-index kernel of the top momentum layer. -/
noncomputable def w (k i : ℕ) (xX xY r dX dY : ℚ) : ℚ :=
  (Nat.choose k i : ℚ) * rf (xY + dY) i * rf (xX + dX) (k - i) * r ^ (k - i)
    / rf (xY + dY) k

/-- The three-term contiguity relation, for every spin index k and every i <= k. -/
theorem contiguity (k i : ℕ) (hi : i ≤ k) (xX xY r : ℚ)
    (hX : xX ≠ 0) (hY : xY ≠ 0) (hYk : rf xY k ≠ 0) (hYk' : xY + k ≠ 0) :
    xX * w k i xX xY r 1 0 + (xY + k) * w k i xX xY r 0 1
      = (xX + xY + k) * w k i xX xY r 0 0 := by
  have hX1 : xX * rf (xX + 1) (k - i) = rf xX (k - i) * (xX + ((k - i : ℕ) : ℚ)) :=
    rf_shift xX (k - i)
  have hYi : xY * rf (xY + 1) i = rf xY i * (xY + i) := rf_shift xY i
  have hYK : xY * rf (xY + 1) k = rf xY k * (xY + k) := rf_shift xY k
  have hY1k : rf (xY + 1) k ≠ 0 := by
    intro h
    have : xY * rf (xY + 1) k = 0 := by rw [h, mul_zero]
    rw [hYK] at this
    exact (mul_ne_zero hYk hYk') this
  have hcast : ((k - i : ℕ) : ℚ) = (k : ℚ) - i := by
    rw [Nat.cast_sub hi]
  -- Replace the shifted factorials by the unshifted ones through the shift identities.
  have eX : rf (xX + 1) (k - i) = rf xX (k - i) * (xX + ((k - i : ℕ) : ℚ)) / xX := by
    rw [eq_div_iff hX, mul_comm, hX1]
  have eYi : rf (xY + 1) i = rf xY i * (xY + i) / xY := by
    rw [eq_div_iff hY, mul_comm, hYi]
  have eYk : rf (xY + 1) k = rf xY k * (xY + k) / xY := by
    rw [eq_div_iff hY, mul_comm, hYK]
  unfold w
  simp only [add_zero]
  rw [eX, eYi, eYk, hcast]
  field_simp
  ring

end IntegrableBoundary.KernelContiguity
