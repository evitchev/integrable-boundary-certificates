/-
  G1 equals the Ward loss-1 layer, for every spin, on all three solutions.

  THE CLAIM (the record's item 311 (3), proved there in sympy).  For each solution,
  Fable's first-quantum-order formula G1 (item 306), applied to the solution's
  Gamma-symbol operator, gives a loss-1 coefficient of X^a Y^(K-1-a).  Normalised
  as the record does,

      R = -(L / top) / CY - (K - a) rho^2,     rho^2 = 2/((t-1)(t+1)),

  it equals Codex's Ward-derived rational function Num(a, K, t)/Den(a, K, t)
  (the record's projection_over_c_a laws, from the audit archive of the Ward
  derivation).

  WHAT IS PROVED, for every spin index K, every 0 <= a <= K-1, and every rational
  t off the degenerate set named by the hypotheses (t != +-1, the solution's
  poles, b != 0, 1 for both exponents, b1 - (K-a) + 1 != 0, top != 0):
  * `g1_ward_sol1`, `g1_ward_sol2`, `g1_ward_sol3`: (R) x Den = Num, with no
    hypothesis on Den;
  * `g1_ward_sol*_ratio`: R = Num/Den where Den != 0.

  HOW.  `Lsum1/2/3` write the G1 loss-1 coefficient as the finite sum of shifted
  generalised binomials C(b - D, j - o) that item 311's term list reduces to. Here
  `gbz` is the generalised binomial with zero for a negative lower index, so the
  rising factorials carry the zeros, and `sh` is the shifted binomial.  Two shift
  lemmas (`sh_le_mul`, `sh_gt_mul`) rewrite every C(b - D, j - o) as C(b, j) times
  a ratio of falling factorials.  The top binomials then cancel, and what remains
  is a rational-function identity in (b0, b1, a, K, t), closed by clearing
  denominators, substituting b0 and b1, and `ring`.  The main theorems use a raised
  heartbeat budget; the file checks in about a minute.

  WHAT THIS FILE DOES NOT COVER, deliberately:
  * The TRANSCRIPTION: `Lsum` is G1 expanded into binomials. The cognitive agent
    generated it mechanically from item 311's own term list and checked it outside
    Lean at 72 exact rational points, (t, K, a) on all three solutions, against
    BOTH Fable's power-series implementation of G1 and Codex's laws: 72/72
    three-way agreement wherever the record's normalisation is defined.  The
    expansion of G1's power-series definition into this sum is not itself
    formalised;
  * G1 itself: Fable's formula is used as defined. Its derivation from the
    operator (item 306) is not re-proved, as in item 311;
  * Codex's laws are INPUTS, copied from the record's audit archive; their Ward
    derivation is not formalised;
  * the operator dictionaries (n, CY, the exponents b0 = nu alpha_0, b1 = nu,
    nu = (2K-1)/n) are inputs from items 298, 303 and 306;
  * nothing at loss >= 2, nothing about rank-drop vanishing beyond what follows
    here, and no spectral statement.

  NEGATIVE CONTROLS, in the file, all at Solution 3, t = 9/4, K = 3, a = 1:
  * the exact values L = 455/88 and top = -1575/21296;
  * NON-VACUITY: every hypothesis of `g1_ward_sol3` is discharged there, and both
    sides equal 21/32;
  * CY is load-bearing: with CY scaled by 11/10 the identity fails there
    (-567/176 against 21/32).

  Contributed by the cognitive agent, 2026-10-02, under the standing
  arrangement, citing the record by description only.
-/
import Mathlib.Tactic

namespace IntegrableBoundary.WardLossOne

open Finset

/-- The falling factorial `x (x-1) ... (x-n+1)`. -/
def ff (x : ℚ) : ℕ → ℚ
  | 0 => 1
  | n + 1 => ff x n * (x - n)

@[simp] lemma ff_zero (x : ℚ) : ff x 0 = 1 := rfl
lemma ff_succ (x : ℚ) (n : ℕ) : ff x (n + 1) = ff x n * (x - n) := rfl
lemma ff_one (x : ℚ) : ff x 1 = x := by simp [ff_succ]
lemma ff_two (x : ℚ) : ff x 2 = x * (x - 1) := by simp [ff_succ]
lemma ff_three (x : ℚ) : ff x 3 = x * (x - 1) * (x - 2) := by
  rw [show (3 : ℕ) = 2 + 1 from rfl, ff_succ, ff_two]; push_cast; ring

lemma ff_add (x : ℚ) (m n : ℕ) : ff x (m + n) = ff x m * ff (x - m) n := by
  induction n with
  | zero => simp
  | succ n ih => rw [← add_assoc, ff_succ, ih, ff_succ]; push_cast; ring

lemma ff_nat_lt (j o : ℕ) (h : j < o) : ff (j : ℚ) o = 0 := by
  induction o with
  | zero => omega
  | succ o ih =>
    rw [ff_succ]
    rcases Nat.lt_succ_iff_lt_or_eq.mp h with h' | h'
    · rw [ih h', zero_mul]
    · rw [h', sub_self, mul_zero]

lemma ff_nat_fact (j o : ℕ) (h : o ≤ j) : ff (j : ℚ) o * (j - o).factorial = j.factorial := by
  induction o with
  | zero => simp
  | succ o ih =>
    rw [ff_succ]
    have h1 : o ≤ j := by omega
    have := ih h1
    have e : j - o = (j - (o + 1)) + 1 := by omega
    rw [e, Nat.factorial_succ] at this
    rw [← this]; push_cast
    rw [show ((j - (o + 1) : ℕ) : ℚ) + 1 = (j : ℚ) - o by
      rw [Nat.cast_sub (by omega)]; push_cast; ring]
    ring

/-- Generalised binomial `C(x, n) = ff x n / n!`. -/
noncomputable def gb (x : ℚ) (n : ℕ) : ℚ := ff x n / n.factorial

/-- `C(x, m)` for an integer lower index, zero when `m < 0`. -/
noncomputable def gbz (x : ℚ) (m : ℤ) : ℚ := if m < 0 then 0 else gb x m.toNat

/-- The shifted binomial `C(x - D, j - o)` (the generator's `sh`). -/
noncomputable def sh (x : ℚ) (D o j : ℕ) : ℚ := gbz (x - D) ((j : ℤ) - o)

lemma gbz_nat (x : ℚ) (j : ℕ) : gbz x j = gb x j := by simp [gbz]

/-- Shift lemma, `o ≤ D`. -/
lemma sh_le_mul (x : ℚ) (D o j : ℕ) (h : o ≤ D) :
    sh x D o j * ff x D = gbz x j * ff (j : ℚ) o * ff (x - j) (D - o) := by
  unfold sh
  by_cases hj : j < o
  · have : (j : ℤ) - o < 0 := by omega
    simp [gbz, this, ff_nat_lt j o hj]
  · push_neg at hj
    have hnn : ¬ ((j : ℤ) - o < 0) := by omega
    have ht : ((j : ℤ) - o).toNat = j - o := by omega
    have hj0 : ¬ ((j : ℤ) < 0) := by omega
    simp only [gbz, hnn, hj0, if_false, ht, Int.toNat_natCast, gb]
    have hf := ff_nat_fact j o hj
    have hadd1 := ff_add x D (j - o)
    have hadd2 := ff_add x j (D - o)
    have e : D + (j - o) = j + (D - o) := by omega
    rw [e] at hadd1
    have key : ff x D * ff (x - D) (j - o) = ff x j * ff (x - j) (D - o) := by rw [← hadd1, hadd2]
    have hjo : ((j - o).factorial : ℚ) ≠ 0 := by positivity
    have hfo : ff (j : ℚ) o ≠ 0 := by
      intro h0; rw [h0, zero_mul] at hf; exact (Nat.factorial_ne_zero j) (by exact_mod_cast hf.symm)
    rw [← hf]; field_simp
    linear_combination key

/-- Shift lemma, `D < o`. -/
lemma sh_gt_mul (x : ℚ) (D o j : ℕ) (h : D < o) :
    sh x D o j * ff x D * ff (x - j + o - D) (o - D) = gbz x j * ff (j : ℚ) o := by
  unfold sh
  by_cases hj : j < o
  · have : (j : ℤ) - o < 0 := by omega
    simp [gbz, this, ff_nat_lt j o hj]
  · push_neg at hj
    have hnn : ¬ ((j : ℤ) - o < 0) := by omega
    have ht : ((j : ℤ) - o).toNat = j - o := by omega
    have hj0 : ¬ ((j : ℤ) < 0) := by omega
    simp only [gbz, hnn, hj0, if_false, ht, Int.toNat_natCast, gb]
    have hf := ff_nat_fact j o hj
    have hadd1 := ff_add x D (j - o)
    have hadd2 := ff_add x (D + (j - o)) (o - D)
    have e : D + (j - o) + (o - D) = j := by omega
    rw [e] at hadd2
    have ex : x - ((D + (j - o) : ℕ) : ℚ) = x - j + o - D := by
      rw [Nat.cast_add, Nat.cast_sub hj]; ring
    rw [ex, hadd1] at hadd2
    have hjo : ((j - o).factorial : ℚ) ≠ 0 := by positivity
    have hfo : ff (j : ℚ) o ≠ 0 := by
      intro h0; rw [h0, zero_mul] at hf; exact (Nat.factorial_ne_zero j) (by exact_mod_cast hf.symm)
    rw [← hf]; field_simp
    linear_combination (-1 : ℚ) * hadd2


lemma sh_le_div (x : ℚ) (D o e j : ℕ) (he : D = o + e) (hx : ff x D ≠ 0) :
    sh x D o j = gbz x j * (ff (j : ℚ) o * ff (x - j) e) / ff x D := by
  have h := sh_le_mul x D o j (by omega)
  rw [show D - o = e by omega] at h
  rw [eq_div_iff hx]; linear_combination h

lemma sh_gt1_div (x : ℚ) (D o j : ℕ) (ho : o = D + 1) (hx : ff x D ≠ 0) (hy : x - (j : ℚ) + 1 ≠ 0) :
    sh x D o j = gbz x j * ff (j : ℚ) o / (ff x D * (x - j + 1)) := by
  have h := sh_gt_mul x D o j (by omega)
  rw [show o - D = 1 by omega, ff_one] at h
  have e : x - (j : ℚ) + (o : ℚ) - (D : ℚ) = x - j + 1 := by rw [ho]; push_cast; ring
  rw [e] at h
  rw [eq_div_iff (mul_ne_zero hx hy)]; linear_combination h

/-- Solution 1: the G1 loss-1 coefficient of `X^a Y^(K-1-a)`, as the finite sum of
    shifted generalised binomials that item 311's term list reduces to (transcribed by
    the generator from that list; see the header). -/
noncomputable def Lsum1 (t : ℚ) (K a : ℕ) (b0 b1 : ℚ) : ℚ :=
  ((2*(K:ℚ) - 1)*(t + 1)^2*(t + 3)/(96*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 (2*a)) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)*(t + 3)/96) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 (2*a)) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)*(t + 3)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 (2*a)) * ((1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)*(t + 3)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 (2*a)) * ((-1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/48) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 2 1 (2*a)) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/(24*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 (2*a)) * ((2:ℚ) * (-1:ℚ)^2 * sh b1 2 2 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/(24*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 (2*a)) * ((-2:ℚ) * (-1:ℚ)^2 * sh b1 2 2 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t - 1)*(t + 1)/96) * ((1:ℚ) * (-1:ℚ)^2 * sh b0 2 2 (2*a)) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/48) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 (2*a)) * ((1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/48) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 (2*a)) * ((-1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/48) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 (2*a)) * ((1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/(24*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 (2*a)) * ((1:ℚ) * (-1:ℚ)^3 * sh b1 2 3 (K - a) + (1:ℚ) * (-1:ℚ)^2 * sh b1 2 2 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/(24*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 (2*a)) * ((1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/48) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 (2*a)) * ((-1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/(24*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 (2*a)) * ((1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/(24*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 (2*a)) * ((1:ℚ) * (-1:ℚ)^3 * sh b1 2 3 (K - a) + (1:ℚ) * (-1:ℚ)^2 * sh b1 2 2 (K - a))
  + (-(2*(K:ℚ) - 1)*(t - 3)*(t + 1)*(t + 3)/(96*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 2 0 (2*a)) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))

/- **Solution 1: G1 = the Ward loss-1 layer.**  `R = -(L/top)/CY - (K-a) rho^2` equals
    Codex's Ward-derived rational function, for every spin index K, every `a < K`, and
    every rational t off the stated degenerate set. -/
set_option maxHeartbeats 20000000 in
theorem g1_ward_sol1 (t : ℚ) (K a : ℕ) (hKa : a + 1 ≤ K) (b0 b1 : ℚ)
    (hb0 : b0 = ((2*(K:ℚ) - 1)*(t - 1)/(t + 3))) (hb1 : b1 = (2*(2*(K:ℚ) - 1)/(t + 3)))
    (ht0 : t - 1 ≠ 0)
    (ht1 : t + 1 ≠ 0)
    (ht2 : t + 3 ≠ 0)
    (hb0z : b0 ≠ 0) (hb0o : b0 - 1 ≠ 0) (hb1z : b1 ≠ 0) (hb1o : b1 - 1 ≠ 0)
    (hb1j : b1 - ((K:ℚ) - a) + 1 ≠ 0)
    (htop : gbz b0 ((2*a : ℕ) : ℤ) * gbz b1 ((K - a : ℕ) : ℤ) ≠ 0) :
    (-(Lsum1 t K a b0 b1 / (gbz b0 ((2*a : ℕ) : ℤ) * gbz b1 ((K - a : ℕ) : ℤ))) / ((t + 1)*(t + 3)^2/(8*(t - 1)))
      - ((K:ℚ) - a) * (2 / ((t - 1) * (t + 1)))) * (6*(t - 1)*(t + 1)) = ((-(K:ℚ) + (a:ℚ))*(-2*(K:ℚ)*t^2 + 2*(K:ℚ) + 2*(a:ℚ)*t^2 + 4*(a:ℚ)*t + 2*(a:ℚ) + t^2 + 11)) := by
  unfold Lsum1
  rw [sh_le_div b0 0 0 0 (2*a) (by norm_num) (by simp)]
  rw [sh_le_div b0 1 1 0 (2*a) (by norm_num) (by rw [ff_one]; exact hb0z)]
  rw [sh_le_div b0 2 0 2 (2*a) (by norm_num) (by rw [ff_two]; exact mul_ne_zero hb0z hb0o)]
  rw [sh_le_div b0 2 1 1 (2*a) (by norm_num) (by rw [ff_two]; exact mul_ne_zero hb0z hb0o)]
  rw [sh_le_div b0 2 2 0 (2*a) (by norm_num) (by rw [ff_two]; exact mul_ne_zero hb0z hb0o)]
  rw [sh_gt1_div b1 0 1 (K - a) (by norm_num) (by simp) (by rw [Nat.cast_sub (by omega)]; exact hb1j)]
  rw [sh_gt1_div b1 1 2 (K - a) (by norm_num) (by rw [ff_one]; exact hb1z) (by rw [Nat.cast_sub (by omega)]; exact hb1j)]
  rw [sh_le_div b1 2 2 0 (K - a) (by norm_num) (by rw [ff_two]; exact mul_ne_zero hb1z hb1o)]
  rw [sh_gt1_div b1 2 3 (K - a) (by norm_num) (by rw [ff_two]; exact mul_ne_zero hb1z hb1o) (by rw [Nat.cast_sub (by omega)]; exact hb1j)]
  have hG0 : gbz b0 ((2*a : ℕ) : ℤ) ≠ 0 := left_ne_zero_of_mul htop
  have hG1 : gbz b1 ((K - a : ℕ) : ℤ) ≠ 0 := right_ne_zero_of_mul htop
  generalize gbz b0 ((2*a : ℕ) : ℤ) = G0 at hG0 ⊢
  generalize gbz b1 ((K - a : ℕ) : ℤ) = G1 at hG1 ⊢
  simp only [ff_zero, ff_one, ff_two, ff_three, Nat.cast_ofNat, Nat.cast_mul, Nat.cast_sub (show a ≤ K by omega)]
  field_simp
  subst hb0 hb1
  field_simp
  ring

/-- Solution 1, in the record's ratio form `R_G1 = Num/Den` (Codex's law). -/
theorem g1_ward_sol1_ratio (t : ℚ) (K a : ℕ) (hKa : a + 1 ≤ K) (b0 b1 : ℚ)
    (hb0 : b0 = ((2*(K:ℚ) - 1)*(t - 1)/(t + 3))) (hb1 : b1 = (2*(2*(K:ℚ) - 1)/(t + 3)))
    (ht0 : t - 1 ≠ 0)
    (ht1 : t + 1 ≠ 0)
    (ht2 : t + 3 ≠ 0)
    (hb0z : b0 ≠ 0) (hb0o : b0 - 1 ≠ 0) (hb1z : b1 ≠ 0) (hb1o : b1 - 1 ≠ 0)
    (hb1j : b1 - ((K:ℚ) - a) + 1 ≠ 0)
    (htop : gbz b0 ((2*a : ℕ) : ℤ) * gbz b1 ((K - a : ℕ) : ℤ) ≠ 0)
    (hden : (6*(t - 1)*(t + 1)) ≠ 0) :
    -(Lsum1 t K a b0 b1 / (gbz b0 ((2*a : ℕ) : ℤ) * gbz b1 ((K - a : ℕ) : ℤ))) / ((t + 1)*(t + 3)^2/(8*(t - 1)))
      - ((K:ℚ) - a) * (2 / ((t - 1) * (t + 1))) = ((-(K:ℚ) + (a:ℚ))*(-2*(K:ℚ)*t^2 + 2*(K:ℚ) + 2*(a:ℚ)*t^2 + 4*(a:ℚ)*t + 2*(a:ℚ) + t^2 + 11)) / (6*(t - 1)*(t + 1)) :=
  (eq_div_iff hden).mpr (g1_ward_sol1 t K a hKa b0 b1 hb0 hb1 ht0 ht1 ht2 hb0z hb0o hb1z hb1o hb1j htop)


/-- Solution 2: the G1 loss-1 coefficient of `X^a Y^(K-1-a)`, as the finite sum of
    shifted generalised binomials that item 311's term list reduces to (transcribed by
    the generator from that list; see the header). -/
noncomputable def Lsum2 (t : ℚ) (K a : ℕ) (b0 b1 : ℚ) : ℚ :=
  ((2*(K:ℚ) - 1)*(t + 1)^2*(t + 5)/(24*(t - 3)^2*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 a) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)*(t + 5)/(24*(t - 3)^2)) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)*(t + 5)/(24*(t - 3)^2)) * ((-1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)*(t + 5)/(48*(t - 3)*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 a) * ((1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)*(t + 5)/(48*(t - 3)*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 a) * ((-1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/(24*(t - 3))) * ((2:ℚ) * (-1:ℚ)^1 * sh b0 2 1 a) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/(24*(t - 3))) * ((-2:ℚ) * (-1:ℚ)^1 * sh b0 2 1 a) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 a) * ((2:ℚ) * (-1:ℚ)^2 * sh b1 2 2 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 a) * ((-2:ℚ) * (-1:ℚ)^2 * sh b1 2 2 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t - 1)*(t + 1)/(12*(t - 3)^2)) * ((1:ℚ) * (-1:ℚ)^2 * sh b0 2 2 a + (1:ℚ) * (-1:ℚ)^1 * sh b0 2 1 a) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t - 1)*(t + 1)/(12*(t - 3)^2)) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/(24*(t - 3))) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/(24*(t - 3))) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((-1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t - 1)*(t + 1)/(12*(t - 3)^2)) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t - 1)*(t + 1)/(12*(t - 3)^2)) * ((1:ℚ) * (-1:ℚ)^2 * sh b0 2 2 a + (1:ℚ) * (-1:ℚ)^1 * sh b0 2 1 a) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/(24*(t - 3))) * ((-1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/(24*(t - 3))) * ((-1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((-1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/(24*(t - 3))) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/(24*(t - 3))) * ((-1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 a) * ((1:ℚ) * (-1:ℚ)^3 * sh b1 2 3 (K - a) + (1:ℚ) * (-1:ℚ)^2 * sh b1 2 2 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 a) * ((1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/(24*(t - 3))) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((-1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/(24*(t - 3))) * ((-1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((-1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 a) * ((1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 a) * ((1:ℚ) * (-1:ℚ)^3 * sh b1 2 3 (K - a) + (1:ℚ) * (-1:ℚ)^2 * sh b1 2 2 (K - a))
  + (-(2*(K:ℚ) - 1)*(t + 1)*(t + 5)*(3*t - 5)/(48*(t - 3)^2*(t - 1))) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 2 1 a + (1:ℚ) * (-1:ℚ)^0 * sh b0 2 0 a) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + (-(2*(K:ℚ) - 1)*(t + 1)*(t + 5)*(3*t - 5)/(48*(t - 3)^2*(t - 1))) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 2 1 a + (1:ℚ) * (-1:ℚ)^0 * sh b0 2 0 a) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))

/- **Solution 2: G1 = the Ward loss-1 layer.**  `R = -(L/top)/CY - (K-a) rho^2` equals
    Codex's Ward-derived rational function, for every spin index K, every `a < K`, and
    every rational t off the stated degenerate set. -/
set_option maxHeartbeats 20000000 in
theorem g1_ward_sol2 (t : ℚ) (K a : ℕ) (hKa : a + 1 ≤ K) (b0 b1 : ℚ)
    (hb0 : b0 = (2*(2*(K:ℚ) - 1)*(t - 1)/(t + 5))) (hb1 : b1 = (-(2*(K:ℚ) - 1)*(t - 3)/(t + 5)))
    (ht0 : t - 1 ≠ 0)
    (ht1 : t + 1 ≠ 0)
    (ht2 : t - 3 ≠ 0)
    (ht3 : t + 5 ≠ 0)
    (hb0z : b0 ≠ 0) (hb0o : b0 - 1 ≠ 0) (hb1z : b1 ≠ 0) (hb1o : b1 - 1 ≠ 0)
    (hb1j : b1 - ((K:ℚ) - a) + 1 ≠ 0)
    (htop : gbz b0 a * gbz b1 ((K - a : ℕ) : ℤ) ≠ 0) :
    (-(Lsum2 t K a b0 b1 / (gbz b0 a * gbz b1 ((K - a : ℕ) : ℤ))) / ((t + 1)*(t + 5)^2/(2*(t - 3)^2*(t - 1)))
      - ((K:ℚ) - a) * (2 / ((t - 1) * (t + 1)))) * (24*(t - 1)*(t + 1)*(-3*(K:ℚ)*t + (K:ℚ) + (a:ℚ)*t + 5*(a:ℚ) + 2*t + 2)) = ((-(K:ℚ) + (a:ℚ))*(12*(K:ℚ)^2*t^3 - 4*(K:ℚ)^2*t^2 - 12*(K:ℚ)^2*t + 4*(K:ℚ)^2 - 8*(K:ℚ)*(a:ℚ)*t^3 - 40*(K:ℚ)*(a:ℚ)*t^2 + 8*(K:ℚ)*(a:ℚ)*t + 40*(K:ℚ)*(a:ℚ) - 10*(K:ℚ)*t^3 - 18*(K:ℚ)*t^2 - 134*(K:ℚ)*t + 66*(K:ℚ) + 2*(a:ℚ)^2*t^3 + 14*(a:ℚ)^2*t^2 + 22*(a:ℚ)^2*t + 10*(a:ℚ)^2 + 5*(a:ℚ)*t^3 + 27*(a:ℚ)*t^2 + 55*(a:ℚ)*t + 225*(a:ℚ) + 2*t^3 + 10*t^2 + 94*t + 86)) := by
  unfold Lsum2
  rw [sh_le_div b0 0 0 0 a (by norm_num) (by simp)]
  rw [sh_le_div b0 1 1 0 a (by norm_num) (by rw [ff_one]; exact hb0z)]
  rw [sh_le_div b0 2 0 2 a (by norm_num) (by rw [ff_two]; exact mul_ne_zero hb0z hb0o)]
  rw [sh_le_div b0 2 1 1 a (by norm_num) (by rw [ff_two]; exact mul_ne_zero hb0z hb0o)]
  rw [sh_le_div b0 2 2 0 a (by norm_num) (by rw [ff_two]; exact mul_ne_zero hb0z hb0o)]
  rw [sh_gt1_div b1 0 1 (K - a) (by norm_num) (by simp) (by rw [Nat.cast_sub (by omega)]; exact hb1j)]
  rw [sh_gt1_div b1 1 2 (K - a) (by norm_num) (by rw [ff_one]; exact hb1z) (by rw [Nat.cast_sub (by omega)]; exact hb1j)]
  rw [sh_le_div b1 2 2 0 (K - a) (by norm_num) (by rw [ff_two]; exact mul_ne_zero hb1z hb1o)]
  rw [sh_gt1_div b1 2 3 (K - a) (by norm_num) (by rw [ff_two]; exact mul_ne_zero hb1z hb1o) (by rw [Nat.cast_sub (by omega)]; exact hb1j)]
  have hG0 : gbz b0 a ≠ 0 := left_ne_zero_of_mul htop
  have hG1 : gbz b1 ((K - a : ℕ) : ℤ) ≠ 0 := right_ne_zero_of_mul htop
  generalize gbz b0 a = G0 at hG0 ⊢
  generalize gbz b1 ((K - a : ℕ) : ℤ) = G1 at hG1 ⊢
  simp only [ff_zero, ff_one, ff_two, ff_three, Nat.cast_ofNat, Nat.cast_mul, Nat.cast_sub (show a ≤ K by omega)]
  field_simp
  subst hb0 hb1
  field_simp
  ring

/-- Solution 2, in the record's ratio form `R_G1 = Num/Den` (Codex's law). -/
theorem g1_ward_sol2_ratio (t : ℚ) (K a : ℕ) (hKa : a + 1 ≤ K) (b0 b1 : ℚ)
    (hb0 : b0 = (2*(2*(K:ℚ) - 1)*(t - 1)/(t + 5))) (hb1 : b1 = (-(2*(K:ℚ) - 1)*(t - 3)/(t + 5)))
    (ht0 : t - 1 ≠ 0)
    (ht1 : t + 1 ≠ 0)
    (ht2 : t - 3 ≠ 0)
    (ht3 : t + 5 ≠ 0)
    (hb0z : b0 ≠ 0) (hb0o : b0 - 1 ≠ 0) (hb1z : b1 ≠ 0) (hb1o : b1 - 1 ≠ 0)
    (hb1j : b1 - ((K:ℚ) - a) + 1 ≠ 0)
    (htop : gbz b0 a * gbz b1 ((K - a : ℕ) : ℤ) ≠ 0)
    (hden : (24*(t - 1)*(t + 1)*(-3*(K:ℚ)*t + (K:ℚ) + (a:ℚ)*t + 5*(a:ℚ) + 2*t + 2)) ≠ 0) :
    -(Lsum2 t K a b0 b1 / (gbz b0 a * gbz b1 ((K - a : ℕ) : ℤ))) / ((t + 1)*(t + 5)^2/(2*(t - 3)^2*(t - 1)))
      - ((K:ℚ) - a) * (2 / ((t - 1) * (t + 1))) = ((-(K:ℚ) + (a:ℚ))*(12*(K:ℚ)^2*t^3 - 4*(K:ℚ)^2*t^2 - 12*(K:ℚ)^2*t + 4*(K:ℚ)^2 - 8*(K:ℚ)*(a:ℚ)*t^3 - 40*(K:ℚ)*(a:ℚ)*t^2 + 8*(K:ℚ)*(a:ℚ)*t + 40*(K:ℚ)*(a:ℚ) - 10*(K:ℚ)*t^3 - 18*(K:ℚ)*t^2 - 134*(K:ℚ)*t + 66*(K:ℚ) + 2*(a:ℚ)^2*t^3 + 14*(a:ℚ)^2*t^2 + 22*(a:ℚ)^2*t + 10*(a:ℚ)^2 + 5*(a:ℚ)*t^3 + 27*(a:ℚ)*t^2 + 55*(a:ℚ)*t + 225*(a:ℚ) + 2*t^3 + 10*t^2 + 94*t + 86)) / (24*(t - 1)*(t + 1)*(-3*(K:ℚ)*t + (K:ℚ) + (a:ℚ)*t + 5*(a:ℚ) + 2*t + 2)) :=
  (eq_div_iff hden).mpr (g1_ward_sol2 t K a hKa b0 b1 hb0 hb1 ht0 ht1 ht2 ht3 hb0z hb0o hb1z hb1o hb1j htop)


/-- Solution 3: the G1 loss-1 coefficient of `X^a Y^(K-1-a)`, as the finite sum of
    shifted generalised binomials that item 311's term list reduces to (transcribed by
    the generator from that list; see the header). -/
noncomputable def Lsum3 (t : ℚ) (K a : ℕ) (b0 b1 : ℚ) : ℚ :=
  ((2*(K:ℚ) - 1)*(t - 7)*(t - 5)*(t + 1)/(24*(t - 3)^2*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 a) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t - 5)*(t + 1)/(24*(t - 3)*(t - 1))) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t - 5)*(t + 1)/(24*(t - 3)*(t - 1))) * ((-1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t - 5)*(t + 1)/(24*(t - 3)*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 a) * ((1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t - 5)*(t + 1)/(24*(t - 3)*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 a) * ((-1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((2:ℚ) * (-1:ℚ)^1 * sh b0 2 1 a) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((-2:ℚ) * (-1:ℚ)^1 * sh b0 2 1 a) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 a) * ((2:ℚ) * (-1:ℚ)^2 * sh b1 2 2 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 a) * ((-2:ℚ) * (-1:ℚ)^2 * sh b1 2 2 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^2 * sh b0 2 2 a + (1:ℚ) * (-1:ℚ)^1 * sh b0 2 1 a) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((-1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^2 * sh b0 2 2 a + (1:ℚ) * (-1:ℚ)^1 * sh b0 2 1 a) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((-1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((-1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((-1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((-1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 a) * ((1:ℚ) * (-1:ℚ)^3 * sh b1 2 3 (K - a) + (1:ℚ) * (-1:ℚ)^2 * sh b1 2 2 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 a) * ((1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((-1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((-1:ℚ) * (-1:ℚ)^1 * sh b0 1 1 a) * ((-1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + (-(2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 a) * ((1:ℚ) * (-1:ℚ)^2 * sh b1 1 2 (K - a))
  + ((2*(K:ℚ) - 1)^2*(t + 1)/(48*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 a) * ((1:ℚ) * (-1:ℚ)^3 * sh b1 2 3 (K - a) + (1:ℚ) * (-1:ℚ)^2 * sh b1 2 2 (K - a))
  + ((2*(K:ℚ) - 1)*(t - 5)*(t + 1)*(3*t - 5)/(24*(t - 3)^2*(t - 1))) * ((1:ℚ) * (-1:ℚ)^0 * sh b0 0 0 a) * ((1:ℚ) * (-1:ℚ)^1 * sh b1 0 1 (K - a))

/- **Solution 3: G1 = the Ward loss-1 layer.**  `R = -(L/top)/CY - (K-a) rho^2` equals
    Codex's Ward-derived rational function, for every spin index K, every `a < K`, and
    every rational t off the stated degenerate set. -/
set_option maxHeartbeats 20000000 in
theorem g1_ward_sol3 (t : ℚ) (K a : ℕ) (hKa : a + 1 ≤ K) (b0 b1 : ℚ)
    (hb0 : b0 = ((2*(K:ℚ) - 1)*(t - 3)/(2*(t - 5)))) (hb1 : b1 = ((2*(K:ℚ) - 1)*(t - 3)/(2*(t - 5))))
    (ht0 : t - 1 ≠ 0)
    (ht1 : t + 1 ≠ 0)
    (ht2 : t - 3 ≠ 0)
    (ht3 : t - 5 ≠ 0)
    (hb0z : b0 ≠ 0) (hb0o : b0 - 1 ≠ 0) (hb1z : b1 ≠ 0) (hb1o : b1 - 1 ≠ 0)
    (hb1j : b1 - ((K:ℚ) - a) + 1 ≠ 0)
    (htop : gbz b0 a * gbz b1 ((K - a : ℕ) : ℤ) ≠ 0) :
    (-(Lsum3 t K a b0 b1 / (gbz b0 a * gbz b1 ((K - a : ℕ) : ℤ))) / (2*(t - 5)^2*(t + 1)/((t - 3)^2*(t - 1)))
      - ((K:ℚ) - a) * (2 / ((t - 1) * (t + 1)))) * (6*(t - 1)*(t + 1)*(4*(K:ℚ) + 2*(a:ℚ)*t - 10*(a:ℚ) + t - 7)) = ((-(K:ℚ) + (a:ℚ))*(-4*(K:ℚ)^2*t^2 + 4*(K:ℚ)^2 - 4*(K:ℚ)*(a:ℚ)*t^3 + 20*(K:ℚ)*(a:ℚ)*t^2 + 4*(K:ℚ)*(a:ℚ)*t - 20*(K:ℚ)*(a:ℚ) - 2*(K:ℚ)*t^3 + 12*(K:ℚ)*t^2 + 2*(K:ℚ)*t + 36*(K:ℚ) + 4*(a:ℚ)^2*t^3 - 20*(a:ℚ)^2*t^2 - 4*(a:ℚ)^2*t + 20*(a:ℚ)^2 + 4*(a:ℚ)*t^3 - 20*(a:ℚ)*t^2 + 20*(a:ℚ)*t - 100*(a:ℚ) + t^3 - 5*t^2 + 11*t - 79)) := by
  unfold Lsum3
  rw [sh_le_div b0 0 0 0 a (by norm_num) (by simp)]
  rw [sh_le_div b0 1 1 0 a (by norm_num) (by rw [ff_one]; exact hb0z)]
  rw [sh_le_div b0 2 1 1 a (by norm_num) (by rw [ff_two]; exact mul_ne_zero hb0z hb0o)]
  rw [sh_le_div b0 2 2 0 a (by norm_num) (by rw [ff_two]; exact mul_ne_zero hb0z hb0o)]
  rw [sh_gt1_div b1 0 1 (K - a) (by norm_num) (by simp) (by rw [Nat.cast_sub (by omega)]; exact hb1j)]
  rw [sh_gt1_div b1 1 2 (K - a) (by norm_num) (by rw [ff_one]; exact hb1z) (by rw [Nat.cast_sub (by omega)]; exact hb1j)]
  rw [sh_le_div b1 2 2 0 (K - a) (by norm_num) (by rw [ff_two]; exact mul_ne_zero hb1z hb1o)]
  rw [sh_gt1_div b1 2 3 (K - a) (by norm_num) (by rw [ff_two]; exact mul_ne_zero hb1z hb1o) (by rw [Nat.cast_sub (by omega)]; exact hb1j)]
  have hG0 : gbz b0 a ≠ 0 := left_ne_zero_of_mul htop
  have hG1 : gbz b1 ((K - a : ℕ) : ℤ) ≠ 0 := right_ne_zero_of_mul htop
  generalize gbz b0 a = G0 at hG0 ⊢
  generalize gbz b1 ((K - a : ℕ) : ℤ) = G1 at hG1 ⊢
  simp only [ff_zero, ff_one, ff_two, ff_three, Nat.cast_ofNat, Nat.cast_mul, Nat.cast_sub (show a ≤ K by omega)]
  field_simp
  subst hb0 hb1
  field_simp
  ring

/-- Solution 3, in the record's ratio form `R_G1 = Num/Den` (Codex's law). -/
theorem g1_ward_sol3_ratio (t : ℚ) (K a : ℕ) (hKa : a + 1 ≤ K) (b0 b1 : ℚ)
    (hb0 : b0 = ((2*(K:ℚ) - 1)*(t - 3)/(2*(t - 5)))) (hb1 : b1 = ((2*(K:ℚ) - 1)*(t - 3)/(2*(t - 5))))
    (ht0 : t - 1 ≠ 0)
    (ht1 : t + 1 ≠ 0)
    (ht2 : t - 3 ≠ 0)
    (ht3 : t - 5 ≠ 0)
    (hb0z : b0 ≠ 0) (hb0o : b0 - 1 ≠ 0) (hb1z : b1 ≠ 0) (hb1o : b1 - 1 ≠ 0)
    (hb1j : b1 - ((K:ℚ) - a) + 1 ≠ 0)
    (htop : gbz b0 a * gbz b1 ((K - a : ℕ) : ℤ) ≠ 0)
    (hden : (6*(t - 1)*(t + 1)*(4*(K:ℚ) + 2*(a:ℚ)*t - 10*(a:ℚ) + t - 7)) ≠ 0) :
    -(Lsum3 t K a b0 b1 / (gbz b0 a * gbz b1 ((K - a : ℕ) : ℤ))) / (2*(t - 5)^2*(t + 1)/((t - 3)^2*(t - 1)))
      - ((K:ℚ) - a) * (2 / ((t - 1) * (t + 1))) = ((-(K:ℚ) + (a:ℚ))*(-4*(K:ℚ)^2*t^2 + 4*(K:ℚ)^2 - 4*(K:ℚ)*(a:ℚ)*t^3 + 20*(K:ℚ)*(a:ℚ)*t^2 + 4*(K:ℚ)*(a:ℚ)*t - 20*(K:ℚ)*(a:ℚ) - 2*(K:ℚ)*t^3 + 12*(K:ℚ)*t^2 + 2*(K:ℚ)*t + 36*(K:ℚ) + 4*(a:ℚ)^2*t^3 - 20*(a:ℚ)^2*t^2 - 4*(a:ℚ)^2*t + 20*(a:ℚ)^2 + 4*(a:ℚ)*t^3 - 20*(a:ℚ)*t^2 + 20*(a:ℚ)*t - 100*(a:ℚ) + t^3 - 5*t^2 + 11*t - 79)) / (6*(t - 1)*(t + 1)*(4*(K:ℚ) + 2*(a:ℚ)*t - 10*(a:ℚ) + t - 7)) :=
  (eq_div_iff hden).mpr (g1_ward_sol3 t K a hKa b0 b1 hb0 hb1 ht0 ht1 ht2 ht3 hb0z hb0o hb1z hb1o hb1j htop)


/-! ## Negative controls and a non-vacuity witness (Solution 3 at t = 9/4, K = 3, a = 1) -/

/-- Exact values at the point: `L = 455/88`, `top = -1575/21296`. -/
theorem control_sol3_values :
    Lsum3 (9/4) 3 1 (15/22) (15/22) = 455/88 ∧
    gbz (15/22) (1 : ℕ) * gbz (15/22) ((3 - 1 : ℕ) : ℤ) = -1575/21296 := by
  constructor
  · norm_num [Lsum3, sh, gbz, gb, ff_succ, ff_one, ff_two, ff_three, Nat.factorial]
  · rw [show ((3 - 1 : ℕ) : ℤ) = ((2 : ℕ) : ℤ) from rfl, gbz_nat, gbz_nat]
    norm_num [gb, ff_one, ff_two, Nat.factorial]

/-- NON-VACUITY: every hypothesis of `g1_ward_sol3` holds at the point, so the theorem
    is not vacuous there (both sides equal 21/32). -/
theorem control_sol3_instance :
    (-(Lsum3 (9/4) 3 1 (15/22) (15/22) / (gbz (15/22) (1 : ℕ) * gbz (15/22) ((3 - 1 : ℕ) : ℤ))) / (2*((9/4:ℚ) - 5)^2*((9/4:ℚ) + 1)/(((9/4:ℚ) - 3)^2*((9/4:ℚ) - 1)))
      - (((3:ℕ):ℚ) - ((1:ℕ):ℚ)) * (2 / (((9/4:ℚ) - 1) * ((9/4:ℚ) + 1)))) * (6*((9/4:ℚ) - 1)*((9/4:ℚ) + 1)*(4*((3:ℕ):ℚ) + 2*((1:ℕ):ℚ)*(9/4:ℚ) - 10*((1:ℕ):ℚ) + (9/4:ℚ) - 7)) = ((-((3:ℕ):ℚ) + ((1:ℕ):ℚ))*(-4*((3:ℕ):ℚ)^2*(9/4:ℚ)^2 + 4*((3:ℕ):ℚ)^2 - 4*((3:ℕ):ℚ)*((1:ℕ):ℚ)*(9/4:ℚ)^3 + 20*((3:ℕ):ℚ)*((1:ℕ):ℚ)*(9/4:ℚ)^2 + 4*((3:ℕ):ℚ)*((1:ℕ):ℚ)*(9/4:ℚ) - 20*((3:ℕ):ℚ)*((1:ℕ):ℚ) - 2*((3:ℕ):ℚ)*(9/4:ℚ)^3 + 12*((3:ℕ):ℚ)*(9/4:ℚ)^2 + 2*((3:ℕ):ℚ)*(9/4:ℚ) + 36*((3:ℕ):ℚ) + 4*((1:ℕ):ℚ)^2*(9/4:ℚ)^3 - 20*((1:ℕ):ℚ)^2*(9/4:ℚ)^2 - 4*((1:ℕ):ℚ)^2*(9/4:ℚ) + 20*((1:ℕ):ℚ)^2 + 4*((1:ℕ):ℚ)*(9/4:ℚ)^3 - 20*((1:ℕ):ℚ)*(9/4:ℚ)^2 + 20*((1:ℕ):ℚ)*(9/4:ℚ) - 100*((1:ℕ):ℚ) + (9/4:ℚ)^3 - 5*(9/4:ℚ)^2 + 11*(9/4:ℚ) - 79)) :=
  g1_ward_sol3 (9/4) 3 1 (by norm_num) (15/22) (15/22) (by norm_num) (by norm_num)
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by rw [control_sol3_values.2]; norm_num)

/-- NEGATIVE CONTROL: the normalisation CY is load-bearing.  With `CY` replaced by
    `(11/10) CY` the cross-multiplied identity fails at the point (-567/176 vs 21/32). -/
theorem control_sol3_CY_load_bearing :
    (-(Lsum3 (9/4) 3 1 (15/22) (15/22) / (gbz (15/22) (1 : ℕ) * gbz (15/22) ((3 - 1 : ℕ) : ℤ)))
        / ((11/10) * (2*((9/4:ℚ) - 5)^2*((9/4:ℚ) + 1)/(((9/4:ℚ) - 3)^2*((9/4:ℚ) - 1))))
      - (((3:ℕ):ℚ) - ((1:ℕ):ℚ)) * (2 / (((9/4:ℚ) - 1) * ((9/4:ℚ) + 1)))) * (6*((9/4:ℚ) - 1)*((9/4:ℚ) + 1)*(4*((3:ℕ):ℚ) + 2*((1:ℕ):ℚ)*(9/4:ℚ) - 10*((1:ℕ):ℚ) + (9/4:ℚ) - 7)) ≠ ((-((3:ℕ):ℚ) + ((1:ℕ):ℚ))*(-4*((3:ℕ):ℚ)^2*(9/4:ℚ)^2 + 4*((3:ℕ):ℚ)^2 - 4*((3:ℕ):ℚ)*((1:ℕ):ℚ)*(9/4:ℚ)^3 + 20*((3:ℕ):ℚ)*((1:ℕ):ℚ)*(9/4:ℚ)^2 + 4*((3:ℕ):ℚ)*((1:ℕ):ℚ)*(9/4:ℚ) - 20*((3:ℕ):ℚ)*((1:ℕ):ℚ) - 2*((3:ℕ):ℚ)*(9/4:ℚ)^3 + 12*((3:ℕ):ℚ)*(9/4:ℚ)^2 + 2*((3:ℕ):ℚ)*(9/4:ℚ) + 36*((3:ℕ):ℚ) + 4*((1:ℕ):ℚ)^2*(9/4:ℚ)^3 - 20*((1:ℕ):ℚ)^2*(9/4:ℚ)^2 - 4*((1:ℕ):ℚ)^2*(9/4:ℚ) + 20*((1:ℕ):ℚ)^2 + 4*((1:ℕ):ℚ)*(9/4:ℚ)^3 - 20*((1:ℕ):ℚ)*(9/4:ℚ)^2 + 20*((1:ℕ):ℚ)*(9/4:ℚ) - 100*((1:ℕ):ℚ) + (9/4:ℚ)^3 - 5*(9/4:ℚ)^2 + 11*(9/4:ℚ) - 79)) := by
  rw [control_sol3_values.1, control_sol3_values.2]; norm_num

end IntegrableBoundary.WardLossOne
