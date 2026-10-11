/-
  Theorem 3 (spectral equivalence at t = 2): the frozen lattice and the indicial
  non-vanishing under the no-resonance hypothesis (H).

  THE CLAIM (the opers paper, Appendix A, preamble and Step 4; the record's items
  on the analytic proof of Theorem 3 and its review).  At t = 2, with s = 5/2,
  P(T) = (T^2 - l0^2)(T^2 - l1^2) and the four exponents
      theta_0 = s + l0, theta_1 = s - l0, theta_2 = s + l1, theta_3 = s - l1,
  * m(nu) = 3^(nu/3) / Gamma(1/3 - nu/3) vanishes exactly on 1 + 3 Z_{>=0};
  * the sixth-order form has indicial polynomial I(mu) = P(s - mu)(mu - 1)(mu - 4);
  * under (H): theta_i - theta_j not in Z for i /= j, and theta_j not in 1 + 3Z,
      I(theta_j + 3k) /= 0 for k >= 1, and I(1 + 3k) /= 0 for k >= 2.

  CONVENTIONS.  Gamma is Mathlib's Complex.Gamma, which is 0 at its poles, so
  (Gamma z)^-1 is the entire function 1/Gamma; 3^z is the principal complex power.
  l0, l1 are arbitrary complex numbers.  Z is the integers, as casts into C.

  WHAT IS PROVED (no sorry, no axioms beyond Lean/Mathlib's):
  * `m_eq_zero_iff`: m(mu) = 0  <->  exists n : N, mu = 1 + 3n;
  * `P_shift_prod`: P(s - mu) = prod_j (mu - theta_j), so the roots of
    P(s - mu) are exactly the four exponents;
  * `indicial_exponent_ne_zero`: under (H), I(theta_j + 3k) /= 0 for EVERY integer
    k /= 0 (the paper needs k >= 1);
  * `indicial_frozen_ne_zero`: under (H), I(1 + 3k) /= 0 for every integer
    k not in {0, 1} (the paper needs k >= 2); this half uses only the second
    clause of (H);
  * `indicial_frozen_zero`: I(1) = I(4) = 0 for all l0, l1, so the exclusion of
    k = 0, 1 is sharp (these are the two frozen exponents);
  * `m_exponent_ne_zero`: under (H), m(theta_j) /= 0 (used in Lemma 3);
  * `lattice_clause_of_difference_clause`: the second clause of (H) FOLLOWS from
    the first at t = 2 (the exponents pair up symmetrically about s = 5/2, and
    theta_j in 1 + 3Z would make 2 theta_j - 5, the difference within its pair,
    an integer).  (H) may be stated as its first clause alone.

  WHAT THIS FILE DOES NOT COVER, deliberately:
  * that I is the indicial polynomial of L6: that is `L6_coeff` in the companion
    file Theorem3Intertwining.lean (the two files are standalone, so the
    polynomial is restated here, with P specialised);
  * existence, convergence and uniqueness of the Frobenius solutions, the
    absence of logarithms at mu = 4 (c_4 = 0 is a choice), and that
    {psi_1, psi_4, psi_j} is a basis: Step 4's analysis;
  * the Hankel formula, Steps 1-3 (companion file Theorem3Intertwining.lean),
    Steps 5-6, omega and (T).

  CONTROLS, in the file:
  * `nonvacuity_H`: (H) holds at l0 = i, l1 = 2i, so the theorems are not vacuous;
    `nonvacuity_H_real`: (H) also holds at real l0 = 1/4, l1 = 1/8;
  * `tamper_lattice`: "m vanishes on all of 1 + 3Z" is FALSE (m(-2) = 3^(-2/3));
  * `tamper_no_difference_clause`: at l0 = 3, l1 = 2i the lattice clause holds,
    the difference clause fails (theta_0 - theta_1 = 6), and I(theta_1 + 6) = 0;
  * `tamper_indicial`: with (mu - 1)(mu - 5) in place of (mu - 1)(mu - 4) the
    frozen exponent 4 is lost (the tampered polynomial is non-zero at 4).
  (There is no control for dropping the lattice clause alone: at t = 2 it
  follows from the difference clause, `lattice_clause_of_difference_clause`.)

  Contributed by the cognitive agent, 2026-10-10, for the lead's commission
  LEAN3, citing the record by description only.
-/
import Mathlib.Analysis.SpecialFunctions.Gamma.Basic
import Mathlib.Analysis.SpecialFunctions.Gamma.Beta
import Mathlib.Analysis.SpecialFunctions.Pow.Complex
import Mathlib.Tactic

namespace IntegrableBoundary.Theorem3FrozenLattice

open Complex

/-! ### The zeros of the multiplier -/

/-- `m(ν) = 3^(ν/3) / Γ(1/3 - ν/3)` (entire in `ν`). -/
noncomputable def m (ν : ℂ) : ℂ := (3 : ℂ) ^ (ν / 3) * (Gamma (1 / 3 - ν / 3))⁻¹

lemma three_pow_ne (z : ℂ) : (3 : ℂ) ^ z ≠ 0 := by
  rw [cpow_def_of_ne_zero (by norm_num : (3 : ℂ) ≠ 0)]
  exact Complex.exp_ne_zero _

/-- The frozen lattice: `m(μ) = 0` iff `μ ∈ 1 + 3ℤ_{≥0}`. -/
theorem m_eq_zero_iff (μ : ℂ) : m μ = 0 ↔ ∃ n : ℕ, μ = 1 + 3 * n := by
  unfold m
  rw [mul_eq_zero, inv_eq_zero, Complex.Gamma_eq_zero_iff]
  simp only [three_pow_ne, false_or]
  constructor
  · rintro ⟨n, hn⟩
    exact ⟨n, by linear_combination -3 * hn⟩
  · rintro ⟨n, hn⟩
    exact ⟨n, by rw [hn]; ring⟩

/-! ### The indicial polynomial -/

/-- `s = 5/2` at `t = 2`. -/
noncomputable def s : ℂ := 5 / 2

/-- `P(T) = (T^2 - l0^2)(T^2 - l1^2)`. -/
def P (l0 l1 T : ℂ) : ℂ := (T ^ 2 - l0 ^ 2) * (T ^ 2 - l1 ^ 2)

/-- The four exponents `s ± l0, s ± l1`. -/
noncomputable def θ (l0 l1 : ℂ) : Fin 4 → ℂ := ![s + l0, s - l0, s + l1, s - l1]

/-- The indicial polynomial of the sixth-order form, `I(μ) = P(s - μ)(μ - 1)(μ - 4)`. -/
noncomputable def Ind (l0 l1 μ : ℂ) : ℂ := P l0 l1 (s - μ) * (μ - 1) * (μ - 4)

/-- The no-resonance hypothesis (H). -/
def H (l0 l1 : ℂ) : Prop :=
  (∀ i j : Fin 4, i ≠ j → ∀ n : ℤ, θ l0 l1 i - θ l0 l1 j ≠ n) ∧
  (∀ j : Fin 4, ∀ n : ℤ, θ l0 l1 j ≠ 1 + 3 * n)

theorem P_shift_prod (l0 l1 μ : ℂ) : P l0 l1 (s - μ) = ∏ j, (μ - θ l0 l1 j) := by
  simp [P, θ, Fin.prod_univ_four]
  ring

/-- Under (H), `I(θ_j + 3k) ≠ 0` for every integer `k ≠ 0`. -/
theorem indicial_exponent_ne_zero {l0 l1 : ℂ} (hH : H l0 l1) (j : Fin 4) (k : ℤ)
    (hk : k ≠ 0) : Ind l0 l1 (θ l0 l1 j + 3 * k) ≠ 0 := by
  obtain ⟨hdiff, hlat⟩ := hH
  unfold Ind
  rw [P_shift_prod]
  refine mul_ne_zero (mul_ne_zero ?_ ?_) ?_
  · rw [Finset.prod_ne_zero_iff]
    intro i _
    by_cases hij : i = j
    · subst hij
      have : (3 : ℂ) * k ≠ 0 := mul_ne_zero (by norm_num) (by exact_mod_cast hk)
      intro h; apply this; linear_combination h
    · intro h
      apply hdiff j i (Ne.symm hij) (-3 * k)
      push_cast
      linear_combination h
  · intro h
    apply hlat j (-k)
    push_cast
    linear_combination h
  · intro h
    apply hlat j (1 - k)
    push_cast
    linear_combination h

/-- Under the second clause of (H), `I(1 + 3k) ≠ 0` for every integer `k ∉ {0, 1}`. -/
theorem indicial_frozen_ne_zero {l0 l1 : ℂ} (hH : H l0 l1) (k : ℤ)
    (hk0 : k ≠ 0) (hk1 : k ≠ 1) : Ind l0 l1 (1 + 3 * k) ≠ 0 := by
  obtain ⟨-, hlat⟩ := hH
  unfold Ind
  rw [P_shift_prod]
  refine mul_ne_zero (mul_ne_zero ?_ ?_) ?_
  · rw [Finset.prod_ne_zero_iff]
    intro i _ h
    exact hlat i k (by linear_combination -h)
  · have : (3 : ℂ) * k ≠ 0 := mul_ne_zero (by norm_num) (by exact_mod_cast hk0)
    intro h; apply this; linear_combination h
  · have hk1' : (k : ℂ) - 1 ≠ 0 := sub_ne_zero.mpr (by exact_mod_cast hk1)
    intro h; apply mul_ne_zero (by norm_num : (3 : ℂ) ≠ 0) hk1'; linear_combination h

/-- The two frozen exponents: `I(1) = I(4) = 0` for all `l0, l1`. -/
theorem indicial_frozen_zero (l0 l1 : ℂ) : Ind l0 l1 1 = 0 ∧ Ind l0 l1 4 = 0 := by
  constructor <;> simp [Ind]

/-- Under (H), `m(θ_j) ≠ 0`. -/
theorem m_exponent_ne_zero {l0 l1 : ℂ} (hH : H l0 l1) (j : Fin 4) : m (θ l0 l1 j) ≠ 0 := by
  rw [Ne, m_eq_zero_iff]
  rintro ⟨n, hn⟩
  exact hH.2 j n (by rw [hn]; push_cast; ring)

/-- At `t = 2` the second clause of (H) follows from the first: the exponents come
in pairs symmetric about `s = 5/2`, so `θ_j ∈ 1 + 3ℤ` would make the difference
within its pair, `2θ_j - 5`, an integer. -/
theorem lattice_clause_of_difference_clause {l0 l1 : ℂ}
    (hd : ∀ i j : Fin 4, i ≠ j → ∀ n : ℤ, θ l0 l1 i - θ l0 l1 j ≠ n) :
    ∀ j : Fin 4, ∀ n : ℤ, θ l0 l1 j ≠ 1 + 3 * n := by
  intro j n hj
  fin_cases j <;> simp [θ, s] at hj
  · exact hd 0 1 (by decide) (6 * n - 3) (by simp [θ, s]; linear_combination 2 * hj)
  · exact hd 0 1 (by decide) (3 - 6 * n) (by simp [θ, s]; linear_combination -2 * hj)
  · exact hd 2 3 (by decide) (6 * n - 3) (by simp [θ, s]; linear_combination 2 * hj)
  · exact hd 2 3 (by decide) (3 - 6 * n) (by simp [θ, s]; linear_combination -2 * hj)

/-! ### Controls -/

/-- NON-VACUITY: (H) holds at `l0 = i`, `l1 = 2i`. -/
theorem nonvacuity_H : H I (2 * I) := by
  have hd : ∀ i j : Fin 4, i ≠ j → ∀ n : ℤ, θ I (2 * I) i - θ I (2 * I) j ≠ n := by
    intro i j hij n h
    have h' := congrArg Complex.im h
    fin_cases i <;> fin_cases j <;> (try simp at hij) <;> norm_num [θ, s] at h'
  exact ⟨hd, lattice_clause_of_difference_clause hd⟩

/-- NON-VACUITY: (H) also holds at real `l0 = 1/4`, `l1 = 1/8`
(exponents 11/4, 9/4, 21/8, 19/8; differences of modulus below 1). -/
theorem nonvacuity_H_real : H (1 / 4) (1 / 8) := by
  have hd : ∀ i j : Fin 4, i ≠ j → ∀ n : ℤ, θ (1 / 4) (1 / 8) i - θ (1 / 4) (1 / 8) j ≠ n := by
    intro i j hij n h
    have h' := congrArg Complex.re h
    fin_cases i <;> fin_cases j <;> (try simp at hij) <;> norm_num [θ, s] at h' <;>
    · have a1 : (n : ℝ) < 1 := by linarith
      have a2 : (-1 : ℝ) < n := by linarith
      have a3 : (n : ℝ) ≠ 0 := by intro h0; rw [h0] at h'; norm_num at h'
      have b1 : n < 1 := by exact_mod_cast a1
      have b2 : -1 < n := by exact_mod_cast a2
      have b3 : n ≠ 0 := by exact_mod_cast a3
      omega
  exact ⟨hd, lattice_clause_of_difference_clause hd⟩

/-- TAMPER: `m` does not vanish on all of `1 + 3ℤ` (`m(-2) = 3^(-2/3) ≠ 0`). -/
theorem tamper_lattice : ¬ ∀ k : ℤ, m (1 + 3 * k) = 0 := by
  intro h
  obtain ⟨n, hn⟩ := (m_eq_zero_iff _).1 (h (-1))
  have := congrArg Complex.re hn
  push_cast at this
  norm_num at this
  linarith [(Nat.cast_nonneg n : (0 : ℝ) ≤ n)]

/-- TAMPER: the difference clause cannot be dropped.  At `l0 = 3`, `l1 = 2i` the
exponents are 11/2, -1/2, 5/2 ± 2i: the lattice clause holds, the difference
clause fails (`θ_0 - θ_1 = 6`), and `I(θ_1 + 3·2) = I(θ_0) = 0`. -/
theorem tamper_no_difference_clause :
    (∀ j : Fin 4, ∀ n : ℤ, θ 3 (2 * I) j ≠ 1 + 3 * n) ∧
    θ 3 (2 * I) 0 - θ 3 (2 * I) 1 = ((6 : ℤ) : ℂ) ∧
    Ind 3 (2 * I) (θ 3 (2 * I) 1 + 3 * ((2 : ℤ) : ℂ)) = 0 := by
  refine ⟨?_, by simp [θ, s]; norm_num, by simp [Ind, P, θ, s]; norm_num⟩
  intro j n h
  fin_cases j
  · have h' := congrArg Complex.re h
    norm_num [θ, s] at h'
    have a1 : (1 : ℝ) < n := by linarith
    have a2 : (n : ℝ) < 2 := by linarith
    have b1 : 1 < n := by exact_mod_cast a1
    have b2 : n < 2 := by exact_mod_cast a2
    omega
  · have h' := congrArg Complex.re h
    norm_num [θ, s] at h'
    have a1 : (-1 : ℝ) < n := by linarith
    have a2 : (n : ℝ) < 0 := by linarith
    have b1 : -1 < n := by exact_mod_cast a1
    have b2 : n < 0 := by exact_mod_cast a2
    omega
  all_goals
    have h' := congrArg Complex.im h
    norm_num [θ, s] at h'

/-- `I` with `(μ - 1)(μ - 5)` in place of `(μ - 1)(μ - 4)`. -/
noncomputable def IndTampered (l0 l1 μ : ℂ) : ℂ := P l0 l1 (s - μ) * (μ - 1) * (μ - 5)

/-- TAMPER: with `(μ - 1)(μ - 5)` the frozen exponent 4 is lost: the tampered
polynomial does not vanish at 4 (at `l0 = i`, `l1 = 2i`), while `I(4) = 0`. -/
theorem tamper_indicial : IndTampered I (2 * I) 4 ≠ 0 ∧ Ind I (2 * I) 4 = 0 := by
  refine ⟨?_, (indicial_frozen_zero _ _).2⟩
  simp only [IndTampered, P, s, mul_pow, I_sq]
  norm_num

end IntegrableBoundary.Theorem3FrozenLattice
