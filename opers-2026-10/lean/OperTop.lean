/-
  The oper top layer: its rank-drop vanishing law, and the Lagrange coefficient
  that gives the F1 top law, for every spin.

  PART 1 -- RANK-DROP VANISHING (the record's item 304, Z0(a)).  The top layer
  of the oper's WKB charge at spin index K is

      (-1)^K sum_a C(nu, a) C(nu, K - a) X^a Y^(K-a),     nu = (2K - 1)/n,

  with C the generalised binomial coefficient (`top`, `gb`).  Proved for every
  K and every rational nu:
  * every coefficient vanishes iff nu = i is a NATURAL number with 2i <= K - 1
    (`top_all_vanish_iff`).  "Integer" in the record means non-negative: at a
    negative integer, as at any non-natural nu, NO coefficient vanishes
    (`top_none_vanish`);
  * at nu = i the number of vanishing coefficients among a = 0..K is
    min(K + 1, 2(K - i)) (`top_vanishing_count`), the record's multiplicity
    law.

  PART 2 -- THE LAGRANGE COEFFICIENT (item 294(a)).  For the curve
  z = z0 (1 - Xz)^e1 (1 - Yz)^e2, the power-form Lagrange-Buermann coefficient
  of (z/z0)^beta is  beta/(k+beta) [z^k] (1 - Xz)^(e1(k+beta)) (1 - Yz)^(e2(k+beta))
  (`lagrangeCoeff`), with (1 - Xz)^r built from Mathlib's binomial series, so
  rational exponents are genuine formal power series.  Proved for every k:
  * its closed form, for any beta (`lagrangeCoeff_eq`);
  * at beta = -1/2, item 294(a)'s law
      -(1/(2k-1)) (1/k!) sum_i C(k,i) (p)_(k-i) (q)_i X^(k-i) Y^i,
      p = -e1 (k - 1/2),  q = -e2 (k - 1/2)        (`lagrange_half`),
    which for X != 0 and (p)_k != 0 is a constant times X^k [W^k] F1(1; p, q;
    p; W, (Y/X)W), the F1 top law of items 168 and 178 (`lagrange_half_eq_F1`);
  * the three solutions' instances, X carrying a_Y and Y carrying a_X as in
    item 178's table (`lagrange_sol1`, `lagrange_sol2`, `lagrange_sol3`).

  WHAT THIS FILE DOES NOT COVER, deliberately:
  * the LAGRANGE-BUERMANN THEOREM ITSELF -- that for z = z0 phi(z), phi(0) = 1,
    the k-th coefficient of (z/z0)^beta in z0 equals beta/(k+beta) [z^k]
    phi(z)^(k+beta) -- is a classical theorem used as an INPUT and NOT
    formalised here (Mathlib has no Lagrange inversion).  What is proved is
    the coefficient side, for every k, which the record had checked by series
    reversion only to k = 6;
  * that the oper's top layer IS the expression of part 1 (item 304 derives it
    from the operator); that the record's exponents are crossed; that the curve
    is the right curve -- all inputs;
  * nothing about the lower layers (loss >= 1), where rank-drop vanishing at
    all layers is still open in the record (item 311 (2));
  * no statement about spectral determinants or excited states.

  NEGATIVE CONTROLS, in the file: at the boundary 2i = K not every
  coefficient vanishes (the middle one is 1; exactly two vanish); off the
  integers a coefficient is nonzero by direct evaluation (-1/8 at nu = 1/2,
  K = 2); and the exponent e (k + beta) is load-bearing (with e k the value is
  1 instead of 1/2).

  Contributed by the cognitive agent, 2026-10-02, under the standing
  arrangement, citing the record by description only.
-/
import Mathlib.Tactic
import Mathlib.RingTheory.PowerSeries.Binomial

namespace IntegrableBoundary.OperTop

open Finset Polynomial

/-! ## Part 1: the top layer's rank-drop vanishing law -/

/-- The falling factorial `ν (ν-1) ... (ν-a+1)`. -/
def ff (ν : ℚ) : ℕ → ℚ
  | 0 => 1
  | a + 1 => ff ν a * (ν - a)

/-- generalized binomial coefficient `C(ν, a)` -/
noncomputable def gb (ν : ℚ) (a : ℕ) : ℚ := ff ν a / a.factorial

lemma ff_eq_descPochhammer (ν : ℚ) (a : ℕ) : ff ν a = (descPochhammer ℚ a).eval ν := by
  induction a with
  | zero => simp [ff]
  | succ a ih => rw [ff, descPochhammer_succ_right, Polynomial.eval_mul, ih]; simp

lemma ff_eq_zero_iff (ν : ℚ) (a : ℕ) : ff ν a = 0 ↔ ∃ i : ℕ, ν = i ∧ i < a := by
  induction a with
  | zero => simp [ff]
  | succ a ih =>
    rw [ff, mul_eq_zero, ih, sub_eq_zero]
    constructor
    · rintro (⟨i, hi, hia⟩ | h)
      · exact ⟨i, hi, by omega⟩
      · exact ⟨a, h, by omega⟩
    · rintro ⟨i, hi, hia⟩
      rcases Nat.lt_succ_iff_lt_or_eq.mp hia with h | h
      · exact Or.inl ⟨i, hi, h⟩
      · exact Or.inr (by rw [hi, h])

lemma gb_eq_zero_iff (ν : ℚ) (a : ℕ) : gb ν a = 0 ↔ ∃ i : ℕ, ν = i ∧ i < a := by
  unfold gb; rw [div_eq_zero_iff, ff_eq_zero_iff]
  have : a.factorial ≠ 0 := Nat.factorial_ne_zero a
  simp [this]

/-- On natural arguments `gb` is `Nat.choose`. -/
lemma gb_nat (n a : ℕ) : gb n a = (n.choose a : ℚ) := by
  unfold gb
  rw [ff_eq_descPochhammer, descPochhammer_eval_eq_descFactorial]
  rw [Nat.descFactorial_eq_factorial_mul_choose]; push_cast
  field_simp

/-- The top layer of item 304 Z0(a): coefficient of `X^a Y^(K-a)` in
    `(-1)^K Σ_a C(ν,a) C(ν,K-a) X^a Y^(K-a)`. -/
noncomputable def top (ν : ℚ) (K a : ℕ) : ℚ := (-1) ^ K * (gb ν a * gb ν (K - a))

lemma top_eq_zero_iff (ν : ℚ) (K a : ℕ) :
    top ν K a = 0 ↔ (∃ i : ℕ, ν = i ∧ i < a) ∨ (∃ i : ℕ, ν = i ∧ i < K - a) := by
  unfold top
  rw [mul_eq_zero, mul_eq_zero, gb_eq_zero_iff, gb_eq_zero_iff]
  simp

/-- **Z0(a), first part.**  Every coefficient of the top layer vanishes iff
    `ν = i` is a natural number with `2i ≤ K - 1`. -/
theorem top_all_vanish_iff (ν : ℚ) (K : ℕ) :
    (∀ a ≤ K, top ν K a = 0) ↔ ∃ i : ℕ, ν = i ∧ 2 * i + 1 ≤ K := by
  constructor
  · intro h
    -- a = 0 gives ν = i with i < K (the a < 0 alternative is impossible)
    have h0 := (top_eq_zero_iff ν K 0).mp (h 0 (Nat.zero_le _))
    obtain ⟨i, hi, hiK⟩ : ∃ i : ℕ, ν = i ∧ i < K := by
      rcases h0 with ⟨i, _, h⟩ | ⟨i, hi, h⟩
      · omega
      · exact ⟨i, hi, by omega⟩
    refine ⟨i, hi, ?_⟩
    by_contra hcon
    -- take a = i (≤ K): neither factor vanishes when K - i ≤ i
    have hai := (top_eq_zero_iff ν K i).mp (h i hiK.le)
    have inj : ∀ j : ℕ, ν = j → j = i := fun j hj => by
      have : (j : ℚ) = i := by rw [← hj, hi]
      exact_mod_cast this
    rcases hai with ⟨j, hj, h⟩ | ⟨j, hj, h⟩
    · have := inj j hj; omega
    · have := inj j hj; omega
  · rintro ⟨i, hi, hK⟩ a ha
    rw [top_eq_zero_iff]
    by_cases hai : i < a
    · exact Or.inl ⟨i, hi, hai⟩
    · exact Or.inr ⟨i, hi, by omega⟩

/-- Off the naturals no coefficient vanishes. -/
theorem top_none_vanish (ν : ℚ) (hν : ∀ i : ℕ, ν ≠ i) (K a : ℕ) : top ν K a ≠ 0 := by
  rw [Ne, top_eq_zero_iff]
  rintro (⟨i, hi, _⟩ | ⟨i, hi, _⟩) <;> exact hν i hi

lemma top_nat_eq_zero_iff (i K a : ℕ) (ha : a ≤ K) :
    top (i : ℚ) K a = 0 ↔ a < K - i ∨ i < a := by
  rw [top_eq_zero_iff]
  have inj : ∀ j : ℕ, (i : ℚ) = j → j = i := fun j hj => by exact_mod_cast hj.symm
  constructor
  · rintro (⟨j, hj, h⟩ | ⟨j, hj, h⟩)
    · have := inj j hj; omega
    · have := inj j hj; omega
  · rintro (h | h)
    · exact Or.inr ⟨i, rfl, by omega⟩
    · exact Or.inl ⟨i, rfl, h⟩

/-- **Z0(a), second part (the multiplicity law).**  At `ν = i` the number of
    vanishing coefficients among `a = 0..K` is `min(K+1, 2(K-i))`. -/
theorem top_vanishing_count (i K : ℕ) :
    ((range (K + 1)).filter (fun a => top (i : ℚ) K a = 0)).card = min (K + 1) (2 * (K - i)) := by
  by_cases h : 2 * i + 1 ≤ K
  · have : (range (K + 1)).filter (fun a => top (i : ℚ) K a = 0) = range (K + 1) := by
      ext a; simp only [mem_filter, mem_range, and_iff_left_iff_imp]
      intro ha; rw [top_nat_eq_zero_iff i K a (by omega)]; omega
    rw [this, card_range]; omega
  · have : (range (K + 1)).filter (fun a => top (i : ℚ) K a = 0)
        = range (K - i) ∪ Ioc i K := by
      ext a; simp only [mem_filter, mem_range, mem_union, mem_Ioc]
      constructor
      · rintro ⟨ha, h0⟩; rw [top_nat_eq_zero_iff i K a (by omega)] at h0; omega
      · intro h'; refine ⟨by omega, ?_⟩
        rw [top_nat_eq_zero_iff i K a (by omega)]; omega
    rw [this, card_union_of_disjoint, card_range, Nat.card_Ioc]
    · omega
    · rw [disjoint_left]; intro a h1 h2; simp only [mem_range, mem_Ioc] at h1 h2; omega


/-! ## Part 2: the Lagrange coefficient, for every k -/

/-- The rising factorial `(x)_n`. -/
def rf (x : ℚ) : ℕ → ℚ
  | 0 => 1
  | n + 1 => rf x n * (x + n)

/-- Mathlib's generalised binomial coefficient on ℚ is `ff r n / n!`. -/
lemma choose_eq_ff (r : ℚ) (n : ℕ) : Ring.choose r n = ff r n / n.factorial := by
  have h := Ring.descPochhammer_eq_factorial_smul_choose r n
  rw [← aeval_eq_smeval, aeval_def, eval₂_eq_eval_map, descPochhammer_map, nsmul_eq_mul,
    ← ff_eq_descPochhammer] at h
  rw [h]; field_simp

lemma ff_neg (A : ℚ) (n : ℕ) : (-1) ^ n * ff A n = rf (-A) n := by
  induction n with
  | zero => simp [ff, rf]
  | succ n ih => rw [ff, rf, ← ih, pow_succ]; ring

/-- `(1 - X z)^r` as a formal power series in z: Mathlib's binomial series
    `(1 + u)^r`, rescaled by `u = -X z`. -/
noncomputable def oneSub (X r : ℚ) : PowerSeries ℚ :=
  PowerSeries.rescale (-X) (PowerSeries.binomialSeries ℚ r)

lemma coeff_oneSub (X r : ℚ) (n : ℕ) :
    PowerSeries.coeff n (oneSub X r) = X ^ n * (rf (-r) n / n.factorial) := by
  simp only [oneSub, PowerSeries.coeff_rescale, PowerSeries.binomialSeries_coeff, smul_eq_mul,
    mul_one, choose_eq_ff, ← ff_neg r n]
  rw [neg_pow]; field_simp

/-- The power-form Lagrange-Buermann coefficient for
    `z = z0 (1 - Xz)^e1 (1 - Yz)^e2` and the power `(z/z0)^beta`:
    `beta/(k+beta) [z^k] (1 - Xz)^(e1 (k+beta)) (1 - Yz)^(e2 (k+beta))`. -/
noncomputable def lagrangeCoeff (β e1 e2 X Y : ℚ) (k : ℕ) : ℚ :=
  β / (k + β) * PowerSeries.coeff k (oneSub X (e1 * (k + β)) * oneSub Y (e2 * (k + β)))

/-- The F1-line numerator in homogeneous form:
    `sum_x C(k, x.2) (p)_(x.1) (q)_(x.2) X^(x.1) Y^(x.2)`. -/
noncomputable def Nh (p q X Y : ℚ) (k : ℕ) : ℚ :=
  ∑ x ∈ antidiagonal k, (k.choose x.2 : ℚ) * rf p x.1 * rf q x.2 * X ^ x.1 * Y ^ x.2

/-- **The Lagrange coefficient in closed form, for every k** (any beta, e1, e2). -/
theorem lagrangeCoeff_eq (β e1 e2 X Y : ℚ) (k : ℕ) :
    lagrangeCoeff β e1 e2 X Y k
      = β / (k + β) / k.factorial * Nh (-(e1 * (k + β))) (-(e2 * (k + β))) X Y k := by
  unfold lagrangeCoeff Nh
  rw [PowerSeries.coeff_mul, mul_sum, mul_sum]
  refine sum_congr rfl (fun x hx => ?_)
  rw [mem_antidiagonal] at hx
  rw [coeff_oneSub, coeff_oneSub]
  have hc := Nat.add_choose_mul_factorial_mul_factorial x.1 x.2
  rw [hx] at hc
  have hc' : (k.choose x.2 : ℚ) * x.1.factorial * x.2.factorial = k.factorial := by exact_mod_cast hc
  have h1 : (x.1.factorial : ℚ) ≠ 0 := by positivity
  have h2 : (x.2.factorial : ℚ) ≠ 0 := by positivity
  have hC : (k.choose x.2 : ℚ) ≠ 0 := by
    have : 0 < k.choose x.2 := Nat.choose_pos (by omega)
    exact_mod_cast this.ne'
  rw [← hc']; field_simp

/-- **Item 294(a)'s law at beta = -1/2**:
    `[z0^k](z/z0)^(-1/2) = -(1/(2k-1)) (1/k!) sum_i C(k,i) (p)_(k-i) (q)_i X^(k-i) Y^i`,
    `p = -e1 (k - 1/2)`, `q = -e2 (k - 1/2)` (given Lagrange-Buermann; see the header). -/
theorem lagrange_half (e1 e2 X Y : ℚ) (k : ℕ) :
    lagrangeCoeff (-1 / 2) e1 e2 X Y k
      = -(1 / (2 * k - 1)) / k.factorial * Nh (-e1 * (k - 1 / 2)) (-e2 * (k - 1 / 2)) X Y k := by
  rw [lagrangeCoeff_eq]
  have e : (-(1 : ℚ) / 2) / (k + -1 / 2) = -(1 / (2 * k - 1)) := by
    rcases Nat.eq_zero_or_pos k with h | h
    · subst h; norm_num
    · have h2 : (2 * (k : ℚ) - 1) ≠ 0 := by
        have : (1 : ℚ) ≤ k := by exact_mod_cast h
        linarith
      rw [show (k : ℚ) + -1 / 2 = (2 * k - 1) / 2 by ring, div_div_eq_mul_div]
      field_simp
  rw [e]; congr 2 <;> ring

/-- The F1 coefficient `[W^k] F1(1; p, q; p; W, zW) = N_k / (p)_k` (as in the
    record's item 178, and the AppellF1Top certificate). -/
noncomputable def F1coef (p q z : ℚ) (k : ℕ) : ℚ :=
  (∑ x ∈ antidiagonal k, (k.choose x.2 : ℚ) * rf p x.1 * rf q x.2 * z ^ x.2) / rf p k

/-- The Lagrange coefficient is the F1 top law: for `X ≠ 0` and `(p)_k ≠ 0`,
    `[z0^k](z/z0)^(-1/2) = -(1/(2k-1)) ((p)_k / k!) X^k F1(1; p, q; p; W, (Y/X)W)_k`.
    (The unnormalised form `lagrange_half` needs neither hypothesis.) -/
theorem lagrange_half_eq_F1 (e1 e2 X Y : ℚ) (hX : X ≠ 0) (k : ℕ)
    (hp : rf (-e1 * (k - 1 / 2)) k ≠ 0) :
    lagrangeCoeff (-1 / 2) e1 e2 X Y k
      = -(1 / (2 * k - 1)) * (rf (-e1 * (k - 1 / 2)) k / k.factorial) * X ^ k
          * F1coef (-e1 * (k - 1 / 2)) (-e2 * (k - 1 / 2)) (Y / X) k := by
  rw [lagrange_half]
  unfold F1coef Nh
  rw [mul_sum, sum_div, mul_sum]
  refine sum_congr rfl (fun x hx => ?_)
  rw [mem_antidiagonal] at hx
  have hk : X ^ k * (Y / X) ^ x.2 = X ^ x.1 * Y ^ x.2 := by
    rw [← hx, pow_add, div_pow]; field_simp
  have hf : (k.factorial : ℚ) ≠ 0 := by positivity
  have hk' : (k.choose x.2 : ℚ) * rf (-e1 * (k - 1 / 2)) x.1 * rf (-e2 * (k - 1 / 2)) x.2 * X ^ x.1 * Y ^ x.2
      = (k.choose x.2 : ℚ) * rf (-e1 * (k - 1 / 2)) x.1 * rf (-e2 * (k - 1 / 2)) x.2 * (X ^ k * (Y / X) ^ x.2) := by
    rw [hk]; ring
  rw [hk']
  generalize rf (-e1 * (k - 1 / 2)) k = R at hp ⊢
  field_simp

/-! The three solutions (item 178's table): X carries a_Y, Y carries a_X. -/

theorem lagrange_sol1 (t X Y : ℚ) (k : ℕ) :
    lagrangeCoeff (-1 / 2) ((t - 1) / (t + 3)) 1 X Y k
      = -(1 / (2 * k - 1)) / k.factorial
          * Nh (-((t - 1) / (t + 3)) * (k - 1 / 2)) (-1 * (k - 1 / 2)) X Y k :=
  lagrange_half _ _ X Y k

theorem lagrange_sol2 (t X Y : ℚ) (k : ℕ) :
    lagrangeCoeff (-1 / 2) (4 * (t - 1) / (t + 5)) (-2 * (t - 3) / (t + 5)) X Y k
      = -(1 / (2 * k - 1)) / k.factorial
          * Nh (-(4 * (t - 1) / (t + 5)) * (k - 1 / 2)) (-(-2 * (t - 3) / (t + 5)) * (k - 1 / 2)) X Y k :=
  lagrange_half _ _ X Y k

theorem lagrange_sol3 (t X Y : ℚ) (k : ℕ) :
    lagrangeCoeff (-1 / 2) ((t - 3) / (t - 5)) ((t - 3) / (t - 5)) X Y k
      = -(1 / (2 * k - 1)) / k.factorial
          * Nh (-((t - 3) / (t - 5)) * (k - 1 / 2)) (-((t - 3) / (t - 5)) * (k - 1 / 2)) X Y k :=
  lagrange_half _ _ X Y k

/-! ## Negative controls -/

/-- Part 1: at the boundary `2i = K` (i = 1, K = 2) NOT every coefficient
    vanishes -- the middle one is `C(1,1)^2 = 1` -- and exactly two do. -/
theorem control_boundary_not_all :
    top (1 : ℕ) 2 1 = 1 ∧
    ((range (2 + 1)).filter (fun a => top ((1 : ℕ) : ℚ) 2 a = 0)).card = 2 := by
  refine ⟨?_, by rw [top_vanishing_count]; decide⟩
  simp only [top]; rw [show (2 : ℕ) - 1 = 1 from rfl, gb_nat]; norm_num

/-- Part 1: off the integers a coefficient is nonzero by direct evaluation:
    `top (1/2) 2 0 = C(1/2,0) C(1/2,2) = -1/8`. -/
theorem control_half_integer : top (1 / 2) 2 0 = -1 / 8 := by
  simp [top, gb, ff, Nat.factorial]; norm_num

/-- Part 2: the exponent `e (k + beta)` is load-bearing.  At k = 1, beta = -1/2,
    e1 = 1, e2 = 0, X = 1 the coefficient is 1/2; with the exponent `e k` instead
    of `e (k + beta)` the same expression gives 1. -/
theorem control_exponent :
    lagrangeCoeff (-1 / 2) 1 0 1 0 1 = 1 / 2 ∧
    (-1 / 2 : ℚ) / ((1 : ℕ) + -1 / 2) * PowerSeries.coeff 1 (oneSub 1 (1 * (1 : ℕ)) * oneSub 0 0) = 1 := by
  constructor
  · rw [lagrangeCoeff_eq]; simp [Nh, rf, Finset.Nat.antidiagonal_succ]; norm_num
  · rw [PowerSeries.coeff_mul]; simp [coeff_oneSub, rf, Finset.Nat.antidiagonal_succ]; norm_num

end IntegrableBoundary.OperTop
