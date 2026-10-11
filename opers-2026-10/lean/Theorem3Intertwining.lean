/-
  Theorem 3 (spectral equivalence at t = 2): the multiplier recurrences and the
  intertwining, at the level of Mellin coefficients.

  THE CLAIM (the opers paper, Appendix A, Steps 1 and 3; the record's items on the
  analytic proof of Theorem 3 and its review).  At t = 2, with s = 5/2,

      K_g x^mu = 3^(mu/3) x^mu / Gamma(g - mu/3),   K = K_{1/3},
      m(mu)   = 3^(mu/3) / Gamma(1/3 - mu/3)        (so K x^mu = m(mu) x^mu),

  and, with T = s - theta, theta = x d/dx,
      L4 = P(T) + x^9 (theta + 2) - E x^6,
      L6 = P(T) (theta - 1)(theta - 4) - x^9 - E x^6,
  the paper's Lemma 2 and Proposition 1 hold:
      (i)   K_{g-1} = K_g (g - 1 - theta/3),
      (ii)  K x^3 = - x^3 K (theta + 2),
      (iii) K = x^6 K x^-6 (theta - 1)(theta - 4),
      Prop. 1:  L4 K = x^6 K x^-6 L6.

  THE MODEL.  A formal sum  sum_nu c(nu) x^nu  is recorded by its coefficient
  function c : C -> C (exponent |-> coefficient; any support, so finite sums,
  monomials and termwise Frobenius series are all instances).  On it
    theta     multiplies the coefficient at nu by nu        (theta x^nu = nu x^nu),
    x^a       shifts:   (x^a c)(nu) = c(nu - a),
    P(T)      multiplies by P(s - nu)                       (T x^nu = (s - nu) x^nu),
    K_g       multiplies by k g nu = 3^(nu/3) / Gamma(g - nu/3).
  Gamma is Mathlib's Complex.Gamma; 3^z is the principal complex power (3 > 0, so
  it is exp(z log 3)).  Mathlib's Gamma is 0 at its poles and 0^-1 = 0, so
  (Gamma z)^-1 is the ENTIRE function 1/Gamma, and every identity below holds at
  every complex exponent, with no "where defined" proviso (`invGamma_sub_one`).

  WHAT IS PROVED (no sorry, no axioms beyond Lean/Mathlib's):
  * `k_sub_one`, `k_add_three`: the two multiplier steps, from 1/Gamma(z) =
    z/Gamma(z+1) only;
  * `m_add_three`: m(mu + 3) = -(mu + 2) m(mu), and `m_sub_six`:
    m(mu) = (mu - 1)(mu - 4) m(mu - 6), for ALL mu in C; the ratio forms of the
    paper, `m_ratio_add_three` and `m_ratio_sub_six`, wherever the denominator
    is non-zero;
  * `lemma2_i`, `lemma2_ii`, `lemma2_iii`, `aux_ii` (K x^3 = 3 x^3 K_{-2/3}),
    `aux_iii` (x^6 K x^-6 = (1/9) K_{7/3}): the operator identities on every
    coefficient function; (ii) and (iii) are derived from (i) and the two
    auxiliary identities, the route of the paper's proof;
  * `prop1`: L4 K c = x^6 K x^-6 L6 c for EVERY function P (not only
    (T^2 - l0^2)(T^2 - l1^2)), every E in C and every coefficient function c;
    `prop1_monomial` is the case c = x^mu;
  * `K_monomial`: K x^mu = m(mu) x^mu in the model;
  * `L6_coeff`, `L4_coeff`: the coefficient of L6 c at x^mu is
    I(mu) c(mu) - c(mu - 9) - E c(mu - 6) with I(mu) = P(s - mu)(mu - 1)(mu - 4),
    and that of L4 a is P(s - mu) a(mu) + (mu - 7) a(mu - 9) - E a(mu - 6): the
    paper's indicial polynomial and both Frobenius recurrences of Step 4.

  WHAT THIS FILE DOES NOT COVER, deliberately:
  * Step 1's Hankel formula (1/(2 pi i)) int_H e^tau tau^-z dtau = 1/Gamma(z): the
    model DEFINES K_g on monomials by its result;
  * Step 2 (growth bounds, convergence, holomorphy, theta K = K theta on
    functions, deformation of the contour) and the integration by parts that
    proves Lemma 2(i) on functions: here (i) is its coefficient form;
  * that K acts termwise on the Frobenius series (Lemma 3's convergence
    argument), Steps 4-6, the constant omega, and (T);
  * the frozen lattice and the indicial non-vanishing: the companion file
    Theorem3FrozenLattice.lean.

  CONTROLS, in the file:
  * `nonvacuity_prop1`: for the paper's P with l0 = i, l1 = 2i, E = 1, both sides
    of Prop. 1 on x^0 at exponent 0 equal P(5/2) m(0), which is NOT zero;
  * `tamper_recurrence`: (mu - 1)(mu - 5) in place of (mu - 1)(mu - 4) is FALSE
    (at mu = 0, because m(-6) is not zero);
  * `tamper_prop1`: x^9 (theta + 3) in place of x^9 (theta + 2) in L4 makes
    Prop. 1 FALSE (on x^0 at exponent 9, because m(0) is not zero).

  Contributed by the cognitive agent, 2026-10-10, for the lead's commission
  LEAN3, citing the record by description only.
-/
import Mathlib.Analysis.SpecialFunctions.Gamma.Basic
import Mathlib.Analysis.SpecialFunctions.Gamma.Beta
import Mathlib.Analysis.SpecialFunctions.Pow.Complex
import Mathlib.Tactic

namespace IntegrableBoundary.Theorem3Intertwining

open Complex

/-! ### The multipliers -/

/-- The multiplier of `K_g` on `x^ν`: `3^(ν/3) / Γ(g - ν/3)` (entire in `ν`). -/
noncomputable def k (g ν : ℂ) : ℂ := (3 : ℂ) ^ (ν / 3) * (Gamma (g - ν / 3))⁻¹

/-- The multiplier of `K = K_{1/3}`: `m(ν) = 3^(ν/3) / Γ(1/3 - ν/3)`. -/
noncomputable def m (ν : ℂ) : ℂ := k (1 / 3) ν

lemma three_ne : (3 : ℂ) ≠ 0 := by norm_num

/-- `1/Γ(z - 1) = (z - 1)/Γ(z)` for every complex `z` (1/Γ is entire). -/
lemma invGamma_sub_one (z : ℂ) : (Gamma (z - 1))⁻¹ = (z - 1) * (Gamma z)⁻¹ := by
  have h := Complex.one_div_Gamma_eq_self_mul_one_div_Gamma_add_one (z - 1)
  rwa [sub_add_cancel] at h

/-- Lemma 2(i) on monomials: `K_{g-1} x^ν = (g - 1 - ν/3) K_g x^ν`. -/
lemma k_sub_one (g ν : ℂ) : k (g - 1) ν = (g - 1 - ν / 3) * k g ν := by
  unfold k
  rw [show g - 1 - ν / 3 = (g - ν / 3) - 1 by ring, invGamma_sub_one]
  ring

/-- Raising the exponent by 3: `K_g x^(ν+3) = 3 x^3 K_{g-1} x^ν`, as multipliers. -/
lemma k_add_three (g ν : ℂ) : k g (ν + 3) = 3 * k (g - 1) ν := by
  unfold k
  rw [show (ν + 3) / 3 = ν / 3 + 1 by ring, cpow_add _ _ three_ne, cpow_one,
    show g - (ν / 3 + 1) = g - 1 - ν / 3 by ring]
  ring

/-- `m(μ + 3) = -(μ + 2) m(μ)` for every complex `μ`. -/
theorem m_add_three (μ : ℂ) : m (μ + 3) = -(μ + 2) * m μ := by
  unfold m
  rw [k_add_three, k_sub_one]
  ring

/-- `m(μ) = (μ - 1)(μ - 4) m(μ - 6)` for every complex `μ`. -/
theorem m_sub_six (μ : ℂ) : m μ = (μ - 1) * (μ - 4) * m (μ - 6) := by
  have h1 := m_add_three (μ - 3)
  have h2 := m_add_three (μ - 6)
  rw [show μ - 3 + 3 = μ by ring] at h1
  rw [show μ - 6 + 3 = μ - 3 by ring] at h2
  rw [h1, h2]
  ring

/-- The paper's ratio form, wherever `m(μ) ≠ 0`. -/
theorem m_ratio_add_three (μ : ℂ) (h : m μ ≠ 0) : m (μ + 3) / m μ = -(μ + 2) := by
  rw [m_add_three, mul_div_assoc, div_self h, mul_one]

/-- The paper's ratio form, wherever `m(μ - 6) ≠ 0`. -/
theorem m_ratio_sub_six (μ : ℂ) (h : m (μ - 6) ≠ 0) :
    m μ / m (μ - 6) = (μ - 1) * (μ - 4) := by
  rw [m_sub_six μ, mul_div_assoc, div_self h, mul_one]

/-! ### The coefficient model -/

/-- A formal sum `∑_ν c ν x^ν`, recorded by its coefficient at each exponent. -/
abbrev Coeff := ℂ → ℂ

/-- `s = 5/2` at `t = 2`. -/
noncomputable def s : ℂ := 5 / 2

/-- `θ = x d/dx`: `θ x^ν = ν x^ν`. -/
def theta (c : Coeff) : Coeff := fun ν => ν * c ν

/-- Multiplication by `x^a`. -/
def xpow (a : ℂ) (c : Coeff) : Coeff := fun ν => c (ν - a)

/-- `K_g`: `K_g x^ν = k g ν x^ν`. -/
noncomputable def Kg (g : ℂ) (c : Coeff) : Coeff := fun ν => k g ν * c ν

/-- `K = K_{1/3}`. -/
noncomputable def K (c : Coeff) : Coeff := Kg (1 / 3) c

/-- `P(T)` with `T = s - θ`: `P(T) x^ν = P(s - ν) x^ν`. -/
noncomputable def PT (P : ℂ → ℂ) (c : Coeff) : Coeff := fun ν => P (s - ν) * c ν

/-- The monomial `x^μ`. -/
noncomputable def mono (μ : ℂ) : Coeff := fun ν => if ν = μ then 1 else 0

/-- `L4 = P(T) + x^9 (θ + 2) - E x^6`. -/
noncomputable def L4 (P : ℂ → ℂ) (E : ℂ) (c : Coeff) : Coeff := fun ν =>
  PT P c ν + xpow 9 (fun μ => theta c μ + 2 * c μ) ν - E * xpow 6 c ν

/-- `L6 = P(T) (θ - 1)(θ - 4) - x^9 - E x^6`. -/
noncomputable def L6 (P : ℂ → ℂ) (E : ℂ) (c : Coeff) : Coeff := fun ν =>
  PT P (fun μ => theta (fun ρ => theta c ρ - 4 * c ρ) μ - (theta c μ - 4 * c μ)) ν
    - xpow 9 c ν - E * xpow 6 c ν

theorem K_monomial (μ : ℂ) : K (mono μ) = fun ν => m μ * mono μ ν := by
  funext ν
  simp only [K, Kg, mono, m]
  split_ifs with h <;> simp [h]

/-- The coefficient of `L6 c` at `x^μ`: `I(μ) c(μ) - c(μ - 9) - E c(μ - 6)` with
`I(μ) = P(s - μ)(μ - 1)(μ - 4)`, the paper's indicial polynomial and recurrence. -/
theorem L6_coeff (P : ℂ → ℂ) (E : ℂ) (c : Coeff) (μ : ℂ) :
    L6 P E c μ = P (s - μ) * (μ - 1) * (μ - 4) * c μ - c (μ - 9) - E * c (μ - 6) := by
  simp only [L6, PT, xpow, theta]
  ring

/-- The coefficient of `L4 a` at `x^μ`: `P(s - μ) a(μ) + (μ - 7) a(μ - 9) - E a(μ - 6)`,
the paper's recurrence for `L4`. -/
theorem L4_coeff (P : ℂ → ℂ) (E : ℂ) (a : Coeff) (μ : ℂ) :
    L4 P E a μ = P (s - μ) * a μ + (μ - 7) * a (μ - 9) - E * a (μ - 6) := by
  simp only [L4, PT, xpow, theta]
  ring

/-! ### Lemma 2 and Proposition 1 -/

/-- Lemma 2(i): `K_{g-1} c = K_g[(g - 1 - θ/3) c]`. -/
theorem lemma2_i (g : ℂ) (c : Coeff) :
    Kg (g - 1) c = Kg g (fun ν => (g - 1) * c ν - theta c ν / 3) := by
  funext ν
  simp only [Kg, theta, k_sub_one]
  ring

/-- The first step of the paper's proof of (ii): `K[x^3 c] = 3 x^3 K_{-2/3} c`. -/
theorem aux_ii (c : Coeff) : K (xpow 3 c) = fun ν => 3 * xpow 3 (Kg (-2 / 3) c) ν := by
  funext ν
  simp only [K, Kg, xpow]
  rw [show k (1 / 3) ν = k (1 / 3) (ν - 3 + 3) by ring_nf, k_add_three,
    show (1 / 3 : ℂ) - 1 = -2 / 3 by norm_num]
  ring

/-- Lemma 2(ii): `K[x^3 c] = -x^3 K[(θ + 2) c]`, from (i) at `g = 1/3` and `aux_ii`. -/
theorem lemma2_ii (c : Coeff) :
    K (xpow 3 c) = fun ν => -xpow 3 (K (fun μ => theta c μ + 2 * c μ)) ν := by
  rw [aux_ii]
  funext ν
  have h := congrFun (lemma2_i (1 / 3) c) (ν - 3)
  rw [show (1 / 3 : ℂ) - 1 = -2 / 3 by norm_num] at h
  simp only [xpow, K, Kg, theta] at h ⊢
  rw [h]
  ring

/-- The first step of the paper's proof of (iii): `x^6 K[x^-6 c] = (1/9) K_{7/3} c`. -/
theorem aux_iii (c : Coeff) :
    xpow 6 (K (xpow (-6) c)) = fun ν => Kg (7 / 3) c ν / 9 := by
  funext ν
  simp only [K, Kg, xpow]
  rw [show k (7 / 3) ν = k (7 / 3) (ν - 6 + 3 + 3) by ring_nf, k_add_three, k_add_three,
    show (7 / 3 : ℂ) - 1 - 1 = 1 / 3 by norm_num, show ν - 6 - -6 = ν by ring]
  ring

/-- Lemma 2(iii): `K c = x^6 K[x^-6 (θ - 1)(θ - 4) c]`, from (i) at `g = 4/3, 7/3`
and `aux_iii`. -/
theorem lemma2_iii (c : Coeff) :
    K c = xpow 6 (K (xpow (-6)
      (fun μ => theta (fun ρ => theta c ρ - 4 * c ρ) μ - (theta c μ - 4 * c μ)))) := by
  rw [aux_iii]
  funext ν
  have h1 := congrFun (lemma2_i (4 / 3) c) ν
  have h2 := congrFun (lemma2_i (7 / 3) (fun ν => (4 / 3 - 1) * c ν - theta c ν / 3)) ν
  rw [show (4 / 3 : ℂ) - 1 = 1 / 3 by norm_num] at h1
  rw [show (7 / 3 : ℂ) - 1 = 4 / 3 by norm_num] at h2
  simp only [K, Kg, theta] at h1 h2 ⊢
  rw [h1]
  linear_combination h2

/-- Proposition 1: `L4 K = x^6 K x^-6 L6`, for every `P`, `E` and coefficient function. -/
theorem prop1 (P : ℂ → ℂ) (E : ℂ) (c : Coeff) :
    L4 P E (K c) = xpow 6 (K (xpow (-6) (L6 P E c))) := by
  funext ν
  have h1 := m_sub_six ν
  have h2 := m_add_three (ν - 9)
  rw [show ν - 9 + 3 = ν - 6 by ring] at h2
  simp only [m] at h1 h2
  simp only [L4, L6, PT, xpow, K, Kg, theta]
  rw [show ν - 6 - -6 = ν by ring]
  linear_combination (P (s - ν) * c ν) * h1 + c (ν - 9) * h2

/-- Proposition 1 on monomials. -/
theorem prop1_monomial (P : ℂ → ℂ) (E μ : ℂ) :
    L4 P E (K (mono μ)) = xpow 6 (K (xpow (-6) (L6 P E (mono μ)))) :=
  prop1 P E (mono μ)

/-! ### Controls -/

/-- `m(ν) ≠ 0` off the poles of `Γ(1/3 - ν/3)`. -/
lemma m_ne_zero (ν : ℂ) (h : ∀ n : ℕ, (1 / 3 : ℂ) - ν / 3 ≠ -n) : m ν ≠ 0 := by
  unfold m k
  refine mul_ne_zero ?_ (inv_ne_zero (Complex.Gamma_ne_zero h))
  rw [cpow_def_of_ne_zero three_ne]
  exact Complex.exp_ne_zero _

lemma m_zero_ne : m 0 ≠ 0 := by
  apply m_ne_zero
  intro n h
  have := congrArg Complex.re h
  simp at this
  linarith [(Nat.cast_nonneg n : (0 : ℝ) ≤ n)]

lemma m_neg_six_ne : m (-6) ≠ 0 := by
  apply m_ne_zero
  intro n h
  have := congrArg Complex.re h
  norm_num at this
  linarith [(Nat.cast_nonneg n : (0 : ℝ) ≤ n)]

/-- The paper's `P(T) = (T^2 - l0^2)(T^2 - l1^2)`. -/
def Ppaper (l0 l1 : ℂ) (T : ℂ) : ℂ := (T ^ 2 - l0 ^ 2) * (T ^ 2 - l1 ^ 2)

/-- NON-VACUITY: on `x^0`, at exponent 0, with the paper's `P` at `l0 = i, l1 = 2i`
and `E = 1`, the left side of Prop. 1 is `P(5/2) m(0)`, which is not zero. -/
theorem nonvacuity_prop1 :
    L4 (Ppaper I (2 * I)) 1 (K (mono 0)) 0 = Ppaper I (2 * I) (5 / 2) * m 0 ∧
    Ppaper I (2 * I) (5 / 2) * m 0 ≠ 0 := by
  refine ⟨?_, mul_ne_zero ?_ m_zero_ne⟩
  · simp [L4, PT, xpow, K, Kg, theta, mono, s, m]
  · simp only [Ppaper, mul_pow, I_sq]
    norm_num

/-- TAMPER: `(μ - 1)(μ - 5)` in place of `(μ - 1)(μ - 4)` is false. -/
theorem tamper_recurrence : ¬ ∀ μ : ℂ, m μ = (μ - 1) * (μ - 5) * m (μ - 6) := by
  intro h
  have h0 := h 0
  rw [m_sub_six 0] at h0
  have : m (0 - 6) = 0 := by linear_combination -h0
  norm_num at this
  exact m_neg_six_ne this

/-- `L4` with `x^9 (θ + 3)` in place of `x^9 (θ + 2)`. -/
noncomputable def L4tampered (P : ℂ → ℂ) (E : ℂ) (c : Coeff) : Coeff := fun ν =>
  PT P c ν + xpow 9 (fun μ => theta c μ + 3 * c μ) ν - E * xpow 6 c ν

/-- TAMPER: with `x^9 (θ + 3)`, Proposition 1 fails (P = 0, E = 0, on `x^0`, at
exponent 9). -/
theorem tamper_prop1 :
    ¬ ∀ c : Coeff, L4tampered (fun _ => 0) 0 (K c) = xpow 6 (K (xpow (-6) (L6 (fun _ => 0) 0 c))) := by
  intro h
  have h9 := congrFun (h (mono 0)) 9
  have h3 := m_add_three 0
  simp [L4tampered, L6, PT, xpow, K, Kg, theta, mono] at h9
  norm_num at h9
  simp only [m] at h3
  norm_num at h3
  apply m_zero_ne
  simp only [m]
  norm_num
  linear_combination h9 - h3

end IntegrableBoundary.Theorem3Intertwining
