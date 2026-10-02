/-
  The three-term ODE and the Gamma-symbol operator: the formal Mellin identities.

  THE CLAIM (the record's item 311 (1), checked there in sympy for deg L = 1, 2, 3).
  With T = s - x d/dx, sigma = n/a, n = a sigma, c = n + sigma, L(T) = prod_r (T - r),
  the Mellin multiplier

      m(nu) = prod_r sigma^(nu/sigma) / Gamma(gamma_r - nu/sigma),
      gamma_r = (s - r)/sigma + (1 - a)/2,

  turns the ODE's Mellin recurrence into the Gamma form's, with
  G_a(T) = prod_r sigma^a Gamma((T - r)/sigma + (1 + a)/2) / Gamma((T - r)/sigma + (1 - a)/2).

  WHAT IS PROVED, for real a, real s, sigma > 0, real nu, and ANY finite list of
  roots r (any degree of L, not only 1-3):
  * (I0) `I0`, `T_monomial`: T(x^(c/2) f) = x^(c/2)((T - c/2) f) at every x > 0
    where f is differentiable, and T x^nu = (s - nu) x^nu;
  * (I1) `I1`: m(mu) = prod_r (sigma gamma_r - mu) m(mu - sigma), UNCONDITIONALLY.
    1/Gamma is entire and Mathlib's Gamma vanishes at the poles, so the
    identity holds as written at every real mu (`inv_Gamma_eq`);
  * (I2) `I2`: m(nu) = G(s - nu) m(nu - n), provided Gamma is finite at
    gamma_r - nu/sigma + a for every root.  These are the poles of the Gamma
    symbol, and the proviso is necessary (control below);
  * (I3) `I3`: L(s - nu + c/2) m(nu - c) = m(nu - n), unconditionally;
  * (I4) `I4`: for ANY P, psi-hat and E, with phi-hat = m psi-hat, the ODE recurrence
    P(s-nu) phi(nu) - L(s-nu+c/2) phi(nu-c) - E phi(nu-n) equals
    m(nu - n) [P G psi(nu) - psi(nu - c) - E psi(nu - n)], under (I2)'s proviso.

  The record states (I1)-(I3) as ratios m(.)/m(.).  Here they are multiplicative,
  which is the form that survives the zeros of m.

  WHAT THIS FILE DOES NOT COVER, deliberately:
  * any CONTOUR or SPECTRAL statement: the passage from these coefficient
    recurrences to solutions, spectral determinants or charges needs the record's
    hypotheses H1-H3 and boundary matching, which item 311 marks CONDITIONAL and
    which are untouched here;
  * the formal-WKB statement that both forms give the same cB x B(A0, B0+1)
    (item 311: NOT proved);
  * the x-space operator L(T) as a differential operator.  (I0) is the
    commutation on functions, and `T_monomial` is its action on monomials; their
    composition over the roots is the record's (I0) corollary and is not stated
    as one theorem;
  * which P, L, a and roots belong to which solution (item 307's table): inputs.

  NEGATIVE CONTROLS, in the file:
  * `control_ordering`: the non-symmetric ordering x^c L(T), i.e. L evaluated at
    s - nu + c, makes (I3) FALSE at an explicit point;
  * `control_pole_proviso`: without (I2)'s proviso, (I2) FAILS at an explicit
    point: the left side is 1/Gamma(-1/2), nonzero, and the right side is 0.

  Contributed by the cognitive agent, 2026-10-02, under the standing
  arrangement, citing the record by description only.
-/
import Mathlib.Analysis.SpecialFunctions.Gamma.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv
import Mathlib.Tactic

namespace IntegrableBoundary.OperMellin

open Real

/-- `1/Γ` is entire: `Γ(x)⁻¹ = x · Γ(x+1)⁻¹` for EVERY real x (at the poles both
    sides are 0 with Mathlib's convention `Γ(-k) = 0`). -/
lemma inv_Gamma_eq (x : ℝ) : (Gamma x)⁻¹ = x * (Gamma (x + 1))⁻¹ := by
  by_cases hx : x = 0
  · subst hx; simp [Real.Gamma_zero]
  · rw [Real.Gamma_add_one hx, mul_inv, ← mul_assoc, mul_inv_cancel₀ hx, one_mul]

variable (a σ s : ℝ)

/-- `γ_r = (s - r)/σ + (1 - a)/2`. -/
noncomputable def gam (r : ℝ) : ℝ := (s - r) / σ + (1 - a) / 2

/-- One root's factor of the Mellin multiplier: `σ^(ν/σ) / Γ(γ_r - ν/σ)`. -/
noncomputable def mf (r ν : ℝ) : ℝ := σ ^ (ν / σ) * (Gamma (gam a σ s r - ν / σ))⁻¹

/-- The multiplier `m(ν) = Π_r σ^(ν/σ) / Γ(γ_r - ν/σ)` over the roots `rs` of `L`. -/
noncomputable def m (rs : List ℝ) (ν : ℝ) : ℝ := (rs.map (fun r => mf a σ s r ν)).prod

/-- `L(T) = Π_r (T - r)`. -/
def L (rs : List ℝ) (T : ℝ) : ℝ := (rs.map (fun r => T - r)).prod

/-- The Gamma symbol `G_a(T) = Π_r σ^a Γ((T-r)/σ + (1+a)/2) / Γ((T-r)/σ + (1-a)/2)`. -/
noncomputable def G (rs : List ℝ) (T : ℝ) : ℝ :=
  (rs.map (fun r => σ ^ a * Gamma ((T - r) / σ + (1 + a) / 2) * (Gamma ((T - r) / σ + (1 - a) / 2))⁻¹)).prod

variable {a σ s}

/-- **(I1)**, multiplicative and unconditional:
    `m(μ) = Π_r (σ γ_r - μ) · m(μ - σ)`. -/
theorem I1 (hσ : 0 < σ) (rs : List ℝ) (μ : ℝ) :
    m a σ s rs μ = (rs.map (fun r => σ * gam a σ s r - μ)).prod * m a σ s rs (μ - σ) := by
  induction rs with
  | nil => simp [m]
  | cons r rs ih =>
    simp only [m, List.map_cons, List.prod_cons] at ih ⊢
    rw [ih]
    have key : mf a σ s r μ = (σ * gam a σ s r - μ) * mf a σ s r (μ - σ) := by
      unfold mf
      have e1 : gam a σ s r - (μ - σ) / σ = (gam a σ s r - μ / σ) + 1 := by field_simp; ring
      rw [e1, inv_Gamma_eq (gam a σ s r - μ / σ)]
      have e2 : σ ^ (μ / σ) = σ ^ ((μ - σ) / σ) * σ := by
        rw [← Real.rpow_add_one hσ.ne']; congr 1; field_simp; ring
      rw [e2]
      have : σ * gam a σ s r - μ = σ * (gam a σ s r - μ / σ) := by field_simp
      rw [this]; ring
    rw [key]; ring

/-- **(I2)**: `m(ν) = G(s - ν) · m(ν - n)`, `n = aσ`, provided Γ is finite at
    `γ_r - ν/σ + a` for every root (the poles of the Gamma symbol). -/
theorem I2 (hσ : 0 < σ) (rs : List ℝ) (ν : ℝ)
    (hpole : ∀ r ∈ rs, Gamma (gam a σ s r - ν / σ + a) ≠ 0) :
    m a σ s rs ν = G a σ rs (s - ν) * m a σ s rs (ν - a * σ) := by
  induction rs with
  | nil => simp [m, G]
  | cons r rs ih =>
    simp only [m, G, List.map_cons, List.prod_cons] at ih ⊢
    rw [ih (fun r' hr' => hpole r' (List.mem_cons_of_mem r hr'))]
    have hp := hpole r (List.mem_cons_self)
    have key : mf a σ s r ν = (σ ^ a * Gamma ((s - ν - r) / σ + (1 + a) / 2)
        * (Gamma ((s - ν - r) / σ + (1 - a) / 2))⁻¹) * mf a σ s r (ν - a * σ) := by
      unfold mf
      have e1 : (s - ν - r) / σ + (1 + a) / 2 = gam a σ s r - ν / σ + a := by unfold gam; field_simp; ring
      have e2 : (s - ν - r) / σ + (1 - a) / 2 = gam a σ s r - ν / σ := by unfold gam; field_simp; ring
      have e3 : gam a σ s r - (ν - a * σ) / σ = gam a σ s r - ν / σ + a := by field_simp; ring
      have e4 : σ ^ (ν / σ) = σ ^ a * σ ^ ((ν - a * σ) / σ) := by
        rw [← Real.rpow_add hσ]; congr 1; field_simp; ring
      rw [e1, e2, e3, e4]
      have hc : Gamma (gam a σ s r - ν / σ + a) * (Gamma (gam a σ s r - ν / σ + a))⁻¹ = 1 :=
        mul_inv_cancel₀ hp
      linear_combination (-(σ ^ a * σ ^ ((ν - a * σ) / σ) * (Gamma (gam a σ s r - ν / σ))⁻¹)) * hc
    rw [key]; ring

/-- **(I3)**, unconditional: `L(s - ν + c/2) · m(ν - c) = m(ν - n)`, with
    `n = aσ`, `c = n + σ`. -/
theorem I3 (hσ : 0 < σ) (rs : List ℝ) (ν : ℝ) :
    L rs (s - ν + (a * σ + σ) / 2) * m a σ s rs (ν - (a * σ + σ)) = m a σ s rs (ν - a * σ) := by
  rw [I1 hσ rs (ν - a * σ)]
  have e : ν - a * σ - σ = ν - (a * σ + σ) := by ring
  rw [e]; congr 1
  unfold L; congr 1; apply List.map_congr_left; intro r _
  unfold gam; field_simp; ring

/-- **(I4)**: with `φ̂ = m ψ̂`, the ODE's Mellin recurrence equals `m(ν - n)` times the
    Gamma form's recurrence `P G ψ̂(ν) - ψ̂(ν - c) - E ψ̂(ν - n)`, for ANY `P`, `ψ̂`, `E`
    (same pole proviso as (I2)). -/
theorem I4 (hσ : 0 < σ) (rs : List ℝ) (P ψ : ℝ → ℝ) (E ν : ℝ)
    (hpole : ∀ r ∈ rs, Gamma (gam a σ s r - ν / σ + a) ≠ 0) :
    P (s - ν) * (m a σ s rs ν * ψ ν)
      - L rs (s - ν + (a * σ + σ) / 2) * (m a σ s rs (ν - (a * σ + σ)) * ψ (ν - (a * σ + σ)))
      - E * (m a σ s rs (ν - a * σ) * ψ (ν - a * σ))
    = m a σ s rs (ν - a * σ) *
        (P (s - ν) * G a σ rs (s - ν) * ψ ν - ψ (ν - (a * σ + σ)) - E * ψ (ν - a * σ)) := by
  have h2 := I2 hσ rs ν hpole
  have h3 := I3 (a := a) (s := s) hσ rs ν
  rw [h2, ← mul_assoc (L rs _), h3]; ring

/-- **(I0)**: `T (x^(c/2) f) = x^(c/2) ((T - c/2) f)` for `T = s - x d/dx`, at any
    `x > 0` where f is differentiable. -/
theorem I0 (c : ℝ) (f : ℝ → ℝ) (x : ℝ) (hx : 0 < x) (f' : ℝ) (hf : HasDerivAt f f' x) :
    ∃ d, HasDerivAt (fun y => y ^ (c / 2) * f y) d x ∧
      s * (x ^ (c / 2) * f x) - x * d = x ^ (c / 2) * ((s - c / 2) * f x - x * f') := by
  have hp := Real.hasDerivAt_rpow_const (p := c / 2) (Or.inl hx.ne')
  refine ⟨_, hp.mul hf, ?_⟩
  have e : x * (c / 2 * x ^ (c / 2 - 1)) = c / 2 * x ^ (c / 2) := by
    rw [Real.rpow_sub_one hx.ne']; field_simp
  linear_combination (-(f x)) * e

/-- `T` acts on monomials by `T x^ν = (s - ν) x^ν`, so
    `x^(c/2) L(T) x^(c/2) x^ν = L(s - ν - c/2) x^(ν + c)` term by term. -/
theorem T_monomial (ν x : ℝ) (hx : 0 < x) :
    HasDerivAt (fun y => y ^ ν) (ν * x ^ (ν - 1)) x ∧ s * x ^ ν - x * (ν * x ^ (ν - 1)) = (s - ν) * x ^ ν := by
  refine ⟨Real.hasDerivAt_rpow_const (Or.inl hx.ne'), ?_⟩
  rw [Real.rpow_sub_one hx.ne']; field_simp

/-! ## Negative controls -/

/-- NEGATIVE CONTROL: the symmetric ordering `x^(c/2) L(T) x^(c/2)` is load-bearing.
    With `x^c L(T)` the polynomial is evaluated at `s - ν + c` instead of `s - ν + c/2`,
    and that version of (I3) is FALSE: one root r = 0, a = σ = 1, s = 0, ν = 1/2
    (there `m(ν - n) = 1/Γ(1/2) ≠ 0`, the correct factor is 1/2 and the wrong one 3/2). -/
theorem control_ordering :
    L [0] ((0:ℝ) - 1/2 + (1 * 1 + 1)) * m 1 1 0 [0] (1/2 - (1 * 1 + 1)) ≠ m 1 1 0 [0] (1/2 - 1 * 1) := by
  intro hbad
  have h := I3 (a := 1) (σ := 1) (s := 0) one_pos [0] (1/2)
  have hm : m 1 1 0 [0] (1/2 - 1 * 1) ≠ 0 := by
    simp only [m, mf, gam, List.map_cons, List.map_nil, List.prod_cons, List.prod_nil]
    norm_num
    exact (Real.Gamma_pos_of_pos (by norm_num)).ne'
  norm_num [L] at h hbad
  apply hm
  linarith

/-- NEGATIVE CONTROL: the pole proviso of (I2) is needed.  At a = 1/2, σ = 1, s = 0,
    one root r = 0, ν = 3/4: the Gamma symbol's numerator sits at the pole Γ(0), so the
    right side of (I2) is 0, while the left side is 1/Γ(-1/2) ≠ 0. -/
theorem control_pole_proviso :
    m (1/2) 1 0 [0] (3/4) ≠ G (1/2) 1 [0] (0 - 3/4) * m (1/2) 1 0 [0] (3/4 - 1/2 * 1) := by
  have hG : G (1/2) 1 [0] (0 - 3/4) = 0 := by
    simp only [G, List.map_cons, List.map_nil, List.prod_cons, List.prod_nil]
    norm_num [Real.Gamma_zero]
  rw [hG, zero_mul]
  simp only [m, mf, gam, List.map_cons, List.map_nil, List.prod_cons, List.prod_nil]
  norm_num
  apply Real.Gamma_ne_zero
  intro k hk
  have h2 : (2 * k : ℝ) = 1 := by linarith
  have : 2 * k = 1 := by exact_mod_cast h2
  omega

end IntegrableBoundary.OperMellin
