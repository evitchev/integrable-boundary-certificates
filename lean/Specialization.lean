/-
Lean v2 certificates: the one-point-vacancy / specialization lemma
and the rational-section vanishing lemma — the abstract skeletons
behind every GF(p)/Q(t) high-spin certificate (code/vacancy_high.py,
code/vacancy_final_modp.py, code/promote_qt.py):

  "Empty kernel at ONE curve point => generically empty
   (semicontinuity) => the commutator, a fixed rational section,
   vanishing generically, vanishes IDENTICALLY."

Division of labor as in Bootstrap.lean: Lean certifies the abstract
implications; the engines certify the instances (the polynomial
ad-matrices after denominator clearing, the single-point kernel
verdicts — over Q directly or rigorously pinned via GF(p) tag/rank
equality — and the degree data). The free-boson realization is not
formalized.
-/
import Mathlib.Algebra.Polynomial.Roots
import Mathlib.Algebra.Polynomial.Div
import Mathlib.Data.Matrix.Mul
import Mathlib.FieldTheory.RatFunc.AsPolynomial

namespace IntegrableBoundary

open Polynomial

/-! ## Rational-section vanishing

A rational function whose numerator degree is bounded by `d`,
evaluated at more than `d` points where its denominator does not
vanish, all evaluations zero, is the zero function — the "vanishing
generically ⟹ vanishing identically" step in the form the engines
use it (denominators cleared, numerator degree bounded). -/

theorem ratfunc_vanishes_identically {K : Type*} [Field K]
    (x : RatFunc K) (d : ℕ) (hdeg : x.num.natDegree ≤ d)
    (s : Finset K) (hcard : d < s.card)
    (heval : ∀ t ∈ s, Polynomial.eval t x.denom ≠ 0 ∧
      RatFunc.eval (RingHom.id K) t x = 0) :
    x = 0 := by
  rw [← RatFunc.num_eq_zero_iff]
  refine Polynomial.eq_zero_of_natDegree_lt_card_of_eval_eq_zero' _ s
    (fun t ht => ?_) (lt_of_le_of_lt hdeg hcard)
  obtain ⟨hden, hzero⟩ := heval t ht
  simp only [RatFunc.eval, Polynomial.eval₂_id] at hzero
  rcases div_eq_zero_iff.mp hzero with hnum | hden'
  · exact hnum
  · exact absurd hden' hden

/-! ## The specialization lemma (one-point vacancy)

Let `M` be a matrix of polynomials and `t₀` a point at which the
entrywise-evaluated matrix `M(t₀)` has trivial kernel. Then `M` has
trivial kernel over the polynomial ring — hence over the rational
function field after clearing denominators, which is how the engines
present kernel vectors.

Proof: infinite descent in divisibility by `(X - C t₀)`. Any
polynomial kernel vector specializes into the kernel of `M(t₀)`, so
it vanishes entrywise at `t₀`, so `(X - C t₀)` divides every entry;
dividing it out yields another kernel vector, so every power divides
every entry, which the degree bound forbids for a nonzero entry. -/

theorem kernel_specialization {R : Type*} [CommRing R] [IsDomain R]
    {m n : Type*} [Fintype n]
    (M : Matrix m n (Polynomial R)) (t₀ : R)
    (hspec : ∀ w : n → R,
      (M.map (Polynomial.eval t₀)).mulVec w = 0 → w = 0)
    (v : n → Polynomial R) (hv : M.mulVec v = 0) :
    v = 0 := by
  have key : ∀ k : ℕ, ∀ u : n → Polynomial R, M.mulVec u = 0 →
      ∀ j, (X - C t₀) ^ k ∣ u j := by
    intro k
    induction k with
    | zero => intro u _ j; simp
    | succ k ih =>
      intro u hu j
      have hdvd : ∀ j, (X - C t₀) ^ k ∣ u j := ih u hu
      choose w hw using hdvd
      have hpow_ne : ((X - C t₀) ^ k : Polynomial R) ≠ 0 :=
        pow_ne_zero _ (X_sub_C_ne_zero t₀)
      -- w is again a kernel vector.
      have hw_ker : M.mulVec w = 0 := by
        funext i
        have hfactor : M.mulVec u i =
            (X - C t₀) ^ k * M.mulVec w i := by
          simp only [Matrix.mulVec, dotProduct]
          rw [Finset.mul_sum]
          exact Finset.sum_congr rfl fun j' _ => by rw [hw j']; ring
        have hzero : (X - C t₀) ^ k * M.mulVec w i = 0 := by
          rw [← hfactor, hu]; rfl
        rcases mul_eq_zero.mp hzero with h | h
        · exact absurd h hpow_ne
        · exact h
      -- The specialization of w lies in the trivial kernel.
      have hw_eval : (fun j => Polynomial.eval t₀ (w j)) = 0 := by
        apply hspec
        funext i
        have hcomp :
            (M.map (Polynomial.eval t₀)).mulVec
              (fun j => Polynomial.eval t₀ (w j)) i =
              Polynomial.eval t₀ (M.mulVec w i) := by
          simpa [Polynomial.coe_evalRingHom, Function.comp] using
            (RingHom.map_mulVec (Polynomial.evalRingHom t₀) M w i).symm
        rw [hcomp, hw_ker]
        simp
      -- Hence (X - C t₀) divides each w j, giving the next power.
      have hroot : (X - C t₀) ∣ w j := by
        rw [dvd_iff_isRoot, IsRoot.def]
        exact congrFun hw_eval j
      obtain ⟨c, hc⟩ := hroot
      exact ⟨c, by rw [hw j, hc, pow_succ]; ring⟩
  funext j
  by_contra hne
  have hdvd := key ((v j).natDegree + 1) v hv j
  have hle := Polynomial.natDegree_le_of_dvd hdvd hne
  rw [Polynomial.natDegree_pow, Polynomial.natDegree_X_sub_C,
    mul_one] at hle
  omega

end IntegrableBoundary
