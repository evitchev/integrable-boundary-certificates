/-
Kernel-checked skeletons for the integrable-boundary draft
(draft/main.tex): Lemma 1 (the bootstrap lemma, Section 2) and the
degree-bound interpolation argument of R3.

Honest scope — division of labor with the exact-arithmetic engines:

* Lean certifies the ABSTRACT LOGICAL SKELETONS: that the bootstrap
  lemma's conclusion follows from its hypotheses in any graded Lie
  algebra, and that a polynomial of bounded degree vanishing at more
  points than its degree is identically zero.
* The sympy engines certify the INSTANCES: that the physical charge
  spaces Q^sigma (local functionals modulo total derivatives, with
  their exact OPE bracket) satisfy bracket-additivity of the grading,
  and that the even-spin vacancy and degree-bound hypotheses hold in
  the verified ranges (code/vacancy.py, code/certifications.json).

Neither replaces the other: the free-boson realization is not
formalized here, and the abstract lemmas prove nothing about it until
the engine-certified hypotheses are supplied.

Checked with Lean 4.24 / Mathlib via the pinned local toolchain
(code/lean_check.py; registered as `lean_certificates` in
code/certifications.json).
-/
import Mathlib.Algebra.Lie.Basic
import Mathlib.Algebra.Polynomial.Roots

namespace IntegrableBoundary

variable {R L : Type*} [CommRing R] [LieRing L] [LieAlgebra R L]

/-- A grading of the Lie algebra `L` by submodules indexed by spin,
additive under the bracket: `⁅G i, G j⁆ ⊆ G (i + j)`. This is the
only structural input the bootstrap lemma uses; in the application
`G σ` is the genuine charge space `Q^σ` and the bracket is the exact
OPE commutator. -/
def BracketGraded (G : ℤ → Submodule R L) : Prop :=
  ∀ ⦃i j : ℤ⦄ ⦃x y : L⦄, x ∈ G i → y ∈ G j → ⁅x, y⁆ ∈ G (i + j)

/-- Even-spin vacancy at one grade: the only element of `G σ` killed
by `ad x` is zero — `ker (ad x) ∩ G σ = 0` in the draft's notation. -/
def VacantAt (G : ℤ → Submodule R L) (x : L) (σ : ℤ) : Prop :=
  ∀ y ∈ G σ, ⁅x, y⁆ = 0 → y = 0

/-- **The bootstrap lemma, core step** (draft Lemma 1, proof body):
if `a` and `b` both commute with `x`, then `⁅a, b⁆` commutes with `x`
by the Jacobi identity; if the grade `i + j` is vacant for `x`, the
commutator therefore vanishes.

Note the hypotheses are exactly the draft's: `x`'s own grade never
enters, and no parity is needed at this level — parity only selects
WHICH grades must be vacant (see `bootstrap_odd`). -/
theorem bootstrap {G : ℤ → Submodule R L} (hG : BracketGraded G)
    (x : L) {i j : ℤ} {a b : L}
    (ha : a ∈ G i) (hb : b ∈ G j)
    (hxa : ⁅x, a⁆ = 0) (hxb : ⁅x, b⁆ = 0)
    (hvac : VacantAt G x (i + j)) :
    ⁅a, b⁆ = 0 := by
  refine hvac _ (hG ha hb) ?_
  rw [leibniz_lie, hxa, hxb, zero_lie, lie_zero, add_zero]

/-- **Draft Lemma 1, as stated**: fix `x` (the distinguished `I₃`) and
suppose even-spin vacancy holds on a range `E`. Then any two odd-spin
charges commuting with `x`, whose spins sum into `E`, commute with
each other. -/
theorem bootstrap_odd {G : ℤ → Submodule R L} (hG : BracketGraded G)
    (x : L) (E : Set ℤ)
    (hvac : ∀ σ ∈ E, Even σ → VacantAt G x σ)
    {σ₁ σ₂ : ℤ} (h₁ : Odd σ₁) (h₂ : Odd σ₂) (hE : σ₁ + σ₂ ∈ E)
    {a b : L} (ha : a ∈ G σ₁) (hb : b ∈ G σ₂)
    (hxa : ⁅x, a⁆ = 0) (hxb : ⁅x, b⁆ = 0) :
    ⁅a, b⁆ = 0 :=
  bootstrap hG x ha hb hxa hxb (hvac _ hE (h₁.add_odd h₂))

/-- **R3's degree-bound interpolation argument, generic form**: over
an integral domain, a polynomial of degree at most `d` vanishing at
every point of a set with more than `d` elements is identically
zero. -/
theorem vanishes_identically {F : Type*} [CommRing F] [IsDomain F]
    (p : Polynomial F) (d : ℕ) (hdeg : p.natDegree ≤ d)
    (s : Finset F) (hcard : d < s.card)
    (heval : ∀ t ∈ s, p.eval t = 0) : p = 0 :=
  Polynomial.eq_zero_of_natDegree_lt_card_of_eval_eq_zero' p s heval
    (lt_of_le_of_lt hdeg hcard)

/-- **R3, as applied**: each entry of the reduced `[I₅, I₇]`
commutator is a polynomial in `n` of degree at most 14 after clearing
denominators; exact vanishing at 20 distinct rational points therefore
forces identical vanishing. (The degree bound and the 20 evaluations
are the engine-certified inputs; see `code/vacancy.py`.) -/
theorem commutator_entry_vanishes (p : Polynomial ℚ)
    (hdeg : p.natDegree ≤ 14) (s : Finset ℚ) (hcard : 20 ≤ s.card)
    (heval : ∀ t ∈ s, p.eval t = 0) : p = 0 :=
  vanishes_identically p 14 hdeg s
    (lt_of_lt_of_le (by norm_num) hcard) heval

end IntegrableBoundary
