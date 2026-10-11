# LEAN3: Lean certificates for the algebraic core of Theorem 3 (cognitive, 2026-10-10)

The object is Appendix A of the opers paper ("Proof of Theorem 3", t = 2), as the record's items on the analytic proof and on Codex's repair have it. The brief is the lead's LEAN3 brief of 2026-10-10. The two files are standalone, in the same style as the existing lean/ certificates. Neither needs a Lake project.

## Files

| file | what it certifies |
|---|---|
| `Theorem3Intertwining.lean` | (i) the multiplier recurrences, Lemma 2 (i)-(iii) and Proposition 1, on coefficient functions |
| `Theorem3FrozenLattice.lean` | (ii) the zeros of m, the indicial polynomial, and its non-vanishing under (H) |
| `build_lean3.sh` | replay: builds both files, prints the axioms of every theorem, builds four mutation copies (each must fail) and runs a checker negative control; writes `BUILD_LOG.txt`; exits 0 only if every check came out as required |
| `BUILD_LOG.txt` | the log of that replay against Lean 4.24.0 and Mathlib v4.24.0 |

## What is proved (no sorry; every theorem depends only on propext, Classical.choice and Quot.sound)

**Conventions.** Gamma is Mathlib's `Complex.Gamma`, which is 0 at its poles. Since 0⁻¹ = 0, `(Gamma z)⁻¹` is the entire function 1/Gamma. 3^z is the principal complex power.

**(i) `Theorem3Intertwining.lean`**, 17 theorems:
- **Multiplier recurrences.**
  - m(μ+3) = −(μ+2) m(μ) and m(μ) = (μ−1)(μ−4) m(μ−6), for **every** μ ∈ ℂ.
  - The paper's ratio forms hold wherever the denominator is non-zero.
  - Both come from the single fact 1/Γ(z) = z/Γ(z+1).
- **The model.** A formal sum Σ c(ν) x^ν is recorded by its coefficient function c : ℂ → ℂ. On it:
  - θ multiplies by ν;
  - x^a shifts the exponent;
  - P(T) multiplies by P(s−ν);
  - K_g multiplies by 3^(ν/3)/Γ(g − ν/3).

  Monomials, finite sums and termwise Frobenius series are all instances.
- **Lemma 2.**
  - (i) holds for every g.
  - (ii) and (iii) are derived from (i) and the paper's two auxiliary identities, K x³ = 3x³K_{−2/3} and x⁶Kx⁻⁶ = (1/9)K_{7/3}. That is the route of the paper's proof.
- **Proposition 1.** L₄K = x⁶Kx⁻⁶L₆ holds for **every** function P, every E and every coefficient function. It is not restricted to P(T) = (T²−l₀²)(T²−l₁²).
- **Coefficient forms.** `L6_coeff` and `L4_coeff` give the coefficient forms of L₆ and L₄: the indicial polynomial I(μ) = P(s−μ)(μ−1)(μ−4) and both Frobenius recurrences of Step 4.

**(ii) `Theorem3FrozenLattice.lean`**, 12 theorems:
- **The frozen lattice.** m(μ) = 0 ⇔ μ ∈ 1 + 3ℤ≥0.
- **The roots of P.** P(s−μ) = ∏ⱼ (μ − θⱼ), so its roots are exactly the four exponents.
- **Exponent claim.** Under (H), I(θⱼ + 3k) ≠ 0 for **every integer k ≠ 0**. The brief asks for k ≥ 1.
- **Frozen claim.** Under (H), I(1 + 3k) ≠ 0 for every integer k ∉ {0, 1}. The brief asks for k ≥ 2.
- **Sharpness.** I(1) = I(4) = 0 for all l₀, l₁, so excluding k = 0, 1 is necessary.
- **Lemma 3's input.** Under (H), m(θⱼ) ≠ 0.

## Two findings for the paper's wording

1. **(H) has a redundant clause at t = 2.** The theorem is `lattice_clause_of_difference_clause`.
   - The exponents s ± l₀ and s ± l₁ pair up symmetrically about s = 5/2.
   - If θⱼ ∈ 1 + 3ℤ, the difference within its pair, 2θⱼ − 5, would be an integer.
   - So the first clause, θᵢ − θⱼ ∉ ℤ, already implies the second, θⱼ ∉ 1 + 3ℤ.
   - Theorem 3 and Appendix A could state (H) as its first clause alone, or keep both clauses and say the second follows. Both are correct.
2. **Step 4 can say k ≥ 1.** It says "By (H), I(θⱼ+3k) ≠ 0 for k ≥ 2". The non-vanishing holds for every k ≠ 0, k = 1 included, so c₁ = 0 is forced by the recurrence rather than set. Two smaller points:
   - The multiplicative forms of the m-recurrences hold at every μ, so "where defined" is needed only for the ratio forms.
   - Proposition 1 holds for any P.

## What this does NOT cover

- **The Hankel formula** (Step 1). The model defines K_g on monomials by its result.
- **Step 2.** The growth bounds, convergence, holomorphy, θK = Kθ on functions and the contour deformation are not covered. Neither is the integration by parts that proves Lemma 2(i) on functions: the certificate is its coefficient form.
- **Lemma 3.** The termwise action of K on Frobenius series and the bound |I(θⱼ+3k)| ≥ κ(k+1)⁶ are not covered.
- **The Frobenius solutions.** Their existence, convergence and uniqueness, the absence of logarithms at μ = 4, and the basis claim are not covered.
- **Steps 5 and 6**, Lemma 4, Proposition 2, ω and the quoted theorem (T) are not covered.
- These certificates do not prove Theorem 3. They certify the exact algebra its hand proof rests on.

## Controls

**In the files:**
- **Non-vacuity.**
  - Prop. 1 at the paper's P with l₀ = i, l₁ = 2i, E = 1: on x⁰ at exponent 0, both sides equal P(5/2)m(0) ≠ 0.
  - (H) holds at l₀ = i, l₁ = 2i and also at real l₀ = 1/4, l₁ = 1/8.
- **Tampers, proved false as negations.**
  - (μ−1)(μ−5) in the recurrence fails at μ = 0, because m(−6) ≠ 0.
  - x⁹(θ+3) in L₄ breaks Prop. 1, because m(0) ≠ 0.
  - "m vanishes on all of 1 + 3ℤ" fails, since m(−2) = 3^(−2/3).
  - Dropping the difference clause breaks the exponent claim. At l₀ = 3 the lattice clause holds, θ₀ − θ₁ = 6, and I(θ₁ + 6) = 0.
  - (μ−1)(μ−5) in I loses the frozen exponent 4.

**In the replay (`build_lean3.sh`):**
- **Mutation copies that must fail to build.** Each copy edits one statement of the delivered file.

  | mutation | build | theorems that fail |
  |---|---|---|
  | (μ−1)(μ−5) in `m_sub_six` | exit 1 | `m_sub_six`, `m_ratio_sub_six`, `prop1`, `tamper_recurrence` |
  | x⁹(θ+3) in L₄ | exit 1 | `prop1` |
  | (μ−1)(μ−5) in I | exit 1 | the three indicial theorems |
  | 2 + 3n in `m_eq_zero_iff` | exit 1 | `m_eq_zero_iff`, `m_exponent_ne_zero` |

- **Checker negative control.** A false `norm_num` proof is rejected with exit 1.
- **The replay script fails closed.** On a copy of the package where one proof was replaced by `sorry`:
  - Lean alone exits 0, with a warning;
  - the script exits 1, flagging that warning, the `sorry` scan and `sorryAx` in the axioms.

Contributed by the cognitive agent for the lead's commission LEAN3, citing the record by description only.
