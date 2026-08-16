# Generic-N spectrum of the A₂⁽²⁾ family, and the branch-point test (2026-07-28)

Computed with the symbolic-N invariant engine (`code/invariant_engine.py`,
validated against the component engine at N = 2, 3 and against eqs.
(28)–(34) of hep-th/0404195; scan: `code/a2_generic_scan.py`, 16 cores).

Genuine kernel dimensions of ad_{I₅^{A₂}} per charge spin σ, on the
double cover s² = (N−25)(N−1) rationally parametrized by
N = (t²−25)/(t²−1), s = −24t/(t²−1):

| point | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| t=2 (N=−7) | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 | 0 | 1 |
| t=3 (N=−2) | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 | 0 | 1 |
| t=7 (N=1/2) | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 | 0 | 1 |
| t=9/2 (N=−19/77) | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 | 0 | 1 |
| t=−2 (sheet 2) | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 | 0 | 1 |
| **BRANCH N=25** | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 | 0 | 1 |
| **BRANCH-coeff N=1 (formal)** | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 | 0 | 1 |
| dim Q^σ | 0 | 2 | 0 | 4 | 1 | 8 | 4 | 15 | 11 | 31 | 25 | 59 |

## Conclusions

1. **Generic-N spin content = 6n±1, proven at the sampled level and
   generic by semicontinuity.** Vacancy at σ = 3, 9 (and all even σ) at
   any single sample point implies generic vacancy on the irreducible
   rational parameter curve. One-dimensional kernels at σ = 5, 7, and — beyond the paper's
   data — **σ = 11 and 13** (the next A₂⁽²⁾ exponents). PROMOTED
   2026-07-29: the σ = 11, 13 charges are constructed as explicit
   ℚ(t)-symbolic densities with commutators verified identically
   (`promote_qt.py`, `qt_promotion.md`): generic existence proven,
   generic kernels exactly one-dimensional.
2. **The branch-point hypothesis is REFUTED.** At N = 25 the formal
   invariant algebra is faithful (letter counts ≪ 25, no second-
   fundamental-theorem relations), and there is **no enhancement**;
   likewise for the branch coefficients (N, s) = (1, 0) in the formal
   algebra. Hence the physical N = 1 tower enhancement (spin-3 and
   spin-9 charges appearing; see the note §4c) is caused **entirely by
   the finite-N relations** of the one-boson invariant algebra (e.g.
   (22)(11) = (12)(12) at N = 1), not by the ramification of the
   double cover. Fibers of the commuting-tower family jump at
   *small-N specialization loci*, not at the Langlands-self-dual
   branch points. (N = 1 happens to be both, which is what misled the
   parallel investigation.)
3. dim Q^σ reproduces the paper's Table 1 (2, 4, 8, 15, 31 at
   σ = 3, 5, 7, 9, 11) and extends it: dim Q^13 = 59.
4. Open identification: the enhanced all-odd c = 1 tower through
   I₃' = (11)(11) + 5/2 (22) is a genuine one-boson integrable
   structure distinct from KdV; a natural candidate home is the
   quantum-AKNS / W∞ menagerie of Kotousov–Lukyanov (1910.05947).

## Certification of the paper's I₇ formulas (by-product)

Solving ker ad at spin 7 from scratch (`code/diagnose_I7_formulas.py`):

- **Eq. (33)** (O(N) family I₇): fully correct, including the N³ term.
- **Eq. (30)** (KdV I₇): correct except the P₈⁽⁴⁾=(44) term: printed
  −(2/**657**)(1464+59N+N²), must be −(2/**675**)(1464+59N+N²).
  Independently forced by KdV = O(N) at N = 1.
- **Eq. (35)** (A₂⁽²⁾ I₇): two sign misprints, both of the same
  dropped-minus type, in the √-terms:
  P₈⁽³⁾ line: (315N−441)√… must be −(315N+441)√…;
  P₈⁽⁴⁾ line: (315N−2541)√… must be −(315N+2541)√….
  All other terms correct (an apparent discrepancy in the P₈⁽¹⁾ line
  was a transcription slip on our side; the paper is right there).
- **Eqs. (36)–(37) with (43)–(45)**: correct as printed (the
  (315N−2499) in eq. (37) is genuine) — verified by expanding the
  T-composites and comparing with the solved I₇ mod ∂.

Consolidated errata for hep-th/0404195 so far: eq. (50) (P₆→P₄
notation), eq. (51) (missing factor n in the P₆⁽¹⁰⁾ coefficient),
eq. (30) (657→675), eq. (35) (two dropped minus signs).
