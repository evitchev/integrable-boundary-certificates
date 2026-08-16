# Spin-content diagnosis of the three Z₂×O(N) cylindrical families (2026-07-28)

Computed with the mixed symbolic-N engine (`code/mixed_engine.py`: one
Z₂-even boson X + N bosons Y via O(N) pairs, index-graph Wick rules;
validated against the component engine at N = 1, 2), scan:
`code/cyl_diagnosis.py`, 16 cores, ~3.5 min.

Genuine kernel dims of ad_{I₃} per charge spin σ, on the double cover
s² = (N−25)(N−1) (rational parametrization by t; both sheets):

| Solution | point | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | all four points | 0 | 1 | 0 | 1 | 0 | 1 | 0 | 1 | 0 | 1 |
| 2 | t=2, t=7, t=−2 | 0 | 1 | 0 | 1 | 0 | 1 | 0 | 1 | 0 | 1 |
| 2 | **t=3 (N=−2)** | 0 | 2 | 0 | 2 | 0 | 2 | 0 | 2 | 0 | 2 |
| 3 | t=2, t=7, t=−2 | 0 | 1 | 0 | 1 | 0 | 1 | 0 | 1 | 0 | 1 |
| 3 | **t=3 (N=−2)** | 0 | 2 | 0 | 2 | 0 | 2 | 0 | 2 | 0 | 2 |

dim Q^σ (Z₂×O(N) sector): 0, 5, 1, 12, 7, 29, 24, 68, 71, 160.

## Conclusions

1. **All three cylindrical families carry all-odd multiplicity-one
   centralizer spectra**: spins
   {1, 3, 5, 7, 9, 11} present with multiplicity one, even spins
   vacant, at generic N on both sheets. Existence and uniqueness of
   I₇, I₉, I₁₁ is new — the paper had only I₃, I₅ ("it seems plausible
   that these are the first members of infinite families").
   PROMOTED 2026-07-29: I₇, I₉, I₁₁ constructed as ℚ(t)-symbolic
   densities, commutators verified identically (`qt_promotion.md`):
   generic existence and exact multiplicity one proven.
2. **Langlands diagnosis**: all-odd exponent content — the pattern of
   A₁⁽¹⁾ (and of the D-type twisted chains in the
   Litvinov–Vilkoviskiy BCD boundary-label framework, 2105.04018) —
   and decisively NOT A₂⁽²⁾-type (which skips 3, 9). This confirms the
   "generalized paperclip" picture for all three cylindrical families
   and sharpens the BCD-identification target: they should sit among
   the D⁽²⁾/A⁽¹⁾-type boundary-label commutants.
3. Each charge commutes with I₃ by construction; higher-pair sums
   exceed the verified vacancy range, so mutual commutativity is
   checked directly: [I₅,I₇] = [I₅,I₉] = 0 at two curve points
   (cyl_mutual_check.py); [I₅,I₁₁] = [I₇,I₉] = 0 at 70–71 rational
   points per solution from the stored ℚ(t) sections
   (mutual_checks.md, 2026-07-29/30). Remaining: [I₇,I₁₁], [I₉,I₁₁]
   (OOM'd; pending the collapsed enumerator).
4. **Jumping fiber at the decoupling locus, fully explained.** At
   t=3 (N=−2, s=−9) the kernel dims double for Solutions 2 and 3 —
   exactly where their X–Y coupling coefficients vanish:
   Sol 2 has (11+N+s)/8 = 0 and Sol 3 has (7−N+s) = 0 there, and
   N=−2 is the unique common zero (subtracting the two conditions
   gives 4+2N=0). With the coupling off, the density decouples into
   X-only + Y-only pieces and the commutant factorizes — an extra
   independent charge per spin. Another instance of the corrected
   picture from `a2_generic_scan.md`: fibers of the commuting-tower
   family jump at *degeneration loci* (small-N relations, decoupling
   points), not at the branch points of the double cover.

## Normalization correction to eqs. (53)–(58)

Reading the printed "4T_Y²" with eq. (52)'s composite
(2T_Y² = (11)(11)_Y − 2(22)_Y) gives densities that (i) fail the
paper's own paperclip-reduction claims at N=1 and (ii) admit NO
commuting I₅. The Y-part must be read as (11)(11)_Y − 2(22)_Y (half).
With that fix, certified: at N=1 (physical, component engine algebra)
the I₃ of Solutions 1 and 3 equals 6·P₄^paperclip(n=−1) and that of
Solution 2 equals 15·P₄^paperclip(n=2), coefficient by coefficient —
precisely the reductions stated below eqs. (54), (56), (58) — and all
three I₃'s then admit unique commuting towers. Practically this is an
erratum/clarification for §3.3 (a factor-of-2 normalization of the
T_Y² composite).
