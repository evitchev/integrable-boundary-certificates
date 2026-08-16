# Super-KdV calibration (audit-corrected)

**Phase 2 of the SUSY-generalization program**
(`[private-note]`).  Engine:
`code/super_engine.py` (one N=1 superfield, c = 3/2); scripts
`code/super_kdv.py` and `code/super_screening_fit.py`; artifacts
`super_kdv_calibration[_naive][_ext].json`,
`super_screening_fit[_naive].json`.

**Correction record (2026-08-04).**  An external audit (OpenAI Codex)
found that the first committed version of this calibration built its
composite pencil with the concatenation product instead of the
point-split normal product :AB:(w) = oint dz (z-w)^-1 A(z)B(w).  The
two differ by non-d-exact terms, so the pencil searched was not the
super-Virasoro composite slice, the charge it selected was NOT the
standard quantum super-KdV I3, and the previously reported screening
quartic (with its "erratum to the phase plan") mislabeled the tower.
Every substantive audit claim was re-verified independently with this
repo's own engine before this rewrite; the phase plan's original
alpha^2 = 1 expectation was correct for the standard hierarchy.  The
first version's tower survives as a distinct exactly-commuting family
(section 4).

## 1. The charge I3 (known) and its selection (the calibration point)

The standard quantum super-KdV I3 is prior art, and this calibration
reproduces it from scratch:

    I3 = oint [ :TT: - (1/4) :G dG: ],

with the fermionic-screening commutant characterization, appears
explicitly in Chistyakova–Litvinov–Orlov, arXiv:2110.05870 ("one can
define the set of local Integrals of Motion I_{2k-1} as a commutant of
a pair of screening fields S± = oint Psi e^{±b Phi} dz", with
I3 = (1/2pi) int (T^2 - (1/4) G G') dx; verified against the paper
text), building on the SUSY N=1 KdV hierarchy construction of
Kulish–Zeitlin — "Superconformal field theory and SUSY N=1 KdV
hierarchy" I (hep-th/0407154) and II (hep-th/0501019), the papers
CLO cite for the hierarchy (their refs [6], [7], verified against
the extracted bibliography; an earlier revision cited the related
hep-th/0312159 instead).  Our point-split density
representative:

    (1/4)(dX)^4 - (3/4)(dX)^2 psi dpsi + (3/8) dX d3X
    - (1/4)(d2X)^2 - (1/8) psi d3psi + (1/2) dpsi d2psi.

Selection inside the pencil :TT: + kappa :G dG: at the free point
c = 3/2 is subtler than expected, and neither of the two "obvious"
principles works:

- ad_G on charges is a differential (ad_G^2 = ad_{oint T} = 0), so
  its kernel contains every G-exact charge and cannot select a ray
  (the ad_G kernel dims are recorded in the JSONs as G-cohomological
  data);
- the dual screening pair psi e^{±X} annihilates EVERY member of the
  point-split pencil (all three pencil conditions carry
  (alpha-1)(alpha+1) factors, and both generators screen exactly at
  alpha = ±1), so screenability is kappa-blind here.

What does select the ray is **G-invariance of the density**:
Res(G, :TT: + kappa :G dG:) is a total derivative iff kappa = -1/4
(unique solution of an exact linear system).

## 2. Centralizer spectrum through sigma = 11 (computed range only)

Full centralizer of I3 in the free super-Fock space: densities of
doubled weight w2 mod total derivatives, sigma = w2/2 - 1 (odd w2 =
fermionic charges, half-integer sigma).  Within the computed range
sigma <= 11 (`_ext.json`):

    dim ker ad_I3 = 1  at sigma in {1/2, 1, 3, 5, 7, 9, 11}
    dim ker ad_I3 = 0  at every other sigma <= 11

(the joint (ad_G, ad_I3) table is identical except sigma = 1/2, whose
kernel charge is oint G itself).  Within this range:

1. the centralizer coincides with the super-KdV tower
   {oint G, oint T, I3, I5, I7, I9, I11} — multiplicity one, and even
   oint psi dpsi is excluded at sigma = 1;
2. even-spin vacancy persists from the bosonic families (sigma = 2,
   4, 6, 8, 10 vacant) — the bootstrap-lemma hypothesis pattern;
3. half-integer spins above 1/2 are vacant (3/2 ... 21/2), so no
   independent fermionic partner charges appear in range:
   supersymmetry acts on the tower by invariance, not pairing.

These are finite-range computations.  No all-spin character,
cohomology, or centralizer theorem is claimed; the graded Jacobi
argument propagates commutativity only across the verified vacancies.

Direct exact verifications (beyond the defining kernel conditions):
[I5,I7], [I5,I9], [I7,I9], [G,I5], [G,I7], [I3,I5], [I3,I7] are all
total derivatives.  I5, I7, I9 densities are in the JSONs; given the
Kulish–Zeitlin/CLO prior art on the hierarchy itself, none of these
are claimed as new charges — the new items here are the explicit
free-field spectral table and the selection analysis above.

## 3. Screening locus alpha^2 = 1 (generic-stratum statement)

For the corrected tower the screening locus — all alpha with
Res(P, psi e^{alpha X}) a twisted total derivative
(D_alpha J = dJ + alpha (dX) J) simultaneously for P = I3, I5, I7 —
is alpha^2 = 1: every screening condition of I5 (7 conditions) and I7
(19 conditions, `--deep`) is divisible by (alpha^2 - 1), and
alpha = 0 is excluded.  This matches the CLO characterization at the
free point: in their conventions c-hat = 1 + 2Q^2 with Q = b + 1/b
(verified against the paper text), so c-hat = 1 forces Q = 0 and
b = -1/b = +-i; with X = i Phi the code exponent is alpha = -ib, and
the dual pair {b, 1/b} = {i, -i} maps to the two screenings
alpha = +-1.  This confirms the phase plan's original expectation.  Stage-0 validation: G and T
impose no alpha-conditions (top components of superfields are
G-invariant for every alpha), so the twisted machinery's Koszul signs
are independently exercised.

Two caveats.  (i) The elimination runs over QQ(alpha) with pivots
that may depend on alpha — a generic localization, so even
substituting the claimed points into its output is not branch-safe;
locus statements are generic-stratum statements, and the claimed
points are certified separately by exact rational-point membership
(integer elimination at fixed alpha): both pencil generators — hence
every kappa by linearity — and the corrected I5, I7 screen exactly at
alpha = ±1, and the bare branch is exact at alpha = 0 (the same
generic-stratum caveat the paperclip screening verification
carries).  (ii) A bare-fermion branch exists away from
the super-KdV ray: oint psi annihilates the ray
:TT: + (1/2) :G dG: (branch {alpha = 0, kappa = 1/2}); "the screening
locus" above refers to nonzero-exponential screenings of the
kappa = -1/4 tower.

## 4. The naive-pencil family (retained as data, with a hypothesis)

The first version's pencil span{TT_concat, G dG_concat} is a
DIFFERENT 2-dimensional slice of the five-dimensional weight-4 charge
space than the composite slice.  Its distinguished ray
(kappa = -1/4 in naive components; density (dX)^4 - 3 (dX)^2 psi dpsi)
generates a finite exactly-commuting set — its bracket verifications
through [I7, I9] are exact theorems, and its centralizer shows the
same multiplicity-one odd pattern through sigma = 11
(`_naive_ext.json`); no infinite-family theorem is claimed — but it
does NOT commute with the point-split I3
(verified directly), so it is not the standard super-KdV tower.
Within the naive slice, nonzero-exponential screenability forces
kappa = -1/4 and the screening locus is the quartic

    2 alpha^4 - 9 alpha^2 + 2 = 0   (alpha^2 pairs multiplying to 1),

i.e. a+ a- = -1 with a+ + a- = ±sqrt(5/2) != 0
(`super_screening_fit_naive.json`).

**Hypothesis CONFIRMED at the I3-ray level (Phase 3,
`results/super_lambda.md`):** the naive I3 ray IS the
lam = +-sqrt(10)/4 member of the standard family —
:T_lam T_lam: - 1/4 :G_lam dG_lam: = (1/4) I3_naive mod Im d
(solution set exact, back-substituted), and at that lam the family's
screening locus (alpha^2 - 2 lam alpha - 1)(alpha^2 + 2 lam alpha - 1)
equals the quartic above exactly.  The "two towers" of the audit
episode are one lambda-family seen through two pencil slices
(concatenation slice: lam^2 = 5/8; point-split slice: lam = 0); in
standard conventions the naive tower sits at c = 3/2 - 12 lam^2 = -6.
Since the ray identity is a charge identity, the two centralizers
coincide, and the multiplicity-one scans pin the stored naive I5, I7,
I9 as the unique in-range centralizer rays; see super_lambda.md
section 3 for the precise scope (no coefficient-level comparison with
published deformed densities, no all-spin claim).

## 5. Caveats and next steps

- The LOCAL polynomial OPE computations are sector-independent; NS/R
  moding matters for oint G as an operator, fermion zero modes,
  half-integer-spin charges, and boundary states, which the
  literature treats per sector.  Nothing here depends on a sector
  choice, but no sector-level statement is made.
- Requires Python >= 3.9 for math.lcm (a gcd-based fallback is now in
  sparse_linalg.py); heavy runs use python3-env.
- Phase 3 DONE (`results/super_lambda.md`): the background-charge
  family, universal kappa = -1/4, locus F- * F+, and the section-4
  identification.  Next (Phase 4): the two-superfield Z2 x Z2
  super-paperclip scan, using pencil-fit screenability with
  POINT-SPLIT composites and the DUAL-PAIR criterion (a single
  screening family is kappa-blind; see super_lambda.md section 2).
