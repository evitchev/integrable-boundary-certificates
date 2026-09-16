# Secondary Opus 5 session: the degree-(k-3) symmetric remainder S_k, and the sealed-prediction record (2026-09-14)

Pinned here because the producing workspace (~/opus5-work) is untracked lab space outside the repository, while the content below is the
durable scientific residue of that session.  Its engines are listed in the handover kept alongside this file; the only implementation of the
canonical end split is that session's end_split.py.  Nothing here was committed by that session; the lead pins it.

## 1. S_k -- the symmetric remainder at degree k-3
S_k = sym(certified - M) - sym(rule E), with M the zero-gauge two-sector model and rule E = d_X d_Y[alpha_k M^(mu) + gamma_k M^(mu+1)],
the parameter-free extension of the degree-(k-2) law.  Symmetric under X <-> Y; the listed entries determine the row.

  spin 5  S[1]       = -(t-3)^2 (8t^2 - 45t + 97) / (378 (t-1)^2 (t+1)^2)
  spin 7  S[X^1]     = -2(t-3)^2 (120t^4 + 41t^3 - 507t^2 + 2167t - 3021) / (135 (t-1)^2 (t+1)^2 (t+9)(3t-1))
  spin 9  S[X^2]     = -(t-3)^2 (275t^4 - t^3 - 731t^2 + 2761t - 3984) / (63 (t-1)^2 (t+1)^3 (t+13))
  spin 9  S[XY]      = -2(t-3)^2 (35t^4 - 325t^3 + 805t^2 + 3373t - 6432) / (21 (t-1)^2 (t+1)^3 (t+13))
  spin 11 S[X^3]     = -22(t-3)^2 (489t^4 - 301t^3 - 911t^2 + 3997t - 5434) / (189 (t-1)^2 (t+1)^2 (t+17)(3t+7))
  spin 11 S[X^2 Y]   = 22(t-3)^2 (855t^5 + 13772t^4 - 55092t^3 - 44942t^2 + 293325t - 228942) / (63 (t-1)^2 (t+1)^2 (t+17)(3t+7)(5t-3))
  spin 13 S[X^4]     = -13(t-3)^2 (762t^4 - 979t^3 - 927t^2 + 5995t - 7491) / (54 (t-1)^2 (t+1)^2 (t+21)(3t+11))
  spin 13 S[X^3 Y]   = 26(t-3)^2 (4290t^5 + 30347t^4 - 105026t^3 - 72860t^2 + 424880t - 306111) / (27 (t-1)^2 (t+1)^2 (t+21)(3t+11)(5t+1))
  spin 13 S[X^2 Y^2] = 13(t-3)^2 (4554t^5 + 6519t^4 - 88918t^3 + 44904t^2 + 321628t - 349455) / (3 (t-1)^2 (t+1)^2 (t+21)(3t+11)(5t+1))

FEATURES: a second-order zero at t = 3 at every spin; the SQUARED curve factor (t-1)^2(t+1)^2 in every denominator (the spin-9 entries carry
(t+1)^3); the spin-dependent poles t + 4k - 7 and, in mixed entries, the twist's Pochhammer factors (2r+1)t + 4k - 7 - 10r; NOT in the twisted
span {P^(mu+1), P^(mu+2)} at spin 13, the only spin where that is testable; and at most three independent symmetric values per spin through
spin 13, which is why the shape cannot be identified from the data we have.  The squared curve factor is the one structural hint: it is what a
PRODUCT of two first-order objects would give.

## 2. Sealed predictions of that session (sha256, what it was sealed before, outcome)
  PREDICTION_pole_residue.md      04f7bd70...  any WKB data at spins 7, 9, 17, 19, 21   ratio form confirmed
  PREDICTION_pole_residue_H2.md   3a6ab565...  any data at spins 23, 25                 absolute law a_k, b_k confirmed on 23 entries
  derive/PREDICTION_D1_G.md       0dec6201...  any spin >= 15 along the c_0 direction    beta_k = 2(2k-1)/(3(k+1)) confirmed
  derive/PREDICTION_deficit_law.md 0757b1f2... the beta hold-out and any spin-15 comparison   the degree-(k-2) law confirmed
  derive/PREDICTION_deficit_k3.md 288da568...  any degree-(k-3) comparison               rule E: antisymmetric part confirmed, symmetric fails

## 3. Its own statement of standing, quoted rather than paraphrased
The t = 3 residue law: preregistered twice and confirmed on withheld spins to 25; every step from the pinned recursion carried out symbolically
at general k with each step checked numerically against an independent route; NOT a written proof -- sympy did the series, limits and
cancellations, while the case analysis (regions in r, the boundary w-degree s+1) was cross-checked numerically for k to 24 rather than argued
on paper.  Its own words: "derived-by-verified-computation; a referee would reasonably want the regional sums and the input-(b) limit written
out by hand."  The same standing applies to the resonance and end-split statements.  The degree-(k-2) deficit law and rule E are EMPIRICAL
closed forms confirmed by sealed hold-outs and by Codex's independent layer, not derivations.  The derivative-form identity is elementary
algebra.  A known limit of its main engine: the removable moment factor at j = s is hard-coded to 2, exact on Solution 3's exponents only.

## 4. Its ranked next steps
  1. Write the residue derivation out by hand (the regional sums and the input-(b) limit).
  2. Identify S_k, once the degree-(k-3) layer exists at spins 15 and 17, sealing hypotheses first -- a third twisted kernel top_(mu+2), or a
     PRODUCT of first-order objects, which is what the squared curve factor suggests.
  3. Only then the rank-2 Bethe construction (T3' step 1).
  4. Derive the amplitude laws from the derivative form, where gamma_k simplifies to
     mu(mu+k)/(mu+1) * A * [6(2k-1)A^2 + (14k-5)A + 2(2k-1)] / (180 (2A+1)(3A+2))   [verified by the lead: rf(mu+1,k)/rf(mu,k) = (mu+k)/mu].
  5. The right-end build for Lambda -> 1.
