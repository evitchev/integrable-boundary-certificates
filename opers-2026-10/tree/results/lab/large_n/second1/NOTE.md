# SECOND1: is Sol 3's generic-t operator a Mellin image of a tensor product of two second-order equations? Opus 5 (2), 2026-10-01.
Seal 3a026fa1 (sent before any comparison; the pre-seal structural computation s0_shifts is disclosed in it). Exact sympy; every run capped (ulimit -v 8-20 GB); full logs.

## Verdict: NO at every t except the paperclip. PROVED for all k within the sealed ansatz.
Ansatz:
- Two second-order theta-equations with a COMMON right side, A_pm: theta^2 chi = u_pm theta chi + (lam_pm + x^a - E x^b) chi. Exponents arbitrary, so the rotation is not imposed.
- Gauge (absorbed in the exponents) and x-scaling (constants cancel) are allowed.
- An ARBITRARY meromorphic Mellin multiplier is allowed, plus an overall factor.
- Both shift assignments are tested: (a, b) = (c, n) and (n, c).
1. Structure (s0_shifts.log, pre-seal): the tensor product has exactly two x-shifts. Its coefficients are
   - x^0: prod (theta - r_i - q_j);
   - x^a: -(2 theta + a - U)(2 theta + 2a - U);
   - x^b: E (2 theta + b - U)(2 theta + 2b - U).
2. Matching Sol 3's Gamma recurrence (item 298) requires one multiplier m solving m(nu)/m(nu - a) = f_a and m(nu)/m(nu - b) = f_b. This forces the COCYCLE identity F(nu) = 1 (s1_cocycle.py), a rational identity, E-free.
3. THEOREM (s2_leading.log): the leading nu-coefficient of the cocycle numerator is -+256 k^7 (k+2)(k+4) for the two assignments.
   - It is INDEPENDENT of all exponents.
   - So closure requires k in {-4, -2, 0}. k = 0 is excluded (sigma = n/k).
   - k = -4 (n = 0, the degenerate N = 0 member of item 299) gives NONE when solved exactly, both assignments.
   - k = -2 (paperclip) CLOSES: it recovers exactly the rotated exponents {1/2 +- l0, 1/2 +- l1} (8 relabellings), reproducing RED2 B2 (item 309) independently.
   - **Hence for every t != the paperclip, Sol 3's operator is NOT a Mellin-multiplier image of a tensor product of two such second-order equations.**
4. Exact solves at the sealed fibres: k = 10/3 (t = 9/4) and k = 2 (t = 2), both assignments: NONE. Sealed prediction NONE 80%: PASS.
   - Scan (s3_scan.log), 18 rational k values in [-6, 4], both assignments: closes ONLY at k = -2, consistent with the theorem.
5. Controls at k = -2:
   - the wrong multiplier shift (string length k + 1/10) gives NONE;
   - the wrong rotation (the unrotated exponent set forced) gives NONE.
   Both fire.
## Extension (post-seal, labelled): a Suzuki-type +-alpha doublet (s4_doublet_shifts.log)
- The tensor product of A_pm with right sides x^a +- alpha x^a' - E x^b has a NON-polynomial relation: denominator 8 alpha x^a' + const, an apparent singularity. It carries eight x-monomials.
- A Mellin multiplier only rescales terms. It cannot remove the extra shifts or the apparent singularity, so it cannot reach Sol 3's three-term recurrence.
- Structural NONE. Not a full classification of three-term second-order pairs.
## Not tested
- Pairs with DIFFERENT right sides (different shifts or spectral dependence per factor).
- Tensor products with an apparent-singularity factor.
- Other non-local transforms (Laplace, Euler).
- Sol 2: per the commission it was conditional on Sol 3 closing, so it was not pursued.
## Reading
- The paperclip's B2 identity is ISOLATED within this ansatz. The leading-coefficient obstruction (k+2)(k+4) is the reason: only there do the shift structures {c, n} = {1, 2} of the tensor product and the Gamma form become compatible.
- RED1's anchor 'product of the two quadratics' is not a tensor product of two equations sharing E; it is a different object.
## Disclosures
- s2_symbolic_k.py (the full joint solve in k) was too slow; I killed it by saved PID (s2_pid.txt; partial log s2_symbolic_k_killed.log). It was replaced by the leading-coefficient argument (s2_leading.py), which suffices.
- The 'wrongrot' control forces one fixed unrotated assignment; it is not an exhaustive scan of non-rotations. The theorem covers all exponents anyway.
