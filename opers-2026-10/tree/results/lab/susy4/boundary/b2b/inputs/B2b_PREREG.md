# PREREGISTRATION B2b (filed by the lead): do Solution 1's charges commute with the N = 2 KdV towers?

Filed 2026-10-13 by the lead, before any B2b computation.  Follows B2 (PREREG c24a7af; PRIMARY 600bfd83; hold-out 7b88494f; register
items 352 and its addendum).  The sha256 of this file and its commit are recorded in the register in a later, separate command.

## 0. Disclosure

No B2b computation has been run.  What is known beforehand comes from B2.  The H-vacuum projection of Z's and K's charges fails to lie
in Solution 1's commutant at weight >= 8.  Z's charges do not commute with e^{-u} (x) 1.  W0 carries a dH~ prefactor.  The seat's
stated prior for a "yes" at spin 3 is about 20%.

## 1. Objects (as in B2, frame eta = (p^2, -p^2), bosonised three-boson algebra (u, v, H~))

- Z: Zeitlin's sl^(1)(2|1) tower (the a = 4 N = 2 KdV; B2 V2).  K: KL's sl(2|2)^(1) tower (a = -2).  Use the charges I_2 (Z only),
  I_3, I_5, as computed and validated in B2, bosonised by B2's validated map.
- L_s: Solution 1's certified charges I_3, I_5, and the X-odd I_2.  They are written as (u,v) densities via B2's isometry
  u = -p phi, v = -2 kappa X, and lifted to the three-boson algebra as H~-free densities.
- Test: the residue Res(L_s(z) I_r(w)) mod d, i.e. [oint L_s, oint I_r] = 0, at the doubled weights of (s, r) in
  {2, 3, 5} x {2, 3, 5}, for Z and K separately.

## 2. Fibres

t = 2, 7/2, 13/5.  Use the Q^6-normalised I_5 wherever the kit normalisation is singular.  No hold-out is needed for a yes/no
commutation test.  If any interpolation in t is attempted afterwards, it requires a new registration.

## 3. Predictions (registered)

- P-b1 (null expectation, following the seat's prior): [L_3, Z_3] != 0 at all three fibres.  Likewise for K.
- P-b2: [L_2, Z_2] = 0.  B2's projection result forces Z_2's vacuum block to be c_2 L_2; full commutation is not forced, so this
  can fail.
- The data are reported for every pair (s, r); the verdict rests on (3, 3) and (5, 5).

## 4. Verdicts

- YES (Z) if [L_s, Z_r] = 0 mod d for all registered pairs at all three fibres.  Statement: "Solution 1 (x) Heis_H~ commutes with
  the a = 4 N = 2 KdV tower through spin 5 at the fibres tested."  Before that is claimed: a negative control (lift with
  kappa -> kappa + 1/7) must fail, and the same test must be run on K.
- NO: a non-zero residue mod d for (3,3) or (5,5) at any fibre.  Certified negative: "Solution 1's charges are not conserved by
  the N = 2 KdV flows (a = 4 or a = -2)."  Together with B2, this closes the coset route for B2 at spins 3-5.
- PARTIAL (discovery-then-check): commutation at some pairs or fibres only.  Reported with the residues, as INCONCLUSIVE, never
  relabelled.

## 5. Rules

Exact arithmetic; ulimit -v 8000000 and timeout on every run; full stdout logs in results/lab/susy4/boundary/b2b/logs/.
Use B2's validated code read-only; new files only.  The lead re-runs every result before it enters the register.
