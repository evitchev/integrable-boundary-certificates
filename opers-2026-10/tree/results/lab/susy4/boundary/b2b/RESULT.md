# B2b RESULT (Fable secondary seat, 2026-10-13): do Solution 1's charges commute with the N = 2 KdV towers?

PREREG results/lab/susy4/boundary/B2b_PREREG.md, commit 3898542, sha256 5feb96f00cdc91162c91f05dddee8c56c08c2ff5a9ae0c08bdb4f118553fc002
(verified on the tree and the blob; item 353).  Run as written: residues [oint L_s, oint I_r] mod d, (s, r) in {2, 3, 5}^2 (I_2 of Z only),
Z and K separately, t = 2, 7/2, 13/5, control kappa -> kappa + 1/7.  Workspace <workspace>/b2b/ (START.time 2026-10-13T14:59:53Z, END.time 2026-10-13T15:02:17Z);
inputs copied unchanged from B2 (inputs/INPUTS.SHA256); B2's validated modules imported read-only; one new file b2b.py; logs/ hold
the full stdout of the six runs; every run under ulimit -v 8000000 and timeout.  Each run takes 12 s.

## VERDICT (PREREG sec. 4): NO, both towers

**[L_s, I_r] != 0 mod d for EVERY registered pair (s, r), at EVERY fibre, for Z and for K.**  In the PREREG's words: *Solution 1's charges
are not conserved by the N = 2 KdV flows (a = 4 or a = -2)*; together with B2 this closes the coset route at spins 2-5.
P-b1 (null expectation): confirmed, (3,3) != 0 at all three fibres, Z and K.  P-b2 ([L_2, Z_2] = 0): **FAILS** at all three fibres: the
vacuum-block proportionality of B2 (Z_2's H~-vacuum block = c_2 L_2) does not extend to full commutation; (2,2) was not forced and
does not hold.  No PARTIAL: no pair commutes anywhere.

## Table (residue mod d: 0 = commute; the number is the monomial count of the raw residue)

| pair (s, r) | Z, t = 2 | Z, 7/2 | Z, 13/5 | K, t = 2 | K, 7/2 | K, 13/5 |
|---|---|---|---|---|---|---|
| (2,2) | NO (24) | NO (24) | NO (24) | - | - | - |
| (2,3) | NO (62) | NO (62) | NO (62) | NO (62) | NO (62) | NO (62) |
| (2,5) | NO (277) | NO (280) | NO (279) | NO (277) | NO (280) | NO (279) |
| (3,2) | NO (47) | NO (47) | NO (47) | - | - | - |
| (3,3) | **NO** (116) | **NO** (116) | **NO** (116) | **NO** (116) | **NO** (116) | **NO** (116) |
| (3,5) | NO (505) | NO (507) | NO (507) | NO (505) | NO (505) | NO (505) |
| (5,2) | NO (159) | NO (159) | NO (157) | - | - | - |
| (5,3) | NO (377) | NO (377) | NO (374) | NO (377) | NO (377) | NO (374) |
| (5,5) | **NO** (1541) | **NO** (1542) | **NO** (1538) | **NO** (1542) | **NO** (1542) | **NO** (1538) |

## Objects and checks

L_2 = Solution 1's X-odd weight-3 (u,v) commutant generator, coeff[(du)^2 dv] = 1 (3 monomials); L_3 = the certified I_3 in (u,v)
(sol1, kit normalisation of cyl_P4; 9 monomials); L_5 = the (u,v) commutant generator at w2 = 12 proportional to the certified I_5
eigenvalue (Q^6-normalised; 16 monomials; at t = 2 the v-even one of two).  Lift: the same monomials in (u, v, H~), metric
(s, 4-s, (s-4)/s).  I_r: B2's tower charges (canonically normalised), psi^+- basis; the charge -+2 parts (I_2 of Z only, single
monomials) ASSERTED total derivatives and dropped; the neutral part bosonised (coset.Bosoniser) and rotated (coset.rotate):
I_2: 9, I_3: 20, I_5: 76 (Z) / 75 (K) monomials.  Total-derivative test: exact sparse elimination on the bosonic monomials of the
three-boson algebra at the residue's doubled weight (up to 22).
Consistency (every fibre, both towers): [L_2, L_3] = [L_2, L_5] = [L_3, L_5] = 0 mod d (the lifted Solution 1 charges commute among
themselves, as FL's AKNS charges must); [I_2, I_3] = [I_2, I_5] = [I_3, I_5] = 0 mod d in the FERMIONIC frame and, independently
computed, after bosonisation (the bosonisation preserves the commutators: the machinery used for the verdict is validated on the
towers' own commutation at the same weights).
Negative control (lift with kappa' = kappa + 1/7, coefficients in Q(kappa) split as A + kappa B, both parts tested): fails for every
pair at every fibre, as required for a YES claim; since the real test also fails everywhere, the control does not discriminate
here (reported, not used).

## Disclosures
- Nothing was fitted or re-normalised; the normalisations only rescale the residues.
- The PREREG's "lifted as H~-free densities" was applied literally; no H~-dependent lift was tried (it would be a new object).
- t = 2: the (u,v) commutant at w2 = 12 is two-dimensional (B2); L_5 is its element proportional to the certified I_5.
- No interpolation in t was attempted.
