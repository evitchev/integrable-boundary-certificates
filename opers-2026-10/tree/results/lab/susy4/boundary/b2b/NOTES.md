# B2b running notes (Fable secondary seat) -- 2026-10-13
PREREG results/lab/susy4/boundary/B2b_PREREG.md, commit 3898542, sha256 5feb96f00cdc91162c91f05dddee8c56c08c2ff5a9ae0c08bdb4f118553fc002
(VERIFIED on the tree and the blob).  Item 353.  Inputs: B2's tower charges (towers_*.json) and certified Solution 1 data (sol1_t*.json)
copied unchanged from <workspace>/b2 (hashes in inputs/INPUTS.SHA256); B2's modules gm.py, coset.py, n2frame.py imported read-only
from <workspace>/b2 (sys.path).  New file: b2b.py.  Every run: ulimit -v 8000000, timeout, full stdout in logs/.
Objects: L_2 = Solution 1's X-odd weight-3 (u,v) generator with coeff[(du)^2 dv] = 1; L_3 = certified I_3 in (u,v) (kit normalisation
of cyl_P4); L_5 = the (u,v) commutant generator at w2 = 12 proportional to the certified I_5 eigenvalue (Q^6-normalised form) -- at
t = 2 the v-even one of the two.  Lift: the same monomials in the three-boson algebra (u, v, H~), metric (s, 4-s, (s-4)/s).
I_r (Z: 2, 3, 5; K: 3, 5): B2's charges, psi^+- basis, charge -+2 parts asserted total derivatives and dropped, neutral part
bosonised (coset.Bosoniser) and rotated (coset.rotate).  Test: Res(L_s(z) I_r(w)) mod d in the three-boson algebra (bosonic
monomials only; d-span by exact sparse elimination).  Control: lift with kappa -> kappa + 1/7 (coefficients in Q(kappa), split
as A + kappa B, both parts must be d-exact for commutation).  Consistency: [L_s, L_r] and [I_r, I_r'] (the latter in the fermionic
frame AND after bosonisation).
