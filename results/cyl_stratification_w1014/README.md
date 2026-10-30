# Claims verified by `cyl_stratification_w1014`

The 54 files in `claims/` are the exact data of the classification over Q-bar (away from t = +-1) of the
rank-drop locus of the family matrix M(t) of ad I_3, at weights 10, 12 and 14 (spins 9, 11, 13) on
sheets 1-3. They were computed on 2026-10-29 by an independent recomputation with open tools (Python,
sympy, numpy, FLINT 3.5.0) and reproduced the same day from byte copies of its scripts. The
certificate `code/cyl_stratification_w1014.py` does not trust them: it rebuilds the nine matrices
from the engine and re-verifies every claim (see its header). `claims/SHA256SUMS` lists every file.
The certificate pins the hash of that list.

Per sheet and weight (`tag` = `sheet{S}_W{W}`):

- `S2_{tag}.json`: the size m x n of M_poly; the scalar cleared from M(t); the generic rank r with
  the prime and the point t* of the lower bound; pivot rows and columns; the free columns.
  - The row and column sets of the two maximal minors (R1 with the pivot columns; R2 with C2), with
    their content, the powers of t -+ 1 removed, the degree bounds, the Hadamard bounds B^2 and the
    62-bit primes, whose product exceeds 2 B 2^64.
  - G = gcd(P1o, P2o), its factors over Q and the certificate prime for the coprime cofactors.
- `S2_P1_{tag}.json`, `S2_P2_{tag}.json`: the exact minors P1, P2 and P1o, P2o (coefficient lists,
  low to high), with the powers of t - 1 and t + 1 divided out.
- `S2_nullvectors_{tag}.json`: n - r exact null vectors of M_poly over Z[t], one per free column.
- `S3_{tag}.json`: the decision for every irreducible factor of G.
  - No drop: a prime p, a root alpha of the factor modulo p, and the rank at alpha.
  - Drop: the rational root, the nullity, the ranks modulo two primes, and the file of exact null
    vectors.
  - Also: the jump set with multiplicities and the deck images.
- `S3_nullvectors_{tag}_t{t0}.json`: the exact null vectors at each rational drop.

The file-name fields inside the claims (`nullvectors_file`, `P1_file`, ...) are relative to the
original working directory; the verifier opens them by base name in this directory. The matrices
are not shipped. The certificate rebuilds them, and their sha256 values are pinned in
`code/strat1014_verify.py` (`EXPECTED`).

Result (t != +-1):

| W  | n x m      | r   | sheet 1 | sheet 2  | sheet 3  |
|----|------------|-----|---------|----------|----------|
| 10 | 68 x 185   | 67  | none    | 3        | 7/5, 3   |
| 12 | 160 x 461  | 159 | none    | 9/7, 3   | 3        |
| 14 | 365 x 1082 | 364 | none    | 3        | 9/7, 3   |

The nullity is 2 at every jump, and beta^2 = (t-1)/(t+1) = 1/2, 1/6 and 1/8 at t = 3, 7/5 and 9/7.
No deck image -t0 of a jump is a jump. These are the expectations preregistered in
`results/cyl_discriminant_preregistration.md` and `results/cyl_discriminant_w14_preregistration.md`.
