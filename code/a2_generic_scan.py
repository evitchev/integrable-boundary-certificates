"""Generic-N spectral scan of the A_2^(2) family, and the branch-point test.

Question: is the N=1 tower enhancement (spin-3 and spin-9 charges
appearing) caused by (a) the branch point of the double cover
s^2 = (N-25)(N-1), or (b) the finite-N relations of the invariant
algebra at N=1?

Design: kernel dims of ad_{I_5^{A2}} for charge spins sigma = 2..13 at
  - generic rational points of the cover (t = 2, 3, 7, 9/2, -2);
  - the branch point (N, s) = (25, 0)  [formal algebra IS faithful at
    N=25 for these letter counts -> clean branch test];
  - the branch-point coefficients (N, s) = (1, 0) in the FORMAL algebra
    [no finite-N relations -> isolates the coefficient degeneration];
compared against the physical N=1 result (component engine, earlier):
enhancement at sigma = 3, 9.

Rank semicontinuity on the irreducible rational curve: full rank
(empty genuine kernel) at ONE sample point proves generic vacancy.

Uses multiprocessing (one worker per sigma; per-process residue caches
are reused across all parameter points).
"""

import sys
from fractions import Fraction as F
from multiprocessing import Pool

from invariant_engine import (gen_inv_basis, genuine_kernel_inv,
                              InvSector, a2_P6, curve_point)

POINTS = [
    ("t=2   (N=-7)", curve_point(2)),
    ("t=3   (N=-2)", curve_point(3)),
    ("t=7   (N=1/2)", curve_point(7)),
    ("t=9/2 (N=-19/77)", curve_point(F(9, 2))),
    ("t=-2  (sheet 2)", curve_point(-2)),
    ("BRANCH N=25", (F(25), F(0))),
    ("BRANCH-coeff N=1 (formal)", (F(1), F(0))),
]

SIGMAS = list(range(2, 14))


def scan_sigma(sigma):
    cand = gen_inv_basis(sigma + 1)
    dims = []
    for label, (Nval, sval) in POINTS:
        P6 = a2_P6(Nval, sval)
        gen = genuine_kernel_inv(P6, cand, Nval)
        dims.append(len(gen))
    qdim = InvSector(sigma + 1).q_dim()
    return sigma, qdim, dims


if __name__ == "__main__":
    with Pool(min(16, len(SIGMAS))) as pool:
        results = pool.map(scan_sigma, SIGMAS)

    w = max(len(lbl) for lbl, _ in POINTS)
    print("ker ad_{I_5^A2} genuine dims per charge spin sigma")
    print("point" + " " * (w - 5) + " | " +
          " ".join(f"{s:2d}" for s in SIGMAS))
    print("-" * (w + 3 + 3 * len(SIGMAS)))
    for i, (lbl, _) in enumerate(POINTS):
        row = " ".join(f"{results[j][2][i]:2d}" for j in range(len(SIGMAS)))
        print(f"{lbl:<{w}} | {row}")
    print(f"{'dim Q^sigma':<{w}} | " +
          " ".join(f"{results[j][1]:2d}" for j in range(len(SIGMAS))))
