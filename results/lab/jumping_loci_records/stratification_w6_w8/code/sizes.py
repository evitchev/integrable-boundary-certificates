"""stratification_w6_w8: sizes of the spin-s monomial bases and their quotients mod total derivatives (engine objects)."""
import os
import sys
import time

# the repository's code/ directory (this file sits in <repository>/results/lab/jumping_loci_records/<name>/code/)
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "..", "code")))   # appended: stdlib and site-packages keep precedence
from mixed_engine import MixedSector, gen_mixed_basis  # noqa: E402

for s in range(1, 12):
    t0 = time.time()
    sec = MixedSector(s)
    print(f"spin {s:2d}: monomials {len(sec.basis):5d}  d-image rank {len(sec._img_pivots):5d}  "
          f"quotient dim {sec.q_dim():5d}  ({time.time() - t0:.1f} s)", flush=True)
