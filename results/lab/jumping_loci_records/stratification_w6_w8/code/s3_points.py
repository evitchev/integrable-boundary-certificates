"""stratification_w6_w8 spot check at the classical points of spins 5 and 7 named in the paper (t = 5/3; t = 11/5, t = 2) and their
deck images: nullity of my M(t0) and len(genuine_kernel_mixed) at t0, for every sheet and W in {6, 8}."""
import os
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import family as fam  # noqa: E402
import linalg as la  # noqa: E402

me = fam.me
PTS = [F(5, 3), F(-5, 3), F(11, 5), F(-11, 5), F(2), F(-2)]
print("sheet W   t0      nullity(M)  len(genuine)")
allok = True
for sheet in (1, 2, 3):
    for W in (6, 8):
        d = fam.load(os.path.join(fam.BASE, "data", f"M_sheet{sheet}_W{W}.json"))
        n = len(d["cols"])
        for t0 in PTS:
            N0 = (t0 * t0 - 25) / (t0 * t0 - 1)
            s0 = F(-24) * t0 / (t0 * t0 - 1)
            nul = n - la.rank(fam.ev_matrix(d["Mpoly"], t0))
            g = len(me.genuine_kernel_mixed(me.cyl_P4(sheet, N0, s0), me.gen_mixed_basis(W), N0))
            allok &= nul == g
            print(f"{sheet:5d} {W:2d}  {str(t0):6s}  {nul:10d}  {g:12d}", flush=True)
print("ALL OK" if allok else "SOME DIFFER")
