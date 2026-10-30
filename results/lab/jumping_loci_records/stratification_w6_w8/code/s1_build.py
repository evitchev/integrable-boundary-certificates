"""stratification_w6_w8 S1: build and save M_poly(t) for sheets 1-3 and W in {6, 8}."""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import family as fam  # noqa: E402

os.makedirs(os.path.join(fam.BASE, "data"), exist_ok=True)
for sheet in (1, 2, 3):
    P4 = fam.p4_family(sheet)
    print(f"sheet {sheet}: I_3 density (engine cyl_P4 evaluated in Q(t)):")
    for m in sorted(P4):
        print(f"    {fam.mono_str(m):16s} {P4[m]}")
    for W in (6, 8):
        t0 = time.time()
        M, rowlab, collab = fam.build_M(P4, W)
        k, c0, Mpoly = fam.clear(M)
        nz = sum(1 for row in Mpoly for e in row if any(e))
        dmax = max((len(e) - 1 for row in Mpoly for e in row if e), default=-1)
        path = os.path.join(fam.BASE, "data", f"M_sheet{sheet}_W{W}.json")
        fam.save(path, sheet, W, k, c0, Mpoly, rowlab, collab)
        print(f"  W={W}: M is {len(rowlab)} x {len(collab)} (target classes x domain classes); "
              f"scalar (t^2-1)^{k} * {c0}; nonzero entries {nz}; max entry degree {dmax}; "
              f"saved {os.path.relpath(path, fam.BASE)} ({time.time() - t0:.1f} s)", flush=True)
        print(f"    domain classes: {' '.join(fam.mono_str(m) for m in collab)}")
