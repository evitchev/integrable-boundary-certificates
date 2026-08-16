"""High-spin even vacancy => mutual commutativity via the bootstrap lemma.

All promoted charges commute with the low charge identically in t, so by
Jacobi every pairwise commutator lies in ker ad_low at an even charge
spin.  Empty kernel at ONE curve point => generically empty
(semicontinuity) => the commutator, a fixed rational section, vanishes
generically, hence IDENTICALLY.  This proves over Q(t):
  cyl:  [I7,I11] via sigma=18, [I9,I11] via sigma=20 (low = I3),
        plus 12, 14, 16 for the record;
  a2:   [I7,I11] via 18, [I7,I13] via 20, [I11,I13] via 24 (low = I5).
Residues involve only the small low-charge density -- no Y-matching or
composition explosion.  Workers hard-capped at 18 GB address space.
"""
import json, sys, resource
from fractions import Fraction as F
from multiprocessing import Pool
import mixed_engine, collapsed
mixed_engine.residue_mixed_mono = collapsed.residue_mixed_mono_smart
from mixed_engine import gen_mixed_basis, genuine_kernel_mixed, cyl_P4
from invariant_engine import (gen_inv_basis, genuine_kernel_inv, a2_P6,
                              curve_point)

T0 = F(2)


def _init():
    lim = 18 * 1024**3
    resource.setrlimit(resource.RLIMIT_AS, (lim, lim))


def task(args):
    fam, sol, sigma = args
    N, s = curve_point(T0)
    try:
        if fam == "cyl":
            g = genuine_kernel_mixed(cyl_P4(sol, N, s),
                                     gen_mixed_basis(sigma + 1), N)
        else:
            g = genuine_kernel_inv(a2_P6(N, s),
                                   gen_inv_basis(sigma + 1), N)
        dim = len(g)
        print(f"{fam}{sol or ''} sigma={sigma}: ker dim {dim}"
              f" {'(VACANT)' if dim == 0 else '(NOT VACANT!)'}", flush=True)
        return (f"{fam}{sol or ''}_{sigma}", dim)
    except Exception as e:
        print(f"{fam}{sol or ''} sigma={sigma}: ERROR {e!r}", flush=True)
        return (f"{fam}{sol or ''}_{sigma}", -1)


if __name__ == "__main__":
    tasks = [("cyl", sol, sig) for sol in (1, 2, 3)
             for sig in (12, 14, 16, 18, 20)] + \
            [("a2", None, sig) for sig in (18, 20, 24)]
    with Pool(min(16, len(tasks)), initializer=_init) as pool:
        results = pool.map(task, tasks)
    ok = all(d == 0 for _, d in results)
    json.dump(dict(results),
              open("../results/vacancy_high.json" if ok
                   else "../results/vacancy_high.partial.json", "w"),
              indent=1)
    print("ALL VACANT: pairwise commutativity follows identically over "
          "Q(t) via Jacobi" if ok else "NON-VACANT OR ERRORS", flush=True)
    sys.exit(0 if ok else 1)
