"""Rerun of the seven open high-spin vacancies with the Phase-2
collapsed engines routed in.  Same Jacobi/semicontinuity logic as
vacancy_high.py; 20 GB per-worker caps; Pool(4) to bound total memory."""
import json, sys, resource
from fractions import Fraction as F
from multiprocessing import Pool
import invariant_engine, mixed_engine, collapsed, collapsed_pairs
invariant_engine.residue_inv_mono = collapsed_pairs.residue_inv_mono_smart
mixed_engine.residue_mixed_mono = collapsed.residue_mixed_mono_smart
from mixed_engine import gen_mixed_basis, genuine_kernel_mixed, cyl_P4
from invariant_engine import (gen_inv_basis, genuine_kernel_inv, a2_P6,
                              curve_point)

T0 = F(2)


def _init():
    lim = 20 * 1024**3
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
    tasks = [("a2", None, 20), ("a2", None, 24),
             ("cyl", 1, 18), ("cyl", 2, 18), ("cyl", 3, 18),
             ("cyl", 1, 20), ("cyl", 2, 20), ("cyl", 3, 20)]
    with Pool(4, initializer=_init) as pool:
        results = pool.map(task, tasks)
    ok = all(d == 0 for _, d in results)
    json.dump(dict(results),
              open("../results/vacancy_rerun.json" if ok
                   else "../results/vacancy_rerun.partial.json", "w"),
              indent=1)
    print("ALL VACANT" if ok else "NON-VACANT OR ERRORS", flush=True)
    sys.exit(0 if ok else 1)
