"""Final vacancy assault: sparse kernels + collapsed enumerators.
Covers the eight remaining high-spin vacancy tasks; each VACANT verdict
closes its mutual-commutativity pair identically over Q(t) via Jacobi
+ semicontinuity (see vacancy_high.py)."""
import json, sys, resource
from fractions import Fraction as F
from multiprocessing import Pool
import invariant_engine, mixed_engine, collapsed, collapsed_pairs
invariant_engine.residue_inv_mono = collapsed_pairs.residue_inv_mono_smart
mixed_engine.residue_mixed_mono = collapsed.residue_mixed_mono_smart
from invariant_engine import (gen_inv_basis, a2_P6, curve_point,
                              d_inv_mono, spin_inv, residue_inv_cur)
from mixed_engine import (gen_mixed_basis, cyl_P4, d_mixed_mono,
                          spin_mixed, residue_mixed_cur)
from sparse_linalg import genuine_kernel_sparse

T0 = F(2)


def _init():
    lim = 20 * 1024**3
    resource.setrlimit(resource.RLIMIT_AS, (lim, lim))


def task(args):
    fam, sol, sigma = args
    N, s = curve_point(T0)
    label = f"{fam}{sol or ''} s{sigma}"
    try:
        if fam == "cyl":
            g = genuine_kernel_sparse(cyl_P4(sol, N, s),
                                      gen_mixed_basis(sigma + 1), N,
                                      residue_mixed_cur, gen_mixed_basis,
                                      d_mixed_mono, spin_mixed,
                                      progress=label)
        else:
            g = genuine_kernel_sparse(a2_P6(N, s),
                                      gen_inv_basis(sigma + 1), N,
                                      residue_inv_cur, gen_inv_basis,
                                      d_inv_mono, spin_inv,
                                      progress=label)
        print(f"{fam}{sol or ''} sigma={sigma}: ker dim {len(g)}"
              f" {'(VACANT)' if not g else '(NOT VACANT!)'}", flush=True)
        return (f"{fam}{sol or ''}_{sigma}", len(g))
    except Exception as e:
        print(f"{fam}{sol or ''} sigma={sigma}: ERROR {e!r}", flush=True)
        return (f"{fam}{sol or ''}_{sigma}", -1)


if __name__ == "__main__":
    tasks = [("a2", None, 24),
             ("cyl", 1, 18), ("cyl", 2, 18), ("cyl", 3, 18),
             ("cyl", 1, 20), ("cyl", 2, 20), ("cyl", 3, 20)]
    with Pool(4, initializer=_init) as pool:
        results = pool.map(task, tasks)
    ok = all(d == 0 for _, d in results)
    json.dump(dict(results),
              open("../results/vacancy_final.json" if ok
                   else "../results/vacancy_final.partial.json", "w"),
              indent=1)
    print("ALL VACANT" if ok else "NON-VACANT OR ERRORS", flush=True)
    sys.exit(0 if ok else 1)
