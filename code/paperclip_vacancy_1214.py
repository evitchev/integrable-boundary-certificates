"""Paperclip even-spin vacancies at sigma = 12, 14.

Extends the sigma <= 10 spectral table of vacancy.py using the sparse
stack (the dense Sector RREF was the wall for these sectors).  One
generic rational n suffices per sigma: kernel dims are upper
semicontinuous, so vacancy at one point implies generic vacancy; a
second point is run as a belt-and-braces cross-check.

Each VACANT verdict extends the bootstrap lemma's even-spin hypothesis
range for the paperclip tower.  Fail-closed JSON artifact; heartbeats
every 10 minutes per worker."""
import json
import resource
import sys
from fractions import Fraction as F
from multiprocessing import Pool

from families import paperclip_P4
from ope_engine import d_mono, gen_basis, residue_cur, spin
from sparse_linalg import genuine_kernel_sparse


def _init():
    lim = 18 * 1024**3
    resource.setrlimit(resource.RLIMIT_AS, (lim, lim))


def task(args):
    sigma, nval = args
    label = f"pc s{sigma} n={nval}"
    try:
        g = genuine_kernel_sparse(
            paperclip_P4(nval),
            gen_basis(sigma + 1, 2, True),
            None,
            lambda P, Q, N: residue_cur(P, Q),
            lambda s: gen_basis(s, 2, True),
            d_mono, spin,
            progress=label)
        print(f"{label}: ker dim {len(g)} "
              f"{'(VACANT)' if not g else '(NOT VACANT!)'}", flush=True)
        return (f"s{sigma}_n{nval}", len(g))
    except Exception as e:
        print(f"{label}: ERROR {e!r}", flush=True)
        return (f"s{sigma}_n{nval}", -1)


if __name__ == "__main__":
    tasks = [(12, F(7, 3)), (14, F(7, 3)),
             (12, F(13, 5)), (14, F(13, 5))]
    with Pool(2, initializer=_init) as pool:
        results = pool.map(task, tasks)
    ok = all(d == 0 for _, d in results)
    json.dump(dict(results),
              open("../results/paperclip_vacancy_1214.json" if ok else
                   "../results/paperclip_vacancy_1214.partial.json",
                   "w"), indent=1)
    print("ALL VACANT" if ok else "NON-VACANT OR ERRORS", flush=True)
    sys.exit(0 if ok else 1)
