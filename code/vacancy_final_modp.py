"""Final vacancy assault, GF(p) edition (replaces the fraction-free
vacancy_final.py run, whose integer-coefficient tail blowup was
measured at weeks-scale on the sigma=18 sectors — see the probe data
of 2026-08-07).

Each task runs genuine_kernel_modp: a tags == d-rank verdict is a
RIGOROUS vacancy certificate over Q (kernels grow, ranks drop mod p,
so equality pins the exact dimensions).  A gap escalates to a second
prime; if the gap persists the task reports the inconclusive upper
bound (NOT a non-vacancy claim) and the batch fails closed.

Each VACANT verdict closes its mutual-commutativity pair identically
over Q(t) via Jacobi + semicontinuity, as before: cyl sigma=18 ->
[I7,I11]; cyl sigma=20 -> [I9,I11]; a2 sigma=24 -> [I11,I13]."""
import json
import os
import resource
import sys
from fractions import Fraction as F
from multiprocessing import Pool

import sympy as sp

import collapsed
import collapsed_pairs
import invariant_engine
import mixed_engine
invariant_engine.residue_inv_mono = collapsed_pairs.residue_inv_mono_smart
mixed_engine.residue_mixed_mono = collapsed.residue_mixed_mono_smart
from invariant_engine import (a2_P6, curve_point, d_inv_mono,
                              gen_inv_basis, residue_inv_cur, spin_inv)
from mixed_engine import (cyl_P4, d_mixed_mono, gen_mixed_basis,
                          residue_mixed_cur, spin_mixed)
from sparse_linalg import MODP_PRIME, genuine_kernel_modp

T0 = F(2)
PRIMES = [MODP_PRIME, int(sp.nextprime(MODP_PRIME))]


HARD_CAP = 16 * 1024**3


def _init():
    # Task-specific soft caps set in task(); the hard cap here bounds
    # them.  16 GB (a2) + 6 x 12 GB (cyl) = 88 GB aggregate, within
    # the machine's ~92 GB available.  A uniform 12 GB was refuted by
    # live data: the a2 sigma=24 worker passed 11.2 GB RSS at 77% of
    # its sweep (audit catch).
    resource.setrlimit(resource.RLIMIT_AS, (HARD_CAP, HARD_CAP))


CKPT_DIR = os.path.expanduser("~/vacancy_logs/ckpts")


def task(args):
    fam, sol, sigma = args
    N, s = curve_point(T0)
    label = f"{fam}{sol or ''} s{sigma}"
    soft = HARD_CAP if fam == "a2" else 12 * 1024**3
    resource.setrlimit(resource.RLIMIT_AS, (soft, HARD_CAP))
    try:
        os.makedirs(CKPT_DIR, exist_ok=True)
        for pi, prime in enumerate(PRIMES):
            ck = os.path.join(
                CKPT_DIR, f"{fam}{sol or ''}_s{sigma}_p{pi}.pkl")
            if fam == "cyl":
                tags, drank, vac = genuine_kernel_modp(
                    cyl_P4(sol, N, s), gen_mixed_basis(sigma + 1), N,
                    residue_mixed_cur, gen_mixed_basis, d_mixed_mono,
                    spin_mixed, prime=prime,
                    progress=f"{label} p{pi}", checkpoint_path=ck)
            else:
                tags, drank, vac = genuine_kernel_modp(
                    a2_P6(N, s), gen_inv_basis(sigma + 1), N,
                    residue_inv_cur, gen_inv_basis, d_inv_mono,
                    spin_inv, prime=prime, progress=f"{label} p{pi}",
                    checkpoint_path=ck)
            if vac:
                print(f"{label}: VACANT (certified over Q, prime "
                      f"#{pi}; tags = d-rank = {tags})", flush=True)
                return (f"{fam}{sol or ''}_{sigma}", 0)
            print(f"{label}: prime #{pi} inconclusive "
                  f"(bound {tags - drank}); escalating", flush=True)
        print(f"{label}: INCONCLUSIVE after {len(PRIMES)} primes "
              f"(genuine_Q <= {tags - drank})", flush=True)
        return (f"{fam}{sol or ''}_{sigma}", tags - drank)
    except Exception as e:
        print(f"{label}: ERROR {e!r}", flush=True)
        return (f"{fam}{sol or ''}_{sigma}", -1)


if __name__ == "__main__":
    # Task selection (argv[1]): the original 7-task batch lost its
    # sigma=20 trio at 90.5% to a session teardown -- reruns select
    # only what is still owed, and production launches are detached
    # (setsid) so no session exit can kill them again.
    FULL = [("a2", None, 24),
            ("cyl", 1, 18), ("cyl", 2, 18), ("cyl", 3, 18),
            ("cyl", 1, 20), ("cyl", 2, 20), ("cyl", 3, 20)]
    SETS = {"full": (FULL, ""),
            "cyl20": ([t for t in FULL if t[2] == 20], "_s20"),
            "smoke": ([("cyl", 1, 12)], "_smoke")}
    which = sys.argv[1] if len(sys.argv) > 1 else "full"
    tasks, suffix = SETS[which]
    # idempotent completion guard (@reboot cron relaunches become
    # no-ops once the batch verdicts are in)
    done_marker = f"../results/vacancy_final_modp{suffix}.json"
    if which != "full" and os.path.exists(done_marker):
        print(f"{done_marker} already exists -- nothing to do")
        sys.exit(0)
    with Pool(len(tasks), initializer=_init) as pool:
        results = pool.map(task, tasks)
    ok = all(d == 0 for _, d in results)
    json.dump(dict(results),
              open(f"../results/vacancy_final_modp{suffix}.json"
                   if ok else
                   f"../results/vacancy_final_modp{suffix}.partial.json",
                   "w"), indent=1)
    print("ALL VACANT (certified)" if ok
          else "INCONCLUSIVE OR ERRORS", flush=True)
    sys.exit(0 if ok else 1)
