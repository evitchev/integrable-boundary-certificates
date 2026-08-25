"""GF(p) vacancy certificates for A2(2) sigma = 18 and 20.

The last two of the 21 ledger instances still resting on the older
fraction-free / sparse+collapsed computations of 2026-08-02/03 (narrative
row in results/vacancy_final.md, no retained per-task artifact): the
A2(2) pairs [I7,I11] at summed spin 18 and [I7,I13] at 20.  The
manuscript's "Every instance is certified by a single-prime GF(p)
certificate" was not literal for these two (Codex finding, 2026-08-24).
This script closes the gap by recomputing both under the same GF(p)
certifier that produced the sigma=24 A2(2) verdict and the cylindrical
12/14/16, 18, 20 verdicts, retaining a per-task artifact.

Unlike the cylindrical 12/14/16 batch, NO prior exact ranks exist for
these two sectors (the old runs retained no rank record), so the tags
are RECORDED here, not preregistered; acceptance is verdict VACANT with
tags == d-rank (the rigorous equality: kernels only grow and ranks only
drop mod p, so equality at one prime pins the exact dimensions over Q).
A gap escalates to a second prime and then fails closed as an
inconclusive bound, never a non-vacancy claim.
"""
import json
import resource
import sys
from fractions import Fraction as F
from multiprocessing import Pool
from pathlib import Path

import sympy as sp

from certificate import stamp
from invariant_engine import (a2_P6, curve_point, d_inv_mono,
                              gen_inv_basis, residue_inv_cur, spin_inv)
from sparse_linalg import MODP_PRIME, genuine_kernel_modp

T0 = F(2)
PRIMES = [MODP_PRIME, int(sp.nextprime(MODP_PRIME))]
CAP = 12 * 1024**3


def _init():
    resource.setrlimit(resource.RLIMIT_AS, (CAP, CAP))


def task(sigma):
    N, s = curve_point(T0)
    label = f"a2 s{sigma}"
    rec = {"family": "a2", "sigma": sigma}
    try:
        for pi, prime in enumerate(PRIMES):
            tags, drank, vac = genuine_kernel_modp(
                a2_P6(N, s), gen_inv_basis(sigma + 1), N,
                residue_inv_cur, gen_inv_basis, d_inv_mono,
                spin_inv, prime=prime, progress=f"{label} p{pi}")
            rec.update(prime=prime, prime_index=pi, tags=tags,
                       derivative_rank=drank, genuine_gap=tags - drank)
            if vac:
                rec["verdict"] = "VACANT"
                print(f"{label}: VACANT (certified over Q, prime #{pi}; "
                      f"tags = d-rank = {tags})", flush=True)
                return rec
            print(f"{label}: prime #{pi} inconclusive "
                  f"(bound {tags - drank}); escalating", flush=True)
        rec["verdict"] = "INCONCLUSIVE"
        print(f"{label}: INCONCLUSIVE after {len(PRIMES)} primes "
              f"(genuine_Q <= {tags - drank})", flush=True)
        return rec
    except Exception as e:                      # fail closed
        rec["verdict"] = "ERROR"
        rec["error"] = repr(e)
        print(f"{label}: ERROR {e!r}", flush=True)
        return rec


if __name__ == "__main__":
    with Pool(2, initializer=_init) as pool:
        records = pool.map(task, [18, 20])
    records.sort(key=lambda r: r["sigma"])

    ok = len(records) == 2 and all(
        r.get("verdict") == "VACANT" and r.get("genuine_gap") == 0
        and r.get("tags") == r.get("derivative_rank")
        for r in records)

    payload = {
        "object": "A2(2) even-spin vacancy at sigma = 18, 20",
        "records": records,
        "ranks_preregistered": False,
        "rank_note": ("no exact ranks were retained by the 2026-08-02/03 "
                      "fraction-free/sparse runs; the tags recorded here "
                      "are established by this computation, with the "
                      "rigorous acceptance tags == d-rank at the prime"),
        "curve_point_t": str(T0),
        "primes_offered": PRIMES,
        "acceptance": ("for each sigma, verdict VACANT with tags == d-rank "
                       "=> genuine vacancy over Q => the A2(2) pairs "
                       "([I7,I11] at 18, [I7,I13] at 20) hold identically "
                       "over Q(t), by Jacobi + the specialization lemma; "
                       "with these two, all 21 ledger instances rest on "
                       "the registered GF(p) certifier"),
        "supersedes": ("the 'fraction-free / sparse+collapsed' row of "
                       "results/vacancy_final.md (2026-08-02/03 runs, no "
                       "retained artifact)"),
    }
    out = Path(__file__).resolve().parent.parent / "results"
    art = out / "vacancy_a2_1820_certificate.json"

    def _canon(d):
        # JSON round-trip comparison (see vacancy_cyl_121416.py: a raw ==
        # against the in-memory payload self-invalidates the pin).
        d = {k: v for k, v in d.items() if k != "provenance"}
        return json.dumps(json.loads(json.dumps(d)), sort_keys=True)

    if art.exists() and _canon(json.loads(art.read_text())) == _canon(payload):
        print("artifact of record unchanged (payload agrees); "
              "provenance preserved", flush=True)
    else:
        art.write_text(json.dumps(stamp(payload), indent=1) + "\n")
        print("artifact written", flush=True)

    print("A2 SIGMA-18/20 VACANCY CERTIFIED" if ok
          else "INCONCLUSIVE OR ERRORS", flush=True)
    sys.exit(0 if ok else 1)
