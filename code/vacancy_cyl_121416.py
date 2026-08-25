"""GF(p) vacancy certificates for cylindrical sigma = 12, 14, 16.

These nine verdicts are the hypothesis under the four LOWEST cylindrical
mutual-commutativity pairs -- [I5,I7] at summed spin 12, [I5,I9] at 14,
[I5,I11] and [I7,I9] at 16 -- and until now they survived only as a
narrative table in results/vacancy_final.md, produced by code/vacancy_high.py,
which is no longer in the working tree (it exists only in git history at
6ebd63c).  The manuscript calls them "certified", which under this project's
own definition means registered and re-runnable; they were not.  This script
closes that gap by recomputing all nine under the same GF(p) certifier that
produced the sigma=18/20/24 verdicts, and retaining a per-task artifact.

Epistemics are those of genuine_kernel_modp: tags == d-rank at the prime is a
RIGOROUS vacancy certificate over Q (kernels only grow and ranks only drop
mod p, so equality pins the exact dimensions).  A gap escalates to a second
prime and, if it persists, is reported as an inconclusive upper bound -- never
as a non-vacancy claim -- and the batch fails closed.
"""
import json
import resource
import sys
from fractions import Fraction as F
from multiprocessing import Pool
from pathlib import Path

import sympy as sp

import collapsed
import mixed_engine
mixed_engine.residue_mixed_mono = collapsed.residue_mixed_mono_smart
from certificate import stamp
from invariant_engine import curve_point
from mixed_engine import (cyl_P4, d_mixed_mono, gen_mixed_basis,
                          residue_mixed_cur, spin_mixed)
from sparse_linalg import MODP_PRIME, genuine_kernel_modp

T0 = F(2)
PRIMES = [MODP_PRIME, int(sp.nextprime(MODP_PRIME))]
CAP = 10 * 1024**3


def _init():
    resource.setrlimit(resource.RLIMIT_AS, (CAP, CAP))


def task(args):
    sol, sigma = args
    N, s = curve_point(T0)
    label = f"cyl{sol} s{sigma}"
    rec = {"solution": sol, "sigma": sigma}
    try:
        for pi, prime in enumerate(PRIMES):
            tags, drank, vac = genuine_kernel_modp(
                cyl_P4(sol, N, s), gen_mixed_basis(sigma + 1), N,
                residue_mixed_cur, gen_mixed_basis, d_mixed_mono,
                spin_mixed, prime=prime, progress=f"{label} p{pi}")
            # Record the prime actually used, not a single top-level one:
            # an escalation to prime #1 must not be reported under #0.
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
    # The exact ranks are what the manuscript prints; success is enforced
    # against them, not merely against a zero gap.
    EXPECTED_TAGS = {12: 379, 14: 929, 16: 2211}

    tasks = [(sol, sigma) for sigma in (12, 14, 16) for sol in (1, 2, 3)]
    with Pool(min(9, len(tasks)), initializer=_init) as pool:
        records = pool.map(task, tasks)
    records.sort(key=lambda r: (r["sigma"], r["solution"]))

    ok = all(
        r.get("verdict") == "VACANT"
        and r.get("genuine_gap") == 0
        and r.get("tags") == EXPECTED_TAGS[r["sigma"]]
        and r.get("derivative_rank") == EXPECTED_TAGS[r["sigma"]]
        for r in records
    ) and len(records) == 9
    for r in records:
        exp = EXPECTED_TAGS[r["sigma"]]
        if r.get("tags") not in (exp, None) :
            print(f"RANK MISMATCH cyl{r['solution']} s{r['sigma']}: "
                  f"tags {r.get('tags')} != expected {exp}", flush=True)

    payload = {
        "object": "cylindrical even-spin vacancy at sigma = 12, 14, 16",
        "records": records,
        "expected_tags": EXPECTED_TAGS,
        "curve_point_t": str(T0),
        "primes_offered": PRIMES,
        "acceptance": ("for each task, verdict VACANT with tags == d-rank == "
                       "the expected rank => genuine vacancy over Q => the "
                       "four lowest cylindrical pairs ([I5,I7] at 12, "
                       "[I5,I9] at 14, [I5,I11] and [I7,I9] at 16) hold "
                       "identically over Q(t) for all three solutions, by "
                       "Jacobi + the specialization lemma"),
        "supersedes": ("the narrative table in results/vacancy_final.md, "
                       "produced by code/vacancy_high.py (git history only, "
                       "6ebd63c), which retained no per-task artifact"),
    }
    out = Path(__file__).resolve().parent.parent / "results"
    out.mkdir(exist_ok=True)
    art = out / "vacancy_cyl_121416_certificate.json"
    # Artifact of record.  The scientific payload is deterministic; only the
    # provenance timestamp moves.  Rewriting on every run would churn the
    # file hash and self-invalidate the registry pin (observed 2026-08-24),
    # so an artifact whose payload agrees is LEFT UNTOUCHED: its provenance
    # records when the verdict was established, not when it was last
    # re-verified.  A payload disagreement is a real event, is written, and
    # breaks the pin -- which is the intended alarm.
    def _canon(d):
        # Compare through a JSON round-trip.  A raw == against the in-memory
        # payload is WRONG: json turns the integer keys of expected_tags into
        # strings, so the comparison would always fail and the artifact would
        # be rewritten on every run -- silently reinstating the very
        # self-invalidation this block exists to prevent (caught 2026-08-24).
        d = {k: v for k, v in d.items() if k != "provenance"}
        return json.dumps(json.loads(json.dumps(d)), sort_keys=True)

    if art.exists() and _canon(json.loads(art.read_text())) == _canon(payload):
        print("artifact of record unchanged (payload agrees); "
              "provenance preserved", flush=True)
    else:
        art.write_text(json.dumps(stamp(payload), indent=1) + "\n")
        print("artifact written", flush=True)

    print("CYL SIGMA-12/14/16 VACANCY CERTIFIED" if ok
          else "INCONCLUSIVE OR ERRORS", flush=True)
    sys.exit(0 if ok else 1)
