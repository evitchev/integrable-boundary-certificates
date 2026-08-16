"""Rerun of the six OOM'd mutual-commutativity tasks using the
multiplicity-class-collapsed enumerator (collapsed.py) routed into the
mixed engine.  Same points, sections, and fail-closed semantics as
mutual_checks.py."""
import json, sys, resource
from multiprocessing import Pool

LIMIT = 18 * 1024**3  # hard per-worker address-space cap: fail, never freeze


def _init():
    resource.setrlimit(resource.RLIMIT_AS, (LIMIT, LIMIT))
import mixed_engine, collapsed
mixed_engine.residue_mixed_mono = collapsed.residue_mixed_mono_smart
import mutual_checks as mc

TASKS = [("cyl", s, a, b) for s in (1, 2, 3) for a, b in ((7, 11), (9, 11))]

if __name__ == "__main__":
    with Pool(len(TASKS), initializer=_init) as pool:
        results = pool.map(mc.task, TASKS)
    ok = all(r[1] for r in results)
    json.dump({k: {"pass": p, "points": n, "skipped": sk}
               for k, p, n, sk in results},
              open("../results/heavy_pairs.json" if ok
                   else "../results/heavy_pairs.partial.json", "w"), indent=1)
    print("ALL PASSED" if ok else "FAILURES", flush=True)
    sys.exit(0 if ok else 1)
