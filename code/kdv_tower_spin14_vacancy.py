"""Spin-14 vacancy of BOTH CENTRALIZERS at the c=1 branch point,
re-derived in-stack -- the manuscript's actual finite-range objects.

History: a first version of this certificate computed the quartet and
pair SCREENING kernels at spin 14.  Review (2026-08-26) correctly
objected that at an exceptional coupling emptiness of a screening
kernel does not prove emptiness of the corresponding centralizer (the
enhancement phenomenon is precisely a centralizer exceeding a screening
kernel's span).  This version computes the two commutator kernels the
manuscript's "spin 14 vacant" claim is actually about, with the SAME
machinery and normalizations as the certified sigma<=9 range
(kdv_tower_classes.py):

  (i)  ker ad_{I5} (the tower side): genuine (mod-d) charges of spin 14
       commuting with the A2(2) defining charge oint P6 at the branch
       point -- expected dimension 0;
  (ii) ker ad_{I'3} (the shifted-KdV side): same for the enhanced
       spin-3 ray -- expected dimension 0.

Both dimensions are enforced separately; 0 == 0 is then the spin-14
instance of the class-by-class equality.  This also closes the vendored
independent program's conditional gap (06_spin14_15_extension.py skips
the second kernel when the first is empty).

Provenance: launch-time state is captured at process start; the
artifact of record is payload-stable (a run whose scientific payload
agrees leaves the file byte-identical, so the registry pin survives
re-runs; disagreement rewrites and breaks the pin -- the intended
alarm).
"""
import hashlib
import json
import subprocess
import sys
from fractions import Fraction as F
from pathlib import Path

from certificate import stamp
from families import L
from ope_engine import Sector, gen_basis
from solve_charges import genuine_kernel

failures = []


def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok:
        failures.append(name)


def acc(pairs):
    out = {}
    for mono, c in pairs:
        out[mono] = out.get(mono, 0) + c
    return {m: c for m, c in out.items() if c}


# the A2(2) I5 density at N=1, exactly as in kdv_tower_classes.py
A2_P6 = acc([(L((3, 0), (3, 0)), F(124, 96)),
             (L((2, 0), (2, 0), (1, 0), (1, 0)), F(-6)),
             (L((1, 0), (2, 0), (1, 0), (2, 0)), -F(-68, 8)),
             (L((1, 0),) * 6, F(1))])

if __name__ == "__main__":
    _root = Path(__file__).resolve().parent.parent
    from launch_provenance import capture_launch
    _launch = capture_launch(_root, Path(__file__).resolve())

    ker3 = genuine_kernel(A2_P6, gen_basis(4, 1, True),
                          Sector(9, 1, True), Sector(4, 1, True))
    check("control: I'3 unique at sigma=3 in ker ad_{I5} (reproduces the "
          "certified low-spin computation)", len(ker3) == 1)
    I3p = ker3[0]

    w = 15
    kerA = genuine_kernel(A2_P6, gen_basis(w, 1, True),
                          Sector(w + 5, 1, True), Sector(w, 1, True))
    dA = len(kerA)
    check(f"tower side: ker ad_{{I5}} at charge spin 14 has genuine class "
          f"dimension {dA} == 0", dA == 0)
    kerK = genuine_kernel(I3p, gen_basis(w, 1, True),
                          Sector(w + 3, 1, True), Sector(w, 1, True))
    dK = len(kerK)
    check(f"shifted-KdV side: ker ad_{{I'3}} at charge spin 14 has genuine "
          f"class dimension {dK} == 0 (the half the vendored program "
          "skips)", dK == 0)
    check("spin-14 class equality: both centralizers empty (0 == 0)",
          dA == 0 and dK == 0)

    print()
    if failures:
        print(f"SPIN-14 CENTRALIZER TEST FAILED ({len(failures)})")
        sys.exit(1)

    payload = {
        "object": ("spin-14 (weight-15) genuine commutator kernels at the "
                   "c=1 branch point: ker ad_{I5} (tower side) and "
                   "ker ad_{I'3} (shifted KdV side), modulo total "
                   "derivatives, one boson"),
        "tower_ker_ad_I5_dim": dA,
        "kdv_ker_ad_I3p_dim": dK,
        "acceptance": "both centralizer dimensions independently zero",
        "launch": _launch,
        "script_sha256": _launch["script_sha256_at_start"],
        "invocation": _launch["invocation"],
    }
    art = _root / "results" / "kdv_tower_spin14_vacancy.json"

    def _canon(d):
        d = {k: v for k, v in d.items()
         if k not in ("provenance", "launch", "script_sha256",
                      "invocation")}   # provenance-class keys: a
    # provenance-only refactor must not fake a scientific change
    # (review, 2026-08-27)
        return json.dumps(json.loads(json.dumps(d)), sort_keys=True)

    if art.exists() and _canon(json.loads(art.read_text())) == _canon(payload):
        print("artifact of record unchanged (payload agrees); "
              "provenance preserved", flush=True)
    else:
        art.write_text(json.dumps(stamp(payload), indent=1) + "\n")
        print("artifact written", flush=True)

    print("SPIN-14 CENTRALIZER VACANCY CERTIFIED IN-STACK (ker ad_{I5} and "
          "ker ad_{I'3} both empty at charge spin 14, computed directly "
          "and unconditionally)")
