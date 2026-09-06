"""Discriminant factorization at weight 14 (spin 13) -- fail-closed
certificate (preregistered: results/cyl_discriminant_w14_preregistration.md).

Third rung of the discriminant program (after weights 10/12,
cyl_discriminant_w1012.py).  Drives wolfram/lab/cyl_discriminant.wls
with weight argument 14 (the certified-w14 lab method: one maximal
minor D(t) is a complete over-approximation of the Q-bar rank-drop
locus away from the parametrization poles t = +-1; every irreducible
factor classified by exact rank, GF(p)-screened for irrational
factors), then pins:

  H_W14  spin-13 jump sets over Q-bar (t != +-1): sol1 {}, sol2 {3},
         sol3 {3, 9/7}; nullity 2 at every jump; every factor
         classified (no unclassified); cleared-denominator factors
         only at t = +-1.
  H_W14  beta^2 images exactly 1/2 (t = 3) and, on sol3, 1/8 (t = 9/7 =
         c_{1,8}: the (1-8x) factor at its second Kac slot, 2q-3 = 13).
  H_W14  no jump anywhere at beta^2 = 1/7, 1/5, 1/3; sheet 2 silent
         beyond c_{1,2} ((3q-2)/2 = 13 has no integer solution) --
         pinned over all of Q-bar minus the poles.

Inputs (three tracked matrices + sidecars) sha-verified; the wls and
inputs are declared sources bound into the launch record.  Writes
results/cyl_discriminant_w14.json (deterministic).  No subprocess
timeout (license hazard; the wls bounds its own Det step; the registry
timeout is the backstop).  Marker: DISCRIMINANT W14 CERTIFIED."""
import hashlib
import json
import os
import subprocess
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.join(HERE, os.pardir, "wolfram", "lab")
OUT = os.path.join(HERE, os.pardir, "results", "cyl_discriminant_w14.json")

fails = []


def check(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg, flush=True)
    if not cond:
        fails.append(msg)


EXPECT_JUMPS = {("14", "1"): set(), ("14", "2"): {"3"}, ("14", "3"): {"3", "9/7"}}
X_OF = {"3": "1/2", "9/7": "1/8"}
for t0, x0 in X_OF.items():
    tq = F(t0)
    if F(x0) != (tq - 1) / (tq + 1):   # fail-closed under -O too (round 15)
        print(f"FAIL  x0 table inconsistent at t0 = {t0}")
        sys.exit(1)

if __name__ == "__main__":
    for w in ("14",):
        for sol in ("1", "2", "3"):
            m = os.path.join(LAB, f"cyl_strata_w{w}_sol{sol}.m")
            side = json.load(open(os.path.join(LAB, f"cyl_strata_w{w}_sol{sol}.json")))
            sha = hashlib.sha256(open(m, "rb").read()).hexdigest()
            check(sha == side["matrix_sha256"],
                  f"w{w} sol{sol}: matrix sha matches the tracked sidecar ({sha[:12]})")
    if fails:
        print(f"DISCRIMINANT: {len(fails)} INPUT FAILURE(S)")
        sys.exit(1)
    wls = os.path.join(LAB, "cyl_discriminant.wls")
    print("running wolframscript (no subprocess timeout; wls bounds its own steps) ...", flush=True)
    proc = subprocess.run(["wolframscript", "-file", wls, "14"], capture_output=True, text=True)
    sys.stdout.write(proc.stdout)
    if proc.returncode != 0 or "DISCRIMINANT-WLS DONE" not in proc.stdout:
        print(f"FAIL  wolframscript exit {proc.returncode} / marker missing; stderr tail: "
              + proc.stderr[-300:])
        sys.exit(1)
    raw = json.load(open(os.path.join(LAB, "cyl_discriminant_w1012_wls.json")))
    pinned = {}
    for (w, sol), exp in sorted(EXPECT_JUMPS.items()):
        e = raw[f"w{w}_sol{sol}"]
        jt = {j["t0"] for j in e["jumps"] if "t0" in j}
        irr = [j["factor"] for j in e["jumps"] if "t0" not in j]
        check(not irr, f"w{w} sol{sol}: no irrational jump factors (H_F1) (got {irr})")
        check(jt == exp, f"w{w} sol{sol}: jump set over Q-bar exactly {sorted(exp)} (got {sorted(jt)})")
        check(all(j["nullity"] == 2 for j in e["jumps"]),
              f"w{w} sol{sol}: nullity 2 at every jump")
        check(not e["unclassified"], f"w{w} sol{sol}: every factor classified (got {e['unclassified']})")
        check(all(p in ("-1 + t", "1 + t") for p in e["poles"]),
              f"w{w} sol{sol}: cleared-denominator factors only at t = +-1 (got {e['poles']})")
        xs = {j.get("x0") for j in e["jumps"]}
        check(xs == {X_OF[t] for t in exp},
              f"w{w} sol{sol}: beta^2 images exactly {sorted(X_OF[t] for t in exp)} (H_F2)")
        pinned[f"w{w}_sol{sol}"] = {k: e[k] for k in
                                    ("rows", "cols", "generic_rank", "generic_multiplicity",
                                     "det_degree", "jumps", "nonjump_factors", "poles")}
    all_jump_x = {j.get("x0") for (w, sol) in EXPECT_JUMPS
                  for j in raw[f"w{w}_sol{sol}"]["jumps"]}
    check(not ({"1/7", "1/5", "1/3"} & all_jump_x),
          "H_W14: no jump anywhere at beta^2 = 1/7, 1/5 or 1/3; sheet 2 silent beyond c_{1,2} "
          "((3q-2)/2 = 13 has no integer q) -- pinned over Q-bar")
    pinned["provenance"] = raw["provenance"]
    print()
    if fails:
        # No artifact on a failed run (Codex-style hygiene: the drill once
        # wrote a partial artifact to the production path via the old
        # write-before-gate order).
        print(f"DISCRIMINANT: {len(fails)} PIN(S) NOT MET -- see FINDING lines")
        sys.exit(1)
    with open(OUT, "w") as fh:
        json.dump(pinned, fh, indent=1, sort_keys=True)
        fh.write("\n")
    print("DISCRIMINANT W14 CERTIFIED (spin-13 jump sets over Q-bar exactly as preregistered: "
          "sheet 1 empty, sheet 2 {3}, sheet 3 {3, 9/7} with the (1-8x) factor; nullity 2 at "
          "every jump; sheet 2 silent beyond c_{1,2} -- (3q-2)/2 = 13 has no integer q -- "
          "pinned over all of Q-bar")
