"""Discriminant factorization at weights 10 and 12 -- fail-closed
certificate (preregistered: results/cyl_discriminant_preregistration.md).

Drives wolfram/lab/cyl_discriminant.wls (the certified-w14 lab method:
one maximal minor D(t) is a complete over-approximation of the Q-bar
rank-drop locus; every irreducible factor classified by exact rank,
GF(p)-screened for irrational factors), then pins:

  H_F1  jump sets over Q-bar: w10 -- sol1 {}, sol2 {3}, sol3 {3, 7/5};
        w12 -- sol1 {}, sol2 {3, 9/7}, sol3 {3}; nullity 2 at every jump;
        every factor classified (no unclassified), poles only t = +-1.
  H_F2  beta^2 images exactly 1/2 (t=3), 1/6 (w10 sol3), 1/8 (w12 sol2).
  H_F3  the odd-q Kac slot is empty ALGEBRAICALLY: no jump anywhere at
        beta^2 = 1/7, 1/5, 1/3 (t = 4/3, 3/2, 2) -- and specifically
        (12, sol3) does not jump at 4/3 although 2q-3 = 11 at q = 7.

Inputs sha-verified against the tracked sidecars' matrix_sha256.
Writes results/cyl_discriminant_w1012.json (deterministic).  Note: the
run needs local Mathematica; per the license hazard the script sets NO
subprocess timeout (the wls bounds its own Det step; the registry
timeout is the backstop).  Marker: DISCRIMINANT W10/W12 CERTIFIED."""
import hashlib
import json
import os
import subprocess
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.join(HERE, os.pardir, "wolfram", "lab")
OUT = os.path.join(HERE, os.pardir, "results", "cyl_discriminant_w1012.json")

fails = []


def check(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg, flush=True)
    if not cond:
        fails.append(msg)


EXPECT_JUMPS = {("10", "1"): set(), ("10", "2"): {"3"}, ("10", "3"): {"3", "7/5"},
                ("12", "1"): set(), ("12", "2"): {"3", "9/7"}, ("12", "3"): {"3"}}
X_OF = {"3": "1/2", "7/5": "1/6", "9/7": "1/8"}
for t0, x0 in X_OF.items():
    tq = F(t0)
    if F(x0) != (tq - 1) / (tq + 1):   # fail-closed under -O too (round 15)
        print(f"FAIL  x0 table inconsistent at t0 = {t0}")
        sys.exit(1)

if __name__ == "__main__":
    for w in ("10", "12"):
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
    proc = subprocess.run(["wolframscript", "-file", wls], capture_output=True, text=True)
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
          "H_F3: no jump anywhere at beta^2 = 1/7, 1/5 or 1/3 -- the odd-q Kac slots "
          "(incl. q = 7 solving 2q-3 = 11) are empty ALGEBRAICALLY")
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
    print("DISCRIMINANT W10/W12 CERTIFIED (jump sets over Q-bar exactly as preregistered on all "
          "three sheets at weights 10 and 12; beta^2 images 1/2, 1/6, 1/8; the odd-q Kac slots "
          "empty algebraically)")
