"""Verification-triangle runner: compare the independent Mathematica
replication (results/wolfram_triangle.json, produced by the clean-
context implementer from wolfram/HANDOFF-SPEC.md)
against the baked targets in code/wolfram_targets.json.

Mirrors lean_check.py: this is an EXPENSIVE, machine-local entry.  The
Wolfram toolchain is proprietary and not part of the portable
sympy-only core; where it is absent the check FAILS LOUDLY, never
silently.  The certificate authority stays in the open stack; this
runner adds a decorrelated check (unrelated implementation, vendor,
and algorithm lineage) and pins the replication output by hash.

Registered only once a kernel exists and results/wolfram_triangle.json
has been produced by an executed script (never register a replication
that has not run -- the fabricated-provenance lesson).
"""
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
TARGETS = ROOT / "code" / "wolfram_targets.json"
OUTPUT = ROOT / "results" / "wolfram_triangle.json"
SCRIPT = ROOT / "wolfram" / "triangle_targets.wls"


def sha256(path):
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ws = shutil.which("wolframscript")
    if ws is None:
        print("WOLFRAM TOOLCHAIN MISSING: wolframscript not on PATH")
        print("(machine-local toolchain; activate against the subscription "
              "or install the free Wolfram Engine; see the hand-off spec)")
        return 1
    if not SCRIPT.exists():
        print(f"REPLICATION SCRIPT MISSING: {SCRIPT.relative_to(ROOT)}")
        return 1
    run = subprocess.run([ws, "-file", str(SCRIPT)], capture_output=True,
                         text=True, timeout=3600)
    if run.returncode != 0:
        print("WOLFRAM RUN FAILED:")
        print(run.stdout[-2000:], run.stderr[-2000:])
        return 1
    if not OUTPUT.exists():
        print(f"REPLICATION OUTPUT MISSING: {OUTPUT.relative_to(ROOT)}")
        return 1
    got = json.loads(OUTPUT.read_text())
    targets = json.loads(TARGETS.read_text())["targets"]
    failures = 0
    for cid, tgt in targets.items():
        entry = got.get(cid)
        if entry is None:
            print(f"FAIL  {cid}: missing from replication output")
            failures += 1
            continue
        val, kind, exp = entry.get("value"), tgt["kind"], tgt["expected"]
        if kind in ("boolean",):
            ok = (val is exp)
        elif kind == "boolean_pair":
            ok = (list(val) == exp) if isinstance(val, list) else False
        elif kind == "list":
            ok = (list(val) == exp)
        elif kind == "set":
            ok = ({sp.nsimplify(str(x)) for x in val}
                  == {sp.nsimplify(x) for x in exp})
        elif kind == "exact":
            ok = sp.simplify(sp.sympify(str(val).replace("Sqrt[", "sqrt(")
                                        .replace("]", ")"))
                             - sp.sympify(exp.replace("Sqrt[", "sqrt(")
                                          .replace("]", ")"))) == 0
        else:
            ok = False
        if entry.get("notes"):
            print(f"      note on {cid}: {entry['notes']}")
        print(f"{'PASS' if ok else 'FAIL'}  {cid}: {tgt['statement']} "
              f"-> {val}")
        failures += 0 if ok else 1
    print(f"replication output pin: {sha256(OUTPUT)}")
    print(f"script pin: {sha256(SCRIPT)}; engine: "
          f"{got.get('T1', {}).get('engine', '?')}")
    if failures:
        print(f"WOLFRAM TRIANGLE FAILED ({failures})")
        return 1
    print("WOLFRAM TRIANGLE VERIFIED (independent Mathematica replication "
          "agrees with every baked target)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
