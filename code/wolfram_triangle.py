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
    # freshness (audit fix): remove any prior artifact so a failed Export
    # cannot leave a stale file in place, and bind the output to THIS
    # script and THIS invocation afterwards
    import time as _time
    t_start = _time.time()
    # the committed file is the artifact of record (the implementer's
    # execution, pinned by the typed-claims mirror); keep it, require a
    # FRESH execution to agree with it value by value, and restore it
    reference = OUTPUT.read_bytes() if OUTPUT.exists() else None
    try:
        return _fresh_run_and_compare(ws, reference, t_start)
    finally:
        # byte-exact restoration on EVERY exit (failure, timeout, missing or
        # malformed output, stale output, hash mismatch, exception): the
        # artifact of record is never left deleted (audit fix, round 4)
        if reference is not None:
            OUTPUT.write_bytes(reference)


def _fresh_run_and_compare(ws, reference, t_start):
    if OUTPUT.exists():
        OUTPUT.unlink()
    run = subprocess.run([ws, "-file", str(SCRIPT)], capture_output=True,
                         text=True, timeout=3600)
    if run.returncode != 0:
        print("WOLFRAM RUN FAILED:")
        print(run.stdout[-2000:], run.stderr[-2000:])
        return 1
    if not OUTPUT.exists():
        print(f"REPLICATION OUTPUT MISSING: {OUTPUT.relative_to(ROOT)}")
        return 1
    if OUTPUT.stat().st_mtime < t_start:
        print("REPLICATION OUTPUT STALE: predates this invocation")
        return 1
    got = json.loads(OUTPUT.read_text())
    script_hash = sha256(SCRIPT)
    recorded = {str(v.get("script_sha256", "")) for v in got.values()
                if isinstance(v, dict) and v.get("script_sha256")}
    if not recorded or any(script_hash.split(":")[1] not in h for h in recorded):
        print(f"SCRIPT HASH MISMATCH: output records {recorded}, current "
              f"script is {script_hash}")
        return 1
    if reference is not None:
        ref = json.loads(reference.decode())
        fresh_vals = {k: v.get("value") for k, v in got.items() if isinstance(v, dict)}
        ref_vals = {k: v.get("value") for k, v in ref.items() if isinstance(v, dict)}
        if fresh_vals != ref_vals:
            print(f"FRESH RUN DISAGREES WITH THE ARTIFACT OF RECORD: {fresh_vals} vs {ref_vals}")
            return 1
        fresh_hash = sha256(OUTPUT)            # (the caller restores the artifact of record)
        print(f"fresh execution agrees with the artifact of record; fresh output "
              f"pin {fresh_hash}")
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
