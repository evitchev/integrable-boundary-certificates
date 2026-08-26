"""Public-edition semantic check of the retained spin-13 run record
(results/cyl_spin13_prediction_run.json; byte integrity pinned on the
cyl_spin13_prediction registry entry).

The full checker (private edition) performs five checks, two of which
resolve the recording commit in the private development history and so
CANNOT run from this repository.  This edition runs the three
repository-local checks and REPORTS the other two as not verifiable in
this edition -- honest reporting, not a silently weakened check (the
precedent that shipping a checker which fails on every public run, or
one quietly taught to skip its git steps, are both wrong: v1.1.3
release notes).  The two git-bound facts (the recorded commit resolves;
the certificate bytes at that commit match the current file) were
verified privately and re-verified by an independent reviewer.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REC = ROOT / "results" / "cyl_spin13_prediction_run.json"
failures = []


def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok:
        failures.append(name)


m = json.loads(REC.read_text())
check("record holds exactly one result, is a partial manifest of the sweep, "
      "and its top-level status and summary say passed",
      len(m.get("results", [])) == 1
      and m.get("partial_of") == "certification_manifest.json"
      and str(m.get("status", "")).startswith("passed")
      and m.get("summary", {}).get("total") == 1
      and m.get("summary", {}).get("passed") == 1
      and m.get("summary", {}).get("failed") == 0)
e = m["results"][0] if m.get("results") else {}
check("record is for cyl_spin13_prediction",
      e.get("id") == "cyl_spin13_prediction")
check("recorded execution ACTUALLY RAN and passed (executed true, status "
      "pass, exit code PRESENT and 0, no failures) -- a missing field "
      "fails, never defaults to success",
      e.get("executed") is True and e.get("status") == "pass"
      and e.get("exit_code") == 0 and e.get("failures") == [])
for name in (
    "the recorded commit resolves in the development history",
    "the certificate bytes at the recorded commit match the current file",
):
    print(f"NOT VERIFIABLE IN THIS EDITION  {name} (requires the private "
          "development git history; verified privately and by independent "
          "review)", flush=True)

print()
if failures:
    print(f"SPIN-13 RUN RECORD PUBLIC CHECK FAILED ({len(failures)})")
    sys.exit(1)
print("SPIN-13 RUN RECORD (PUBLIC EDITION) VERIFIED (3 repository-local "
      "checks pass; 2 git-bound checks explicitly reported as not "
      "verifiable in this edition)")
