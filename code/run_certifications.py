"""Unified certification runner (peer-review audit defects 7-8).

Every numbered claim of draft/main.tex (R1-R8) maps, via
certifications.json, to retained automated checks and/or hash-pinned
result artifacts. One command re-runs everything cheap and verifies
every pinned artifact; --full also re-executes the expensive
regenerations (compute_p8, promote_qt, the multiprocessing scans).

    python3 run_certifications.py            # fast checks + hash pins
    python3 run_certifications.py --full     # + expensive regenerations
    python3 run_certifications.py --only ID  # a single entry

Output: results/certification_manifest.json — per-entry status, exit
codes, runtimes, marker hits, artifact hashes, environment capture,
and a claim-coverage map. Exit code is nonzero if ANY executed check
fails, ANY required marker is missing, or ANY pinned hash mismatches;
the manifest is still written (recording failure is the point).
Documented coverage gaps (claims certified narratively, or data-backed
conjecture) are listed explicitly in the registry and echoed in the
manifest rather than silently omitted.

Design notes: entries with mode "run" execute by default; mode
"artifact" entries only verify pinned hashes unless --full, because
their regeneration costs minutes to hours (promote_qt: ~160 min on 11
workers). Regeneration writes the scripts' fixed output paths in
results/; a post-regeneration hash mismatch is reported as a failure —
if a regeneration is legitimately nondeterministic, that fact should
surface here, not be absorbed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from datetime import date, datetime, timezone
from pathlib import Path

CODE_DIR = Path(__file__).resolve().parent
ROOT = CODE_DIR.parent
REGISTRY_PATH = CODE_DIR / "certifications.json"
MANIFEST_PATH = ROOT / "results" / "certification_manifest.json"
FORBIDDEN_DEFAULT = ("Traceback", "FAIL", "MISMATCH")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return "sha256:" + digest.hexdigest()


def git_commit() -> str:
    result = subprocess.run(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD"],
        capture_output=True, text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else "unknown"


def environment() -> dict:
    try:
        import sympy
        sympy_version = sympy.__version__
    except ImportError:
        sympy_version = "MISSING"
    return {
        "python": sys.version.split()[0],
        "python_executable": sys.executable,
        "sympy": sympy_version,
        "platform": platform.platform(),
        "machine": platform.machine(),
    }


def worktree_state() -> dict:
    """Dirty-worktree recording (certifier upgrade, Codex round 4).

    A manifest generated from a dirty tree says so: `git_commit` alone
    over-promises reproducibility when uncommitted edits were present.
    """

    status = subprocess.run(
        ["git", "-C", str(ROOT), "status", "--porcelain"],
        capture_output=True, text=True,
    )
    if status.returncode != 0:
        # provenance must fail UNKNOWN, never clean (audit, round 7)
        return {
            "worktree_dirty": None,
            "provenance_error": (f"git status exited {status.returncode}: "
                                 f"{status.stderr.strip()[:200]}"),
        }
    paths = [line[3:] for line in status.stdout.splitlines() if line]
    tracked = [line[3:] for line in status.stdout.splitlines()
               if line and not line.startswith("??")]
    untracked = [line[3:] for line in status.stdout.splitlines()
                 if line.startswith("??")]
    return {
        "worktree_dirty": bool(paths),
        "dirty_path_count": len(paths),
        "dirty_paths": paths[:40],
        # tracked-dirty means the committed sources themselves differ from
        # `git_commit`; untracked files can still affect execution (modules,
        # configuration, artifacts), so both are recorded (rounds 5, 7)
        "tracked_dirty": bool(tracked),
        "tracked_dirty_paths": tracked[:40],
        "untracked_paths": untracked[:40],
    }


# Every key the runner actually consults.  A key outside this set is a
# typo, and a typo in a pin key is silent: the entry would appear to be
# pinned while nothing was ever hashed (observed 2026-08-24, where an
# entry carried "raw_certificate_pins" -- a phrase from another entry's
# prose -- instead of "artifacts").  Fail closed instead.
KNOWN_ENTRY_KEYS = {
    "id", "claims", "description", "mode", "command", "cwd", "timeout_s",
    "expensive", "required_markers", "forbidden_markers", "artifacts",
}


def run_entry(entry: dict, *, execute: bool) -> dict:
    """Execute (or hash-verify) one certification entry."""

    unknown = sorted(set(entry) - KNOWN_ENTRY_KEYS)
    if unknown:
        return {
            "id": entry.get("id", "<no id>"),
            "claims": entry.get("claims", []),
            "description": entry.get("description", ""),
            "mode": entry.get("mode", "<none>"),
            "status": "failed",
            "artifact_hashes": {},
            "failures": [f"unknown_entry_key:{k}" for k in unknown],
        }

    result = {
        "id": entry["id"],
        "claims": entry.get("claims", []),
        "description": entry.get("description", ""),
        "mode": entry["mode"],
        "executed": False,
        "status": "pass",
        "failures": [],
        "runtime_s": None,
        "artifact_hashes": {},
    }

    command = entry.get("command")
    if execute and command:
        argv = [sys.executable if token == "PYTHON" else token
                for token in command]
        workdir = (ROOT / entry["cwd"]) if entry.get("cwd") else CODE_DIR
        started = time.monotonic()
        try:
            proc = subprocess.run(
                argv, cwd=str(workdir), capture_output=True, text=True,
                timeout=entry.get("timeout_s", 1800),
            )
            rc, output = proc.returncode, proc.stdout + proc.stderr
        except subprocess.TimeoutExpired:
            rc, output = -1, "<timeout>"
        result["runtime_s"] = round(time.monotonic() - started, 1)
        result["executed"] = True
        result["exit_code"] = rc
        if rc != 0:
            result["failures"].append(f"exit_code_{rc}")
        for marker in entry.get("required_markers", []):
            if marker not in output:
                result["failures"].append(f"missing_marker:{marker}")
        for marker in entry.get("forbidden_markers", FORBIDDEN_DEFAULT):
            if marker in output:
                result["failures"].append(f"forbidden_marker:{marker}")
    elif command:
        result["status"] = "skipped_expensive"

    for artifact in entry.get("artifacts", []):
        path = ROOT / artifact["path"]
        if not path.exists():
            result["failures"].append(f"artifact_missing:{artifact['path']}")
            continue
        actual = sha256_file(path)
        result["artifact_hashes"][artifact["path"]] = actual
        if artifact.get("sha256") and actual != artifact["sha256"]:
            result["failures"].append(f"hash_mismatch:{artifact['path']}")

    if (
        result["status"] == "skipped_expensive"
        and entry.get("artifacts")
        and not result["failures"]
    ):
        # The pins WERE verified; only the regeneration was deferred.
        result["status"] = "pass_pins_verified"
    if result["failures"]:
        result["status"] = "fail"
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--full", action="store_true",
                        help="also run expensive regenerations")
    parser.add_argument("--only", default=None,
                        help="run a single certification id")
    args = parser.parse_args()

    registry = json.loads(REGISTRY_PATH.read_text(encoding="utf-8"))
    entries = registry["certifications"]
    if args.only:
        entries = [e for e in entries if e["id"] == args.only]
        if not entries:
            print(f"no certification id {args.only!r}", file=sys.stderr)
            return 2

    results = []
    for entry in entries:
        execute = (
            (entry["mode"] == "run" and not entry.get("expensive"))
            or args.full
            or bool(args.only)
        )
        label = "RUN " if (execute and entry.get("command")) else "PIN "
        print(f"[{label}] {entry['id']} ...", flush=True)
        outcome = run_entry(entry, execute=execute)
        status = outcome["status"].upper()
        runtime = (f" ({outcome['runtime_s']}s)"
                   if outcome["runtime_s"] is not None else "")
        print(f"        {status}{runtime}"
              + (f"  {outcome['failures']}" if outcome["failures"] else ""),
              flush=True)
        results.append(outcome)

    coverage: dict[str, list[str]] = {}
    for entry in registry["certifications"]:
        for claim in entry.get("claims", []):
            coverage.setdefault(claim, []).append(entry["id"])

    failed = [r for r in results if r["status"] == "fail"]
    manifest = {
        "schema_version": 1,
        "generated_on": date.today().isoformat(),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(
            timespec="seconds"),
        "git_commit": git_commit(),
        "worktree": worktree_state(),
        "environment": environment(),
        "invocation": {
            "full": args.full,
            "only": args.only,
        },
        "summary": {
            "total": len(results),
            "passed": sum(
                1 for r in results
                if r["status"] in ("pass", "pass_pins_verified")),
            "failed": len(failed),
            "skipped_expensive": sum(
                1 for r in results if r["status"] == "skipped_expensive"),
            "pins_verified": sum(
                1 for r in results if r["status"] == "pass_pins_verified"),
        },
        "status": "failed" if failed else (
            "passed_fast_checks_and_pins" if not args.full
            else "passed_full"),
        "claim_coverage": dict(sorted(coverage.items())),
        "coverage_gaps": registry.get("coverage_gaps", []),
        "results": results,
    }
    # --only runs a subset: write a partial manifest alongside, never
    # replacing the general one (certifier upgrade, Codex round 4).
    manifest_path = (
        MANIFEST_PATH.with_name("certification_manifest.partial.json")
        if args.only else MANIFEST_PATH
    )
    if args.only:
        manifest["partial_of"] = MANIFEST_PATH.name
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    staged = manifest_path.with_suffix(".json.tmp")
    staged.write_text(json.dumps(manifest, indent=2) + "\n",
                      encoding="utf-8")
    staged.replace(manifest_path)
    print(f"\nmanifest: {manifest_path}")
    s = manifest["summary"]
    print(f"status: {manifest['status']} "
          f"({s['passed'] - s['pins_verified']} executed passes, "
          f"{s['pins_verified']} pins-verified, "
          f"{s['skipped_expensive']} skipped-expensive, "
          f"of {s['total']})")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
