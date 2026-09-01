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

import os
import sys

# ARCHIVE TRUSTED BOOTSTRAP (Codex rounds 2-3, 2026-08-28, each reproduced).
# In a git-free tree NOTHING from the tree may execute before the whole
# tree is verified: a sourceless code/fractions.pyc shadows the stdlib, a
# stamp-matched __pycache__ pyc is imported over the verified source, and
# an unmanifested code/argparse.py ran before the previous gate.  So:
#   1. re-exec in isolated mode (-I): the script directory is NOT on
#      sys.path and PYTHON* environment variables are ignored;
#   2. import only trusted stdlib modules (os, sys, hashlib, json);
#   3. verify EXPORT_RECORD.json + every file against RELEASE_MANIFEST.json
#      inline, refusing bytecode and unmanifested files;
#   4. only then expose code/ on sys.path (lazily, in main) and run
#      certificates with a sanitized environment and -B -s.
_CODE_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT_STR = os.path.dirname(_CODE_DIR)
_ARCHIVE = not os.path.exists(os.path.join(_ROOT_STR, ".git"))   # a linked worktree has a .git FILE (Codex round 4)
if _ARCHIVE:
    if not sys.flags.isolated:
        os.execv(sys.executable, [sys.executable, "-I", "-B",
                                  os.path.abspath(__file__)] + sys.argv[1:])
    sys.dont_write_bytecode = True
    import hashlib
    import json as _json

    def _archive_bootstrap():
        gen = {"RELEASE_MANIFEST.json", "results/certification_manifest.json",
               "results/certification_manifest.partial.json"}
        man_path = os.path.join(_ROOT_STR, "RELEASE_MANIFEST.json")
        rec_path = os.path.join(_ROOT_STR, "EXPORT_RECORD.json")
        if not (os.path.isfile(man_path) and os.path.isfile(rec_path)):
            sys.exit("archive gate REFUSED: no RELEASE_MANIFEST.json + "
                     "EXPORT_RECORD.json at the tree root (not a verified "
                     "release archive, not a repository)")
        with open(man_path, "rb") as fh:
            man_bytes = fh.read()
        try:
            files = _json.loads(man_bytes.decode("utf-8")).get("files", {})
            with open(rec_path, "rb") as fh:
                rec_bytes = fh.read()
            freeze = _json.loads(rec_bytes.decode("utf-8")).get(
                "private_freeze_commit", "")
        except (ValueError, UnicodeDecodeError) as exc:
            sys.exit(f"archive gate REFUSED: unreadable release record ({exc})")
        if not (isinstance(freeze, str) and len(freeze) == 40
                and all(c in "0123456789abcdef" for c in freeze)):
            sys.exit("archive gate REFUSED: release freeze commit missing/malformed")
        rec_pin = files.get("EXPORT_RECORD.json")
        if rec_pin != "sha256:" + hashlib.sha256(rec_bytes).hexdigest():
            sys.exit("archive gate REFUSED: EXPORT_RECORD.json does not match "
                     "its RELEASE_MANIFEST pin")
        problems = []
        for rel, pinned in sorted(files.items()):
            p = os.path.join(_ROOT_STR, rel)
            if not os.path.isfile(p):
                problems.append("missing: " + rel)
                continue
            with open(p, "rb") as fh:
                if "sha256:" + hashlib.sha256(fh.read()).hexdigest() != pinned:
                    problems.append("modified: " + rel)
        # directories that hold manifested files (all ancestors)
        expected_dirs = {""}
        for rel in files:
            parts = rel.split("/")[:-1]
            for i in range(1, len(parts) + 1):
                expected_dirs.add("/".join(parts[:i]))
        for d, dirs, fnames in os.walk(_ROOT_STR):
            reld = os.path.relpath(d, _ROOT_STR).replace(os.sep, "/")
            reld = "" if reld == "." else reld
            if os.path.basename(d) == "__pycache__":
                problems.append("bytecode: " + reld + "/")
                dirs[:] = []
                continue
            if reld and reld not in expected_dirs:
                problems.append("unmanifested directory: " + reld + "/")
            keep = []
            for x in dirs:
                if x == ".git":
                    continue
                if os.path.islink(os.path.join(d, x)):
                    problems.append("symlink: " + (reld + "/" if reld else "") + x + "/")
                    continue
                keep.append(x)
            dirs[:] = keep
            for f in fnames:
                full = os.path.join(d, f)
                rel = (reld + "/" if reld else "") + f
                if os.path.islink(full):
                    problems.append("symlink: " + rel)
                elif not os.path.isfile(full):
                    problems.append("non-regular: " + rel)
                elif f.endswith(".pyc"):
                    problems.append("bytecode: " + rel)
                elif f.endswith(".tmp") or rel in gen or rel in files:
                    continue
                else:
                    problems.append("unmanifested: " + rel)
        if problems:
            sys.exit("archive gate REFUSED: release tree does not match "
                     "RELEASE_MANIFEST.json: " + "; ".join(sorted(problems)[:10]))
        for k in [k for k in os.environ if k.startswith("PYTHON")]:
            del os.environ[k]                      # children: no injected paths
        os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
        # coordination flag for capture_launch()/stamp(): value = manifest
        # hash, re-checked there.  NOT an attestation (any process can set
        # it); its role is to make direct execution -- where verification
        # would be late -- refuse.  The runner's own manifest is the record.
        os.environ["IB_TRUSTED_BOOTSTRAP"] = hashlib.sha256(man_bytes).hexdigest()
        print(f"archive bootstrap: isolated mode, {len(files)} files verified, "
              f"freeze {freeze[:12]}")

    _archive_bootstrap()

import argparse
import hashlib
import json
import platform
import posixpath
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
    try:
        result = subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"],
            capture_output=True, text=True,
        )
    except FileNotFoundError:
        return "unavailable (no git executable)"
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

    try:
        status = subprocess.run(
            ["git", "-C", str(ROOT), "status", "--porcelain"],
            capture_output=True, text=True,
        )
    except FileNotFoundError:
        # provenance unknown, never clean (archive execution context)
        return {"worktree_dirty": None,
                "provenance_error": "no git executable"}
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
    "sources",            # extra source/input files bound into the launch record
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
    # LAUNCH-TIME capture (Codex round 7): the certificate source AND its
    # local import closure are hashed BEFORE the subprocess starts, with
    # the commit and tracked-worktree state; a multi-hour run can outlive
    # edits, so the record must describe what Python actually loaded.
    script = (command[1] if command and command[0] == "PYTHON"
              and len(command) > 1 else None)
    closure = []
    strict_run = execute and any(a.get("kind") == "runner_record"
                                 for a in entry.get("artifacts", []))
    if script and (CODE_DIR / script).is_file():
        closure, unresolved = _local_import_closure(script)
        for extra in entry.get("sources", []):
            if not (CODE_DIR / extra).is_file():
                result["failures"].append(f"declared_source_missing:{extra}")
            elif extra not in closure:
                closure.append(extra)
        closure = sorted(closure)
        result["script_sha256"] = sha256_file(CODE_DIR / script)
        result["sources_sha256"] = {m: sha256_file(CODE_DIR / m) for m in closure}
        result["runner_sha256"] = sha256_file(Path(__file__).resolve())
        result["closure_unresolved"] = unresolved
        launch = {"commit": git_commit(), "worktree": worktree_state()}
        result["launch"] = launch
        if strict_run:
            # Codex round 8: a certified run must be reconstructible -- a
            # resolvable commit, a clean TRACKED worktree, a fully resolved
            # closure, and no cached bytecode anywhere under code/ (Python
            # loads a stamp-matched .pyc in preference to the hashed .py).
            if unresolved:
                result["failures"].append("closure_unresolved:" + ",".join(unresolved))
            if not _ARCHIVE:
                if len(launch["commit"]) != 40:
                    result["failures"].append("launch_commit_unresolvable")
                if launch["worktree"].get("tracked_dirty") is not False:
                    result["failures"].append("launch_worktree_tracked_dirty")
                for pyc in list(CODE_DIR.rglob("*.pyc")):
                    pyc.unlink()                         # regenerable cache only
                for cache in list(CODE_DIR.rglob("__pycache__")):
                    if cache.is_dir() and not any(cache.iterdir()):
                        cache.rmdir()
            leftover = [str(x.relative_to(CODE_DIR)) for x in CODE_DIR.rglob("*.pyc")]
            launch["bytecode_clean"] = not leftover
            if leftover:
                result["failures"].append("bytecode_present_at_launch")
            if result["failures"]:
                result["status"] = "fail"
                return result                            # refuse to execute
    if execute and command:
        # archive: certificates run bytecode-free, without user site-packages
        argv = [sys.executable if token == "PYTHON" else token
                for token in command]
        if argv and argv[0] == sys.executable:
            argv[1:1] = ["-B", "-s"] if _ARCHIVE else ["-B"]   # never write bytecode
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
        if closure:
            after = {m: sha256_file(CODE_DIR / m) for m in closure
                     if (CODE_DIR / m).is_file()}
            if after != result["sources_sha256"]:
                result["failures"].append("source_changed_during_run")
    elif command:
        result["status"] = "skipped_expensive"

    for artifact in entry.get("artifacts", []):
        path = ROOT / artifact["path"]
        kind = artifact.get("kind")
        if execute and kind == "runner_record":
            # being regenerated by this very execution: the replacement
            # record is written natively afterwards (Codex round 7: the old
            # record must not block its own successor)
            continue
        if not path.exists():
            result["failures"].append(f"artifact_missing:{artifact['path']}")
            continue
        actual = sha256_file(path)
        result["artifact_hashes"][artifact["path"]] = actual
        if artifact.get("sha256") and actual != artifact["sha256"]:
            result["failures"].append(f"hash_mismatch:{artifact['path']}")
        pin = str(artifact.get("sha256", ""))
        if kind == "runner_record" and not (
                pin.startswith("sha256:") and len(pin) == 71
                and all(c in "0123456789abcdef" for c in pin[7:])):
            # Codex round 10: a predeclared record without a pin must not
            # turn the registry green merely by existing
            result["failures"].append(f"runner_record_unpinned:{artifact['path']}")
        # Pinned RUNNER RECORDS are checked semantically, not only bytewise
        # (Codex round 6: a record generated before a claim was rescoped
        # still matched its pin and passed green).
        result["failures"].extend(
            f"stale_run_record:{why}:{artifact['path']}"
            for why in _check_run_record(path, entry, script, result,
                                         strict=(kind == "runner_record")))

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


def _local_import_closure(script: str) -> list:
    """Files under code/ reachable from the certificate through local
    imports (static scan of `import X` / `from X import`), including the
    script itself; sorted, code/-relative (Codex round 7: an unchanged
    wrapper over a changed engine must not pass source agreement)."""
    import ast
    seen, todo, unresolved = [], [script], []
    while todo:
        name = todo.pop()
        if name in seen or not (CODE_DIR / name).is_file():
            continue
        seen.append(name)
        try:
            tree = ast.parse((CODE_DIR / name).read_text(encoding="utf-8"))
        except (SyntaxError, UnicodeDecodeError):
            unresolved.append(f"unparsable:{name}")
            continue
        for node in ast.walk(tree):
            mods = []
            if isinstance(node, ast.Import):
                mods = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                if node.level:                         # relative import
                    unresolved.append(f"relative_import:{name}")
                    continue
                mods = [node.module or ""]
            elif isinstance(node, ast.Call):
                fn = node.func
                fname = (fn.id if isinstance(fn, ast.Name)
                         else fn.attr if isinstance(fn, ast.Attribute) else "")
                if fname in ("__import__", "import_module", "exec_module",
                             "spec_from_file_location", "run_path", "load_source"):
                    unresolved.append(f"dynamic_import:{name}:{fname}")
            for mod in mods:
                top = mod.split(".")[0]
                if top in ("importlib", "runpy"):
                    unresolved.append(f"dynamic_import_module:{name}:{top}")
                if (CODE_DIR / top).is_dir():              # local package
                    unresolved.append(f"package_import:{name}:{top}")
                cand = top + ".py"
                if (CODE_DIR / cand).is_file() and cand not in seen:
                    todo.append(cand)
    return sorted(seen), sorted(set(unresolved))


def _check_run_record(path: Path, entry: dict, script, result: dict,
                      strict: bool = False) -> list:
    """Semantic checks on a pinned runner record (a partial manifest):
    the record must be for this entry, a clean pass, carry the CURRENT
    registry claims, and describe the CURRENT certificate source -- via
    its recorded script hash, or for legacy records via the committed
    bytes at the recorded commit.  Returns a list of problems."""
    try:
        rec = json.loads(path.read_text(encoding="utf-8"))
    except (ValueError, UnicodeDecodeError):
        return ["malformed_runner_record"] if strict else []
    if not (isinstance(rec, dict) and "partial_of" in rec and "results" in rec):
        return ["malformed_runner_record"] if strict else []
    problems = []
    rr = rec.get("results")
    if not (isinstance(rr, list) and len(rr) == 1 and isinstance(rr[0], dict)
            and rr[0].get("id") == entry["id"]):
        return ["record_id"]
    r0 = rr[0]
    summ = rec.get("summary")
    inv = rec.get("invocation")
    if not (isinstance(summ, dict) and isinstance(inv, dict)):
        return ["malformed_runner_record"]
    if (inv.get("only") != entry["id"] or summ.get("total") != 1
            or summ.get("passed") != 1 or summ.get("failed") != 0
            or summ.get("skipped_expensive") != 0 or summ.get("pins_verified", 0) != 0
            or not str(rec.get("status", "")).startswith("passed")):
        problems.append("record_summary_inconsistent")
    if (r0.get("status") != "pass" or r0.get("exit_code") != 0
            or r0.get("executed") is not True or r0.get("failures")):
        problems.append("record_not_a_clean_pass")
    if r0.get("claims") != entry.get("claims", []):
        problems.append("claims_differ_from_registry")
    if r0.get("description") != entry.get("description", ""):
        result.setdefault("notes", []).append(
            f"record_description_drift:{path.name}")   # prose: noted, not fatal
    if script:
        current = result.get("sources_sha256", {})   # launch-time closure hashes
        recorded = r0.get("sources_sha256")
        launch = r0.get("launch")
        if isinstance(recorded, dict) and isinstance(launch, dict):
            if set(recorded) != set(current):
                problems.append("source_closure_differs")
            for mod, h in current.items():
                if recorded.get(mod) != h:
                    problems.append(f"source_changed_since_run:{mod}")
                    break
            # launch attestation (Codex round 8): reconstructible run
            commit = str(launch.get("commit", ""))
            wt = launch.get("worktree") if isinstance(launch.get("worktree"), dict) else {}
            if len(commit) != 40:
                problems.append("launch_commit_unresolvable")
            if wt.get("tracked_dirty") is not False:
                problems.append("launch_worktree_not_clean")
            if launch.get("bytecode_clean") is not True:
                problems.append("launch_bytecode_unattested")
            if r0.get("closure_unresolved"):
                problems.append("launch_closure_unresolved")
            if _ARCHIVE:
                result.setdefault("notes", []).append(
                    f"blob_agreement_unverifiable_in_archive:{path.name}")
            elif len(commit) == 40 and subprocess.run(
                    ["git", "-C", str(ROOT), "cat-file", "-e", commit + "^{commit}"],
                    capture_output=True).returncode != 0:
                # public/foreign-history context: the launch commit lives in
                # the private repository (same epistemic situation as the
                # archive branch above); current-tree hash agreement is
                # still verified below via sources_changed checks.
                result.setdefault("notes", []).append(
                    f"blob_agreement_unverifiable_foreign_history:{path.name}")
            elif len(commit) == 40:
                to_check = dict(recorded)
                if r0.get("runner_sha256"):
                    to_check["run_certifications.py"] = r0["runner_sha256"]
                else:
                    problems.append("launch_runner_unattested")
                for mod, h in sorted(to_check.items()):
                    # declared sources may be CODE_DIR-relative with "..":
                    # git show does not normalize pathspecs (Codex round 15)
                    spec = posixpath.normpath("code/" + mod)
                    shown = subprocess.run(
                        ["git", "-C", str(ROOT), "show", f"{commit}:{spec}"],
                        capture_output=True)
                    if shown.returncode != 0:
                        problems.append(f"launch_source_not_at_commit:{mod}")
                        break
                    if "sha256:" + hashlib.sha256(shown.stdout).hexdigest() != h:
                        problems.append(f"launch_source_differs_from_commit:{mod}")
                        break
        elif recorded:
            problems.append("launch_block_missing")
        else:
            result.setdefault("notes", []).append(
                f"legacy_launch_unverified:{path.name}")
            # legacy record (no launch attestation): the current closure's
            # files must equal their committed bytes at the COMPLETION commit
            # recorded in the manifest -- weaker than launch attestation
            commit = str(rec.get("git_commit", ""))
            if _ARCHIVE or len(commit) != 40:
                result.setdefault("notes", []).append(
                    f"legacy_record_source_unverified:{path.name}")
            elif subprocess.run(
                    ["git", "-C", str(ROOT), "cat-file", "-e", commit + "^{commit}"],
                    capture_output=True).returncode != 0:
                # public/foreign-history context (see the launch branch above)
                result.setdefault("notes", []).append(
                    f"legacy_source_unverifiable_foreign_history:{path.name}")
            else:
                for mod, h in current.items():
                    shown = subprocess.run(
                        ["git", "-C", str(ROOT), "show", f"{commit}:code/{mod}"],
                        capture_output=True)
                    if shown.returncode != 0:
                        problems.append(f"source_not_at_recorded_commit:{mod}")
                        break
                    if "sha256:" + hashlib.sha256(shown.stdout).hexdigest() != h:
                        problems.append(f"source_changed_since_run:{mod}")
                        break
    return problems


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--full", action="store_true",
                        help="also run expensive regenerations")
    parser.add_argument("--only", default=None,
                        help="run a single certification id")
    args = parser.parse_args()

    if not (ROOT / ".git").exists():
        # Archive gate (Codex round, 2026-08-28): in a git-free tree the
        # runner executes nothing unless the WHOLE release tree matches
        # RELEASE_MANIFEST.json.  Before this gate a stubbed certificate
        # PASSed here in 0.1 s; public_lint.py caught it, the runner did not.
        if str(CODE_DIR) not in sys.path:
            sys.path.insert(0, str(CODE_DIR))   # tree verified above; safe now
        from launch_provenance import capture_launch
        try:
            gate = capture_launch(ROOT, Path(__file__).resolve())
        except SystemExit as exc:
            print(f"archive gate REFUSED: {exc}", file=sys.stderr)
            return 2
        print(f"archive gate: {gate['context']} "
              f"({gate.get('release_tree_verified_files')} files verified)")

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
        # native replacement record (Codex round 7): a successful --only
        # execution of an entry with a runner_record artifact writes the
        # new record itself; the registry pin is then updated by hand.
        [entry] = entries
        [res] = results
        for artifact in entry.get("artifacts", []):
            if artifact.get("kind") == "runner_record" and res["status"] == "pass":
                target = ROOT / artifact["path"]
                target.parent.mkdir(parents=True, exist_ok=True)
                tmp = target.with_suffix(".json.tmp")
                tmp.write_text(json.dumps(manifest, indent=2) + "\n",
                               encoding="utf-8")
                tmp.replace(target)
                print(f"runner record written: {artifact['path']} "
                      f"{sha256_file(target)} (update the registry pin)")
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
