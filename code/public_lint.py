"""Public-repository provenance lint (Decision 2's two-tier model).

In the public certificate repository, private-development commit SHAs
inside artifacts are OPAQUE ARCHIVAL REFERENCES: never mutated, never
publicly resolvable. Verification is therefore content-based:

1. every file listed in RELEASE_MANIFEST.json exists with the exact
   SHA-256 the manifest records (and no manifest entry is missing);
2. artifact relationships hold: runner_launch_commit fields agree
   with their referenced launch declarations, and merged-record
   source_launch_commits agree with every evidence file they cite —
   the same semantic checks as the private lint, minus git
   resolution;
3. the privacy screen passes: none of the deny patterns recorded in
   EXPORT_RECORD.json appear anywhere — a reader can re-verify the
   export's privacy attestation, not merely trust it.

Exit 0 with "PUBLIC PROVENANCE LINT PASSED" on success; nonzero with
diagnostics otherwise. This file ships in the public repository and
replaces certificate_lint.py's registry entry there (that lint
requires the private history and is intentionally not runnable in a
fresh clone).
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
_SHA_RE = re.compile(r"^[0-9a-f]{7,40}$")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return "sha256:" + digest.hexdigest()


def check_manifest(failures: list[str]) -> int:
    manifest_path = ROOT / "RELEASE_MANIFEST.json"
    if not manifest_path.exists():
        failures.append("RELEASE_MANIFEST.json missing")
        return 0
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    checked = 0
    for rel, pinned in sorted(manifest.get("files", {}).items()):
        path = ROOT / rel
        if not path.exists():
            failures.append(f"manifest file missing: {rel}")
            continue
        checked += 1
        actual = sha256_file(path)
        if actual != pinned:
            failures.append(
                f"hash mismatch: {rel} is {actual[:20]}…, manifest "
                f"says {pinned[:20]}…")
    listed = set(manifest.get("files", {}))
    # Root-RELATIVE paths: a basename match would let an unmanifested
    # file hide anywhere under a whitelisted name (v1.1.1 audit fix;
    # demonstrated with evil/RELEASE_MANIFEST.json against v1.1.0).
    generated = {"RELEASE_MANIFEST.json",
                 "results/certification_manifest.json",
                 "results/certification_manifest.partial.json"}
    for path in sorted(ROOT.rglob("*")):
        if (not path.is_file() or ".git" in path.parts
                or "__pycache__" in path.parts or path.suffix == ".pyc"):
            continue
        rel = path.relative_to(ROOT).as_posix()
        if rel in generated:
            continue
        if rel not in listed:
            failures.append(f"unmanifested file in release: {rel}")
    return checked


def _semantic_checks(payload, name: str, failures: list[str],
                     agreements: list[int]) -> None:
    prov = payload.get("provenance", {})
    rlc = prov.get("runner_launch_commit")
    decl_ref = payload.get("launch_declaration")
    if "launch_declaration" in payload and (
            not isinstance(decl_ref, str) or not decl_ref):
        failures.append(f"{name}: launch_declaration present but "
                        "empty or null (fail-closed)")
    if decl_ref and (not isinstance(rlc, str) or not rlc):
        failures.append(f"{name}: launch_declaration present but "
                        "runner_launch_commit missing or null "
                        "(fail-closed)")
    if rlc and decl_ref:
        decl_path = ROOT / decl_ref
        if not decl_path.exists():
            failures.append(f"{name}: launch_declaration {decl_ref} missing")
        else:
            decl = json.loads(decl_path.read_text(encoding="utf-8"))
            want = decl.get("provenance", {}).get("launch_commit")
            if not (isinstance(want, str) and _SHA_RE.match(want)):
                failures.append(
                    f"{name}: declaration {decl_ref} lacks a valid "
                    "launch_commit (fail-closed)")
            elif rlc != want:
                failures.append(
                    f"{name}: runner_launch_commit {rlc[:12]} != "
                    f"declaration {want[:12]}")
            else:
                agreements[0] += 1
    slc = payload.get("_provenance", {}).get("source_launch_commits")
    if slc is None:
        return
    if not isinstance(slc, list) or not slc:
        failures.append(
            f"{name}: source_launch_commits not a nonempty list "
            "(fail-closed)")
        return
    for entry in slc:
        if not isinstance(entry, dict):
            failures.append(
                f"{name}: source_launch_commits entry not a mapping")
            continue
        lc = entry.get("launch_commit", "")
        refs = re.findall(r"results/[\w.]+\.json", str(entry.get("see", "")))
        if not (isinstance(lc, str) and _SHA_RE.match(lc)) or not refs:
            failures.append(
                f"{name}: malformed source_launch_commits entry")
            continue
        for ref in refs:
            ref_path = ROOT / ref
            if not ref_path.exists():
                failures.append(f"{name}: evidence {ref} missing")
                continue
            evidence = json.loads(ref_path.read_text(encoding="utf-8"))
            candidate = (evidence.get("launch_commit")
                         or evidence.get("provenance", {})
                         .get("launch_commit"))
            if candidate != lc:
                failures.append(
                    f"{name}: source {lc[:12]} disagrees with {ref} "
                    f"({str(candidate)[:12]})")
            else:
                agreements[0] += 1


def check_privacy(failures: list[str]) -> int:
    record_path = ROOT / "EXPORT_RECORD.json"
    if not record_path.exists():
        failures.append("EXPORT_RECORD.json missing")
        return 0
    record = json.loads(record_path.read_text(encoding="utf-8"))
    patterns = [re.compile(p) for p in
                record.get("privacy_screen", {}).get("deny_patterns", [])]
    if not patterns:
        failures.append(
            "EXPORT_RECORD.json carries no deny patterns (fail-closed)")
        return 0
    screened = 0
    skip = {"EXPORT_RECORD.json", "certification_manifest.json"}
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or ".git" in path.parts:
            continue
        # Run-generated files (bytecode, the per-run manifest) embed the
        # local user's absolute paths by nature; the screen certifies
        # shipped content, not run residue.
        if ("__pycache__" in path.parts or path.suffix == ".pyc"
                or path.name in skip):
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        screened += 1
        for pattern in patterns:
            if pattern.search(text):
                failures.append(
                    f"privacy screen hit: {path.relative_to(ROOT)} "
                    f"matches {pattern.pattern!r}")
    return screened


def main() -> int:
    failures: list[str] = []
    agreements = [0]
    hashed = check_manifest(failures)
    semantic_files = 0
    for path in sorted((ROOT / "results").glob("*.json")):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError) as exc:
            failures.append(f"{path.name}: unreadable ({exc})")
            continue
        semantic_files += 1
        _semantic_checks(payload, path.name, failures, agreements)
    screened = check_privacy(failures)
    if failures:
        print()
        for failure in failures:
            print(f"  FAIL  {failure}")
        print(f"\n{len(failures)} failure(s).")
        return 1
    print(f"PUBLIC PROVENANCE LINT PASSED ({hashed} file hashes, "
          f"{agreements[0]} semantic agreement(s) across "
          f"{semantic_files} artifacts, {screened} files "
          "privacy-screened)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
