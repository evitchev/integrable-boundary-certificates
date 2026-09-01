"""Launch-time provenance capture, shared by the certificate scripts.

Context classification is decided by the presence of .git AT THE
SUPPLIED ROOT (parent repositories never influence it), and every
degraded state aborts rather than recording something hollow
(review, 2026-08-27):

- DEVELOPMENT CLONE (root/.git exists): git MUST work -- a missing git
  executable or any git error aborts; the repository toplevel must be
  the root itself.  Records the launch commit and exact dirty-path
  lists (tracked/untracked split).
- RELEASE ARCHIVE (no root/.git): requires a parseable
  EXPORT_RECORD.json with a 40-hex freeze commit AND agreement between
  the export record's bytes and its hash pinned in
  RELEASE_MANIFEST.json.  Anything else aborts: an unpacked tree that
  is neither a git clone nor a verifiable release is not a supported
  execution context.

The invocation is recorded literally (basenames): actual executable,
argv, and sys.flags.optimize (catches -O, which argv does not show).
The shared helper's own hash is recorded alongside the script's.
"""
import hashlib
import os
import json
import re
import subprocess
import sys
from pathlib import Path


def _invocation():
    return {
        "executable": Path(sys.executable).name,
        "argv": [Path(sys.argv[0]).name] + list(sys.argv[1:]),
        "optimize": sys.flags.optimize,
    }


_GENERATED = {"RELEASE_MANIFEST.json",
              "results/certification_manifest.json",
              "results/certification_manifest.partial.json"}


def verify_release_tree(root, manifest):
    """Every file pinned in RELEASE_MANIFEST.json must exist and match its
    hash, no unmanifested file or directory may be present, symlinks and
    other non-regular entries are refused, and any bytecode (__pycache__/,
    *.pyc) is REFUSED (generated manifests and *.tmp excluded).  Returns a
    list of problems; empty
    means the whole tree is the released one.  Codex round 2026-08-28:
    checking only EXPORT_RECORD.json let a stubbed certificate PASS in
    0.1 s under 'archive-verified-release'."""
    root = Path(root)
    files = manifest.get("files", {})
    problems = []
    for rel, pinned in sorted(files.items()):
        p = root / rel
        if not p.is_file():
            problems.append(f"missing: {rel}")
            continue
        actual = "sha256:" + hashlib.sha256(p.read_bytes()).hexdigest()
        if actual != pinned:
            problems.append(f"modified: {rel}")
    expected_dirs = {""}
    for rel in files:
        parts = rel.split("/")[:-1]
        for i in range(1, len(parts) + 1):
            expected_dirs.add("/".join(parts[:i]))
    for p in sorted(root.rglob("*")):
        if ".git" in p.parts:
            continue
        rel = p.relative_to(root).as_posix()
        if p.is_symlink():                       # Codex round 4: whole-tree
            problems.append(f"symlink: {rel}")   # means no non-regular entries
            continue
        if p.is_dir():
            if p.name == "__pycache__":
                problems.append(f"bytecode: {rel}/")
            elif rel not in expected_dirs:
                problems.append(f"unmanifested directory: {rel}/")
            continue
        if not p.is_file():
            problems.append(f"non-regular: {rel}")
            continue
        if "__pycache__" in p.parts or p.suffix == ".pyc":
            # Codex round 2 (2026-08-28), reproduced: a sourceless .pyc in
            # code/ shadows even stdlib modules, and a __pycache__ pyc with
            # a matching source stamp is imported OVER the verified .py.
            # Bytecode is never part of a release: refuse it outright.
            problems.append(f"bytecode: {rel}")
            continue
        if p.suffix == ".tmp":
            continue
        if rel in _GENERATED or rel in files:
            continue
        problems.append(f"unmanifested: {rel}")
    return problems


def capture_launch(root, script_path):
    root = Path(root).resolve()
    base = {
        "script_sha256_at_start": hashlib.sha256(
            Path(script_path).read_bytes()).hexdigest(),
        "helper_sha256": hashlib.sha256(
            Path(__file__).read_bytes()).hexdigest(),
        "invocation": _invocation(),
    }
    if (root / ".git").exists():
        try:
            rp = subprocess.run(
                ["git", "-C", str(root), "rev-parse", "HEAD"],
                capture_output=True, text=True)
            tl = subprocess.run(
                ["git", "-C", str(root), "rev-parse", "--show-toplevel"],
                capture_output=True, text=True)
            st = subprocess.run(
                ["git", "-C", str(root), "status", "--porcelain"],
                capture_output=True, text=True)
        except FileNotFoundError:
            raise SystemExit(
                "launch provenance: repository metadata present but no "
                "git executable; refusing to run without provenance")
        if rp.returncode != 0 or st.returncode != 0 or tl.returncode != 0:
            raise SystemExit(
                "launch provenance capture failed (git error in a "
                "repository); refusing to run without provenance")
        if Path(tl.stdout.strip()).resolve() != root:
            raise SystemExit(
                "launch provenance: git resolves a DIFFERENT toplevel "
                f"({tl.stdout.strip()}); refusing a nested/ambiguous tree")
        lines = [l for l in st.stdout.splitlines() if l]
        base.update({
            "context": "development-clone",
            "launch_commit_at_start": rp.stdout.strip(),
            "worktree_dirty_at_start": bool(lines),
            "tracked_dirty_at_start": [l[3:] for l in lines
                                       if not l.startswith("??")],
            "untracked_at_start": [l[3:] for l in lines
                                   if l.startswith("??")],
        })
        return base
    # archive candidate: STRICT validation, fail-closed
    rec = root / "EXPORT_RECORD.json"
    man = root / "RELEASE_MANIFEST.json"
    if not (rec.exists() and man.exists()):
        raise SystemExit(
            "launch provenance: no repository metadata and no verifiable "
            "release records (EXPORT_RECORD.json + RELEASE_MANIFEST.json "
            "required); refusing an unclassifiable tree")
    try:
        rec_bytes = rec.read_bytes()
        rdata = json.loads(rec_bytes)
        mdata = json.loads(man.read_text())
    except Exception as exc:
        raise SystemExit(f"launch provenance: unreadable release "
                         f"records ({exc!r}); refusing")
    freeze = rdata.get("private_freeze_commit", rdata.get("freeze_commit"))
    if not (isinstance(freeze, str) and re.fullmatch(r"[0-9a-f]{40}", freeze)):
        raise SystemExit("launch provenance: release freeze commit "
                         f"missing or malformed ({freeze!r}); refusing")
    pinned = mdata.get("files", {}).get("EXPORT_RECORD.json")
    actual = "sha256:" + hashlib.sha256(rec_bytes).hexdigest()
    if pinned != actual:
        raise SystemExit(
            "launch provenance: EXPORT_RECORD.json does not match its "
            f"RELEASE_MANIFEST pin ({pinned} != {actual}); refusing a "
            "tampered or partial archive")
    # Codex round 3 (2026-08-28, reproduced): verification here is LATE --
    # the calling script and its imports have already executed.  In a
    # git-free tree the only fail-closed entry point is the runner's
    # isolated-mode bootstrap, which verifies the tree before any import
    # and leaves this marker (the manifest hash, re-checked here).
    # The marker is a COORDINATION FLAG, not an attestation: any process
    # can set it (Codex round 4).  Its job is to make direct execution --
    # where this verification is necessarily late -- refuse; the tree is
    # re-verified here regardless, and the runner's manifest is the record
    # of what ran.  The launch block says so explicitly.
    marker = os.environ.get("IB_TRUSTED_BOOTSTRAP", "")
    if marker != hashlib.sha256(man.read_bytes()).hexdigest():
        raise SystemExit(
            "launch provenance: git-free tree without the trusted bootstrap "
            "marker; run through run_certifications.py (isolated mode, tree "
            "verified before any import) -- direct execution cannot be "
            "fail-closed against modules that were already imported")
    problems = verify_release_tree(root, mdata)
    if problems:
        shown = "\n  ".join(problems[:12])
        more = "" if len(problems) <= 12 else f"\n  ... {len(problems) - 12} more"
        raise SystemExit(
            "launch provenance: release tree does not match "
            "RELEASE_MANIFEST.json; refusing to execute in a tampered or "
            f"partial archive:\n  {shown}{more}")
    base.update({"context": "archive-verified-release",
                 "release_freeze_commit": freeze,
                 "release_tree_verified_files": len(mdata.get("files", {})),
                 "bootstrap_marker": "present (coordination flag set by "
                 "run_certifications after its isolated-mode verification; "
                 "not an attestation -- the tree was re-verified here)"})
    return base
