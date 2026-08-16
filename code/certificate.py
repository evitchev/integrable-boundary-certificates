"""Programmatic certificate provenance — never transcribe, always resolve.

Born from the 2026-08-08 review catch (commit b1b6b9f): a certificate's
launch_commit had been hand-typed in a heredoc with a fabricated tail.
The recorded lesson — "provenance fields must be resolved
programmatically, never transcribed" — is this module. Certificate
writers import `stamp()` (or run the CLI) instead of typing SHAs,
dates, or environment strings by hand; `certificate_lint.py` checks
the results.

Library:

    from certificate import stamp, declare_launch
    payload = stamp({...})            # adds the provenance block
    declare_launch("../results/foo_launch.json",
                   tasks=[...], acceptance="all kernels empty")

CLI (stamps an existing JSON file in place, refusing to overwrite an
existing provenance block unless --force):

    python3 certificate.py stamp ../results/foo_certificate.json

The provenance block:

    "provenance": {
      "launch_commit": <full HEAD SHA, resolved via git rev-parse>,
      "worktree_dirty": <bool>,
      "dirty_paths": [...],          # capped, so a dirty launch is visible
      "stamped_at_utc": <ISO 8601>,
      "python": ..., "sympy": ..., "platform": ...
    }

A dirty worktree does not block stamping (mid-session launches are
real), but it is recorded — a certificate claiming to be reproducible
from `launch_commit` while the tree was dirty says so on its face.
"""

from __future__ import annotations

import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
_DIRTY_PATH_CAP = 40


def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(ROOT), *args],
                          capture_output=True, text=True)


def resolve_head_commit() -> str:
    """The full current commit SHA, verified to exist — never typed."""

    result = _git("rev-parse", "--verify", "HEAD^{commit}")
    if result.returncode != 0:
        raise RuntimeError(f"cannot resolve HEAD: {result.stderr.strip()}")
    sha = result.stdout.strip()
    if len(sha) != 40:
        raise RuntimeError(f"unexpected rev-parse output: {sha!r}")
    return sha


def worktree_state() -> dict:
    status = _git("status", "--porcelain")
    paths = [line[3:] for line in status.stdout.splitlines() if line]
    return {
        "worktree_dirty": bool(paths),
        "dirty_paths": paths[:_DIRTY_PATH_CAP],
        "dirty_path_count": len(paths),
    }


def provenance_block() -> dict:
    try:
        import sympy
        sympy_version = sympy.__version__
    except ImportError:
        sympy_version = "MISSING"
    block = {
        "launch_commit": resolve_head_commit(),
        "stamped_at_utc": datetime.now(timezone.utc).isoformat(
            timespec="seconds"),
        "python": sys.version.split()[0],
        "sympy": sympy_version,
        "platform": platform.platform(),
    }
    block.update(worktree_state())
    return block


def stamp(payload: dict, *, force: bool = False) -> dict:
    """Return payload with a resolved provenance block attached."""

    if "provenance" in payload and not force:
        raise ValueError(
            "payload already carries a provenance block; pass force=True "
            "to restamp (the old block is replaced, not merged)"
        )
    stamped = dict(payload)
    stamped["provenance"] = provenance_block()
    return stamped


def declare_launch(path: str | Path, *, tasks, acceptance: str,
                   notes: str = "") -> Path:
    """Write a pre-run launch declaration for a detached computation.

    Declaring the task list and acceptance criterion BEFORE the run
    (and stamping the declaration) is the lightweight analogue of a
    preregistration: the eventual certificate can point at the
    declaration instead of describing its own goalposts after the
    fact.
    """

    declaration = stamp({
        "kind": "launch_declaration",
        "tasks": list(tasks),
        "acceptance": acceptance,
        "notes": notes,
    })
    dest = Path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    staged = dest.with_name(dest.name + ".tmp")
    staged.write_text(json.dumps(declaration, indent=2, sort_keys=True)
                      + "\n", encoding="utf-8")
    staged.replace(dest)
    return dest


def _cli() -> int:
    args = sys.argv[1:]
    force = "--force" in args
    args = [a for a in args if a != "--force"]
    if len(args) != 2 or args[0] != "stamp":
        print(__doc__.split("CLI", 1)[1].split("\n\n")[0], file=sys.stderr)
        print("usage: certificate.py stamp [--force] <file.json>",
              file=sys.stderr)
        return 2
    target = Path(args[1])
    payload = json.loads(target.read_text(encoding="utf-8"))
    try:
        stamped = stamp(payload, force=force)
    except ValueError as exc:
        print(f"refusing: {exc}", file=sys.stderr)
        return 1
    staged = target.with_name(target.name + ".tmp")
    staged.write_text(json.dumps(stamped, indent=2, sort_keys=True) + "\n",
                      encoding="utf-8")
    staged.replace(target)
    print(f"stamped {target} at "
          f"{stamped['provenance']['launch_commit'][:12]} "
          f"(dirty={stamped['provenance']['worktree_dirty']})")
    return 0


if __name__ == "__main__":
    sys.exit(_cli())
