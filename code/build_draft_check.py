"""Draft build-and-hash check (certifier upgrade, Codex round 4).

Compiles every draft/*.tex with pdflatex (two passes each, for
references) into a scratch directory — tracked PDFs are never
touched — and prints per-file source hashes plus pass/fail markers:

    <name>.tex sha256:<hash>
    <name> build OK: <pages> pages
    ALL DRAFTS BUILD OK

Fails (nonzero exit) on any LaTeX error, missing pdflatex, or missing
output. Registered as `draft_builds_and_hash`; the certifier records
the printed source hash in its manifest via the run log, deliberately
unpinned — the draft is expected to churn, and the check certifies
"the committed draft compiles", not "the draft never changes".
"""

from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DRAFT_DIR = ROOT / "draft"


def build_one(tex: Path) -> bool:
    digest = hashlib.sha256(tex.read_bytes()).hexdigest()
    print(f"{tex.name} sha256:{digest}")
    with tempfile.TemporaryDirectory(prefix="draft-build-") as scratch:
        # compile a COPY inside the scratch dir with cwd=scratch:
        # TeX's input path falls back to the cwd, so stale .aux
        # residue in draft/ (e.g. from an in-place non-hyperref run)
        # would otherwise poison pass 1 (observed: hyperref's
        # five-field \newlabel vs a stale two-field aux)
        shutil.copy(tex, Path(scratch) / tex.name)
        for attempt in (1, 2):
            proc = subprocess.run(
                ["pdflatex", "-interaction=nonstopmode",
                 "-halt-on-error", tex.name],
                cwd=scratch, capture_output=True, text=True,
                timeout=300,
            )
            if proc.returncode != 0 or "!" == proc.stdout.lstrip()[:1]:
                tail = "\n".join(proc.stdout.splitlines()[-15:])
                print(f"{tex.stem} build FAILED (pass {attempt}):\n{tail}")
                return False
        if not (Path(scratch) / f"{tex.stem}.pdf").exists():
            print(f"{tex.stem} build FAILED: no PDF produced")
            return False
        match = re.search(r"Output written on .*?\((\d+) pages?",
                          proc.stdout)
        pages = match.group(1) if match else "?"
        print(f"{tex.stem} build OK: {pages} pages")
    return True


def main() -> int:
    if shutil.which("pdflatex") is None:
        print("draft build SKIP-FAIL: pdflatex not installed")
        return 1
    sources = sorted(DRAFT_DIR.glob("*.tex"))
    if not sources:
        print("draft build FAILED: no .tex sources found")
        return 1
    if not all([build_one(tex) for tex in sources]):
        return 1
    print(f"ALL DRAFTS BUILD OK ({len(sources)} file(s))")
    return 0


if __name__ == "__main__":
    sys.exit(main())
