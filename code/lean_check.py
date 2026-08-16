"""Kernel-check the Lean certificates in lean/ (with negative control).

Runs the pinned Lean 4.24 + Mathlib toolchain (the lean-mathlib
wrapper with baked LEAN_PATH) over every lean/*.lean file, then — the
negative-control discipline — verifies the checker itself rejects a
deliberately false proof, so a silently broken toolchain cannot report
green. Registered as `lean_certificates` in certifications.json
(expensive: the toolchain is machine-local, not part of the sympy-only
portable core; the entry fails loudly where the toolchain is absent).

Exit 0 with "ALL LEAN CERTIFICATES KERNEL-CHECKED" on success.
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LEAN_DIR = ROOT / "lean"
WRAPPER = Path.home() / "projects" / "lean-benchmarks" / "lean-mathlib"

_NEGATIVE = (
    "import Mathlib.Algebra.Lie.Basic\n"
    "theorem bogus (L : Type) [LieRing L] (x y : L) : "
    "⁅x, y⁆ = 0 := rfl\n"
)


def main() -> int:
    if not WRAPPER.exists():
        print(f"LEAN TOOLCHAIN MISSING: {WRAPPER}")
        print("(machine-local toolchain; see lean/README.md)")
        return 1
    files = sorted(LEAN_DIR.glob("*.lean"))
    if not files:
        print("no lean/*.lean files found")
        return 1
    failures = 0
    for path in files:
        proc = subprocess.run(
            [str(WRAPPER), str(path)], capture_output=True, text=True,
            timeout=1800,
        )
        ok = proc.returncode == 0 and "error" not in proc.stdout.lower()
        print(f"  {'ok  ' if ok else 'FAIL'} {path.relative_to(ROOT)}")
        if not ok:
            failures += 1
            print("\n".join(proc.stdout.splitlines()[-12:]))
    with tempfile.NamedTemporaryFile(
        "w", suffix=".lean", delete=False) as handle:
        handle.write(_NEGATIVE)
        negative_path = handle.name
    try:
        proc = subprocess.run(
            [str(WRAPPER), negative_path], capture_output=True, text=True,
            timeout=1800,
        )
    finally:
        Path(negative_path).unlink(missing_ok=True)
    if proc.returncode == 0:
        print("  FAIL negative control: the checker accepted a false proof")
        failures += 1
    else:
        print("  ok   negative control: false proof rejected")
    if failures:
        print(f"\n{failures} failure(s).")
        return 1
    print(f"\nALL LEAN CERTIFICATES KERNEL-CHECKED ({len(files)} file(s), "
          "negative control passed)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
