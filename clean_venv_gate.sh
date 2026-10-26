#!/usr/bin/env bash
# Clean-environment reproduction gate (new in v1.2.7).
#
# Builds a fresh virtualenv containing ONLY the packages declared in the root
# requirements.txt (sympy>=1.12 and mpmath>=1.3;
# numpy and python-flint are commented as optional), runs every quick
# certification of the public repository, the full opers quick suite with its
# integrity gate and drill rounds, and the super-KdV checks, and exits 0 only
# if every line passes.
#
# Usage: bash clean_venv_gate.sh
#   (run from anywhere; paths are resolved relative to this file's directory)
set -eu
ROOT="$(cd "$(dirname "$0")" && pwd)"
PY="${PYTHON:-python3}"
LOG="$ROOT/logs/clean_venv_gate.log"
mkdir -p "$ROOT/logs"
{
echo "=== clean_venv_gate $(date -Iseconds) ==="
echo "--- declared packages ---"
cat "$ROOT/requirements.txt"
# Two environments, exactly as the root requirements.txt declares them (P4):
#   core  = sympy + mpmath              (required; everything except the two vir1r checks)
#   opers = core + numpy                (numpy optional at the root; the two vir1r
#                                        explicit-component rank checks need it)
# The opers quick suite and its drills must run on 'opers' (sympy+mpmath+numpy) or the
# two vir1r checks fail; the root certifications and the super-KdV checks run on 'core'.
VENV_CORE="$ROOT/.venv-clean-gate-core"
VENV_OPERS="$ROOT/.venv-clean-gate-opers"
rm -rf "$VENV_CORE" "$VENV_OPERS"
mk_venv() { # mk_venv <dir> [extra packages]
    local d="$1"; shift
    echo "--- creating fresh venv: $d ---"
    "$PY" -m venv "$d" || exit 1
    "$d/bin/pip" install --quiet --no-input -r "$ROOT/requirements.txt" "$@"
    "$d/bin/pip" list | grep -Ei '^(sympy|mpmath|numpy)'
}
echo "--- environment 1: core (declared required packages only) ---"
mk_venv "$VENV_CORE"
P=0

echo "--- root certifications (default quick set), core env ---"
(cd "$ROOT/code" && "$VENV_CORE/bin/python" -B run_certifications.py) || P=1

echo "--- super-KdV checks, core env ---"
(cd "$ROOT/super-kdv-2026-10" && "$VENV_CORE/bin/python" -I -B run_checks.py) || P=1

echo "--- environment 2: opers (core + numpy, the declared optional dependency of the two vir1r checks) ---"
mk_venv "$VENV_OPERS" numpy

echo "--- opers quick suite (integrity gate included), opers env ---"
(cd "$ROOT/opers-2026-10" && "$VENV_OPERS/bin/python" -I -B run_all.py --quick) || P=1

echo "--- opers drills: every changed/new check of v1.2.7, opers env ---"
(cd "$ROOT/opers-2026-10" && "$VENV_OPERS/bin/python" -I -B run_all.py \
    --drill --only kit_integrity_manifest,b2_level2_w8_agreement_z,b2_level2_w8_agreement_k) || P=1

if [ "$P" -eq 0 ]; then
    echo "CLEAN-VENV GATE: ALL PASS ($(date -Iseconds))"
else
    echo "CLEAN-VENV GATE: FAILURES (see lines above)"
fi
exit $P
} 2>&1 | tee "$LOG"
