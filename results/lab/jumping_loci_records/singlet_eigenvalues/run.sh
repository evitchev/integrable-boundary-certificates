#!/usr/bin/env bash
# Re-runs the two scripts of this folder (code/b5_fusion.py, code/b5_vacuity.py) and writes their logs
# into a NEW directory given as $1; the shipped logs in logs/ are never touched.
#
# Usage, from the repository root:
#   bash results/lab/jumping_loci_records/singlet_eigenvalues/run.sh <new-output-dir>
# then
#   bash results/lab/jumping_loci_records/singlet_eigenvalues/compare.sh <new-output-dir>
#
# Interpreter: $PYTHON (default python3), with sympy and mpmath (the repository's requirements.txt).
# Each step runs under ulimit -v 8000000 (KiB) and timeout 900 (seconds), as the archived run did.
# Each log is the step's stdout and stderr followed by a line "exit <status>", the format of the shipped logs.
# The output directory must not exist yet, and must lie outside this folder and outside the repository.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../../.." && pwd)"   # the repository root when this folder is results/lab/jumping_loci_records/singlet_eigenvalues/
PYTHON="${PYTHON:-python3}"
MEM_KB=8000000
TIME_S=900
STEPS="b5_fusion b5_vacuity"

if [ "$#" -ne 1 ] || [ -z "$1" ]; then
    echo "usage: $0 <new-output-dir>" >&2; exit 2
fi
if [ -e "$1" ]; then
    echo "refusing: $1 already exists; give a new directory" >&2; exit 2
fi
PARENT="$(dirname "$1")"
if [ ! -d "$PARENT" ]; then
    echo "refusing: the parent directory $PARENT does not exist" >&2; exit 2
fi
OUT="$(cd "$PARENT" && pwd)/$(basename "$1")"
case "$OUT" in
    "$HERE"|"$HERE"/*) echo "refusing: the output directory must lie outside this folder: $OUT" >&2; exit 2 ;;
esac
if [ -f "$ROOT/RELEASE_MANIFEST.json" ]; then
    case "$OUT" in
        "$ROOT"|"$ROOT"/*) echo "refusing: the output directory must lie outside the repository: $OUT" >&2; exit 2 ;;
    esac
fi
command -v timeout >/dev/null 2>&1 || { echo "timeout (GNU coreutils) is required" >&2; exit 2; }
command -v "$PYTHON" >/dev/null 2>&1 || { echo "interpreter not found: $PYTHON" >&2; exit 2; }

mkdir "$OUT"
cd "$OUT"   # nothing is written next to the scripts
echo "interpreter: $("$PYTHON" -I -B -c 'import sys; print(sys.version.split()[0])'); output: $OUT"

STATUS=0
for name in $STEPS; do
    rc=0
    start=$SECONDS
    ( ulimit -v "$MEM_KB"; exec timeout "$TIME_S" "$PYTHON" -I -B "$HERE/code/$name.py" ) > "$OUT/$name.log" 2>&1 || rc=$?
    echo "exit $rc" >> "$OUT/$name.log"
    echo "$name: exit $rc, $((SECONDS - start)) s, log $OUT/$name.log"
    [ "$rc" -eq 0 ] || STATUS=1
done
exit "$STATUS"
