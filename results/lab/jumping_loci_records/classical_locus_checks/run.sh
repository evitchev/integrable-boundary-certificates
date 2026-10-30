#!/bin/bash
# run.sh OUTPUT_DIR -- re-runs the five scripts of this folder (code/run_a2.py, run_a3.py, run_a4.py, run_a5.py and
# open_tool_crosscheck.py) from the release tree. OUTPUT_DIR must not exist yet; it receives logs/*.log (stdout and stderr of
# each script), LEDGER.json and LEDGER.md (written by run_a5.py), and run_info.txt (times and exit codes; not compared). The records shipped in this
# folder are never written. Interpreter: $PYTHON (default python3), run as "$PYTHON -I -B". Each step runs under
# ulimit -v 8000000 (kB) and its own timeout. This folder must sit at results/lab/jumping_loci_records/classical_locus_checks/
# in a clone of the repository (the repository root is ../../../.. from here); the scripts read repository files from there.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
PYTHON="${PYTHON:-python3}"
MEM_KB=8000000

if [ "$#" -ne 1 ] || [ -z "$1" ]; then
  echo "usage: run.sh OUTPUT_DIR (a directory that does not exist yet)" >&2; exit 2
fi
if [ -e "$1" ]; then
  echo "run.sh: $1 exists; give a new output directory (the shipped records are never overwritten)" >&2; exit 2
fi
for f in code/mixed_engine.py code/invariant_engine.py code/cyl_mo2_profiles.py code/certifications.json; do
  if [ ! -f "$REPO/$f" ]; then
    echo "run.sh: $f not found at the repository root ../../../..; place this folder at" \
         "results/lab/jumping_loci_records/classical_locus_checks/ in a clone of the repository" >&2; exit 2
  fi
done

mkdir -p "$1/logs"
OUT="$(cd "$1" && pwd)"
INFO="$OUT/run_info.txt"
echo "interpreter: $("$PYTHON" -I -B -c 'import sys, sympy, mpmath; print("Python", sys.version.split()[0], "sympy", sympy.__version__, "mpmath", mpmath.__version__)')" | tee "$INFO"

step() {   # step NAME TIMEOUT_S [ARGS...]: run code/NAME.py, stdout and stderr to OUTPUT_DIR/logs/NAME.log
  local name="$1" to="$2"; shift 2
  local t0 rc=0
  t0=$(date +%s)
  echo "$name: start $(date -Is), timeout ${to} s, ulimit -v $MEM_KB" | tee -a "$INFO"
  ( ulimit -v "$MEM_KB"; timeout "$to" "$PYTHON" -I -B "$HERE/code/$name.py" "$@" ) > "$OUT/logs/$name.log" 2>&1 || rc=$?
  echo "$name: exit $rc after $(( $(date +%s) - t0 )) s" | tee -a "$INFO"
  if [ "$rc" -ne 0 ]; then
    echo "run.sh: $name failed (exit $rc); see $OUT/logs/$name.log" >&2; exit "$rc"
  fi
}

step run_a2 1800
step run_a3 1800
step run_a4 900
step run_a5 1800 "$OUT"
step open_tool_crosscheck 600 "$OUT/LEDGER.json"
echo "run.sh: all five scripts exit 0; compare with: bash compare.sh $OUT" | tee -a "$INFO"
