#!/bin/bash
# run.sh NEWDIR [STEP ...]: re-runs the computations of this folder and writes their outputs into NEWDIR, a directory that
# must not exist yet (the shipped records in out/ and logs/ are never written).  Then compare.sh NEWDIR compares them.
#
# Steps, by default all, in this order (the order of the original runs):
#   g1          code/run_g1.py: the classical multiplicities at weight 22            ulimit -v 8000000,  timeout 600  (seconds)
#   g2          the positive control, sheet 3 at t = 13/11: fast_kernel3.py 22 3 -119/2 -143/2
#                                                                                      ulimit -v 32000000, timeout 5400 (~23 min)
#   compare_g2  code/compare_g2.sh: g2's new record against the repository's results/lab/fast_kernel3/w22_t1311_sheet3.json
#               (needs g2's record in NEWDIR)                                          ulimit -v 8000000,  timeout 600  (seconds)
#   s3          the null cell on sheet 3 at t = 11/5: fast_kernel3.py 22 3 -21/4 -55/4  ulimit -v 32000000, timeout 5400 (~23 min)
#   s2          the null cell on sheet 2 at t = 2: fast_kernel3.py 22 2 -7 -16          ulimit -v 32000000, timeout 5400 (~23 min)
# Each kernel run uses --primes=12 --workers=4 --solvers=4, as originally run: at most 4 busy processes at a time, and about
# 60 GiB of memory at its peak (4 parallel eliminations of about 15 GiB address space each).  Before each kernel run the script
# stops if MemAvailable is below MIN_MEM_GIB (default 80, the floor of the original runs; MIN_MEM_GIB=0 skips the check).
#
# The interpreter is $PYTHON (default python3), run as $PYTHON -I -B; it needs sympy and python-flint.  compare_g2 needs jq.
# This folder must sit at results/lab/jumping_loci_records/spin21_null_cells/ of the repository (the root is ../../../..).
set -euo pipefail
PYTHON=${PYTHON:-python3}
MIN_MEM_GIB=${MIN_MEM_GIB:-80}
HERE=$(cd "$(dirname "$0")" && pwd -P)
ROOT=$(cd "$HERE/../../../.." && pwd)
[ $# -ge 1 ] && [ -n "$1" ] || { echo "usage: run.sh NEWDIR [g1 g2 compare_g2 s3 s2]"; exit 2; }
OUT=$1; shift
[ -e "$OUT" ] && { echo "run.sh: $OUT exists; give a new directory: nothing run"; exit 2; }
case "$(realpath -m -- "$OUT")/" in
  "$HERE/out/"*|"$HERE/logs/"*|"$HERE/code/"*) echo "run.sh: $OUT is inside the shipped out/, logs/ or code/: nothing run"; exit 2;;
esac
STEPS=("$@"); [ ${#STEPS[@]} -gt 0 ] || STEPS=(g1 g2 compare_g2 s3 s2)
for st in "${STEPS[@]}"; do
  case "$st" in g1|g2|compare_g2|s3|s2) ;; *) echo "run.sh: unknown step $st: nothing run"; exit 2;; esac
done

# the code files must be the shipped ones, and the repository files the scripts read must be there
grep -v '^#' "$HERE/code/SOURCE_SHA256.txt" | awk '{print $1 "  " $2}' | (cd "$HERE/code" && sha256sum --check --strict --quiet -) \
  || { echo "run.sh: a code file differs from code/SOURCE_SHA256.txt: nothing run"; exit 2; }
for f in code/mixed_engine.py code/invariant_engine.py results/lab/fast_kernel3/w22_t1311_sheet3.json; do
  [ -f "$ROOT/$f" ] || { echo "run.sh: $ROOT/$f is missing (is this folder at results/lab/jumping_loci_records/spin21_null_cells/?): nothing run"; exit 2; }
done

mkdir -p "$OUT/out" "$OUT/logs"
OUT=$(cd "$OUT" && pwd)
META="$OUT/run_meta.txt"
{
  echo "run.sh (sha256 $(sha256sum "$HERE/run.sh" | cut -d' ' -f1)) started $(date -u +%Y-%m-%dT%H:%M:%SZ); steps: ${STEPS[*]}"
  echo "uname: $(uname -srm)"
  echo "python: $PYTHON = $(command -v "$PYTHON")"
  "$PYTHON" -I -B -c 'import sys, sympy, mpmath, flint; print("Python", sys.version.split()[0], "| sympy", sympy.__version__, "| mpmath", mpmath.__version__, "| python-flint", flint.__version__)'
  echo "jq: $(jq --version 2>&1 || true)"
} | tee "$META"

memcheck() {
  local mg
  mg=$(awk '/^MemAvailable/ {printf "%d", $2/1048576}' /proc/meminfo 2>/dev/null || true)
  echo "  MemAvailable ${mg:-unknown} GiB (floor $MIN_MEM_GIB GiB)" | tee -a "$META"
  [ "$MIN_MEM_GIB" = 0 ] || [ "${mg:-0}" -ge "$MIN_MEM_GIB" ] || { echo "run.sh: MemAvailable below $MIN_MEM_GIB GiB: stopped" | tee -a "$META"; exit 3; }
}

kernel() {  # kernel NAME SOL N S: one exact kernel computation at weight 22; record out/w22_NAME.json, log logs/w22_NAME.log
  local name=$1 sol=$2 n=$3 s=$4 rc=0
  local rec="$OUT/out/w22_$name.json" log="$OUT/logs/w22_$name.log"
  echo "[$name] $(date -u +%Y-%m-%dT%H:%M:%SZ) fast_kernel3.py 22 $sol $n $s --primes=12 --workers=4 --solvers=4 (ulimit -v 32000000, timeout 5400)" | tee -a "$META"
  memcheck
  ( set -o pipefail
    ulimit -v 32000000 && cd "$HERE" && FASTK_RECORD="$rec" IB_CODE="$HERE/code/vendor/code" timeout 5400 \
      "$PYTHON" -I -B "$HERE/code/vendor/lab/fast_kernel3.py" 22 "$sol" "$n" "$s" --primes=12 --workers=4 --solvers=4 2>&1 | tee "$log" ) || rc=$?
  echo "[$name] $(date -u +%Y-%m-%dT%H:%M:%SZ) exit status $rc" | tee -a "$META"
  [ $rc = 0 ] || exit $rc
  [ -f "$rec" ] || { echo "[$name] no record written: stopped" | tee -a "$META"; exit 4; }
  echo "[$name] record sha256 $(sha256sum "$rec" | cut -d' ' -f1), $(stat -c %s "$rec") bytes" | tee -a "$META"
}

small() {  # small NAME COMMAND...: a light step under ulimit -v 8000000 and timeout 600; log logs/NAME.log
  local name=$1 rc=0; shift
  echo "[$name] $(date -u +%Y-%m-%dT%H:%M:%SZ) $* (ulimit -v 8000000, timeout 600)" | tee -a "$META"
  ( set -o pipefail
    ulimit -v 8000000 && cd "$HERE" && timeout 600 "$@" 2>&1 | tee "$OUT/logs/$name.log" ) || rc=$?
  echo "[$name] $(date -u +%Y-%m-%dT%H:%M:%SZ) exit status $rc" | tee -a "$META"
  [ $rc = 0 ] || exit $rc
}

for st in "${STEPS[@]}"; do
  case "$st" in
    g1)         small run_g1 "$PYTHON" -I -B "$HERE/code/run_g1.py";;
    g2)         kernel t1311_sheet3 3 -119/2 -143/2;;
    compare_g2) [ -f "$OUT/out/w22_t1311_sheet3.json" ] || { echo "run.sh: compare_g2 needs step g2's record in $OUT/out/: stopped" | tee -a "$META"; exit 2; }
                small compare_g2 bash "$HERE/code/compare_g2.sh" "$OUT/out/w22_t1311_sheet3.json";;
    s3)         kernel t115_sheet3 3 -21/4 -55/4;;
    s2)         kernel t2_sheet2 2 -7 -16;;
  esac
done
echo "run.sh done $(date -u +%Y-%m-%dT%H:%M:%SZ); next: $HERE/compare.sh $OUT" | tee -a "$META"
