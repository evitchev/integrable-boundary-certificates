#!/usr/bin/env bash
# Re-run the multiplicity stratification at weights 6 and 8 from the release tree.
#
# Usage, from the repository root:
#   bash results/lab/jumping_loci_records/stratification_w6_w8/run.sh OUTDIR
# OUTDIR must not exist; every output (data/, results/, logs/) is written there, never into this folder.
# Then compare with the shipped records:
#   bash results/lab/jumping_loci_records/stratification_w6_w8/compare.sh OUTDIR
#
# Environment: PYTHON (default python3; Python >= 3.9 with sympy), MEM_KB (ulimit -v per step, default 8000000),
# TIMEOUT_S (timeout per step in seconds, default 900). SYMPY_GROUND_TYPES defaults to "python", the setting of the
# archived run.
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "usage: run.sh OUTDIR   (OUTDIR must not exist yet)" >&2
  exit 2
fi
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
PYTHON="${PYTHON:-python3}"
MEM_KB="${MEM_KB:-8000000}"
TIMEOUT_S="${TIMEOUT_S:-900}"
export SYMPY_GROUND_TYPES="${SYMPY_GROUND_TYPES:-python}"
export PYTHONDONTWRITEBYTECODE=1

for f in code/mixed_engine.py code/invariant_engine.py; do
  if [ ! -f "$REPO/$f" ]; then
    echo "run.sh: $f not found under the repository root $REPO" >&2
    exit 2
  fi
done
if [ -e "$1" ]; then
  echo "run.sh: $1 exists; give a new output directory" >&2
  exit 2
fi
mkdir "$1"                    # its parent directory must exist
OUT="$(cd "$1" && pwd)"
case "$OUT/" in
  "$HERE"/*) rmdir "$OUT"; echo "run.sh: the output directory must lie outside $HERE" >&2; exit 2 ;;
esac
cd "$OUT"                     # the scripts write data/ and results/ under the current directory
mkdir -p data results logs

# capped CMD...: run CMD under ulimit -v and timeout
capped() {
  ( ulimit -v "$MEM_KB"; exec timeout "$TIMEOUT_S" "$@" )
}
# step LOG CMD...: run CMD capped, appending its standard output to logs/LOG
step() {
  local log="$1"
  shift
  echo "  $(date '+%H:%M:%S') $* -> logs/$log"
  capped "$@" >> "logs/$log"
}
# verdict LOG: the last line of logs/LOG must be "ALL OK"
verdict() {
  if [ "$(tail -n 1 "logs/$1")" != "ALL OK" ]; then
    echo "run.sh: logs/$1 does not end with ALL OK" >&2
    exit 1
  fi
}
PY=("$PYTHON" -B)
C="$HERE/code"

capped "${PY[@]}" -c 'import sys, sympy
from sympy.external.gmpy import GROUND_TYPES
print("python", sys.version.split()[0]); print("sympy", sympy.__version__, "ground types", GROUND_TYPES)' > versions.txt
echo "run.sh: output in $OUT; $(tr '\n' ' ' < versions.txt)"

# basis sizes (engine objects) at weights 1-11 (the engine calls the weight "spin")
step 03_sizes.txt "${PY[@]}" "$C/sizes.py"
# S1: the family matrices M_poly(t) and their validation against the engine at 42 points
step 10_s1_build.txt "${PY[@]}" "$C/s1_build.py"
step 11_s1_validate.txt "${PY[@]}" "$C/s1_validate.py"
verdict 11_s1_validate.txt
# S2 + S3: rank over Q(t), the gcd G of the computed maximal minors, exact decisions (seed 68100 + 10*sheet + W)
for W in 6 8; do
  if [ "$W" = 6 ]; then log="20_s2s3_W6.txt"; else log="21_s2s3_W8.txt"; fi
  for sheet in 1 2 3; do
    step "$log" "${PY[@]}" "$C/s2s3.py" "data/M_sheet${sheet}_W${W}.json" --seed "$((68100 + 10 * sheet + W))" \
      --out "results/S23_sheet${sheet}_W${W}.json"
  done
done
# negative control (sheet 3, W = 8: one I_3 coefficient multiplied by 11/10), run when it is shipped
if [ -f "$C/s4_tamper.py" ]; then
  step 30_s4_tamper.txt "${PY[@]}" "$C/s4_tamper.py"
  step 30_s4_tamper.txt "${PY[@]}" "$C/s2s3.py" data/M_sheet3_W8_tamper.json --seed 69138 \
    --out results/S23_sheet3_W8_tamper.json
  mkdir control
  # (a) S2 + S3 on the tampered matrix give another generic rank or another G
  capped "${PY[@]}" - <<'EOF' | tee control/s2s3_check.txt
import json
a = json.load(open("results/S23_sheet3_W8.json"))
b = json.load(open("results/S23_sheet3_W8_tamper.json"))
new = sorted({d["factor"] for d in b["decisions"]} - {d["factor"] for d in a["decisions"]})
print(f"negative control (a), sheet 3, W = 8: untampered r = {a['r']}, deg G = {a['G_degree']}; "
      f"tampered r = {b['r']}, deg G = {b['G_degree']}, new factors of G: {new}")
ok = a["r"] != b["r"] or a["G"] != b["G"]
print("TAMPER DETECTED" if ok else "TAMPER NOT DETECTED")
raise SystemExit(0 if ok else 1)
EOF
  # (b) the S1 validation against the engine, run on the six matrices with sheet 3, W = 8 replaced by the tampered
  #     one (in control/), fails at that matrix and only there
  mkdir control/data
  cp data/M_sheet[123]_W[68].json control/data/
  cp data/M_sheet3_W8_tamper.json control/data/M_sheet3_W8.json
  ( cd control && capped "${PY[@]}" "$C/s1_validate.py" ) > control/11_s1_validate_tampered.txt
  awk '$NF == "FAIL" { n++; if (!($1 == 3 && $2 == 8)) bad++ } $NF == "OK" && $1 == 3 && $2 == 8 { bad++ }
       END { printf "negative control (b), S1 validation with the tampered matrix: %d FAIL lines, %d unexpected lines\n%s\n",
             n, bad, (n == 7 && bad == 0) ? "TAMPER DETECTED" : "TAMPER NOT DETECTED"; exit !(n == 7 && bad == 0) }' \
    control/11_s1_validate_tampered.txt | tee control/s1_validate_check.txt
fi
# S3 checks: direct nullities and the factors of one minor; cross-check over number fields; engine at named points
step 31_s3_extra.txt "${PY[@]}" "$C/s3_extra.py"
step 33_s3_crosscheck.txt "${PY[@]}" "$C/s3_crosscheck.py"
verdict 33_s3_crosscheck.txt
step 34_s3_points.txt "${PY[@]}" "$C/s3_points.py"
verdict 34_s3_points.txt
echo "run.sh: done; now run compare.sh $1"
