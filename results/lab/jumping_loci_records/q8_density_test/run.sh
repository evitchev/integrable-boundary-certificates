#!/usr/bin/env bash
# Re-runs the q = 8 density-level test from the release tree, writing into a NEW output directory.
#
#   bash run.sh OUTDIR        (OUTDIR must not exist yet and must lie outside this folder; it is created)
#
# Interpreter: $PYTHON (default python3), called as "$PYTHON -I -B". Each step runs under ulimit -v (default 8000000 KiB,
# override with Q8_ULIMIT_KB) and a timeout. This folder must sit at results/lab/jumping_loci_records/q8_density_test/
# inside a clone of the repository (the repository root is ../../../.. from here). The five steps run one after another.
set -euo pipefail

if [ "$#" -ne 1 ]; then echo "usage: bash run.sh OUTDIR   (a new directory)" >&2; exit 2; fi
HERE=$(cd "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
ROOT=$(cd "$HERE/../../../.." && pwd -P)
PYTHON=${PYTHON:-python3}
MEM_KB=${Q8_ULIMIT_KB:-8000000}

# the output directory: new, and never inside this folder (the shipped records stay untouched)
if [ -e "$1" ]; then echo "run.sh: $1 already exists; give a new directory" >&2; exit 2; fi
PARENT=$(cd "$(dirname -- "$1")" 2>/dev/null && pwd -P) || { echo "run.sh: the parent directory of $1 does not exist" >&2; exit 2; }
OUT="$PARENT/$(basename -- "$1")"
case "$OUT/" in
  "$HERE/"*) echo "run.sh: the output directory must lie outside $HERE" >&2; exit 2 ;;
esac

# the inputs: three files of the repository, the vendored engine, the shipped cross-check input
echo "run.sh: checking the inputs (sha256)"
( cd "$ROOT" && sha256sum -c --quiet - ) <<'EOF'
6dd97090651f54574ee541339ffcdbfd184756ddd5ce78b840439411e6c0caf0  results/lab/first_order/vev_w4_points.json
51eecbbcfb14e32c67e7f031bb9bba65d3e0bbfe232a552ddf0917559f5e1f35  results/lab/vev/vev_sol3_w10.json
be12400d2584e9f27b57e702e5fb6aa352b5f574d2bb5a010cf6cc3a3f63dbbb  results/lab/anchor13/vev_profile_sol3_w14.json
EOF
( cd "$HERE/code/vendor" && sha256sum -c --quiet SOURCE_SHA256.txt )
( cd "$HERE/inputs" && sha256sum -c --quiet SOURCE_SHA256.txt )

mkdir "$OUT"
echo "run.sh: writing to $OUT"

FAILED=0
step() {   # step NAME TIMEOUT_SECONDS EXPECTED_EXIT [ARGS...]
  local name=$1 secs=$2 want=$3 rc=0
  shift 3
  echo "run.sh: $name (recorded exit $want, cap ${secs}s)"
  ( ulimit -v "$MEM_KB"; exec timeout --signal=TERM --kill-after=60s "$secs" "$PYTHON" -I -B "$HERE/code/$name.py" "$@" ) \
    > "$OUT/$name.log" 2>&1 || rc=$?
  echo "EXIT $rc" >> "$OUT/$name.log"
  if [ "$rc" -ne "$want" ]; then echo "run.sh: $name exited $rc, the record says $want" >&2; FAILED=1; fi
}

step s1_frame2   1800 0
step w0_engine   1800 0
step w1_gates    1800 0
step w2_bilinear 3600 1           # exit 1 as recorded: control C4 does not fire at the eigenvalue level (see README)
step w3_kernel   7200 0 "$OUT"    # also writes $OUT/w3_kernel.json

if [ "$FAILED" -ne 0 ]; then echo "run.sh: at least one step did not exit as recorded" >&2; exit 1; fi
echo "run.sh: every step exited as recorded; next: bash $HERE/compare.sh $OUT"
