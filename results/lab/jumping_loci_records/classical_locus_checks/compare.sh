#!/bin/bash
# compare.sh OUTPUT_DIR -- compares the outputs of "run.sh OUTPUT_DIR" with the records shipped in this folder, byte for
# byte (cmp). Prints one line per file (IDENTICAL / DIFFERS) and a final line "ALL IDENTICAL" or "DIFFERENCES".
# No mask is applied: the scripts print no times, host names or absolute paths (run.sh keeps its times in
# OUTPUT_DIR/run_info.txt, which is not compared).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ "$#" -ne 1 ] || [ ! -d "$1" ]; then
  echo "usage: compare.sh OUTPUT_DIR (the directory given to run.sh)" >&2; exit 2
fi
OUT="$1"
ndiff=0
for f in logs/run_a2.log logs/run_a3.log logs/run_a4.log logs/run_a5.log logs/open_tool_crosscheck.log LEDGER.json LEDGER.md; do
  if [ -f "$OUT/$f" ] && cmp -s "$HERE/$f" "$OUT/$f"; then
    echo "IDENTICAL  $f"
  elif [ -f "$OUT/$f" ]; then
    echo "DIFFERS    $f"; ndiff=$((ndiff + 1))
  else
    echo "DIFFERS    $f (not in $OUT)"; ndiff=$((ndiff + 1))
  fi
done
if [ "$ndiff" -eq 0 ]; then
  echo "ALL IDENTICAL"
else
  echo "DIFFERENCES"; exit 1
fi
