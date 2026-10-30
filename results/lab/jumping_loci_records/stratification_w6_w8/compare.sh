#!/usr/bin/env bash
# Compare the outputs of run.sh (directory OUTDIR) with the records shipped in this folder.
#
# Usage, from the repository root:
#   bash results/lab/jumping_loci_records/stratification_w6_w8/compare.sh OUTDIR
# One line per shipped record (IDENTICAL / DIFFERS), then "ALL IDENTICAL" or "DIFFERENCES".
# data/*.json: compared byte for byte.
# results/*.json: compared after the mask  "seconds": <number>  ->  "seconds": MASKED  (the run time of each run).
# logs/*.txt: compared after the mask  (<number> s)  ->  (MASKED s)  (run times printed by the scripts).
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "usage: compare.sh OUTDIR" >&2
  exit 2
fi
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ ! -d "$1" ]; then
  echo "compare.sh: $1 is not a directory" >&2
  exit 2
fi
NEW="$(cd "$1" && pwd)"
cd "$HERE"                    # the shipped records are listed relative to this folder

mask_results() { sed -E 's/"seconds": [0-9]+(\.[0-9]+)?/"seconds": MASKED/' "$1"; }
mask_logs() { sed -E 's/\([0-9]+(\.[0-9]+)? s\)/(MASKED s)/g' "$1"; }

ndiff=0
nfiles=0
for rel in data/*.json results/*.json logs/*.txt; do
  ref="$rel"
  [ -f "$ref" ] || continue
  nfiles=$((nfiles + 1))
  new="$NEW/$rel"
  same=0
  if [ -f "$new" ]; then
    case "$rel" in
      data/*) cmp -s "$ref" "$new" && same=1 ;;
      results/*) cmp -s <(mask_results "$ref") <(mask_results "$new") && same=1 ;;
      logs/*) cmp -s <(mask_logs "$ref") <(mask_logs "$new") && same=1 ;;
    esac
  fi
  if [ "$same" = 1 ]; then
    echo "IDENTICAL $rel"
  else
    [ -f "$new" ] || rel="$rel (missing in the output directory)"
    echo "DIFFERS   $rel"
    ndiff=$((ndiff + 1))
  fi
done
if [ "$nfiles" -gt 0 ] && [ "$ndiff" = 0 ]; then
  echo "ALL IDENTICAL"
else
  echo "DIFFERENCES"
  exit 1
fi
