#!/usr/bin/env bash
# Compares the outputs of run.sh with the shipped records in records/.
#
#   bash compare.sh OUTDIR
#
# w3_kernel.json is compared byte for byte. Each .log is compared after the same mask on both sides:
#   (1) line 1, "<script>.py sha256 <64 hex digits>": the hash is replaced by <sha256>;
#   (2) run times: a number followed by "s", after "(" or a space and before ")", ";" or the end of the line, becomes <N>s.
# Prints one line per file (IDENTICAL / DIFFERS) and a final line ALL IDENTICAL or DIFFERENCES (exit status 0 or 1).
set -euo pipefail

if [ "$#" -ne 1 ]; then echo "usage: bash compare.sh OUTDIR" >&2; exit 2; fi
HERE=$(cd "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
NEW=$1

mask() {
  sed -E -e '1s/^([A-Za-z0-9_]+\.py sha256) [0-9a-f]{64}$/\1 <sha256>/' \
         -e 's/([( ])[0-9]+s([);]|$)/\1<N>s\2/g' "$1"
}

bad=0
for f in s1_frame2.log w0_engine.log w1_gates.log w2_bilinear.log w3_kernel.log w3_kernel.json; do
  if [ ! -f "$NEW/$f" ]; then echo "DIFFERS    $f (missing)"; bad=1; continue; fi
  case "$f" in
    *.log)
      if cmp -s <(mask "$HERE/records/$f") <(mask "$NEW/$f"); then echo "IDENTICAL  $f (masked)"; else echo "DIFFERS    $f (masked)"; bad=1; fi ;;
    *)
      if cmp -s "$HERE/records/$f" "$NEW/$f"; then echo "IDENTICAL  $f"; else echo "DIFFERS    $f"; bad=1; fi ;;
  esac
done
if [ "$bad" -eq 0 ]; then echo "ALL IDENTICAL"; else echo "DIFFERENCES"; exit 1; fi
