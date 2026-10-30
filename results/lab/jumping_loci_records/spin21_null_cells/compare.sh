#!/bin/bash
# compare.sh NEWDIR: compares the outputs of run.sh in NEWDIR with the shipped records (out/) and logs (logs/) of this folder.
# One line per shipped file, IDENTICAL or DIFFERS, then ALL IDENTICAL or DIFFERENCES (exit status 0 or 1).
#
# What is compared (the mask; the same for old and new files):
#   out/*.json   canonical JSON (jq -S -c) without the field "seconds" (the run time); the field "tool_sha256" is replaced by
#                <TOOL> when it is one of the two known hashes of the kernel program: 07c3e3f2..., the program that wrote the
#                shipped records, and the sha256 of the shipped code/vendor/lab/fast_kernel3.py (the same program with four
#                comment lines changed).  Any other value is kept, and so differs.  Everything else, the kernel vectors
#                included, must be equal.
#   logs/*.log   line by line after these substitutions:
#                - every duration, a number directly followed by "s" (e.g. "in 126s", "total 1353s", "1.8s") -> <T>s;
#                - "record written: <path>/w22_X.json (tool sha256 <16 hex>)": the directory -> <OUT>, and the hash -> <TOOL>
#                  if it is the first 16 hex digits of one of the two known hashes above;
#                - "code/jl_lib.py is the pinned copy of the shared library (sha256 <16 hex>)": the hash -> <JL_LIB> if it is
#                  9ff3e9d0f544908b (the library of the original run) or the first 16 hex digits of the shipped code/jl_lib.py;
#                - in compare_g2.log: the directory before "w22_t1311_sheet3.json: sha256" -> <DIR>; on line 1 (the new
#                  record) its sha256 and size -> <SHA256>, <SIZE>; the new record's tool_sha256 -> <TOOL> as above; its
#                  seconds -> <T>.
# The masked copies of both sides are written to NEWDIR/compare_masked/ (<file>.shipped, <file>.new) for inspection.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd -P)
[ $# = 1 ] && [ -d "$1" ] || { echo "usage: compare.sh NEWDIR (the directory written by run.sh)"; exit 2; }
NEW=$(cd "$1" && pwd -P)
[ "$NEW" != "$HERE" ] || { echo "compare.sh: NEWDIR is this folder"; exit 2; }
command -v jq > /dev/null || { echo "compare.sh: jq is required"; exit 2; }

ORIG_TOOL=07c3e3f227a719703e8a2e6108083871f14f91d997aa4201bcb03f0d30d4846e
SHIP_TOOL=$(sha256sum "$HERE/code/vendor/lab/fast_kernel3.py" | cut -d' ' -f1)
ORIG_JL=9ff3e9d0f544908b
SHIP_JL=$(sha256sum "$HERE/code/jl_lib.py" | cut -c1-16)
[ ${#SHIP_TOOL} = 64 ] && [ ${#SHIP_JL} = 16 ] || { echo "compare.sh: cannot hash the shipped code"; exit 2; }

mask_json() {  # canonical JSON without "seconds", known tool hashes -> <TOOL>
  jq -S -c --arg o "$ORIG_TOOL" --arg s "$SHIP_TOOL" \
    'del(.seconds) | if (.tool_sha256 == $o or .tool_sha256 == $s) then .tool_sha256 = "<TOOL>" else . end' "$1"
}
mask_log() {
  sed -E \
    -e 's/\b[0-9]+(\.[0-9]+)?s\b/<T>s/g' \
    -e "s#^record written: (.*/)?(w22_[A-Za-z0-9_]+\\.json) \\(tool sha256 (${ORIG_TOOL:0:16}|${SHIP_TOOL:0:16})\\)\$#record written: <OUT>/\\2 (tool sha256 <TOOL>)#" \
    -e "s#^(PASS  code/jl_lib\\.py is the pinned copy of the shared library \\(sha256 )($ORIG_JL|$SHIP_JL)\\)\$#\\1<JL_LIB>)#" \
    -e 's#^.*/(w22_t1311_sheet3\.json): sha256 #<DIR>/\1: sha256 #' \
    -e '1s#^(<DIR>/w22_t1311_sheet3\.json): sha256 [0-9a-f]{64}, [0-9]+ bytes$#\1: sha256 <SHA256>, <SIZE> bytes#' \
    -e "s#^(  tool_sha256: mine )\"($ORIG_TOOL|$SHIP_TOOL)\" \\|#\\1\"<TOOL>\" |#" \
    -e 's#^(  seconds: mine )[0-9]+ \|#\1<T> |#' \
    "$1"
}

bad=0
report() {  # report STATUS FILE [NOTE]
  echo "$1  $2${3:+  ($3)}"
  [ "$1" = IDENTICAL ] || bad=1
}
M="$NEW/compare_masked"   # the masked copies: <file>.shipped and <file>.new
mkdir -p "$M"

for f in out/w22_t1311_sheet3.json out/w22_t115_sheet3.json out/w22_t2_sheet2.json \
         logs/run_g1.log logs/w22_t1311_sheet3.log logs/compare_g2.log logs/w22_t115_sheet3.log logs/w22_t2_sheet2.log; do
  [ -f "$NEW/$f" ] || { report DIFFERS "$f" "not in $NEW"; continue; }
  old="$M/$(basename "$f").shipped"; new="$M/$(basename "$f").new"
  case "$f" in
    *.json) mask_json "$HERE/$f" > "$old" || { report DIFFERS "$f" "the shipped file does not parse"; continue; }
            mask_json "$NEW/$f" > "$new" 2> /dev/null || { report DIFFERS "$f" "not valid JSON"; continue; };;
    *.log)  mask_log "$HERE/$f" > "$old"; mask_log "$NEW/$f" > "$new";;
  esac
  if cmp -s "$old" "$new"; then report IDENTICAL "$f"; else report DIFFERS "$f"; fi
done
if [ $bad = 0 ]; then echo "ALL IDENTICAL"; else echo "DIFFERENCES"; fi
exit $bad
