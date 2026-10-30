#!/usr/bin/env bash
# Compares the logs written by run.sh (in the directory given as $1) with the shipped logs in logs/.
# Prints one line per log (IDENTICAL or DIFFERS) and a final line "ALL IDENTICAL" or "DIFFERENCES".
#
# Mask: the first line of every log is "script sha256 <64 hex digits>", the hash of the script file that wrote it.
# The shipped scripts differ from the ones that wrote the shipped logs in their first docstring line only, so this
# line is replaced by "script sha256 <masked>" in both logs before they are compared. In its place, the first line
# of each new log must equal "script sha256 " followed by the sha256 of the shipped script in code/.
# Every other line, including the final "exit 0", must match byte for byte.
# For a log that differs, the masked diff is written next to the new log as <name>.masked.diff.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STEPS="b5_fusion b5_vacuity"

if [ "$#" -ne 1 ] || [ ! -d "$1" ]; then
    echo "usage: $0 <output-dir-of-run.sh>" >&2; exit 2
fi
NEW="$(cd "$1" && pwd)"

sha() {
    if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1
    else shasum -a 256 "$1" | cut -d' ' -f1; fi
}
mask() { sed -E '1s/^script sha256 [0-9a-f]{64}$/script sha256 <masked>/' "$1"; }

ALL=1
for name in $STEPS; do
    shipped="$HERE/logs/$name.log"
    new="$NEW/$name.log"
    verdict=DIFFERS
    note=""
    if [ ! -f "$new" ]; then
        note=" (no new log)"
    elif [ "$(head -n 1 "$new")" != "script sha256 $(sha "$HERE/code/$name.py")" ]; then
        note=" (first line is not the sha256 of code/$name.py)"
    elif cmp -s <(mask "$shipped") <(mask "$new"); then
        verdict=IDENTICAL
    else
        diff <(mask "$shipped") <(mask "$new") > "$NEW/$name.masked.diff" || true
        note=" (masked diff: $NEW/$name.masked.diff)"
    fi
    echo "$verdict $name.log$note"
    [ "$verdict" = IDENTICAL ] || ALL=0
done
if [ "$ALL" -eq 1 ]; then
    echo "ALL IDENTICAL"
else
    echo "DIFFERENCES"; exit 1
fi
