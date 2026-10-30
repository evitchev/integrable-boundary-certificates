#!/bin/bash
# compare_g2.sh [RECORD]: RECORD (default G2's record out/w22_t1311_sheet3.json; another path only for drills) against the
# repository's results/lab/fast_kernel3/w22_t1311_sheet3.json.
# Fields side by side (tool, tool_sha256, weight, solution, N, s, genuine_dim, kc, primes_used, seconds, status); the witnesses
# compared in canonical form (jq -S -c '.witnesses': sorted keys), by sha256 and size; a control that the comparison can fail
# (the repository record's two witnesses, compared with each other, must differ).  No coefficient is printed.  Exit 0 iff
# every jq call succeeds; the verdict lines are WITNESSES EQUAL / WITNESSES DIFFER and FIELD DIFFERS: <key>.
set -u; set -o pipefail
W=$(cd "$(dirname "$0")/.." && pwd) || exit 1                 # this folder (code/..)
A=${1:-$W/out/w22_t1311_sheet3.json}
B=$(cd "$W/../../../.." && pwd)/results/lab/fast_kernel3/w22_t1311_sheet3.json || exit 1   # the repository root is ../../../.. from this folder
EMPTY=e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
canon() { local h; h=$(jq -S -c "$2" "$1" | sha256sum | cut -d' ' -f1) || { echo "jq failed: $1 $2" >&2; exit 1; }
          [ "$h" != "$EMPTY" ] || { echo "empty output: $1 $2" >&2; exit 1; }; echo "$h"; }
for f in "$A" "$B"; do echo "$f: sha256 $(sha256sum "$f" | cut -d' ' -f1), $(stat -c %s "$f") bytes"; done
for k in tool tool_sha256 weight solution N s genuine_dim kc primes_used seconds status; do
  echo "  $k: mine $(jq -c ".$k" "$A") | repository $(jq -c ".$k" "$B")" || exit 1
done
for k in weight solution N s genuine_dim kc primes_used; do
  [ "$(jq -c ".$k" "$A")" = "$(jq -c ".$k" "$B")" ] || echo "FIELD DIFFERS: $k"
done
echo "  witness counts: mine $(jq '.witnesses | length' "$A") | repository $(jq '.witnesses | length' "$B")"
ha=$(canon "$A" '.witnesses') || exit 1; hb=$(canon "$B" '.witnesses') || exit 1
echo "  witnesses, canonical form: mine sha256 $ha ($(jq -S -c '.witnesses' "$A" | wc -c) bytes) | repository sha256 $hb ($(jq -S -c '.witnesses' "$B" | wc -c) bytes)"
h0=$(canon "$B" '.witnesses[0]') || exit 1; h1=$(canon "$B" '.witnesses[1]') || exit 1
[ "$h0" != "$h1" ] && echo "  control: the repository record's witnesses 0 and 1 differ (${h0:0:16} vs ${h1:0:16}), so the comparison can fail" \
                   || { echo "  control FAILED: witnesses 0 and 1 compare equal"; exit 1; }
[ "$ha" = "$hb" ] && echo "WITNESSES EQUAL" || echo "WITNESSES DIFFER"
