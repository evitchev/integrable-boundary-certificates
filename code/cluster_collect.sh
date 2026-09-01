#!/usr/bin/env bash
# Print every worker's ~/ib-logs/*.out (job results) and running kernel jobs.
set -uo pipefail
HOSTS_FILE=${IB_HOSTS:-$HOME/ib-cluster/hosts}
for h in $(grep -v '^\s*\(#\|$\)' "$HOSTS_FILE"); do
  echo "== $h: $(timeout 20 ssh -n -o BatchMode=yes "$h" 'ps -o etime,args -C python3 | grep -c kdim' 2>/dev/null) kdim job(s) running"
  timeout 20 ssh -n -o BatchMode=yes "$h" 'for f in ~/ib-logs/*.out; do [ -f "$f" ] && echo "  $(basename $f): $(tail -1 $f)"; done' 2>/dev/null
done
