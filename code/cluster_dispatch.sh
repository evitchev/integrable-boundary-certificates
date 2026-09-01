#!/usr/bin/env bash
# Dispatch single-point kernel jobs to worker nodes, round-robin.
#   cluster_dispatch.sh JOBFILE        lines: TOOL SOL N S  (blank/# ignored)
# Hosts: one per line in $IB_HOSTS (default ~/ib-cluster/hosts) -- an
# untracked file; hostnames never enter the repository.  Workers run a
# git-HEAD snapshot of code/ at ~/ib-worker/code with the python312 stack
# at the same path as on the main machine; logs at ~/ib-logs/<tag>.out.
set -euo pipefail
JOBFILE=${1:?jobfile}
HOSTS_FILE=${IB_HOSTS:-$HOME/ib-cluster/hosts}
mapfile -t HOSTS < <(grep -v '^\s*\(#\|$\)' "$HOSTS_FILE")
[ ${#HOSTS[@]} -gt 0 ] || { echo "no hosts in $HOSTS_FILE" >&2; exit 1; }
PY=${IB_PY:-$HOME/python312/bin/python3}
i=0
grep -v '^\s*\(#\|$\)' "$JOBFILE" | while read -r tool sol N S; do
  h=${HOSTS[$((i % ${#HOSTS[@]}))]}; i=$((i+1))
  tag="${tool%.py}_${sol}_$(echo "$N" | tr '/-' 'd_m')_$(echo "$S" | tr '/-' 'd_m')"
  # NB: the launch must be its own statement -- 'test && cmd &' would
  # background the whole and-list as a subshell holding the ssh pipe open
  # (this hung the first dispatch until the kernel finished).
  ssh -n -o BatchMode=yes "$h" "cd ~/ib-worker/code || exit 1; test -f $tool || exit 1
    setsid nohup $PY $tool $sol $N $S > ~/ib-logs/$tag.out 2>&1 < /dev/null &
    sleep 1; pgrep -f '$tool $sol $N $S' >/dev/null && echo 'launched'" \
    | sed "s/^/$h [$tag]: /"
done
