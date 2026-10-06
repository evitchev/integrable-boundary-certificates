#!/bin/bash
# EXC3 batch: one momentum point; u_one.sh <mode> <t> <a0> <a1>
export PYTHONDONTWRITEBYTECODE=1
TAG="t$(echo $2 | tr '/' '_')_$(echo $3 | tr '/' '-')_$(echo $4 | tr '/' '-')"
( ulimit -v 16000000; timeout 3000 <home>/python312-lab/bin/python3 h_two_sol2.py $1 $2 $3 $4 > u_one_${1}_$TAG.log 2>&1; echo "exit $?" >> u_one_${1}_$TAG.log )
grep -A1 "^vdim" u_one_${1}_$TAG.log | tail -1 | sed "s/^/$TAG vdim /"
