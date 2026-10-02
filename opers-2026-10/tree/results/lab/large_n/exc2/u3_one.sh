#!/bin/bash
# EXC2 U3: one momentum point, self-adjoint two-point solutions.  u3_one.sh <t> <a0> <a1> <b>
export PYTHONDONTWRITEBYTECODE=1
PY=<home>/python312-lab/bin/python3
T=$1; A0=$2; A1=$3; B=$4
TAG="t$(echo $T | tr '/' '_')_$(echo $A0 | tr '/' '-')_$(echo $A1 | tr '/' '-')"
{ $PY h_two_v2.py $T $A0 $A1 && $PY h_reduce.py $TAG $B && $PY h_reduce2.py $TAG ; } > u3one_$TAG.log 2>&1 || { echo "$TAG PIPELINE FAILED"; exit 1; }
$PY u3_point.py $T $A0 $A1 2>&1 | tee -a u3one_$TAG.log | cut -c1-200
