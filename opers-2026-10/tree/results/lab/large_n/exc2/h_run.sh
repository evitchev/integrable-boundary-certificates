#!/bin/bash
# EXC2 two-point pipeline: h_run.sh <t> <a0> <a1> <b>
set -e
export PYTHONDONTWRITEBYTECODE=1
PY=<home>/python312-lab/bin/python3
T=$1; A0=$2; A1=$3; B=$4
TAG="t$(echo $T | tr '/' '_')_$(echo $A0 | tr '/' '-')_$(echo $A1 | tr '/' '-')"
$PY h_two_v2.py $T $A0 $A1 | tail -2
$PY h_reduce.py $TAG $B | tail -3
$PY h_reduce2.py $TAG | head -8
timeout 1500 Singular -q < h_red2_$TAG.sing > h_red2_$TAG.out 2>&1
head -c 1200 h_red2_$TAG.out
