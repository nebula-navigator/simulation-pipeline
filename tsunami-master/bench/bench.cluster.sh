#!/bin/bash

source lib/benchfunc.sh

line="./tsunami.x $input/cluster_N200.dat -N 200 -L au -ft 1.316 -dt 1e2"
message="Cluster N=200 benchmark"
key="cluster"
clean=${1:-1}
updatebs=${2:-0}

tsunami_benchmark  "$line" "$key" "$message" "$clean" "$updatebs"