#!/bin/bash

source lib/benchfunc.sh

line="./tsunami.x $input/eccentric_binary.dat -N 2 -PN -L au -ft 1.445e3 -dt 1e4"
message="Inspiraling eccentric binary benchmark"
key="binarypn"
clean=${1:-1}
updatebs=${2:-0}

tsunami_benchmark  "$line" "$key" "$message" "$clean" "$updatebs"
