#!/bin/bash

source lib/benchfunc.sh

line="./tsunami.x $input/solar_system.dat -N 10 -L au -ft 1.14e3 -dt 6.28e10"
message="Solar system benchmark"
key="solar"
clean=${1:-1}
updatebs=${2:-0}

tsunami_benchmark  "$line" "$key" "$message" "$clean" "$updatebs"