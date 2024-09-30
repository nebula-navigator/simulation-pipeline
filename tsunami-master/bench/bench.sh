#!/bin/bash
#
#    Benchmark script
#    Inspired by NEMO bench5
#    https://teuben.github.io/nemo/man_html/bench.5.html
#

# Default arguments
clean=${1:-1}
updatebs=${2:-0}

echo Running all benchmarks
echo
source bench.solar.sh $clean $updatebs
echo
source bench.cluster.sh $clean $updatebs
echo
source bench.binarypn.sh $clean $updatebs
