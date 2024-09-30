#!/bin/bash
#
#    Benchmark script
#    Inspired by NEMO bench5
#    https://teuben.github.io/nemo/man_html/bench.5.html
#
input=input
data=data
bin=../bin
time=/usr/bin/time

tsunami_benchmark () {
  line=$1
  key=$2
  message=$3
  clean=${4:-1}
  updatebs=${5:-0}

  # Copy exec files
  cp -r $bin/tsunami.x $bin/$input -t .

  [ ! -e $data/baseline.current.$key.log ] && cp $data/baseline.$key.log $data/baseline.current.$key.log

  echo $message
  tmp=bench.$key.$$
  mkdir $tmp

  echo $line
  echo "Running 5 times, ~5sec CPU for each on a Intel Core i7-7700HQ CPU @ 2.80GHz"

  for i in {0..4}
  do
      echo -n $i; $time $line > $tmp/bench.$key$i.log 2>&1
  done

  echo

  sed -n 's/^\(.*\)user.*system.*:\([0-9\.]*\)elapsed.*$/\1 \2/p' $tmp/bench.$key*.log > $tmp/bench.$key.cpu_times.log
  cpu_mean=$(awk 'BEGIN{s=0;}{s=s+$2;}END{print s/NR;}' $tmp/bench.$key.cpu_times.log)
  cat $tmp/bench.$key.cpu_times.log
  echo mean: $cpu_mean sec

  score=$(awk -v m=$cpu_mean 'BEGIN { print (5000 / m) }')

  echo "score          :"  $score

  if [ $updatebs == 1 ]; then
  	echo "Updating baseline files with current score (use updatebs=0 to prevent this)"
    echo $score > $data/baseline.current.$key.log
  else
    baseline=$(awk 'NR==1 {print; exit}' $data/baseline.current.$key.log)
    echo "baseline       :"  $baseline
    echo "(use updatebs=1 to update the baseline score in $data/baseline.$key.log)"
  fi

  if [ $clean == 1 ]; then
  	echo "Cleaning $tmp (use clean=0 to prevent this)"
  	rm -rf $tmp
  else
  	echo
  fi
}

