#!/bin/bash
rm out.txt
rm ../test_output -r


t=5
u=(10 20 30 40 50 60 70 80 90 100)
c=2
mkdir ../test_output
echo "running tests for $t tasks"
mkdir ../test_output/t"$t"

echo "  using $c cores"
mkdir ../test_output/t"$t"/c"$c"

mkdir ../test_output/t"$t"/c"$c"/raw
echo "cores, utilization, tasks, correct, true postive, false postive, false negative, true avg state reductio, avg state reduction for true positives, true avg cpu reduction, avg cpu reduction for true positives, true mem reduction, avg mem reduction for true positives, total successfull reductions, total failed reductions, avg jobs in por, timeout normal, timeout por, sched % normal, sched % por" >> ../test_output/t"$t"/c"$c"/raw/out.txt

for utilization in "${u[@]}"; do
    echo "      with a utilization of $utilization"
    ./compare-script.sh -t $t -c $c -u $(($utilization*$c)) -b 0 -j 1 -s 100

done