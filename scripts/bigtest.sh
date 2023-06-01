#!/bin/bash
c=4
rm -r ../../experiment_results/experiment1
mkdir ../../experiment_results/experiment1


for ((i = 10; i <= 80; i += 10)); do
    echo $(($i*$c))
    mkdir ../../experiment_results/experiment1/$i
    ./compare-script.sh -c $c -t $((2*$c)) -u $(($i*$c)) -j 1 -b 0 -s 10
    mv ../../real-time-task-generators-main/tasksets ../../experiment_results/experiment1/$i/
    mv ../../real-time-task-generators-main/jobsets ../../experiment_results/experiment1/$i/
    mv result.txt ../../experiment_results/experiment1/$i/
    mkdir ../../experiment_results/experiment1/$i/RTA
    mv ../../experiment_results/experiment1/$i/jobsets/*.rta.csv ../../experiment_results/experiment1/$i/RTA
done