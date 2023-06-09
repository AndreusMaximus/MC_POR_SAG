#!/bin/bash

cores="1"

util="90"

set_size="100"

generate_jobsets=0

set_tasks="1"

build=0
ex_name="temp"
while getopts "b:c:s:u:j:t:e:" opt; do
    case ${opt} in
        b )
        build="${OPTARG}"
        ;;
        j)
        generate_jobsets="${OPTARG}"
        ;;
        c)
        cores="${OPTARG}"
        ;;
        u)
        util="${OPTARG}"
        ;;
        s)
        set_size="${OPTARG}"
        ;;
        t)
        set_tasks="${OPTARG}"
        ;;
        e)
        ex_name="${OPTARG}"
        ;;

    esac
done

if [ $generate_jobsets -eq 1 ]; then
    echo "(1) Making test sets"
    cd ~/Downloads/real-time-task-generators-main/
    mkdir jobsets
    cd jobsets
    rm *
    cd ..
    mkdir segmentsets
    cd segmentsets
    rm *
    cd ..
    mkdir tasksets
    cd tasksets
    rm *
    cd ..
    ./group_generator.sh -s $set_size -c $cores -u $util -p 2 -g 2 -n $set_tasks
fi

if [ $build -eq 1 ]; then
    echo "(2) BUILDING PHASE"
    cd ~/Downloads/np-schedulability-analysis-partial_order_reduction/build
    cmake --fresh ..
    make clean
    make -j
 fi

echo "(3) RUNNING COMPARISONS"
cd ~/Downloads/np-schedulability-analysis-partial_order_reduction/scripts
truncate -s 0 result.txt
for i in $(seq 0 $(expr $set_size - 1 ));
do
    lc=$(wc -l < ~/Downloads/real-time-task-generators-main/jobsets/jobset-taskset-$i.csv)
    if [ "$lc" -lt "100000" ]; then
        ../build/nptest -r -m $cores ~/Downloads/real-time-task-generators-main/jobsets/jobset-taskset-$i.csv | tee -a result.txt 
        ../build/nptest -r -m $cores ~/Downloads/real-time-task-generators-main/jobsets/jobset-taskset-$i.csv --por=priority  | tee -a result.txt 
        ../build/nptest -r -m $cores ~/Downloads/real-time-task-generators-main/jobsets/jobset-taskset-$i.csv --por=priority --interfering=all | tee -a result.txt 

    fi
done

#echo "(4) PARSING RESULTS"
#python3 parse-results.py $cores $util
#echo "|| For $cores cores with a utilization of $util"