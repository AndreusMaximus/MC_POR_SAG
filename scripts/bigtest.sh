#!/bin/bash
ex_name="$1"
c="4"
s=10
normal=1
por=1
interfering=0
task_scaler=8
make_set=1


while getopts "n:p:i:e:j:" opt; do
    case ${opt} in
        n )
        normal="${OPTARG}"
        ;;
        p )
        por="${OPTARG}"
        ;;
        i )
        interfering="${OPTARG}"
        ;;
        e )
        ex_name="${OPTARG}"
        ;;
        j )
        make_set="${OPTARG}"
        ;;

    esac
done


rm -r ../../experiment_results/"$ex_name"
mkdir ../../experiment_results/"$ex_name"
mkdir ../../experiment_results/"$ex_name"/task-sets
if [ $normal -eq 1 ]; then
    mkdir ../../experiment_results/"$ex_name"/sag/
fi
if [ $por -eq 1 ]; then
    mkdir ../../experiment_results/"$ex_name"/sag-por/
fi
if [ $interfering -eq 1 ]; then
    mkdir ../../experiment_results/"$ex_name"/sag-por-bulk/
fi


for ((u = 10; u <= 100; u += 10)); do

    if [ $normal -eq 1 ]; then
        mkdir ../../experiment_results/"$ex_name"/sag/$u
        mkdir ../../experiment_results/"$ex_name"/sag/$u/rta
        fi
    if [ $por -eq 1 ]; then
        mkdir ../../experiment_results/"$ex_name"/sag-por/$u
        mkdir ../../experiment_results/"$ex_name"/sag-por/$u/rta
        fi
    if [ $interfering -eq 1 ]; then
        mkdir ../../experiment_results/"$ex_name"/sag-por-bulk/$u
        mkdir ../../experiment_results/"$ex_name"/sag-por-bulk/$u/rta
    fi

    mkdir ../../experiment_results/$ex_name/task-sets/$u
    if [ $make_set -eq 1 ]; then
        ./create_sets.sh -c $c -t $(($task_scaler)) -u $(($u*$c)) -s $s -e ../experiment_results/"$ex_name" -g 0 -p 2
    fi

    cd ~/Downloads/np-schedulability-analysis-partial_order_reduction/scripts
    truncate -s 0 result.txt
    for i in $(seq 0 $(expr $s - 1 ));
    do
        lc=$(wc -l < ../../experiment_results/$ex_name/task-sets/$u/jobsets/jobset-taskset-$i.csv )
        if [ "$lc" -lt "100000" ]; then
            if (($c > 1)); then

                if [ $normal -eq 1 ]; then
                    ../build/nptest -r -m $c ../../experiment_results/$ex_name/task-sets/$u/jobsets/jobset-taskset-$i.csv | tee -a result.txt ../../experiment_results/$ex_name/sag/$u/results.csv
                    mv ../../experiment_results/$ex_name/task-sets/$u/jobsets/jobset-taskset-$i.rta.csv ../../experiment_results/$ex_name/sag/$u/rta
                fi
                if [ $por -eq 1 ]; then
                    ../build/nptest -r -m $c ../../experiment_results/$ex_name/task-sets/$u/jobsets/jobset-taskset-$i.csv  --por=priority  | tee -a result.txt ../../experiment_results/$ex_name/sag-por/$u/results.csv
                    mv ../../experiment_results/$ex_name/task-sets/$u/jobsets/jobset-taskset-$i.rta.csv ../../experiment_results/$ex_name/sag-por/$u/rta
                fi

                if [ $interfering -eq 1 ]; then
                    ../build/nptest -r -m $c ../../experiment_results/$ex_name/task-sets/$u/jobsets/jobset-taskset-$i.csv --por=priority  --interfering=all | tee -a result.txt ../../experiment_results/$ex_name/sag-por-bulk/$u/results.csv
                    mv ../../experiment_results/$ex_name/task-sets/$u/jobsets/jobset-taskset-$i.rta.csv ../../experiment_results/$ex_name/sag-por-bulk/$u/rta
                fi

            
            fi
        fi
    done

done