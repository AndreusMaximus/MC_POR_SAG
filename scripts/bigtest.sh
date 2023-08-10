#!/bin/bash
ex_name="$1"
c="4"
s=25
sag=1
base_por=1
base_interfering=1
limit_por=1
limit_interfering=1
task_scaler=8
make_set=0


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


mkdir ../../SAG_datasets/"$ex_name"
if [ $make_set -eq 1 ]; then
    rm -r ../../SAG_datasets/"$ex_name"/task-sets
fi
mkdir ../../SAG_datasets/"$ex_name"/task-sets
if [ $sag -eq 1 ]; then
    rm -r ../../SAG_datasets/"$ex_name"/sag/
    mkdir ../../SAG_datasets/"$ex_name"/sag/
fi
if [ $base_por -eq 1 ]; then
    rm -r ../../SAG_datasets/"$ex_name"/sag-por-base/
    mkdir ../../SAG_datasets/"$ex_name"/sag-por-base/
fi
if [ $base_interfering -eq 1 ]; then
    rm -r ../../SAG_datasets/"$ex_name"/sag-por-interfering/
    mkdir ../../SAG_datasets/"$ex_name"/sag-por-interfering/
fi
if [ $limit_por -eq 1 ]; then
    rm -r ../../SAG_datasets/"$ex_name"/sag-por-limit/
    mkdir ../../SAG_datasets/"$ex_name"/sag-por-limit/
fi
if [ $limit_interfering -eq 1 ]; then
    rm -r ../../SAG_datasets/"$ex_name"/sag-por-interfering-limit/
    mkdir ../../SAG_datasets/"$ex_name"/sag-por-interfering-limit/
fi


for ((u = 10; u <= 70; u += 10)); do

    if [ $sag -eq 1 ]; then
        mkdir ../../SAG_datasets/"$ex_name"/sag/$u
        mkdir ../../SAG_datasets/"$ex_name"/sag/$u/rta
        fi
    if [ $base_por -eq 1 ]; then
        mkdir ../../SAG_datasets/"$ex_name"/sag-por-base/$u
        mkdir ../../SAG_datasets/"$ex_name"/sag-por-base/$u/rta
        fi
    if [ $base_interfering -eq 1 ]; then
        mkdir ../../SAG_datasets/"$ex_name"/sag-por-interfering/$u
        mkdir ../../SAG_datasets/"$ex_name"/sag-por-interfering/$u/rta
    fi
    if [ $limit_por -eq 1 ]; then
        mkdir ../../SAG_datasets/"$ex_name"/sag-por-limit/$u
        mkdir ../../SAG_datasets/"$ex_name"/sag-por-limit/$u/rta
        fi
    if [ $limit_interfering -eq 1 ]; then
        mkdir ../../SAG_datasets/"$ex_name"/sag-por-interfering-limit/$u
        mkdir ../../SAG_datasets/"$ex_name"/sag-por-interfering-limit/$u/rta
    fi

    
    if [ $make_set -eq 1 ]; then
        mkdir ../../SAG_datasets/$ex_name/task-sets/$u
        ./create_sets.sh -c $c -t $(($task_scaler)) -u $(($u*$c)) -s $s -e ../SAG_datasets/"$ex_name" -g 3 -p 2
    fi

    cd ~/Downloads/np-schedulability-analysis-partial_order_reduction/scripts
    truncate -s 0 result.txt
    for i in $(seq 0 $(expr $s - 1 ));
    do
        lc=$(wc -l < ../../SAG_datasets/$ex_name/task-sets/$u/jobsets/jobset-taskset-$i.csv )
        if [ "$lc" -lt "100000" ]; then
            if (($c > 1)); then

                if [ $sag -eq 1 ]; then
                    ../build/nptest -r -m $c ../../SAG_datasets/$ex_name/task-sets/$u/jobsets/jobset-taskset-$i.csv | tee -a result.txt ../../SAG_datasets/$ex_name/sag/$u/results.csv
                    mv ../../SAG_datasets/$ex_name/task-sets/$u/jobsets/jobset-taskset-$i.rta.csv ../../SAG_datasets/$ex_name/sag/$u/rta
                fi
                if [ $base_por -eq 1 ]; then
                    ../build/nptest -r -m $c ../../SAG_datasets/$ex_name/task-sets/$u/jobsets/jobset-taskset-$i.csv  --por=priority  | tee -a result.txt ../../SAG_datasets/$ex_name/sag-por-base/$u/results.csv
                    mv ../../SAG_datasets/$ex_name/task-sets/$u/jobsets/jobset-taskset-$i.rta.csv ../../SAG_datasets/$ex_name/sag-por-base/$u/rta
                fi

                if [ $base_interfering -eq 1 ]; then
                    ../build/nptest -r -m $c ../../SAG_datasets/$ex_name/task-sets/$u/jobsets/jobset-taskset-$i.csv --por=priority  --interfering=all | tee -a result.txt ../../SAG_datasets/$ex_name/sag-por-interfering/$u/results.csv
                    mv ../../SAG_datasets/$ex_name/task-sets/$u/jobsets/jobset-taskset-$i.rta.csv ../../SAG_datasets/$ex_name/sag-por-interfering/$u/rta
                fi
                if [ $limit_por -eq 1 ]; then
                    ../build/nptest -r -m $c ../../SAG_datasets/$ex_name/task-sets/$u/jobsets/jobset-taskset-$i.csv  --por=priority --limit_fail=yes | tee -a result.txt ../../SAG_datasets/$ex_name/sag-por-limit/$u/results.csv
                    mv ../../SAG_datasets/$ex_name/task-sets/$u/jobsets/jobset-taskset-$i.rta.csv ../../SAG_datasets/$ex_name/sag-por-limit/$u/rta
                fi

                if [ $limit_interfering -eq 1 ]; then
                    ../build/nptest -r -m $c ../../SAG_datasets/$ex_name/task-sets/$u/jobsets/jobset-taskset-$i.csv --por=priority  --interfering=all --limit_fail=yes| tee -a result.txt ../../SAG_datasets/$ex_name/sag-por-interfering-limit/$u/results.csv
                    mv ../../SAG_datasets/$ex_name/task-sets/$u/jobsets/jobset-taskset-$i.rta.csv ../../SAG_datasets/$ex_name/sag-por-interfering-limit/$u/rta
                fi
            fi
        fi
    done
done