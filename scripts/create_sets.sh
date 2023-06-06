#!/bin/bash

cores="1"

util="90"

set_size="100"

set_tasks="1"

ex_name="temp"
while getopts "c:s:u:t:e:" opt; do
    case ${opt} in
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
./group_generator.sh -s $set_size -c $cores -u $util -p 2 -g 0 -n $set_tasks


mv tasksets $ex_name/task-sets/$(($util/$cores))
mv jobsets $ex_name/task-sets/$(($util/$cores))