#!/bin/bash
truncate -s 0 result.txt
for i in {0..999}
do
    
    ../build/nptest -r -m 4 ~/Downloads/real-time-task-generators-main/jobsets/jobset-taskset-$i.csv >> result.txt
    ../build/nptest -r -m 4 ~/Downloads/real-time-task-generators-main/jobsets/jobset-taskset-$i.csv --por=priority  >> result.txt
done