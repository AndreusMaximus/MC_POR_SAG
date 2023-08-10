#!/bin/bash

for ((u = 10; u <= 80; u += 10)); do
echo "running $u"
    for ((j = 0; j <= 199; j += 1)); do
        ./build/nptest -m 1 ../SC/task/task\ sets/$u/jobsets/jobset-taskset-$j.csv --por=priority  >> ../SC/task/task\ sets/$u/result-POR.txt
        ./build/nptest -m 1 ../SC/task/task\ sets/$u/jobsets/jobset-taskset-$j.csv --por=priority --limit_fail=yes  >> ../SC/task/task\ sets/$u/result-limit.txt
        ./build/nptest -m 1 ../SC/task/task\ sets/$u/jobsets/jobset-taskset-$j.csv --por=priority  --interfering=all >> ../SC/task/task\ sets/$u/result-interfering.txt
        ./build/nptest -m 1 ../SC/task/task\ sets/$u/jobsets/jobset-taskset-$j.csv --por=priority --limit_fail=yes --interfering=all >> ../SC/task/task\ sets/$u/result-limit-interfering.txt
        ./build/nptest -m 1 ../SC/task/task\ sets/$u/jobsets/jobset-taskset-$j.csv  >> ../SC/task/task\ sets/$u/result-original.txt
   done
done