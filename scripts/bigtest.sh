#!/bin/bash
rm out.txt
for ((i = 10; i <= 100; i += 10)); do
    ./compare-script.sh -b 0 -j 1 -c 4 -s 10000 -t 10 -u "$i"
done