#!/bin/bash
rm out.txt
for ((i = 20; i <= 100; i += 10)); do
    ./compare-script.sh -b 0 -j 1 -c 4 -s 1000 -t 10 -u "$i"
done