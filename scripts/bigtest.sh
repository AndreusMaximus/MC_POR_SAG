#!/bin/bash
rm out.txt
for ((i = 10; i <= 100; i += 10)); do
    ./compare-script.sh -b 0 -j 1 -c 8 -s 100 -t 15 -u "$i"
done