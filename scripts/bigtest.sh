#!/bin/bash
rm out.txt
for ((i = 10; i <= 100; i += 10)); do
    echo $(($i*4))
done