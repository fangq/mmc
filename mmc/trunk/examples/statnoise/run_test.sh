#!/bin/bash

for i in {1..10}
do
    echo "run #$i"
    ../../src/bin/mmc -f vartest.json -s v_$i -n 10000000 -b 0 -l -D TP -E 9876543$i
done
