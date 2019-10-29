#!/bin/sh

time ../../src/bin/mmcl -f mesh0.inp -s mesh0 -n 30000000 -b 0 -l -D TP
time ../../src/bin/mmcl -f mesh1.inp -s mesh1 -n 30000000 -b 0 -l -D TP
time ../../src/bin/mmcl -f mesh2.inp -s mesh2 -n 30000000 -b 0 -l -D TP

