#!/bin/sh

time ../../bin/mmc -f mesh0.inp -s mesh0 -n 30000000 -b 0 -l -D TP
time ../../bin/mmc -f mesh1.inp -s mesh1 -n 30000000 -b 0 -l -D TP
time ../../bin/mmc -f mesh2.inp -s mesh2 -n 30000000 -b 0 -l -D TP

