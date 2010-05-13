#!/bin/sh

time ../../src/bin/mmc -f sph1.inp -s sph1 -D 128 -n 30000000 -l
time ../../src/bin/mmc -f sph2.inp -s sph2 -D 128 -n 30000000 -l
time ../../src/bin/mmc -f sph3.inp -s sph3 -D 128 -n 30000000 -l

