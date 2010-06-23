#!/bin/sh

time ../../src/bin/mmc -f sph1.inp -s sph1 -n 30000000 -b 0 -l
time ../../src/bin/mmc -f sph2.inp -s sph2 -n 30000000 -b 0 -l
time ../../src/bin/mmc -f sph3.inp -s sph3 -n 30000000 -b 0 -l

