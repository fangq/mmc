#!/bin/sh

# run with -x 2 to record detected photons in the form of time-resolved images
../../src/bin/mmc -f widedet.inp -s widedet -n 1000000 -x 2 -b 1 -D TP
