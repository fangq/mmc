#!/bin/sh

time ../../src/bin/mmcl -f brain.inp -s brain -n 1e6 -b 0 -D TP -F bin $@
