#!/bin/sh

time ../../src/bin/mmcl -f cube2.inp -s cube2 -n 1e8 -b 0 -D TP -M G -F bin $@
