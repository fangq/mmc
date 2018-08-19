#!/bin/sh

time ../../src/bin/mmc -f cube.inp -s cube -n 30000000 -b 0 -D TP -M G -F bin $@
