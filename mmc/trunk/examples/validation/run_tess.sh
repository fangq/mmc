#!/bin/sh

time ../../src/bin/mmc -f cube.inp -s cube -n 30000000 -b 0 -D TP --outputdomain 1 -M S --atomic 0 -F bin
