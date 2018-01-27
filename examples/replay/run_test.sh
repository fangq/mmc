#!/bin/sh

# first run to get the seeds for each detected photon (-q 1)
../../src/bin/mmc -f replaytest.inp -s step1 -n 3e6 -b 0 -q 1 -x 1 -d 1 -D TP

# replay the detected photons (-E *.mch) detected by detector#1 (-P 1) without normalization (-U 0)
../../src/bin/mmc -f replaytest.inp -E step1.mch -s step2 -P 1 -b 0 -q 1 -x 1 -d 1 -O L -D TP

# replay the detected photons (-E *.mch) detected by detector#1 (-P 1) without normalization (-U 0)
../../src/bin/mmc -f replaytest.inp -E step1.mch -s step3 -P 1 -b 0 -q 1 -x 1 -d 1 -O P -D TP

