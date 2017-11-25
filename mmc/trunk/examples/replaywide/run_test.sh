#!/bin/sh

# first run to get the seeds for each detected photon (-q 1)
# ../../src/bin/mmc -f init_MC.inp -s init -n 3e7 -b 1 -q 1 -x 1 -M H -D TP

# replay the detected photons (-E *.mch) for full illumination & half detection (-P 0) under wl mode
../../src/bin/mmc -f replay1.inp -E init.mch -s replay1 -P 0 -b 1 -O L -M H -D TP

# replay the detected photons (-E *.mch) for half illumination & full detection (-P 0) under wp mode
../../src/bin/mmc -f replay2.inp -E init.mch -s replay2 -P 0 -b 1 -O L -M H -D TP

# replay the detected photons (-E *.mch) for half illumination & half detection (-P 0) under wl mode
../../src/bin/mmc -f replay3.inp -E init.mch -s replay3 -P 0 -b 1 -O L -M H -D TP

# replay the two detection patterns together
../../src/bin/mmc -f replay4.inp -E init.mch -s replay4 -P 0 -b 1 -O L -M H -D TP