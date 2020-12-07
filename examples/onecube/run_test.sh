#!/bin/sh

../../src/bin/mmc -n 20 -f onecube.inp -s onecube -G -1 -D M  | grep '^M' | sed -e 's/^M/1/g' -e 's/^B/0/g' -e 's/P/2/g'| sed '$d' > mov.txt
../../src/bin/mmc -n 20 -f onecube.inp -s onecube -G -1 -D MA | grep '^[A-Z]' | sed -e 's/^[A-Z] //g' | sed '$d' > ad.txt
