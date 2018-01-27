#!/bin/sh

../../src/bin/mmc -n 100000 -f planar.inp -s planar -b 1 -x 1 -D E > debug.txt
grep '^x ' debug.txt | sed -e 's/^x //g'> exitpos.txt
grep '^e ' debug.txt | sed -e 's/^e //g'> enterpos.txt

#../../src/bin/mmc -n 1 -f planar.inp -s planar -b 1 -x 1 -D M  | sed -e 's/^M/1/g' -e 's/^B/0/g' -e 's/P/2/g'| sed '$d' > mov.txt
#../../src/bin/mmc -n 1 -f planar.inp -s planar -b 1 -x 1 -D MA | sed -e 's/^[A-Z] //g' | sed '$d' > ad.txt

#../../src/bin/mmc -n 1000000 -f planar.inp -s planar -b 1 -x 1 -D TP
