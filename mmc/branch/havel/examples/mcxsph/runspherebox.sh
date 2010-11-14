#!/bin/sh

# you need to fist add the path to mcx to PATH environment variable

#time mcx -t 5120 -T 256 -g 50 -m 1000000 -f spherebox.inp -s spherebox -r 6 -a 0 -b 0
#time mcx -t 40320 -T 576 -g 50 -m 1429954 -f spherebox.inp -s spherebox -a 0 -b 0 -G 1
#time mcx -t 4032 -T 576 -g 50 -m  100000 -f spherebox.inp -s spherebox -a 0 -b 0 -G 1

# please use MCX v0.4.9 or later
time mcx -A -g 50 -n 3e8 -f spherebox.inp -s spherebox -a 0 -b 0 -G 1

