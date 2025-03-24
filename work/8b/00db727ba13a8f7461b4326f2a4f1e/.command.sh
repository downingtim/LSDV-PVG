#!/bin/bash -ue
# visualise output, reading in heaps file, 3rd column only of interest
for N in {1..11}
do
   odgi heaps -i out.og -n 1000 -d $N -t 12 | sort -nk 2 | tail -n 1
done > heaps.txt
