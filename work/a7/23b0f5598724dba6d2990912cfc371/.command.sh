#!/bin/bash -ue
odgi flatten -i out.og -b out.flatten.tsv -t 22

# convert to BED
sed '1d' out.flatten.tsv | awk -v OFS='	' '{print($4,$2,$3,"step.rank_"$6,".",$5)}' > out.flatten.bed

# get BED file
odgi pav -i out.og  -b out.flatten.bed > out.flatten.pavs.tsv
grep -c "" out.flatten.pavs.tsv  > flatten.pavs.count.txt # gives all PAVs

# make a plot
