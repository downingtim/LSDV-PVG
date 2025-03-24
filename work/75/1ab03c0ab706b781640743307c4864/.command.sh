#!/bin/bash -ue
# get general info on gfa
gfastats pggb.gfa > gfa.stats.txt

# get lengths of genomes etc
 lengths.py genomes.fasta
