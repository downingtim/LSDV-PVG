#!/bin/bash -ue
mafft  --thread 8 --auto genomes.fasta > genomes.aln 
raxml-ng  --all --msa genomes.aln --model GTR+G4 --prefix T14 --seed 21231 --bs-metric fbp,tbe --redo
