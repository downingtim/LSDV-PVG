#!/bin/bash -ue
# you need to install panancus
# mamba install -c conda-forge -c bioconda panacus
# get haplotypes
grep '^P' pggb.gfa | cut -f2 > haplotypes.txt

# run panacus to get data - just use defaults
RUST_LOG=info panacus histgrowth -t4 -l 1,2,1,1,1 -q 0,0,1,0.5,0.1 -S -s haplotypes.txt pggb.gfa > histgrowth.node.tsv

# visualise plot as PDF
panacus-visualize -e histgrowth.node.tsv > histgrowth.node.pdf
