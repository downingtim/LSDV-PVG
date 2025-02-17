#!/bin/bash -ue
# need to split files into individual files in folder SEQS
# wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faSplit and chmod it
mkdir SEQS -p
faSplit byname test_genomes.FMDV_WRL.A.fa SEQS/
for FILE in SEQS/*.fa 
do
# Extract sample name from file name
sample_name=$(basename "$FILE" | cut -f 1 -d '.')
# Execute ~/bin/fastix with sample name and append output to "COMMUNITIES/genomes.fasta"
fastix -p "${sample_name}#1#" "$FILE" >> communities.genomes.fasta
done

# run pangrowth
pangrowth hist -k 17 -t 12 SEQS/*a > hist.txt

# plot AFS for single dataset - not very useful! 
plot_hist.py hist.txt pangrowth.pdf

# model growth - prepare dataset 1 - not useful
pangrowth growth -h hist.txt > temp1
plot_growth.py temp1 growth.pdf

# do core - this does not converge to a solution for n<10 samples, causes error
pangrowth core -h hist.txt -q 0.5 > p_core

if [ 9 -gt 10 ]; then # if the core is too small, this crashes
   plot_core.py p_core p_core.pdf     # so need at least 10 genomes
fi
