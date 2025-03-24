#!/bin/bash -ue
#Splir refFasta by name i to different files
mkdir SEQS -p
awk '{
    if (substr($1, 1, 1) == ">") {
        NAME = substr($1, 2);
        split(NAME, arr, " ");
        NAME = arr[1];
        print ">" NAME > "SEQS/" NAME ".fa"
    } else {
        if ($1) {
            print $1 >> "SEQS/" NAME ".fa"
        }
    }
}' genomes.fasta

for FILE in SEQS/*.fa 
do
  # Extract sample name from file name
  sample_name=$(basename "$FILE" .fasta)
  # Execute ~/bin/fastix with sample name and append output
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

# if [ 11 -gt 8 ]; then # if the core is too small, this crashes
plot_core.py p_core p_core.pdf     # so need at least 10 genomes
# fi
