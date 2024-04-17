#!/usr/bin/bash

rm -rf ../../../COMMUNITIES/genomes.fasta
rm -rf ../../../CURRENT/COMMUNITIES/genomes.fasta

# Iterate over files ending with ".fa" in the specified directory
ls ../../../CURRENT/PANGROWTH/SEQS/*.fa | while read -r f; do
    # Extract sample name from file name
    sample_name=$(basename "$f" | cut -f 1 -d '.')
    # Execute ~/bin/fastix with sample name and append output to "COMMUNITIES/genomes.fasta"
    ~/bin/fastix -p "${sample_name}#1#" "$f" >> ../../../CURRENT/COMMUNITIES/genomes.fasta
done

