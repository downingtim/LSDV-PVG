#!/bin/bash -ue
# create Blast indexes for self-self blast
# makeblastdb -dbtype nucl -in genomes.fasta -out /mnt/lustre/RDS-live/downing/LSDV-PVG/BLAST/fasta.db 
# blastall -p blastn -d /mnt/lustre/RDS-live/downing/LSDV-PVG/BLAST/fasta.db -i genomes.fasta -e 1.0 -outfmt 6 > self-blast.out

# SAMtools index
bgzip genomes.fasta
samtools faidx genomes.fasta.gz > genomes.fasta.gz.fai

# run PGGB - you need to specify the number of haplotypes in number
pggb -i genomes.fasta.gz -m -S -o . -t 20 -p 90 -s 1k -n 11
mv *.gfa pggb.gfa
# if issues, remove CURRENT/ in /mnt/lustre/RDS-live/downing/LSDV-PVG  before re-running 

# run pangenome in nf-core - don't work as reliably as PGGB, so best avoided
# nextflow run nf-core/pangenome --input genomes.fasta.gz --n_haplotypes number --outdir     #   /mnt/lustre/RDS-live/downing/LSDV-PVG/CURRENT  --communities  --wfmash_segment_length 1000
#odgi stats -m -i /mnt/lustre/RDS-live/downing/LSDV-PVG/CURRENT/FINAL_ODGI/genomes.fasta.gz.squeeze.og -S > /mnt/lustre/RDS-live/downing/LSDV-PVG/odgi.stats.txt
