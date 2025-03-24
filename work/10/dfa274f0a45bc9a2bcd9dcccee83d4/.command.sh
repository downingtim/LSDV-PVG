#!/bin/bash -ue
# script to create initial input to visualise with R
Ref=$(cat path1)
for ((N=300; N<=7300; N+=2)) 
# genome ends skipped, even numbers only
do
    gfautil --quiet -t 35 -i pggb.gfa snps --ref $Ref  --snps $N
done| sort -nk 3 | grep -v \: |grep -v path > variation_map.txt

plot_variation_map.R || true # plot image of variation map in PDF

# plot SNP density across genome
plot_SNP_density.R gfavariants.vcf || true

# plot AFS (allele freq spectrum) - input VCF and number of samples
afs.py out 11
