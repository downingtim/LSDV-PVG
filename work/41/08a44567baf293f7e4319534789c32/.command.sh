#!/bin/bash -ue
esearch -db nucleotide -query "human coxsackievirus a1 [organism] AND complete genome [title]"|efetch -format genbank > gb.1 
esearch -db nucleotide -query "human coxsackievirus a1 [organism] AND genomic sequence [title]"|efetch -format genbank > gb.2
cat gb.1 gb.2 > genbank.gb 
rm gb.1 gb.2
#parseGB.py genbank.gb|sed -e "s%complete genome\|human coxsackievirus a1\| human coxsackievirus\| coxsackievirus\|isolate\|genomic sequence\|_NULL\|strain%%g"|tr -d "',)(:;\""|sed -e "s%/\| \|-%_%g" |sed -e "s%[_]+%_%g"|tr -s _|sed -e "s/_$//"  > genomes.fasta
parseGB.py genbank.gb > genomes.fasta
ls -lt genomes.fasta genbank.gb
echo "File created with size: $(ls -la genomes.fasta)"
