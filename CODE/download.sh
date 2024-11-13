esearch -db nucleotide -query "lumpy skin disease virus [organism] AND complete genome [title]"|efetch -format genbank >t1.gb
esearch -db nucleotide -query "lumpy skin disease virus [organism] AND genomic sequence [title]"|efetch -format genbank >t2.gb

esearch -db nucleotide -query "lumpy skin disease virus [organism] AND complete genome [title]"|efetch -format fasta >t1.fa
esearch -db nucleotide -query "lumpy skin disease virus [organism] AND genomic sequence [title]"|efetch -format fasta >t2.fa
#seqret -sequence t.gb -outseq stdout|awk '{$2="";print $0}' >tt
#seqret -sequence t.gb -outseq stdout >tt
#
python parseGB.py t.gb|sed -e "s%complete genome\|Lumpy skin disease virus \|Lumpy_skin_disease_virus\|LSDV\|isolate\|genomic sequence\|_NULL\|strain\|LSD%%g"|tr -d "',)(:;\""|sed -e "s%/\| \|-%_%g" |sed -e "s%[_]+%_%g"|tr -s _|sed -e "s/_$//" >t
