awk '{if(substr($1,0,1)==">"){h=$0}else{if(length($0) >= 13850){print h"\n"substr($0, 1, 13850)} }}' bin/lumpy_2_skin_2_disease_2_virus_2_.fasta >t.5end
awk '{if(substr($1,0,1)==">"){h=$0}else{if(length($0) >= 106910){print h"\n"substr($0,13851,106910-13851+1)} }}' bin/lumpy_2_skin_2_disease_2_virus_2_.fasta >t.core
awk '{if(substr($1,0,1)==">"){h=$0}else{if(length($0) >= 106911){print h"\n"substr($0,106911)} }}' bin/lumpy_2_skin_2_disease_2_virus_2_.fasta >t.3end
